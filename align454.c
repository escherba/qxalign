/*
 * =====================================================================================
 *
 *       Filename:  align454.c
 *
 *    Description:  Collection of routines for quality-aware alignment of Roche/454
 *                  reads. Functions with asw_* prefix implement asymmetric Smith-
 *                  Waterman-like algorithm with inverse scores (URL:
 *                  http://dx.doi.org/10.1101/gr.6468307)
 *
 *        Version:  1.0
 *        Created:  04/05/2011 22:05:57
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Eugene Scherba (es), escherba@bu.edu
 *        Company:  Laboratory for Biocomputing and Informatics, Boston University
 *
 * =====================================================================================
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>

#include "align454.h"

#define AMBIGUOUS_BASE 'N'

#define min(a,b) ((b)<(a)?(b):(a))
#define max(a,b) ((b)>(a)?(b):(a))

#define IS_MATCH(a,b) ((a) == (b) || (b) == AMBIGUOUS_BASE)

/**
 * Describing how CIGAR operation/length is packed in a 32-bit integer.
 */
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)

/*
  CIGAR operations.
 */
/*! @abstract CIGAR: match */
#define BAM_CMATCH      0
/*! @abstract CIGAR: insertion in the read/donor, deletion in reference */
#define BAM_CINS        1
/*! @abstract CIGAR: deletion in the read/donor, insertion in reference */
#define BAM_CDEL        2
/*! @abstract CIGAR: skip on the reference (e.g. spliced alignment) */
#define BAM_CREF_SKIP   3
/*! @abstract CIGAR: clip on the read with clipped sequence present in qseq */
#define BAM_CSOFT_CLIP  4
/*! @abstract CIGAR: clip on the read with clipped sequence trimmed off */
#define BAM_CHARD_CLIP  5
/*! @abstract CIGAR: padding */
#define BAM_CPAD        6
#define BAM_CSEQ_MATCH        7
#define BAM_CSEQ_MISMATCH 8
/*
 * CIGAR operation codes (http://samtools.sourceforge.net/SAM-1.3.pdf):
 *
 * M - match or mismatch
 * I - insertion in the read/donor, deletion in reference
 * D - deletion in the read/donor, insertion in reference
 * N - skipped region from reference
 * S - soft-clip in the read
 * H - hard-clip in the read
 * P - padding (silent deletion from padded reference)
 * = - match
 * X - mismatch
 *                           0    1    2    3    4    5    6    7    8 */
const char cigar_chars[] = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'};

/*
 * calculate number of digits (plus sign) needed to show a number in decimal format
 */
int ndigits(int i);

int ndigits(int i) {
    int n = i < 0 ? 2 : 1;
    while (i > 9) {
        n++;
        i /= 10;
    }
    return n;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_new
 *  Description:  Allocate and initialize an Alignment struct
 * =====================================================================================
 */
Alignment_ASW* asw_new(int MATCH_PEN, int MISMATCH_PEN, int GAP_OPEN_EXTEND, int GAP_EXTEND)
{
        /* use stdlib functions */
        return asw_init(
                asw_alloc(malloc, realloc, free),
                MATCH_PEN,
                MISMATCH_PEN,
                GAP_OPEN_EXTEND,
                GAP_EXTEND);
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_alloc
 *  Description:  Allocate an Alignment struct and set its members to zero
 * =====================================================================================
 */
Alignment_ASW* asw_alloc(
        void *(*p_malloc)(size_t size),
        void *(*p_realloc)(void * ptr, size_t size),
        void (p_free)(void * ptr))
{
        /* step 1: allocate new object */

        Alignment_ASW *al;
        if ((al = (Alignment_ASW*)p_malloc(sizeof(Alignment_ASW))) == NULL)
                goto cleanup;

        al->p_malloc = p_malloc;
        al->p_realloc = p_realloc;
        al->p_free = p_free;

        if ((al->match_penalty = (int*)p_malloc(sizeof(int)*PHRED_RANGE)) == NULL)
                goto cleanup;
        if ((al->mismatch_penalty = (int*)p_malloc(sizeof(int)*PHRED_RANGE)) == NULL)
                goto cleanup;
        if ((al->gopen_penalty = (int*)p_malloc(sizeof(int)*PHRED_RANGE)) == NULL)
                goto cleanup;
        if ((al->gext_penalty = (int*)p_malloc(sizeof(int)*PHRED_RANGE)) == NULL)
                goto cleanup;

        if ((al->matTra = (cigar_t**)p_malloc(sizeof(cigar_t*) * 1u)) == NULL)
                goto cleanup;

#ifdef DEBUG
        if ((al->matPen = (int**)p_malloc(sizeof(int*) * 1u)) == NULL)
                goto cleanup;
        if ((al->matIns = (int**)p_malloc(sizeof(int*) * 1u)) == NULL)
                goto cleanup;
        if ((al->matDel = (int**)p_malloc(sizeof(int*) * 1u)) == NULL)
                goto cleanup;
#endif

        /* step 2: set object members to zero */

        al->phred_offset = 0;

        al->db = NULL;
        al->subdb = NULL;
        al->query = NULL;
        al->subquery = NULL;
        al->qual = NULL;
        al->subqual = NULL;

        al->vecPen_m1_act = NULL;
        al->vecPen_m_act = NULL;
        al->vecIns_m1_act = NULL;
        al->vecIns_m_act = NULL;
        al->I_ext_m_act = NULL;
        al->I_ext_m1_act = NULL;

        al->db_len = 0u;
        al->subdb_len = 0u;
        al->query_len = 0u;
        al->subquery_len = 0u;
        al->offset = 0u;

        al->matTra[0] = NULL;
        al->rcigar = NULL;

#ifdef DEBUG
        al->matPen[0] = NULL;
        al->matIns[0] = NULL;
        al->matDel[0] = NULL;
#endif
        return al;
cleanup:
        asw_free(al);
        return NULL;
}


/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_set_phoffset
 *  Description:  Set  PHRED offset in the ASCII encoding
 * =====================================================================================
 */
void asw_set_phoffset(Alignment_ASW* al, int phred_offset)
{
        al->phred_offset = phred_offset;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_init
 *  Description:  Initialize an Alignment struct using provided penalty scores. The
 *                function is written in such a way as to allow multiple calls on a
 *                single object instance.
 * =====================================================================================
 */
Alignment_ASW* asw_init(Alignment_ASW* al, int MATCH_PEN, int MISMATCH_PEN, int GAP_OPEN_EXTEND, int GAP_EXTEND) {

        al->GAP_OPEN_EXTEND = GAP_OPEN_EXTEND;
        al->GAP_EXTEND = GAP_EXTEND;

        /* Initialize penalty look-up vectors */
        const double qN = -10.0 * log10(0.75); /* P(error | N) = 0.75 */

        unsigned int i;
        for (i = 0u; i < PHRED_RANGE; ++i) {
                double weight = 1.0 - pow(10.0, -((double)i + qN)/10.0);
                al->match_penalty[i] = 10 + round(weight * (double)MATCH_PEN);
                al->mismatch_penalty[i] = 10 + round(weight * (double)MISMATCH_PEN);
                al->gopen_penalty[i] = 10 + round(weight * (double)GAP_OPEN_EXTEND);
                al->gext_penalty[i] = 10 + round(weight * (double)GAP_EXTEND);
        }
        return al;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_free
 *  Description:  Free an Alignment struct
 * =====================================================================================
 */
void asw_free(Alignment_ASW *al) {
        if (al == NULL) return;

        if (al->match_penalty != NULL) al->p_free(al->match_penalty);
        if (al->mismatch_penalty != NULL) al->p_free(al->mismatch_penalty);
        if (al->gopen_penalty != NULL) al->p_free(al->gopen_penalty);
        if (al->gext_penalty != NULL) al->p_free(al->gext_penalty);

        if (al->vecPen_m_act != NULL) al->p_free(al->vecPen_m_act);
        if (al->vecPen_m1_act != NULL) al->p_free(al->vecPen_m1_act);
        if (al->vecIns_m_act != NULL) al->p_free(al->vecIns_m_act);
        if (al->vecIns_m1_act != NULL) al->p_free(al->vecIns_m1_act);

        if (al->I_ext_m_act != NULL) al->p_free(al->I_ext_m_act);
        if (al->I_ext_m1_act != NULL) al->p_free(al->I_ext_m1_act);

        if (al->rcigar != NULL) al->p_free(al->rcigar);

        if (al->matTra != NULL) {
                cigar_t ** matTra_p = al->matTra;
                cigar_t ** matTra_end = al->matTra + (al->subquery_len + 1u);
                for (; matTra_p < matTra_end; ++matTra_p) {
                        if (*matTra_p != NULL) al->p_free(*matTra_p);
                }
                al->p_free(al->matTra);
        }

#ifdef DEBUG
        if (al->matPen != NULL) {
                int ** matPen_p = al->matPen;
                int ** matPen_end = al->matPen + (al->subquery_len + 1u);
                for (; matPen_p < matPen_end; ++matPen_p) {
                        if (*matPen_p != NULL) al->p_free(*matPen_p);
                }
                al->p_free(al->matPen);
        }
        if (al->matIns != NULL) {
                int ** matIns_p = al->matIns;
                int ** matIns_end = al->matIns + (al->subquery_len + 1u);
                for (; matIns_p < matIns_end; ++matIns_p) {
                        if (*matIns_p != NULL) al->p_free(*matIns_p);
                }
                al->p_free(al->matIns);
        }
        if (al->matDel != NULL) {
                int ** matDel_p = al->matDel;
                int ** matDel_end = al->matDel + (al->subquery_len + 1u);
                for (; matDel_p < matDel_end; ++matDel_p) {
                        if (*matDel_p != NULL) al->p_free(*matDel_p);
                }
                al->p_free(al->matDel);
        }
#endif
        al->p_free(al);
}


/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_print_cigar
 *  Description:  Print CIGAR traceback to specified file or stream
 * =====================================================================================
 */
void asw_print_cigar(const Alignment_ASW* al, FILE *fp)
{
        const cigar_t * cigar_p = al->cigar_begin;
        const cigar_t * cigar_end = al->cigar_end;
        for (; cigar_p < cigar_end; ++cigar_p) {
                cigar_t cigar = *cigar_p;
                fprintf(fp, "%d%c ",
                        cigar >> BAM_CIGAR_SHIFT,
                        cigar_chars[cigar & BAM_CIGAR_MASK]);
        }
        fputc('\n', fp);
        fflush(fp);
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_show_cigar
 *  Description:  Print CIGAR traceback to a (string) buffer
 * =====================================================================================
 */

const char* asw_show_cigar(const Alignment_ASW* al)
{
        size_t len = 1;
        const cigar_t * cigar_p = al->cigar_begin;
        const cigar_t * cigar_end = al->cigar_end;
        for (; cigar_p < cigar_end; ++cigar_p) {
            cigar_t cigar = *cigar_p;
            len += ndigits(cigar >> BAM_CIGAR_SHIFT) + 2;
        }
        if (len > 1) { len -= 1; }  // subtract last separator
        char *buf = (char*)calloc(len, sizeof(char));
        cigar_p = al->cigar_begin;
        cigar_end = al->cigar_end;
        for (char *p = buf; cigar_p < cigar_end; ++cigar_p) {
            cigar_t cigar = *cigar_p;
            int num = cigar >> BAM_CIGAR_SHIFT;
            sprintf(p, "%d%c ", num, cigar_chars[cigar & BAM_CIGAR_MASK]);
            p += sizeof(char) * ndigits(num) + 2;
        }
        buf[len - 1] = '\0';
        return buf;
}

#ifdef DEBUG
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_print_matrix1
 *  Description:  Prints a 2D array together with query sequence, database sequence,
 *                and query-associated quality vector
 * =====================================================================================
 */
void asw_print_matrix1(const Alignment_ASW *al, FILE *fp)
{
        int ** mat = al->matPen;
        const char *db = al->subdb;
        const char *query = al->subquery;
        const uint8_t *qual = al->subqual;
        size_t x_len = al->subdb_len;
        size_t y_len = al->subquery_len;

        const char format_c[] = "%5c";
        const char format_d[] = "%5d";
        const char format_u[] = "%5u";
        const char format_ucc[] = "%5u%5c%5c";
        const char format_uuc[] = "%5u%5u%5c";

        unsigned int i;

        /* Row 1: horizontal index */
        fprintf(fp, format_c, ' ');
        fprintf(fp, format_c, ' ');
        fprintf(fp, format_c, ' ');
        for (i = 0u; i <= x_len; ++i) {
                fprintf(fp, format_u, i);
        }
        fputc('\n', fp); fflush(fp);

        /* Row 2: database sequence */
        fprintf(fp, format_c, ' ');
        fprintf(fp, format_c, ' ');
        fprintf(fp, format_c, ' ');
        fprintf(fp, format_c, '-');
        const char *db_end = db + x_len;
        for (; db < db_end; ++db) {
                fprintf(fp, format_c, *db);
        }
        fputc('\n', fp); fflush(fp);

        /* Row 3: first matrix row */
        fprintf(fp, format_ucc, 0u, '-', '-');
        const int *cell_p = *mat,
                *cell_end = *mat + x_len + 1u;
        for (; cell_p < cell_end; ++cell_p) {
                fprintf(fp, format_d, *cell_p);
        }
        fputc('\n', fp); fflush(fp);
        ++mat;

        /* The rest of the matrix (note: **mat has already
         * been incremented by one) */
        int **row_end = mat + y_len;
        for (i = 1u; mat < row_end; ++mat, ++query, ++qual, ++i) {
                fprintf(fp, format_uuc, i, (unsigned int)*qual, *query);
                const int *cell_p = *mat,
                        *cell_end = *mat + x_len + 1u;
                for (; cell_p < cell_end; ++cell_p) {
                        fprintf(fp, format_d, *cell_p);
                }
                fputc('\n', fp); fflush(fp);
        }
        fputc('\n', fp); fflush(fp);
}
#endif

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  resize_matrix1
 *  Description:  Reallocates a 2D array such that the amount of work is minimized if
 *                the new dimensions are similar to the old
 * =====================================================================================
 */
cigar_t ** resize_matrix1(Alignment_ASW* al,
                          cigar_t **matTra,
                          size_t new_x_len,
                          size_t new_y_len)
{
        size_t old_x_len = al->subdb_len,
               old_y_len = al->subquery_len;

        size_t cigar_hor = sizeof(cigar_t) * (new_x_len + 1u);
        if (old_y_len != new_y_len) {

                /* first free the bottom half */
                cigar_t ** matTra_p = matTra + (new_y_len + 1u);
                cigar_t ** matTra_end = matTra + (old_y_len + 1u);
                for (; matTra_p < matTra_end; ++matTra_p) {
                        al->p_free(*matTra_p);
                }
                cigar_t ** tmp_matTra
                      = (cigar_t**)al->p_realloc(matTra, sizeof(cigar_t*) * (new_y_len + 1u));
                if (tmp_matTra != NULL) matTra = tmp_matTra; else goto error;

                /* fill the bottom half */
                matTra_p = matTra + (old_y_len + 1u);
                matTra_end = matTra + (new_y_len + 1u);
                for (; matTra_p < matTra_end; ++matTra_p) {
                        if ((*matTra_p = (cigar_t*)al->p_malloc(cigar_hor)) == NULL)
                                goto error;
                }
        }
        if (old_x_len != new_x_len) {

                /* resize the top half */
                cigar_t ** matTra_p = matTra;
                cigar_t ** matTra_end = matTra + (min(old_y_len, new_y_len) + 1u);
                for (; matTra_p < matTra_end; ++matTra_p) {
                        cigar_t *tmp = (cigar_t*)al->p_realloc(*matTra_p, cigar_hor);
                        if (tmp != NULL) *matTra_p = tmp; else goto error;
                }
        }
        return matTra;
error:
        return NULL;
}

#ifdef DEBUG
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  resize_matrix2
 *  Description:  Reallocates a 2D array such that the amount of work is minimized if
 *                the new dimensions are similar to the old
 * =====================================================================================
 */
int ** resize_matrix2(Alignment_ASW* al,
                      int **matTra,
                      size_t new_x_len,
                      size_t new_y_len)
{
        size_t old_x_len = al->subdb_len,
               old_y_len = al->subquery_len;

        size_t cigar_hor = sizeof(int) * (new_x_len + 1u);
        if (old_y_len != new_y_len) {

                /* first free the bottom half */
                int ** matTra_p = matTra + (new_y_len + 1u);
                int ** matTra_end = matTra + (old_y_len + 1u);
                for (; matTra_p < matTra_end; ++matTra_p) {
                        al->p_free(*matTra_p);
                }
                int ** tmp_matTra
                        = (int**)al->p_realloc(matTra, sizeof(int*) * (new_y_len + 1u));
                if (tmp_matTra != NULL) matTra = tmp_matTra; else goto error;

                /* fill the bottom half */
                matTra_p = matTra + (old_y_len + 1u);
                matTra_end = matTra + (new_y_len + 1u);
                for (; matTra_p < matTra_end; ++matTra_p) {
                        if ((*matTra_p = (int*)al->p_malloc(cigar_hor)) == NULL)
                                goto error;
                }
        }
        if (old_x_len != new_x_len) {

                /* resize the top half */
                int ** matTra_p = matTra;
                int ** matTra_end = matTra + (min(old_y_len, new_y_len) + 1u);
                for (; matTra_p < matTra_end; ++matTra_p) {
                        int *tmp = (int*)al->p_realloc(*matTra_p, cigar_hor);
                        if (tmp != NULL) *matTra_p = tmp; else goto error;
                }
        }
        return matTra;
error:
        return NULL;
}
#endif

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_prepare_query
 *  Description:  Assign data fields to prepare for the alignment
 * =====================================================================================
 */
int asw_prepare_query(Alignment_ASW *al,
                 const char* m_query,
                 const uint8_t* m_qual,
                 size_t m_query_len,
                 uint32_t clip_head,
                 uint32_t clip_tail)
{
        size_t m_subquery_len = m_query_len - clip_head - clip_tail;

        al->query_len = m_query_len;
        al->query = m_query;
        al->subquery = m_query + clip_head;
        al->qual = m_qual;
        al->subqual = m_qual + clip_head;

        /* Resize matrices and vectors
         *
         * TODO: prevent buffers from shrinking as input sequence sizes
         * decrease in size to reduce the amount of memory allocations
         * performed on a stream of many short sequences
         */
        cigar_t ** tmp1;
        tmp1 = resize_matrix1(al, al->matTra,
                                    al->subdb_len,
                                    m_subquery_len);
        if (tmp1 != NULL) al->matTra = tmp1; else goto error;
#ifdef DEBUG
        int ** tmp2;
        tmp2 = al->matPen = resize_matrix2(al, al->matPen,
                                    al->subdb_len,
                                    m_subquery_len);
        if (tmp2 != NULL) al->matPen = tmp2; else goto error;
        tmp2 = resize_matrix2(al, al->matIns,
                                    al->subdb_len,
                                    m_subquery_len);
        if (tmp2 != NULL) al->matIns = tmp2; else goto error;
        tmp2 = resize_matrix2(al, al->matDel,
                                    al->subdb_len,
                                    m_subquery_len);
        if (tmp2 != NULL) al->matDel = tmp2; else goto error;
#endif
        al->subquery_len = m_subquery_len;

        return 0;
error:
        //asw_free(al);
        return -1;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_prepare_db
 *  Description:  Assign data fields to prepare for the alignment
 * =====================================================================================
 */
int asw_prepare_db(Alignment_ASW *al,
                 const char* m_db,
                 size_t m_db_len,
                 uint32_t clip_head,
                 uint32_t clip_tail)
{
        size_t m_subdb_len = m_db_len - clip_head - clip_tail;

        al->db = m_db;
        al->db_len = m_db_len;
        al->subdb = m_db + clip_head;

        /* Resize matrices and vectors
         *
         * TODO: prevent buffers from shrinking as input sequence sizes
         * decrease in size to reduce the amount of memory allocations
         * performed on a stream of many short sequences
         */
        cigar_t ** tmp1;
        tmp1 = resize_matrix1(al, al->matTra,
                                    m_subdb_len,
                                    al->subquery_len);
        if (tmp1 != NULL) al->matTra = tmp1; else goto error;
#ifdef DEBUG
        int ** tmp2;
        tmp2 = al->matPen = resize_matrix2(al, al->matPen,
                                    m_subdb_len,
                                    al->subquery_len);
        if (tmp2 != NULL) al->matPen = tmp2; else goto error;
        tmp2 = resize_matrix2(al, al->matIns,
                                    m_subdb_len,
                                    al->subquery_len);
        if (tmp2 != NULL) al->matIns = tmp2; else goto error;
        tmp2 = resize_matrix2(al, al->matDel,
                                    m_subdb_len,
                                    al->subquery_len);
        if (tmp2 != NULL) al->matDel = tmp2; else goto error;
#endif
        if (al->subdb_len != m_subdb_len) {
                int *tmp;

                tmp = (int*)al->p_realloc(al->vecPen_m1_act,
                                    sizeof(int) * (m_subdb_len + 1u));
                if (tmp != NULL) al->vecPen_m1_act = tmp; else goto error;

                tmp = (int*)al->p_realloc(al->vecPen_m_act,
                                    sizeof(int) * (m_subdb_len + 1u));
                if (tmp != NULL) al->vecPen_m_act = tmp; else goto error;

                tmp = (int*)al->p_realloc(al->vecIns_m1_act,
                                    sizeof(int) * (m_subdb_len + 1u));
                if (tmp != NULL) al->vecIns_m1_act = tmp; else goto error;

                tmp = (int*)al->p_realloc(al->vecIns_m_act,
                                    sizeof(int) * (m_subdb_len + 1u));
                if (tmp != NULL) al->vecIns_m_act = tmp; else goto error;

                uint32_t *utmp;

                utmp = (uint32_t*)al->p_realloc(al->I_ext_m_act,
                                          sizeof(uint32_t) * (m_subdb_len + 1u));
                if (utmp != NULL) al->I_ext_m_act = utmp; else goto error;

                utmp = (uint32_t*)al->p_realloc(al->I_ext_m1_act,
                                          sizeof(uint32_t) * (m_subdb_len + 1u));
                if (utmp != NULL) al->I_ext_m1_act = utmp; else goto error;

                al->subdb_len = m_subdb_len;
        }
        return 0;
error:
        //asw_free(al);
        return -1;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_prepare
 *  Description:  Assign data fields to prepare for the alignment
 * =====================================================================================
 */
int asw_prepare(Alignment_ASW *al,
                 const char* m_db,
                 size_t m_db_len,
                 const char* m_query,
                 const uint8_t* m_qual,
                 size_t m_query_len,
                 uint32_t clip_head,
                 uint32_t clip_tail)
{
        size_t m_subdb_len = m_db_len - clip_head - clip_tail,
               m_subquery_len = m_query_len - clip_head - clip_tail;

        al->db = m_db;
        al->db_len = m_db_len;
        al->query_len = m_query_len;
        al->subdb = m_db + clip_head;
        al->query = m_query;
        al->subquery = m_query + clip_head;
        al->qual = m_qual;
        al->subqual = m_qual + clip_head;

        /* Resize matrices and vectors
         *
         * TODO: prevent buffers from shrinking as input sequence sizes
         * decrease in size to reduce the amount of memory allocations
         * performed on a stream of many short sequences
         */
        cigar_t ** tmp1;
        tmp1 = resize_matrix1(al, al->matTra,
                                    m_subdb_len,
                                    m_subquery_len);
        if (tmp1 != NULL) al->matTra = tmp1; else goto error;
#ifdef DEBUG
        int ** tmp2;
        tmp2 = al->matPen = resize_matrix2(al, al->matPen,
                                    m_subdb_len,
                                    m_subquery_len);
        if (tmp2 != NULL) al->matPen = tmp2; else goto error;
        tmp2 = resize_matrix2(al, al->matIns,
                                    m_subdb_len,
                                    m_subquery_len);
        if (tmp2 != NULL) al->matIns = tmp2; else goto error;
        tmp2 = resize_matrix2(al, al->matDel,
                                    m_subdb_len,
                                    m_subquery_len);
        if (tmp2 != NULL) al->matDel = tmp2; else goto error;
#endif
        al->subquery_len = m_subquery_len;
        if (al->subdb_len != m_subdb_len) {
                int *tmp;

                tmp = (int*)al->p_realloc(al->vecPen_m1_act,
                                    sizeof(int) * (m_subdb_len + 1u));
                if (tmp != NULL) al->vecPen_m1_act = tmp; else goto error;

                tmp = (int*)al->p_realloc(al->vecPen_m_act,
                                    sizeof(int) * (m_subdb_len + 1u));
                if (tmp != NULL) al->vecPen_m_act = tmp; else goto error;

                tmp = (int*)al->p_realloc(al->vecIns_m1_act,
                                    sizeof(int) * (m_subdb_len + 1u));
                if (tmp != NULL) al->vecIns_m1_act = tmp; else goto error;

                tmp = (int*)al->p_realloc(al->vecIns_m_act,
                                    sizeof(int) * (m_subdb_len + 1u));
                if (tmp != NULL) al->vecIns_m_act = tmp; else goto error;

                uint32_t *utmp;

                utmp = (uint32_t*)al->p_realloc(al->I_ext_m_act,
                                          sizeof(uint32_t) * (m_subdb_len + 1u));
                if (utmp != NULL) al->I_ext_m_act = utmp; else goto error;

                utmp = (uint32_t*)al->p_realloc(al->I_ext_m1_act,
                                          sizeof(uint32_t) * (m_subdb_len + 1u));
                if (utmp != NULL) al->I_ext_m1_act = utmp; else goto error;

                al->subdb_len = m_subdb_len;
        }
        return 0;
error:
        //asw_free(al);
        return -1;
}


/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_align_init_semi
 *  Description:  Fill out top row in the alignment matrix (semiglobal)
 * =====================================================================================
 */
void asw_align_init_semi(Alignment_ASW *al)
{
        const uint8_t* m_subqual = al->subqual;

        size_t m_subdb_len = al->subdb_len;

        int *gopen_penalty = al->gopen_penalty - al->phred_offset,
            *gext_penalty = al->gext_penalty - al->phred_offset;

        int *vecPen_m = al->vecPen_m_act,
            *vecIns_m = al->vecIns_m_act;
        uint32_t *I_ext_m = al->I_ext_m_act;

        cigar_t ** matTra = al->matTra;
#ifdef DEBUG
        int ** matPen = al->matPen;
        int ** matIns = al->matIns;
        int ** matDel = al->matDel;

        int GAP_OPEN_EXTEND = al->GAP_OPEN_EXTEND,
            GAP_EXTEND = al->GAP_EXTEND;
#endif
        /* qq - quality in the query at position m */
        unsigned int qq = (unsigned int)m_subqual[0];
        int gopen_true_pen = gopen_penalty[qq] - gext_penalty[qq];

        /* Top-left cell */
        vecPen_m[0] = 0;
        vecIns_m[0] = vecPen_m[0] + gopen_true_pen;
        I_ext_m[0] = 0;

#ifdef DEBUG
        int storedDel_score = vecPen_m[0] + (GAP_OPEN_EXTEND - GAP_EXTEND);
        matPen[0][0] = vecPen_m[0];
        matDel[0][0] = storedDel_score;
        matIns[0][0] = vecIns_m[0];
#endif
        cigar_t *rowTra = matTra[0];
        rowTra[0] = (0 << BAM_CIGAR_SHIFT) | BAM_CSEQ_MATCH;

        size_t n, n1;
        for (n = 0u, n1 = 1u; n < m_subdb_len; ++n, ++n1) {
                vecPen_m[n1] = 0;
                vecIns_m[n1] = vecPen_m[n1] + gopen_true_pen;
                /* topmost row consists of only horizontal moves (deletions) */
                rowTra[n1] = (n1 << BAM_CIGAR_SHIFT) | BAM_CDEL;
                I_ext_m[n1] = 0;
#ifdef DEBUG
                matPen[0][n1] = vecPen_m[n1];
                matDel[0][n1] = storedDel_score;
                matIns[0][n1] = vecIns_m[n1];
#endif
        }
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_align_init
 *  Description:  Fill out top row in the alignment matrix (global)
 * =====================================================================================
 */
void asw_align_init(Alignment_ASW *al)
{
        const uint8_t* m_subqual = al->subqual;

        size_t m_subdb_len = al->subdb_len;

        int *gopen_penalty = al->gopen_penalty - al->phred_offset,
            *gext_penalty = al->gext_penalty - al->phred_offset;

        int GAP_OPEN_EXTEND = al->GAP_OPEN_EXTEND,
            GAP_EXTEND = al->GAP_EXTEND;

        int *vecPen_m = al->vecPen_m_act,
            *vecIns_m = al->vecIns_m_act;
        uint32_t *I_ext_m = al->I_ext_m_act;

        cigar_t ** matTra = al->matTra;
#ifdef DEBUG
        int ** matPen = al->matPen;
        int ** matIns = al->matIns;
        int ** matDel = al->matDel;
#endif
        /* qq - quality in the query at position m */
        unsigned int qq = (unsigned int)m_subqual[0];
        int gopen_true_pen = gopen_penalty[qq] - gext_penalty[qq];

        /* Top-left cell */
        vecPen_m[0] = 0;
        vecIns_m[0] = vecPen_m[0] + gopen_true_pen;
        I_ext_m[0] = 0;
        int storedDel_score = vecPen_m[0] + (GAP_OPEN_EXTEND - GAP_EXTEND);

#ifdef DEBUG
        matPen[0][0] = vecPen_m[0];
        matDel[0][0] = storedDel_score;
        matIns[0][0] = vecIns_m[0];
#endif
        cigar_t *rowTra = matTra[0];
        rowTra[0] = (0 << BAM_CIGAR_SHIFT) | BAM_CSEQ_MATCH;

        size_t n, n1;
        for (n = 0u, n1 = 1u; n < m_subdb_len; ++n, ++n1) {
                storedDel_score += GAP_EXTEND;
                vecPen_m[n1] = storedDel_score;
                vecIns_m[n1] = vecPen_m[n1] + gopen_true_pen;
                /* topmost row consists of only horizontal moves (deletions) */
                rowTra[n1] = (n1 << BAM_CIGAR_SHIFT) | BAM_CDEL;
                I_ext_m[n1] = 0;
#ifdef DEBUG
                matPen[0][n1] = vecPen_m[n1];
                matDel[0][n1] = storedDel_score;
                matIns[0][n1] = vecIns_m[n1];
#endif
        }
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_align
 *  Description:  Performs global affine-gap alignment according to Gotoh algorithm
 *                using asymmetric quality-weighted scoring
 * =====================================================================================
 */
void asw_align(Alignment_ASW *al)
{
        const char *m_subdb = al->subdb,
                   *m_subquery = al->subquery;
        const uint8_t* m_subqual = al->subqual;
        size_t m_subdb_len = al->subdb_len,
               m_subquery_len = al->subquery_len;

        int GAP_OPEN_EXTEND = al->GAP_OPEN_EXTEND,
            GAP_EXTEND = al->GAP_EXTEND;

        int *match_penalty = al->match_penalty - al->phred_offset,
            *mismatch_penalty = al->mismatch_penalty - al->phred_offset;

        int *gopen_penalty = al->gopen_penalty - al->phred_offset,
            *gext_penalty = al->gext_penalty - al->phred_offset;

        cigar_t ** matTra = al->matTra;
#ifdef DEBUG
        int ** matPen = al->matPen;
        int ** matIns = al->matIns;
        int ** matDel = al->matDel;
#endif

        int *vecPen_m = al->vecPen_m_act,
            *vecIns_m = al->vecIns_m_act;
        uint32_t *I_ext_m = al->I_ext_m_act;

        int *vecIns_m1 = al->vecIns_m1_act,
            *vecPen_m1 = al->vecPen_m1_act;
        uint32_t *I_ext_m1 = al->I_ext_m1_act;

        /* Initialize first row */

        size_t m, m1;

        /* Fill out the rest of the matrix */

        for (m = 0u, m1 = 1u; m < m_subquery_len; ++m, ++m1) {
                /* cq - character in the query at position m */
                char cq = m_subquery[m];
                /* qq - quality in the query at position m */
                unsigned int qq = (unsigned int)m_subqual[m];
                int match_pen = match_penalty[qq],
                    mismatch_pen = mismatch_penalty[qq],
                    gopen_pen = gopen_penalty[qq],
                    gext_pen = gext_penalty[qq];
                cigar_t *rowTra = matTra[m1];
                uint32_t cD = 0u;

                /* leftmost column consists of only vertical moves (insertions) */
                uint32_t cI;
                int wI_extend = vecIns_m[0] + gext_pen;
                vecIns_m1[0] = wI_extend;
                I_ext_m1[0] = cI = I_ext_m[0] + 1u;
                rowTra[0] = (cI << BAM_CIGAR_SHIFT) | BAM_CINS;

                vecPen_m1[0] = wI_extend;
                int storedDel_score = vecPen_m1[0] + (GAP_OPEN_EXTEND - GAP_EXTEND);

#ifdef DEBUG
                int * rowPen = matPen[m1];
                int * rowDel = matDel[m1];
                int * rowIns = matIns[m1];
                rowPen[0] = vecPen_m1[0];
                rowDel[0] = storedDel_score;
                rowIns[0] = wI_extend;
#endif
                size_t n, n1;

                for (n = 0u, n1 = 1u; n < m_subdb_len; ++n, ++n1) {
                        int wD, wI, wM;
                        uint32_t cI;
                        int is_seq_match = IS_MATCH(m_subdb[n], cq);

                        /* deletion: horizontal move */
                        int wD_open = vecPen_m1[n] + GAP_OPEN_EXTEND;
                        int wD_extend = storedDel_score + GAP_EXTEND;

                        /* insertion: vertical move */
                        int wI_open = vecPen_m[n1] + gopen_pen;
                        int wI_extend = vecIns_m[n1] + gext_pen;

                        /* given equal scores, prefer extending
                         * existing gaps to opening new ones */
                        if (wD_open < wD_extend) {
                                storedDel_score = wD = wD_open;
                                cD = 1u;
                        } else {
                                storedDel_score = wD = wD_extend;
                                ++cD;
                        }
                        if (wI_open < wI_extend) {
                                vecIns_m1[n1] = wI = wI_open;
                                I_ext_m1[n1] = cI = 1u;
                        } else {
                                vecIns_m1[n1] = wI = wI_extend;
                                I_ext_m1[n1] = cI = I_ext_m[n1] + 1u;
                        }

                        int mstate;
                        if (is_seq_match) {
                                wM = vecPen_m[n] + match_pen;
                                mstate = BAM_CSEQ_MATCH;
                        } else {
                                wM = vecPen_m[n] + mismatch_pen;
                                mstate = BAM_CSEQ_MISMATCH;
                        }

                        /* Order of preference: M, I, D */
                        if (wI < wM) {
                                /* either insertion or deletion */
                                if (wD < wI) {
                                        /* deletion */
                                        rowTra[n1] = (cD << BAM_CIGAR_SHIFT) | BAM_CDEL;
                                        vecPen_m1[n1] = wD;
                                } else {
                                        /* insertion */
                                        rowTra[n1] = (cI << BAM_CIGAR_SHIFT) | BAM_CINS;
                                        vecPen_m1[n1] = wI;
                                }
                        } else if (wD < wM) {
                                /* deletion */
                                rowTra[n1] = (cD << BAM_CIGAR_SHIFT) | BAM_CDEL;
                                vecPen_m1[n1] = wD;
                        } else {
                                /* either match or mismatch */
                                rowTra[n1] = (1u << BAM_CIGAR_SHIFT) | mstate;
                                vecPen_m1[n1] = wM;
                        }
#ifdef DEBUG
                        rowDel[n1] = wD;
                        rowIns[n1] = wI;
                        rowPen[n1] = vecPen_m1[n1];
#endif
                }
                int* tmp;
                /* Swap vecIns_m1 and vecIns_m */
                tmp = vecIns_m1, vecIns_m1 = vecIns_m, vecIns_m = tmp;
                /* Swap vecPen_m1 and vecPen_m */
                tmp = vecPen_m1, vecPen_m1 = vecPen_m, vecPen_m = tmp;

                unsigned int* utmp;
                /* Swap I_ext_m1 and I_ext_m */
                utmp = I_ext_m1, I_ext_m1 = I_ext_m, I_ext_m = utmp;
        }
        /* At this point, vecIns_m and vecPen_m point to their corresponding
         * last rows, and vecIns_m1 and vecPen_m1 point to penultimate rows */

        al->vecPen_lastRow = vecPen_m;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_locate_minscore
 *  Description:  Find minimum alignment score (in last row)
 * =====================================================================================
 */

int asw_locate_minscore(Alignment_ASW* al)
{
        int * vecPen_m = al->vecPen_lastRow;
        int opt_score = vecPen_m[0];
        size_t opt_score_col = 0u;

        size_t n1;
        size_t m_subdb_len = al->subdb_len;
        for (n1 = 1u; n1 <= m_subdb_len; ++n1) {
                if (vecPen_m[n1] < opt_score) {
                        opt_score = vecPen_m[n1];
                        opt_score_col = n1;
                }
        }
        al->opt_score = opt_score;
        al->opt_score_col = opt_score_col;
        return opt_score;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_trace
 *  Description:  Produce a traceback given an Alignment_ASW struct
 *
 *     Modifies:  al->offset
 *                al->rcigar
 *                al->cigar_begin
 *                al->cigar_end
 * =====================================================================================
 */
int asw_trace(Alignment_ASW* al)
{
        assert(al->query_len >= al->subquery_len);

        /* resize cigar string to query length */
        cigar_t *rcigar = (cigar_t*)al->p_realloc(al->rcigar, sizeof(cigar_t) * (al->subquery_len + 4u));
        if (rcigar != NULL) al->rcigar = rcigar; else goto error;

        /* fill out cigar string */
        int m1 = (int)al->subquery_len,
            n1 = (int)al->opt_score_col;

        cigar_t ** matTra = al->matTra;

        cigar_t cigar = matTra[m1][n1];
        /* z - length of CIGAR operation */
        uint32_t z = cigar >> BAM_CIGAR_SHIFT;
        /* state - CIGAR operation */
        uint32_t state = cigar & BAM_CIGAR_MASK;

        /* rc - emulates a reverse iterator (except two elements at the end are omitted for padding) */
        cigar_t * rc = rcigar + al->subquery_len + 1u;
        unsigned int num_matches = 0u;

        while (m1 > 0) {
        switch (state) {
        case BAM_CSEQ_MATCH:
                num_matches = 0u;
                do {
                        /* simply accumulate num_matches */
                        num_matches += z;
                        m1 -= z, n1 -= z;
                        cigar = matTra[m1][n1];
                        z = cigar >> BAM_CIGAR_SHIFT;
                        state = cigar & BAM_CIGAR_MASK;
                } while (state == BAM_CSEQ_MATCH && m1 > 0);
                *rc = (num_matches << BAM_CIGAR_SHIFT) | BAM_CSEQ_MATCH;
                --rc;
                break;
        case BAM_CSEQ_MISMATCH:
                num_matches = 0;
                do {
                        /* simply accumulate num_matches */
                        num_matches += z;
                        m1 -= z, n1 -= z;
                        cigar = matTra[m1][n1];
                        z = cigar >> BAM_CIGAR_SHIFT;
                        state = cigar & BAM_CIGAR_MASK;
                } while (state == BAM_CSEQ_MISMATCH && m1 > 0);
                *rc = (num_matches << BAM_CIGAR_SHIFT) | BAM_CSEQ_MISMATCH;
                --rc;
                break;
        case BAM_CDEL:
                *rc = cigar;
                --rc;
                n1 -= z;
                cigar = matTra[m1][n1];
                z = cigar >> BAM_CIGAR_SHIFT;
                state = cigar & BAM_CIGAR_MASK;
                break;
        case BAM_CINS:
                *rc = cigar;
                --rc;
                m1 -= z;
                cigar = matTra[m1][n1];
                z = cigar >> BAM_CIGAR_SHIFT;
                state = cigar & BAM_CIGAR_MASK;
                break;
        default:
                fprintf(stderr, "ERROR: unknown CIGAR operation %u\n", state);
                goto error;
        }}

        cigar_t * fc5p = rc + 1u,
                * fc3p = rcigar + al->subquery_len + 2u;

        al->offset = n1;
        al->cigar_begin = fc5p;
        al->cigar_end = fc3p;
        return 0;
error:
        return -1;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_append_softclip
 *  Description:  extend CIGAR trace outer boundaries to previous clipping
 *
 *     Modifies:  al->rcigar
 *                al->cigar_begin
 *                al->cigar_end
 *                al->offset
 * =====================================================================================
 */
void asw_append_softclip(Alignment_ASW* al)
{
        assert(al->subquery >= al->query);
        assert(al->subdb >= al->db);
        uint32_t clip_head = al->subquery - al->query;
        if (clip_head > 0u) {
                /* clipped beginnning */
                cigar_t cigar = *al->cigar_begin;
                uint32_t state = cigar & BAM_CIGAR_MASK,
                         z = cigar >> BAM_CIGAR_SHIFT;
                if (state == BAM_CSOFT_CLIP) {
                        /* extend clipping */
                        *al->cigar_begin = ((clip_head + z) << BAM_CIGAR_SHIFT)
                                       | BAM_CSOFT_CLIP;
                } else if (state == BAM_CSEQ_MATCH || state == BAM_CMATCH) {
                        /* try to contract clipping */
                        uint32_t match_add = 0u;
                        const char *subquery = al->subquery,
                                   *subdb = al->subdb + al->offset;
                        while (clip_head > 0u && *--subquery == *--subdb) {
                                ++match_add;
                                --clip_head;
                        }
                        if (match_add > 0u) {
                                *al->cigar_begin = ((z + match_add) << BAM_CIGAR_SHIFT)
                                               | state;
                                al->offset -= match_add;
                        }
                        if (clip_head > 0u) {
                                *(--al->cigar_begin) = (clip_head << BAM_CIGAR_SHIFT)
                                                   | BAM_CSOFT_CLIP;
                        }
                } else {
                        /* add clipping */
                        *(--al->cigar_begin) = (clip_head << BAM_CIGAR_SHIFT)
                                               | BAM_CSOFT_CLIP;
                }
        }
        uint32_t clip_tail = (al->query + al->query_len) - (al->subquery + al->subquery_len);

        if (clip_tail > 0u) {
                /* clipped end */
                cigar_t cigar = *(al->cigar_end - 1u);
                uint32_t state = cigar & BAM_CIGAR_MASK,
                         z = cigar >> BAM_CIGAR_SHIFT;
                if (state == BAM_CSOFT_CLIP) {
                        /* extend clipping */
                        *(al->cigar_end - 1u) = ((clip_tail + z) << BAM_CIGAR_SHIFT)
                                    | BAM_CSOFT_CLIP;
                } else if (state == BAM_CSEQ_MATCH || state == BAM_CMATCH) {
                        /* try to contract clipping */

                        uint32_t match_add = 0u;
                        const char *subquery = al->subquery + al->subquery_len - 1u,
                                   *subdb = al->subdb + al->offset + al->subdb_len - 1u;

                        while (clip_tail > 0u && *++subquery == *++subdb) {
                                ++match_add;
                                --clip_tail;
                        }
                        if (match_add > 0u) {
                                *(al->cigar_end - 1u)
                                        = ((z + match_add) << BAM_CIGAR_SHIFT) | state;
                        }
                        if (clip_tail > 0u) {
                                *al->cigar_end = (clip_tail << BAM_CIGAR_SHIFT)
                                                   | BAM_CSOFT_CLIP;
                                ++al->cigar_end;
                        }
                } else {
                        /* add clipping */
                        *al->cigar_end = (clip_tail << BAM_CIGAR_SHIFT)
                                               | BAM_CSOFT_CLIP;
                        ++al->cigar_end;
                }
        }
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_append_hardclip
 *  Description:  extend CIGAR trace outer boundaries to previous clipping
 *
 *     Modifies:  al->rcigar
 *                al->cigar_begin
 *                al->cigar_end
 * =====================================================================================
 */
void asw_append_hardclip(Alignment_ASW* al, uint32_t clip_head, uint32_t clip_tail)
{
        if (clip_head > 0u) {
                /* clipped beginnning */
                cigar_t cigar = *al->cigar_begin;
                uint32_t state = cigar & BAM_CIGAR_MASK,
                         z = cigar >> BAM_CIGAR_SHIFT;
                if (state == BAM_CHARD_CLIP) {
                        /* extend clipping */
                        *al->cigar_begin = ((clip_head + z) << BAM_CIGAR_SHIFT)
                                       | BAM_CHARD_CLIP;
                } else {
                        /* add clipping */
                        *(--al->cigar_begin) = (clip_head << BAM_CIGAR_SHIFT)
                                               | BAM_CHARD_CLIP;
                }
        }
        if (clip_tail > 0u) {
                /* clipped end */
                cigar_t cigar = *(al->cigar_end - 1u);
                uint32_t state = cigar & BAM_CIGAR_MASK,
                         z = cigar >> BAM_CIGAR_SHIFT;
                if (state == BAM_CHARD_CLIP) {
                        /* extend clipping */
                        *(al->cigar_end - 1u) = ((clip_tail + z) << BAM_CIGAR_SHIFT)
                                    | BAM_CHARD_CLIP;
                } else {
                        /* add clipping */
                        *al->cigar_end = (clip_tail << BAM_CIGAR_SHIFT)
                                               | BAM_CHARD_CLIP;
                        ++al->cigar_end;
                }
        }
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_softclip_trace
 *  Description:  replace edits at the end of the alignment that are not exact matches
 *                with soft clipping
 *
 *     Modifies:  al->rcigar
 *                al->cigar_begin
 *                al->cigar_end
 *                al->offset
 * =====================================================================================
 */
void asw_softclip_trace(Alignment_ASW* al)
{
        /* scan CIGAR from the tail backwards until the last match:
         *                   |<-----
         * 5= 1X 2D 20= 1I 30= 3I 1X
         */
        unsigned int soft_clip_3p = 0;
        cigar_t * rc3p = al->cigar_end - 1u;
        cigar_t * rc = al->cigar_begin - 1u;
        for (; rc3p != rc; --rc3p) {
                cigar_t cigar = *rc3p;
                uint32_t state = cigar & BAM_CIGAR_MASK;
                if (state == BAM_CSEQ_MATCH) {
                        break;
                } else if (state != BAM_CDEL && state != BAM_CHARD_CLIP) {
                        /* BAM_CSEQ_MISMATCH or BAM_CINS */
                        soft_clip_3p += (cigar >> BAM_CIGAR_SHIFT);
                }
        }
        if (soft_clip_3p > 0u) {
                ++rc3p;
                *rc3p = (soft_clip_3p << BAM_CIGAR_SHIFT) | BAM_CSOFT_CLIP;
        }
        cigar_t * fc3p = rc3p + 1u;

        /* scan CIGAR from forward until the first match:
         * --->|
         * 1X 5= 2D 20= 1I 30= 3I 1X
         */
        size_t offset = al->offset;
        unsigned int soft_clip_5p = 0u;
        cigar_t * fc5p = al->cigar_begin;

        for (; fc5p != fc3p; ++fc5p) {
                cigar_t cigar = *fc5p;
                uint32_t state = cigar & BAM_CIGAR_MASK;
                if (state == BAM_CSEQ_MATCH) {
                        break;
                } else if (state != BAM_CHARD_CLIP) {
                        int op_len = cigar >> BAM_CIGAR_SHIFT;
                        if (state != BAM_CDEL) {
                                /* BAM_CSEQ_MISMATCH or BAM_CINS */
                                soft_clip_5p += op_len;
                        }
                        if (state != BAM_CINS) {
                                /* BAM_CDEL, BAM_CSEQ_MISMATCH */
                                offset += op_len;
                        }
                }
        }
        if (soft_clip_5p > 0) {
                --fc5p;
                *fc5p = (soft_clip_5p << BAM_CIGAR_SHIFT) | BAM_CSOFT_CLIP;
        }
        rc = fc5p - 1u;

        assert(fc3p >= fc5p);
        al->offset = offset;
        al->cigar_begin = fc5p;
        al->cigar_end = fc3p;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_compact_trace
 *  Description:  collapse the CIGAR string by treating matches and mismatches as the
 *                same state
 *
 *     Modifies:  al->rcigar
 *                al->cigar_begin
 * =====================================================================================
 */
void asw_compact_trace(Alignment_ASW* al)
{
        cigar_t *rbucket, *start_riter, *rc;
        rbucket = start_riter
                = al->cigar_end - 1u;
        rc = al->cigar_begin - 1u;

        while (rbucket != rc) {
                cigar_t cigar;
                uint32_t num_matches = 0u;
                while (1) {
                        cigar = *rbucket;
                        uint32_t state = cigar & BAM_CIGAR_MASK;
                        --rbucket;
                        if (state == BAM_CSEQ_MATCH || state == BAM_CSEQ_MISMATCH) {
                                num_matches += (cigar >> BAM_CIGAR_SHIFT);
                                if (rbucket == rc) {
                                        cigar = (num_matches << BAM_CIGAR_SHIFT) | \
                                                BAM_CMATCH;
                                        break;
                                }
                        } else if (num_matches > 0u) {
                                *start_riter = (num_matches << BAM_CIGAR_SHIFT) | \
                                               BAM_CMATCH;
                                --start_riter;
                                break;
                        } else {
                                break;
                        }
                }
                *start_riter = cigar;
                --start_riter;
        }
        cigar_t * fc5p = start_riter + 1u;

        assert(al->cigar_end >= fc5p);
        al->cigar_begin = fc5p;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  getAlignmentStart
 *  Description:  Get starting position of alignment in the genome given offset
 *                    obtained from realignment.
 * =====================================================================================
 */
int32_t asw_getAlignmentStart(const Alignment_ASW* al, int alstart) {
        assert(al->subdb >= al->db);
        return max(0, alstart) + al->offset + (al->subdb - al->db);
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_getBasicAlignPair
 *  Description:  Allocates and initializes a BASICALIGNPAIR struct for use in
 *                composition alignment from Alignment_ASW struct
 * =====================================================================================
 */
BASICALIGNPAIR* asw_getBasicAlignPair(const Alignment_ASW* al) {

        /* In BASICALIGNPAIR, the length of the alignment is the total number of
         * traversal instructions in traceback */
        size_t len;
        cigar_t *cigar_p = al->cigar_begin;
        cigar_t *cigar_end = al->cigar_end;
        for (len = 0u; cigar_p < cigar_end; ++cigar_p) {
                len += (*cigar_p >> BAM_CIGAR_SHIFT);
        }
        char * seq1 = (char*)al->p_malloc(len + 1u);
        char * seq2 = (char*)al->p_malloc(len + 1u);
        char * seq1_p = seq1;
        char * seq2_p = seq2;
        const char * seq1source_p = al->subdb + al->offset;
        const char * seq2source_p = al->subquery;
        for (cigar_p = al->cigar_begin; cigar_p < cigar_end; ++cigar_p) {
                cigar_t cigar = *cigar_p;
                uint32_t i = 0u;
                uint32_t op = cigar & BAM_CIGAR_MASK;
                uint32_t op_len = cigar >> BAM_CIGAR_SHIFT;
                switch (op) {
                case BAM_CHARD_CLIP:
                        /* do nothing */
                        break;
                case BAM_CSOFT_CLIP:
                        /* diagonal move: either match or mismatch */
                        for (; i < op_len; ++i) {
                                ++seq1_p, ++seq1source_p;
                                ++seq2_p, ++seq2source_p;
                        }
                        break;
                case BAM_CMATCH:
                case BAM_CSEQ_MATCH:
                case BAM_CSEQ_MISMATCH:
                        /* diagonal move: either match or mismatch */
                        for (; i < op_len; ++i) {
                                *seq1_p = *seq1source_p;
                                ++seq1_p, ++seq1source_p;
                                *seq2_p = *seq2source_p;
                                ++seq2_p, ++seq2source_p;
                        }
                        break;
                case BAM_CINS:
                        /* vertical move: letters in query but not in the reference */
                        for (; i < op_len; ++i) {
                                *seq1_p = '-';
                                ++seq1_p;
                                *seq2_p = *seq2source_p;
                                ++seq2_p, ++seq2source_p;
                        }
                        break;
                case BAM_CDEL:
                        /* horizontal move: letters in reference but not in the query */
                        for (; i < op_len; ++i) {
                                *seq1_p = *seq1source_p;
                                ++seq1_p, ++seq1source_p;
                                *seq2_p = '-';
                                ++seq2_p;
                        }
                        break;
                default:
                        fprintf(stderr, "ERROR (getBasicAlignPair): unknown CIGAR operation %u\n", op);
                        return NULL;
                }
        }
        assert(seq1_p - seq1 == seq2_p - seq2);
        *seq1_p = '\0';
        *seq2_p = '\0';

        /* fprintf(stdout, "1: %s\n", seq1);
         * fprintf(stdout, "2: %s\n", seq2);
         */

        BASICALIGNPAIR* bap = (BASICALIGNPAIR*)malloc(sizeof(BASICALIGNPAIR));
        bap->sequence1side = seq1;
        bap->sequence2side = seq2;
        bap->sequence1start = al->offset;
        bap->sequence1end = al->opt_score_col - 1;
        bap->sequence2start = 0;
        bap->sequence2end = al->subquery_len - 1;
        bap->score = al->opt_score;
        bap->length = (int)len;
        return bap;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  freeBasicAlignPair
 *  Description:  Frees a BASICALIGNPAIR struct
 * =====================================================================================
 */
void freeBasicAlignPair(BASICALIGNPAIR* ap)
{
        free(ap->sequence1side);
        free(ap->sequence2side);
        free(ap);
        return;
}
