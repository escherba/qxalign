/*
 * =====================================================================================
 *
 *       Filename:  align454.h
 *
 *    Description:  Collection of routines for quality-aware alignment of Roche/454
 *                  reads. Functions with asw_* prefix implement asymmetric Smith-
 *                  Waterman-like algorithm with inverse scores (URL:
 *                  http://dx.doi.org/10.1101/gr.6468307)
 *
 *        Version:  1.0
 *        Created:  04/05/2011 10:05:57
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Eugene Scherba (es), escherba@bu.edu
 *        Company:  Laboratory for Biocomputing and Informatics, Boston University
 *
 * =====================================================================================
 */

#ifndef NDEBUG
#define DEBUG
#endif

// Sanger PHRED scores range from 0 to 93
#define PHRED_RANGE 94

#ifdef __cplusplus
extern "C" {
#endif

typedef uint32_t cigar_t;

struct Alignment_ASW {

        /* PHRED offset in the ASCII encoding: 33 for Sanger format */
        int phred_offset;

        /* look-up tables for quality-based scoring */
        int *match_penalty,
            *mismatch_penalty,
            *gopen_penalty,
            *gext_penalty;

        int GAP_OPEN_EXTEND,
            GAP_EXTEND;

        /* pointers to external arrays */
        const char *db,
                   *subdb,
                   *query,
                   *subquery;
        const uint8_t *qual,
                      *subqual;

        /* The following must always hold:
         *
         * subdb >= db &&
         * subquery >= query &&
         * subqual >= qual &&
         * subdb_len <= db_len &&
         * subquery_len <= query_len &&
         * subquery - query == subqual - qual
         */

        /* lengths of these external arrays */
        size_t db_len,
               subdb_len,
               query_len,
               subquery_len;

        int *vecPen_m_act,      /* "previous" row in matPen */
            *vecPen_m1_act,     /* "current" row in matPen */
            *vecIns_m_act,      /* "previous" row in matIns */
            *vecIns_m1_act;     /* "current" row in matIns */

        uint32_t *I_ext_m_act,  /* "previous" vector with insertion lengths */
                 *I_ext_m1_act; /* "current" vector with insertion lengths */

        int *vecPen_lastRow;    /* vector corresponding to last row in matPen, always
                                 * equal to either vecPen_m_act or vecPen_m1_act */

        cigar_t **matTra;       /* trace matrix */

#ifdef DEBUG
        int **matPen;           /* scoring matrix */
        int **matIns;           /* scoring matrix */
        int **matDel;           /* scoring matrix */
#endif

        int opt_score;          /* minimum score in the last row */
        size_t opt_score_col;     /* column containing cell with minimum score */
        /* size_t opt_score_row; */  /* row containing cell with minimum score */

        cigar_t *rcigar,        /* buffer containing the CIGAR segment */
                *cigar_begin,   /* pointer to the start of the CIGAR segment */
                *cigar_end;     /* pointer to the end of the CIGAR segment */

        size_t offset;          /* position in reference where to start the alignment */

        /* "virtual table" */

        void *(*p_malloc)(size_t size);
        void *(*p_realloc)(void * ptr, size_t size);
        void (*p_free)(void * ptr);

};                             /* ----------  end of struct Alignment_ASW  ---------- */

typedef struct Alignment_ASW Alignment_ASW;

typedef struct
{
        char*   sequence1side;
        char*   sequence2side;
        int     sequence1start;
        int     sequence1end;
        int     sequence2start;
        int     sequence2end;
        int     score;
        int     length;
}   BASICALIGNPAIR;


/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_new
 *  Description:  Allocate and initialize an Alignment_ASW struct
 * =====================================================================================
 */
Alignment_ASW* asw_new(int MATCH_PEN,
                       int MISMATCH_PEN,
                       int GAP_OPEN_EXTEND,
                       int GAP_EXTEND);
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_alloc
 *  Description:  Allocate an Alignment_ASW struct and set its members to zero
 * =====================================================================================
 */
Alignment_ASW* asw_alloc(
        void *(*p_malloc)(size_t size),
        void *(*p_realloc)(void * ptr, size_t size),
        void (p_free)(void * ptr));

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_init
 *  Description:  Initialize an Alignment_ASW struct using provided penalty scores
 * =====================================================================================
 */
Alignment_ASW* asw_init(Alignment_ASW *al,
                        int MATCH_PEN,
                        int MISMATCH_PEN,
                        int GAP_OPEN_EXTEND,
                        int GAP_EXTEND);
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_free
 *  Description:  Free an Alignment_ASW struct
 * =====================================================================================
 */
void asw_free(Alignment_ASW *al);

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_set_phoffset
 *  Description:  Set  PHRED offset in the ASCII encoding
 * =====================================================================================
 */
void asw_set_phoffset(Alignment_ASW* al, int phred_offset);

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
                 uint32_t clip_tail);
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
                 uint32_t clip_tail);

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_prepare
 *  Description:  Assign data fields and prepare for alignment
 * =====================================================================================
 */
int asw_prepare(Alignment_ASW *al,
                 const char* m_db,
                 size_t ref_string_len,
                 const char* m_query,
                 const uint8_t* m_qual,
                 size_t query_string_len,
                 uint32_t clip_head,
                 uint32_t clip_tail);
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_align_init
 *  Description:  Fill out top row in the alignment matrix (global)
 * =====================================================================================
 */
void asw_align_init(Alignment_ASW *al);

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_align_init_semi
 *  Description:  Fill out top row in the alignment matrix (semiglobal)
 * =====================================================================================
 */
void asw_align_init_semi(Alignment_ASW *al);


/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_align
 *  Description:  Performs an global affine-gap alignment according to Gotoh algorithm
 *                using asymmetric quality-weighted scoring
 * =====================================================================================
 */
void asw_align(Alignment_ASW *al);

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_locate_minscore
 *  Description:  Find minimum alignment score (in last row)
 * =====================================================================================
 */
int asw_locate_minscore(Alignment_ASW* al);

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_trace
 *  Description:  Produce a traceback given an Alignment_ASW struct
 * =====================================================================================
 */
int asw_trace(Alignment_ASW* al);

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
void asw_append_softclip(Alignment_ASW* al);

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
void asw_append_hardclip(Alignment_ASW* al, uint32_t clip_head, uint32_t clip_tail);

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_softclip_trace
 *  Description:  Replace non-matching edits at the ends with soft clipping
 * =====================================================================================
 */
void asw_softclip_trace(Alignment_ASW* al);

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_compact_trace
 *  Description:  Replace BAM_CSEQ_MATCH and BAM_CSEQ_MISMATCH with BAM_CMATCH
 * =====================================================================================
 */
void asw_compact_trace(Alignment_ASW* al);

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_getAlignmentStart
 *  Description:  Get starting position of alignment in the genome given offset
 *                obtained from realignment.
 * =====================================================================================
 */
int32_t asw_getAlignmentStart(const Alignment_ASW* al, int alstart);

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_getBasicAlignPair
 *  Description:  Allocates and initializes a BASICALIGNMENTPAIR struct for use in
 *                composition alignment from Alignment_ASW struct
 * =====================================================================================
 */
BASICALIGNPAIR* asw_getBasicAlignPair(const Alignment_ASW* al);

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  freeBasicAlignPair
 *  Description:  Frees a BASICALIGNPAIR struct
 * =====================================================================================
 */
void freeBasicAlignPair(BASICALIGNPAIR* ap);

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_print_cigar
 *  Description:  Print CIGAR traceback to specified file or stream
 * =====================================================================================
 */
void asw_print_cigar(const Alignment_ASW* al, FILE *fp);

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_show_cigar
 *  Description:  Print CIGAR traceback to specified file or stream
 * =====================================================================================
 */
const char* asw_show_cigar(const Alignment_ASW* al);

#ifdef DEBUG
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  asw_print_matrix1
 *  Description:  Prints a 2D array together with query sequence, database sequence,
 *                and query-associated quality vector
 * =====================================================================================
 */
void asw_print_matrix1(const Alignment_ASW *al, FILE *fp);
#endif

#ifdef __cplusplus
}
#endif
