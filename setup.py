#!/usr/bin/env python3

from distutils.core import setup, Extension

setup(
        ext_modules=[Extension("qxalign",
                sources=["qxalign.c", "align454.c"])
            ],
        name="qxalign",
        author="Eugene Scherba, Boston University",
        author_email="escherba@bu.edu",
        url="http://tandem.bu.edu/",
        license="BSD",
        version="0.0.1",
        description="Module for quality-aware realignment",
        platforms = ["Linux", "Windows", "Mac OS X"],
        long_description="""\
This is a collection of routines for quality-aware alignment of Roche/454 reads,
for use either from C or from Python 3 (including a Python 3 C extension).
Functions with asw_* prefix implement asymmetric Smith- Waterman-like algorithm
with inverse scores (URL: http://dx.doi.org/10.1101/gr.6468307)
 
"""
    )
