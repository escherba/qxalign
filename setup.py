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
        long_description="qxalign is..."
    )
