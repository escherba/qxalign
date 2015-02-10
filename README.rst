Qxalign
=======

This is a collection of routines for quality-aware alignment of Roche/454 reads,
for use either from C or from Python 3 (including a Python 3 C extension).
Functions with asw_* prefix implement asymmetric Smith- Waterman-like algorithm
with inverse scores (URL: http://dx.doi.org/10.1101/gr.6468307)

.. code-block:: python

    >>> from qxalign import Qxalign
    ...
    >>> q = Qxalign()
    >>> q.prepare("AAAACGT", "TGCA", b"!!!!")
    >>> q.align()
    60
    >>> q.trace()
    >>> q.show_trace()
    '3I 1='
