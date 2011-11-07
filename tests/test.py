#!/usr/bin/env python3

from qxalign import Qxalign

q = Qxalign()

q.prepare("AAAACGT","TGCA",b"!!!!")
print(q.align())
q.trace()
q.print_trace()

q.prepare_query(query_seq="CAAC")
print(q.align(semi=True))
q.trace()
q.print_trace()

q.prepare("", "", "")
try:
    print(q.align())
except IndexError:
    print("caught qxalign exception")
    raise
