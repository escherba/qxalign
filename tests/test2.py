#!/usr/bin/env python3

import qxalign

q = qxalign.Qxalign()

q.prepare("AAAACGT","TGCA","!!!!!!!!!!!")
print(q.align())
