#!/usr/bin/env python
import re
from numpy import asarray, float32
from numpy.linalg import norm
import sys

fname = sys.argv[1]

a = [re.findall(r'^.*MC.*$', line) for line in open(fname)]
c = [re.findall("[-+]?\d+[\.]?\d*", s[0]) for s in a if s != []]
cc = asarray(list(map(float32, c)))
adj_msc = cc[:, 3]

print(norm(adj_msc))
