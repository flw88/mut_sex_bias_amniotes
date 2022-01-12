#!/usr/bin/env python

import sys
import numpy as np

reptimings = []
lengths = []

for line in sys.stdin:
    line = line.strip()
    fields = line.split()
    l = int(fields[2])-int(fields[1])
    lengths.append(l)
    reptimings.append(float(fields[-1]))

if len(lengths)>0:
    sys.stdout.write("{:.5f}\n".format(np.average(reptimings,weights=lengths)))
else:
    sys.stdout.write("{}\n".format(np.nan))
    
