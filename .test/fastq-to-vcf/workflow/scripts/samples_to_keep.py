#!/usr/bin/env python3
import sys

infile = sys.argv[1]
DEPTH_CUTOFF = 2
samps = []
with open(infile, "r") as f:
    next(f)
    for line in f:
        line = line.strip().split()
        if float(line[2]) >= DEPTH_CUTOFF:
            samps.append(line[0])

for s in samps:
    print(s)
        