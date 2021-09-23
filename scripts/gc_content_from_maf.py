#!/usr/bin/env python

import sys
from optparse import OptionParser
import numpy as np

def safediv(n,d):
    return np.nan if d==0 else n/float(d)

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-s", "--start", action="store", type="str", dest="species_list",help="species list to compute gc% on")
(options, args) = parser.parse_args()

# Useful variables
gc = {}
gc_letters = ["G","g","C","c"]

# Read species
with open(options.species_list, "r") as spfile:
    species = [line.strip().split(",") for line in spfile][0]

# Read from stdin
for line in sys.stdin:

    line = line.strip()
    fields = line.split()

    # If line with sequences, process
    if line.startswith("s"):

        sp = fields[1].split(".")[0]  # Assumes there are not "." in species names
        
        # Skip species that are not in the list
        if sp not in species:
            continue

        # Initiate key in dictionary if not present yet
        if sp not in gc:
            gc[sp] = {"GC":0, "total":0}

        # Count G, g, C and c
        seq = fields[-1]
        gc_count = sum([1 for nuc in seq if nuc in gc_letters])
        gc[sp]["GC"] += gc_count
        gc[sp]["total"] += len(seq)
        
# Print result
# If no data
if len(gc)==0:
    print("\t".join(["nan" for sp in species]), end="")
# If data
else:
    for i,sp in enumerate(species):
        gc_frac = safediv(gc[sp]["GC"], float(gc[sp]["total"])) if sp in gc else np.nan
        if len(species)-1==i:
            print("{:.3f}".format(gc_frac), end="")
        else:
            print("{:.3f}\t".format(gc_frac), end="")
print("")
