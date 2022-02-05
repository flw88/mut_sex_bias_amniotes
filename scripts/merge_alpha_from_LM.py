#!/usr/bin/env python

import sys
from optparse import OptionParser
import numpy as np

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-o", "--outgroups", action="store", type="str", dest="outgroup_file",help="Tab separated file with name of group and species that is an outgroup")
parser.add_option("-g", "--group", action="store", type="str", dest="group_name",help="Name of the group")
parser.add_option("-s", "--suffix", action="store", type="str", dest="suffix",help="Suffix (e.g., LM.tsv)")
(options, args) = parser.parse_args()

# Read outgroups
with open(options.outgroup_file, "r") as out_f:
    group2out = {line.split()[0]:line.split()[1] for line in out_f}

# Read alpha from LM
group = options.group_name
suffix = options.suffix
directory="alphas/"
with open("{}/{}.{}".format(directory, group, suffix),"r") as alpha_results:
    for i,line in enumerate(alpha_results):
        line = line.strip()
        fields = line.split()

        if i==0:
            sp_col = fields.index("species")
            continue

        if group not in group2out:
            print(line)
        elif fields[sp_col]!=group2out[group]:
            print(line)

