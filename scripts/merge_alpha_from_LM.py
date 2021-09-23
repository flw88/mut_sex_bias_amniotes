#!/usr/bin/env python

import sys
from optparse import OptionParser
import numpy as np

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-o", "--outgroups", action="store", type="str", dest="outgroup_file",help="Tab separated file with name of group and species that is an outgroup")
parser.add_option("-g", "--group", action="store", type="str", dest="group_name",help="Name of the group")
(options, args) = parser.parse_args()

# Read outgroups
with open(options.outgroup_file, "r") as out_f:
    group2out = {line.split()[0]:line.split()[1] for line in out_f}

# Read alpha from LM
group = options.group_name
directory="alphas/"
with open("{}/{}.LM.tsv".format(directory, group),"r") as alpha_results:
    for i,line in enumerate(alpha_results):
        if i==0:
            continue
        line = line.strip()
        fields = line.split()
        if group not in group2out:
            print(line)
        elif fields[-2]!=group2out[group]:
            print(line)

