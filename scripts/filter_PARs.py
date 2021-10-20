#!/usr/bin/env python                                                                                                                                                                                      

import sys
from optparse import OptionParser

def read_species(filename):
    '''Returns list of species from comma separated file'''
    with open(filename,"r") as fh:
        return [line.strip().split(",") for line in fh][0]

def read_pars(filename, species_list):
    '''Returns 2 dictionaries, relating selected species to their sex 
    chromosomes and PAR coordinates (if available)'''
    pars = {}
    sex_chroms = {}
    with open(filename,"r") as fh:
        for line in fh:
            line = line.strip()
            fields = line.split()
            species, chrom, coordinates = fields
            if species not in species_list:
                continue
            elif coordinates in ["TBD","None","Unk"]:
                continue
            else:
                sex_chroms[species] = chrom
                pars[species] = [list(map(int,p.split("-"))) for p in coordinates.split(",")]

    return sex_chroms, pars

def process_block_line(seq):
    '''Returns species, chrom, position, strand and DNA seq in block line'''
    sp = seq[1].split(".")[0]
    chrom = ".".join(seq[1].split(".")[1:])
    position = int(seq[2])
    strand = seq[4]
    dna = seq[-1]
    return sp, chrom, position, strand, dna

def distance_bw_intervals(r1, r2):
    '''Returns distance between two genomic intervals'''
    x, y = sorted((r1, r2))
    if x[0] <= x[1] < y[0] and all( y[0] <= y[1] for y in (r1,r2)):
        return y[0] - x[1]
    return 0

def intersect_with_par(position, strand, dna_l, par_coordinates):
    '''Return true if there's overlap with current block and species' PAR, 
    else False'''
    start = position if strand=="+" else position-dna_l
    end = position+dna_l if strand=="+" else position
    interval = [start, end]
    for coor in par_coordinates:
        d = distance_bw_intervals(interval, coor)
        if d==0:
            return True
    return False
    

def main(stdin, sp2sexchrom, sp2par):

    block = []
    masked_bps = {sp:0 for sp in sp2par}

    for line in stdin:
        line = line.strip()
        # Alignment block
        if line.startswith("a"):
            # Jump the first block, no information yet
            if len(block)==0:
                continue
                
            masked_block = []
            for seq in block:
                sp, chrom, position, strand, dna = process_block_line(seq)
                # Only mask PAR if species has information
                if sp in sp2par:
                    if chrom!=sp2sexchrom[sp]:
                        masked_block.append(seq)
                        continue
                    bool_par = intersect_with_par(position, strand, len(dna), sp2par[sp])
                    if bool_par:
                        seq[-1] = "N"*len(dna)
                        masked_bps[sp] += len(dna)
                masked_block.append(seq)

            # Print PAR masked block
            sys.stdout.write("a\n")
            sys.stdout.write("\n".join(["\t".join(b) for b in masked_block]) + "\n")
            sys.stdout.write("\n")
            block = []

        # If not in start of block, just keep information and process it when hitting the next block (marked by string 'a')
        if line.startswith("s"):
            block.append(line.split())
    # Print stats
    for sp in masked_bps:
        sys.stderr.write("Masked {} PAR base-pairs in {}\n".format(masked_bps[sp], sp))
    sys.stdout.write("a\n")

if __name__ == "__main__":
    
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", "--sp_list", action="store", type="str", dest="species_list",help="Comma separated file with species names")
    parser.add_option("-p", "--par_file", action="store", type="str", dest="par_file",help="Tab separated file with species name, chrX PAR coordinates")
    (options, args) = parser.parse_args()

    species = read_species(options.species_list)  
    sp2sexchrom, sp2par = read_pars(options.par_file, species)

    main(sys.stdin, sp2sexchrom, sp2par)

