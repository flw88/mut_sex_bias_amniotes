#!/usr/bin/env python                                                                                                                                                                                      

import sys
from optparse import OptionParser

def read_species(filename):
    '''Returns list of species from comma separated file'''
    with open(filename,"r") as fh:
        return [line.strip().split(",") for line in fh][0]

def parse_species_in_block(block):
    '''Returns list of species in a given block'''
    return [seq[1].split(".")[0] for seq in block]

def parse_speciesChroms_in_block(block):
    '''Returns dictionary of species to chroms in a given block'''
    return {seq[1].split(".")[0]:".".join(seq[1].split(".")[1:]) for seq in block}

def divide_block(block):
    '''Divides alignment block into contiguous segments without gaps'''
    ref_position = int(block[0][2])
    ref_seq = block[0][6]
    all_seqs = [b[6] for b in block]
    filon,non_gapped_regions = [],[]
    all_nucleotides = ["A","G","T","C"]

    # Identify positions in which refseq carries a nucleotide (ie. no-gap no-N)
    for block_position,nucleotide in enumerate(ref_seq):
        if all(n.upper() in all_nucleotides for n in [seq[block_position] for seq in all_seqs]):
            if filon==[]:
                filon.append(block_position)
            else:
                if block_position-1==filon[-1]:
                    filon.append(block_position)
                else:
                    non_gapped_regions.append(filon)
                    filon = [block_position]
    if filon!=[]:
        non_gapped_regions.append(filon)

    # Stop if there are no free gap/N regions
    if len(non_gapped_regions)==0:
        return []

    # Iterate over all free gap/N regions
    total_blocks = []
    for region in non_gapped_regions:
        new_block = []
        for species in block:
            # Original values
            seq = species[-1]
            position = int(species[2])
            block_length = int(species[3])
            strand = species[4]
            move_pos = len(seq[:region[0]].replace("-",""))

            # New values
            new_seq = "".join([seq[i] for i in region])
            new_block_length = sum([new_seq.upper().count(n) for n in all_nucleotides])
            new_position = position+move_pos if strand=="+" else position-move_pos
            new_species = ['s', species[1], str(new_position), str(new_block_length), strand, species[5], new_seq]
            new_block.append(new_species)
        total_blocks.append(new_block)

    return total_blocks

def read_chroms(chrom_file, selected_species):
    '''Returns X/Y strings in chromosome-level assemblies'''

    species2chrom = {}
    
    # Read file and build dictionary
    with open(chrom_file, "r") as fh:
        for line in fh:
            line = line.strip()
            fields = line.split()
            species = fields[0]
            chrx = fields[1]
            chry = fields[2]
            if species in selected_species:
                species2chrom[species] = [chrx, chry]

    return species2chrom

def check_if_mapping_XYA(chrom2species, species2chrom_in_block, X, Y):
    '''Checks if all species' sequences are within the same chromosome type'''

    # Return jump=True if anything is mapping where it shouldn't
    jump = False

    # Note that we iterate over species in chrom2species. chrom2species should only contain
    # species for which at least chrX is known in the assembly.
    for sp in chrom2species:
        true_x = chrom2species[sp][0]
        true_y = chrom2species[sp][1]
        chrom_block = species2chrom_in_block[sp]

        # If in the X, jump if anything doesn't map to the chrX
        if X:
            if true_x!=chrom_block:
                jump = True
                return jump
    
        # If in the Y, jump if anything doesn't map to the chrY
        if Y:
            if true_y!=chrom_block:
                jump = True
                return jump

        # If in the autosomes, jump if anything maps to chrX or chrY
        if not X and not Y:
            if true_x==chrom_block or true_y==chrom_block:
                jump = True
                return jump

    return jump

def main(stdin, selected_species, chrom2species, X, Y, bedoutput, homology):

    bedlist = []
    block = []
    included_bps = 0
    total_bps = 0

    for line in stdin:

        line = line.strip()

        # Alignment block
        if line.startswith("a"):

            # Jump the first block, no information yet
            if len(block)==0:
                continue

            block_size = int(block[0][3])
            total_bps += block_size

            # Print if all species are present in block (in case of -h, else do not check)
            species2chrom_in_block = parse_speciesChroms_in_block(block)
            if not homology:
                for sp in selected_species:
                    if sp not in species2chrom_in_block:
                        species2chrom_in_block[sp] = "NaN"

            if all(sp in species2chrom_in_block.keys() for sp in selected_species):
                
                # Mapping to XYA when it shouldn't? Jump if that's the case
                jump = check_if_mapping_XYA(chrom2species, species2chrom_in_block, X, Y)

                if jump:
                    block = []
                    continue

                # Remove positions in which there's a non-nucleotide character & gaps
                new_blocks = divide_block(block)                

                # Filter out blocks where no data is left
                if len(new_blocks)==0:
                    block = []
                    continue
                    
                # Print all the divided blocks
                for nb in new_blocks:
                    block_size = int(nb[0][3])
                    included_bps += block_size
                    start = int(nb[0][2])
                    end = start+block_size
                    chrom = nb[0][1]
                    bedlist.append([chrom,str(start+1),str(end)])
                    sys.stdout.write("a\n"+"\n".join(["\t".join(b) for b in nb if b[1].split(".")[0] in selected_species])+"\n\n")
            
            block = []

        # If not in start of block, just keep information and process it when hitting the next block (marked by string 'a')
        if line.startswith("s"):
            block.append(line.split())


    # Write the bed file
    with open(bedoutput,"w") as bf:
        for region in bedlist:
            bf.write("\t".join(region) + "\n")

    # Print stats
    sys.stderr.write("Total bps = {}\n".format(total_bps))
    sys.stderr.write("Included bps = {}\n".format(included_bps))
    sys.stderr.write("Fraction of included bps = {:.3f}\n".format(included_bps/float(total_bps)))
        

if __name__ == "__main__":
    
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", "--list", action="store", type="str", dest="species_list",help="Comma separated file")
    parser.add_option("-c", "--chrom", action="store", type="str", dest="chrom_species",help="Tab separated file with species name, chrX and chrY")
    parser.add_option("-X", "--x_chrom", action="store_true", default=False, dest="chrx",help="Is the chunk in the X of the reference sequence?")
    parser.add_option("-Y", "--y_chrom", action="store_true", default=False, dest="chry",help="Is the chunk in the Y of the reference sequence?")
    parser.add_option("-a", "--all", action="store_true", default=False, dest="homology",help="All species must align")
    parser.add_option("-b", "--bed", action="store", type="str", dest="bedoutput",help="bed file output")
    (options, args) = parser.parse_args()

    species = read_species(options.species_list)  
    chrom2species = read_chroms(options.chrom_species, species)

    main(sys.stdin, species, chrom2species, options.chrx, options.chry, options.bedoutput, options.homology)
