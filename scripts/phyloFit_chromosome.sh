#!/usr/bin/env bash

# Arguments for script
while getopts w:l:r:c:o:x:y:n:m flag
do
    case "${flag}" in
        w) region_file=${OPTARG};;  # Bed files with coordinates of windows to extract
        l) species_file=${OPTARG};;  # Comma separated list of genomes to keep when filtering
        r) ref_name=${OPTARG};;  # Genome to be used as reference sequence
        c) spname=${OPTARG};;  # Alias/code name for experiment, whatever string one likes/prefers
        o) out_spname_list=${OPTARG};;  # Alias/code name for subexperiment, whatever string one likes/prefers
        x) chrX=${OPTARG};;  # String for chrX in reference sequence
        y) chrY=${OPTARG};;  # String for chrY in reference sequence
        n) newick_file=${OPTARG};; # Newick file for phylofit runs
        m) cpg_mask='true';; # Filter positions that are CpG sites in any species
    esac
done

# Important directories
maf_path="/moto/palab/projects/male_mutation_bias_XA/MAFs_2021";
region_path="/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes/regions/";
maf_output_path="/moto/palab/projects/male_mutation_bias_XA/MAFs_2021";

# Make new directory if maf_output_path!=maf_path
mkdir -p $maf_output_path/$ref_name;

# Make sure to not repeatedly write on top of already existing command files
find qu/ -name filter.*.$ref_name.$spname.sh -delete;

# Get largest autosomal chromosome name
chr_sizes=$( dirname $region_file )/$ref_name.chrom_sizes.tab

chrom=$( awk -v chrX=$chrX -v chrY=$chrY '{if($1 != chrX && $1 != chrY){print $0}}' $chr_sizes | sort -nk 2 | tail -n 1 | cut -f 1 )

if [ -z $chrom ] ; then
    echo "ERROR: No chromosome found?"
    exit 1
fi

cat_maf=$maf_output_path/$ref_name/$chrom/$chrom.all_windows.$spname.filtered.maf;
cat_bed=$maf_output_path/$ref_name/$chrom/$chrom.all_windows.$spname.filtered.bed;
for fn in $cat_maf $cat_bed; do # Remove if pre-existing
    if [ -e $fn ]; then rm -f $fn; fi
done


# Iterate over windows and cat the desired MAFs
while read chunk;do 
    
    read chrom start end <<<$(echo $chunk);
    mkdir -p $maf_output_path/$ref_name/$chrom;

    cur_maf=$maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$spname.filtered.maf;
    cur_bed=$maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$spname.filtered.bed;
    
    cat $cur_maf >> $cat_maf
    cat $cur_bed >> $cat_bed
done < <( fgrep $chrom$'\t' $region_file );


# 6 replicates of phylofit on cat maf
for replicate in {1..6}; do 
    qu_fn="qu/phylo_chrom.$chrom.$ref_name.$spname.$replicate.sh"
    printf "phyloFit -r --msa-format MAF --tree $newick_file --subst-mod UNREST -Z $cat_maf -e $maf_output_path/$ref_name/$chrom/$chrom.all_windows.$spname.$replicate.filtered.err_est.txt -o $maf_output_path/$ref_name/$chrom/$chrom.all_windows.$spname.$replicate.filtered --EM --precision MED;\n" > $qu_fn;

    #Send jobs
    sbatch --account="palab" -t 00:30:00 --mem-per-cpu=1G -c 1 -J $ref_name.$chrom.$spname.phylo_chrom -o out/$ref_name.$chrom.$spname.$replicate.phylo_chrom -e out/$ref_name.$chrom.$spname.$replicate.phylo_chrom.err --wrap="bash $qu_fn";
done;
