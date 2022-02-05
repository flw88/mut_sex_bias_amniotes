#!/usr/bin/env bash

# Arguments for script
while getopts w:l:r:c:o:x:y:n:s:b flag
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
        s) suffix=${OPTARG};;
        b) biggest='true';; # Only analyze largest chromosome
    esac
done

# Important directories
maf_path="/moto/palab/projects/male_mutation_bias_XA/MAFs_2021";
region_path="/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes/regions/";
maf_output_path="/moto/palab/projects/male_mutation_bias_XA/MAFs_2021";

# Make new directory if maf_output_path!=maf_path
mkdir -p $maf_output_path/$ref_name;

# Make sure to not repeatedly write on top of already existing command files
# find qu/ -name phylo_window.*.$ref_name.$spname.sh -delete;


while read chunk; do 
    read chrom start end <<<$(echo $chunk);
    
    cur_maf=$maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$spname.${suffix}.maf;


    # 6 replicates of phylofit on cat maf
    qu_fn="qu/phylo_window.$chrom.$start.$end.$ref_name.$spname.sh"
    for replicate in {1..6}; do 
        printf "phyloFit -r --msa-format MAF --tree $newick_file --subst-mod UNREST -Z $cur_maf -e $maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$spname.$replicate.${suffix}.err_est.txt -o $maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$spname.$replicate.${suffix} --EM --precision MED;\n" >> $qu_fn;
        printf "phyloFit -r --msa-format MAF --tree $newick_file --subst-mod U2S -Z $cur_maf -e $maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$spname.$replicate.${suffix}.err_est.txt -o $maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$spname.$replicate.${suffix}.U2S --EM --precision MED;\n" >> $qu_fn;

    done;
    #Send jobs
    # sbatch --account="palab" -t 4:00:00 --mem-per-cpu=1G -c 1 -J $ref_name.$chrom.$start.$end.$spname.phylo_window -o out/$ref_name.$chrom.$start.$end.$spname.phylo_window -e out/$ref_name.$chrom.$start.$end.$spname.phylo_window.err --wrap="bash $qu_fn";
done < <( fgrep -e chr1$'\t' -e chr2$'\t' $region_file );

