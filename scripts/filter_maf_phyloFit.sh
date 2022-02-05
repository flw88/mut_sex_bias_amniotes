#!/usr/bin/env bash

# Arguments for script
while getopts w:l:r:c:o:x:y:n:s:ma flag
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
    s) region_suffix=${OPTARG};; # Suffix for filtering file
    m) cpg_mask='true';; # Filter positions that are CpG sites in any species
    a) all_align='true';; # All species must align
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

# Iterate over windows and send prepare jobs for each chromosome
while read chunk;do 
    
    read chrom start end <<<$(echo $chunk);
    mkdir -p $maf_output_path/$ref_name/$chrom;
    
    # If current chromosome is the X/Z
    if [ "$chrom" == "$chrX" ];then 
        lex="-X";
    elif [ "$chrom" == "$chrY" ]; then
        lex="-Y";
    else 
        lex="";
    fi;  

    # If requiring that all species align to keep sequence
    if [ "$all_align" == true ]; then argA="-a"; else argA=""; fi;

    # If masking CpGs
    if [ "$cpg_mask" == true ]; then cpgmask="-c"; else cpgmask=""; fi;
    
    # Useful variables for the command
    regions_to_mask=$region_path/$ref_name/$ref_name.$chrom.$region_suffix.bed;
    input_maf=$maf_path/$ref_name/$chrom/$chrom.$start.$end.$spname.maf;
    output_maf=$maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$spname.filtered.maf;
    output_bed=$maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$spname.filtered.bed;
    add_a_string="cat - <(printf 'a Just so the python filtering script processes the last alignment block')" # Python script relies on 'a' at the start of line to process an alignment block, so we add one at the end of the file to process the last block

    # Filtering
    printf "maf_parse --features $regions_to_mask -M $ref_name $input_maf --seqs $(cat $species_file) | $add_a_string | python filter_PARs_micros_CpGs.py -l $species_file -p ../data/Species_to_PARs.tsv -m ../data/Species_to_micros.tsv $cpgmask | python keep_species_XYA-synteny.py -l $species_file -b $output_bed -c ../data/Species_to_chromosomes.txt $lex $argA > $output_maf \n" >> qu/filter.$chrom.$ref_name.$spname.sh;
       
    # 6 replicates of phylofit, only if filtered sequence >= 10kb
    for out_spname in $(echo $out_spname_list | tr "," "\n");do 
	newick_file="../trees/$out_spname.nwk";
	printf "if [ -s $output_bed ];then filt_size=\$(awk '{size+=\$3-\$2+1}END{print size}' $output_bed); else filt_size=0; fi; if [ \$filt_size -gt 9999  ]; then for replicate in \$(seq 1 6);do phyloFit -r --msa-format MAF --tree $newick_file --subst-mod UNREST -Z $output_maf -e $maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$out_spname.\$replicate.filtered.err_est.txt -o $maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$out_spname.\$replicate.filtered --EM --precision MED;done;fi \n" >> qu/filter.$chrom.$ref_name.$spname.sh;
    done;
done < $region_file;

#Send jobs
for file in $(ls qu/filter.*.$ref_name.$spname.sh);do
    chrom=$(echo $file | rev | cut -d"/" -f1 | rev | cut -d"." -f2);
    sbatch --account="palab" -t 11:50:00 -J $ref_name.$chrom.$out_spname.filter -o out/$ref_name.$chrom.$out_spname.filter --wrap="bash $file";
done




