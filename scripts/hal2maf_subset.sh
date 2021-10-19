#!/usr/bin/env bash

# Arguments for script
while getopts w:l:r:c:g:a: flag
do
    case "${flag}" in
        w) region_file=${OPTARG};;  # Bed files with coordinates of windows to extract
        l) species_file=${OPTARG};;  # Comma separated list of genomes to extract
        r) ref_name=${OPTARG};;  # Genome to be used as reference sequence
	c) spname=${OPTARG};;  # Alias/code name for subexperiment, whatever string one likes/prefers
	g) job_grouping=${OPTARG};;  # Number of processes to be run in a single job, pseudo-array 
	a) hal_file=${OPTARG};;  # HAL file from where to convert to MAF
    esac
done

hal_dir="/moto/palab/users/md3914/Male_mutation_bias/ALIGNMENTS/";
maf_path="/moto/palab/projects/male_mutation_bias_XA/MAFs_2021";
mkdir -p $maf_path $maf_path/$ref_name;

# Iterate over windows and send hal2maf job for grouped regions
counter=0;
cmd="";
while read chunk;do 
    read chrom start end <<<$(echo $chunk);
    mkdir -p $maf_path/$ref_name/$chrom;
    l=$(echo "${end}-${start}" | bc);
    cmd=$(echo $cmd"hal2maf $hal_dir/$hal_file $maf_path/$ref_name/$chrom/$chrom.$start.$end.$spname.maf --targetGenomes $(cat $species_file) --refGenome $ref_name --refSequence $chrom --onlyOrthologs --noDupes --start $start --length $l --noAncestors;>&2 echo '$chrom.$start.$end.$spname';");
    let "counter++";
    if (( $counter % $job_grouping == 0 ));then
	sbatch --account="palab" -t 11:55:00 -c 1 -J $ref_name.$counter.$spname -o out/$ref_name.$chrom.$counter.$spname --wrap="$cmd";
	cmd="";
    fi;
done < $region_file;

# Also last commands if necessary
if (( $counter % $job_grouping != 0 ));then
    sbatch --account="palab" -t 11:55:00 -c 1 -J $ref_name.$counter.$spname -o out/$ref_name.$chrom.$counter.$spname --wrap="$cmd";
fi;

