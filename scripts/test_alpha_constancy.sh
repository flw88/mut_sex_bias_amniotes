#!/usr/bin/env bash

# Arguments for script
while getopts r:c:o: flag
do
	case "${flag}" in
		r) ref_name=${OPTARG};;  # Genome to be used as reference sequence
		c) spname=${OPTARG};;  # Alias/code name for experiment
		o) out_spname_list=${OPTARG};;  # Alias/code name for subexperiment, whatever string one likes/prefers
	esac
done

qu_fn="qu/test_alpha.$ref_name.$spname.sh"

# Make sure to not repeatedly write on top of already existing command files
find qu/ -name $( basename ${qu_fn} ) -delete;

echo "./test_alpha_constancy.R -r $ref_name -c $spname -o $out_spname_list" > $qu_fn

# Send job
sbatch --account="palab" -t 00:30:00 -c 1 --mem-per-cpu=6G -J $ref_name.$spname.test_alpha -o out/$ref_name.$spname.test_alpha --wrap="bash $qu_fn";
