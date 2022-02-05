#!/usr/bin/env bash

# Arguments for script
######
# First two arguments should be reference (ref_name) and 
# experiment name (spname) *without* flags.
# All other arguments should be specified afterwards with flags.
argArray=( "$@" )

ref_name=${argArray[0]}
spname=${argArray[1]}

# while getopts r:c:s:f:o:zm12 flag
# do
# 	case "${flag}" in
# 		r) ref_name=${OPTARG};;  # Genome to be used as reference sequence
# 		c) spname=${OPTARG};;  # Alias/code name for experiment
# 		s) seed=${OPTARG};;  # RNG seed
# 		o) out_prefix=${OPTARG};;  # Output prefix
# 		f) features=${OPTARG};;  # Feature(s) to regress on
# 		z) is_zw='true' ;;       # Is ZW system
# 		m) report_mean='true' ;; # Also report alpha estimated from mean across windows (no regression)
# 		1) odd_chrom='true' ;; # Only analyze odd windows
# 		2) even_chrom='true' ;; # Only analyze even windows
# 	esac
# done

qu_fn="qu/a_unrest.$ref_name.$spname.sh"

# if [ "$is_zw" == true ]; then
#     zflag="-z"
# else
# 	zflag=""
# fi;

# if [ "$report_mean" == true ]; then
#     mflag="-m"
# else
# 	mflag=""
# fi;

# if [ "$odd_chrom" == true ]; then
#     oddflag="-m"
# else
# 	mflag=""
# fi;

# Make sure to not repeatedly write on top of already existing command files
find qu/ -name $( basename ${qu_fn} ) -delete;

echo "./alpha_from_subrates.regression.R -r $ref_name -c $spname ${argArray[@]:2} -p -g -u" > $qu_fn

# Send job
sbatch --account="palab" -t 01:30:00 -c 1 --mem-per-cpu=1G -J $ref_name.$spname.a_unrest -o out/$ref_name.$spname.a_unrest --wrap="bash $qu_fn";
