#!/usr/bin/env bash

# Arguments for script
while getopts r:c:s:f:z flag
do
	case "${flag}" in
		r) ref_name=${OPTARG};;  # Genome to be used as reference sequence
		c) spname=${OPTARG};;  # Alias/code name for experiment
		s) seed=${OPTARG};;  # RNG seed
		f) features=${OPTARG};;  # Feature(s) to regress on
		z) is_zw='true' ;;       # Is ZW system
	esac
done

qu_fn="qu/a_unrest.$ref_name.$spname.sh"

if [ "$is_zw" == true ]; then
    zflag="-z"
else
	zflag=""
fi;

# Make sure to not repeatedly write on top of already existing command files
find qu/ -name $( basename ${qu_fn} ) -delete;

echo "./alpha_from_subrates.regression.R -r $ref_name -c $spname -f $features -s $seed -p -g -u ${zflag} -o $spname.$ref_name" > $qu_fn

# Send job
sbatch --account="palab" -t 01:30:00 -c 1 --mem-per-cpu=1G -J $ref_name.$spname.a_unrest -o out/$ref_name.$spname.a_unrest --wrap="bash $qu_fn";
