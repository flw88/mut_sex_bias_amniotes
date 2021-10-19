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

ucsc_dir="/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes/ucsc_data/multiz_100way/maf/";
maf_path="/moto/palab/projects/male_mutation_bias_XA/MAFs_2021";
mkdir -p $maf_path $maf_path/$ref_name;

# Get genome names; create sed command string
seq_str=""
sed_cmd="sed "
for sp in $(cat $species_file | tr ',' ' '); do
    gnm=$(grep -m1 $sp ../ucsc_data/genome_names.txt | cut -f 2)
    seq_str="$seq_str,$gnm"

    chrX=$(grep -m1 "^$sp"$'\t' ../data/Species_to_chromosomes.txt | cut -f 2)
    chrY=$(grep -m1 "^$sp"$'\t' ../data/Species_to_chromosomes.txt | cut -f 3)

    if [[ $sp != "Homo_sapiens" && $chrX != "NaN" ]]; then
        sed_cmd=$( echo "$sed_cmd -e 's/s $gnm\.chrX/s $sp.$chrX/' " )
    fi

    if [[ $sp != "Homo_sapiens" && $chrY != "NaN" ]]; then
        sed_cmd=$( echo "$sed_cmd -e 's/s $gnm\.chrY/s $sp.$chrY/' " )
    fi

    sed_cmd=$( echo "$sed_cmd -e 's/s $gnm\./s $sp./' " )
done

# Iterate over windows and send maf jobs for grouped regions
counter=0;
while read chunk;do 
    read chrom start end <<<$(echo $chunk);
    mkdir -p $maf_path/$ref_name/$chrom;
    l=$(echo "${end}-${start}" | bc);
    echo ">&2 echo '$chrom.$start.$end.$spname';" > qu/ucsc2maf/$ref_name.$counter.$spname.sh
    echo "zcat $ucsc_dir/$chrom.maf.gz | maf_parse --seqs $seq_str --start $(expr $start + 1) --end $end -I -E /dev/stdin | \\" >> qu/ucsc2maf/$ref_name.$counter.$spname.sh
    echo "$sed_cmd > $maf_path/$ref_name/$chrom/$chrom.$start.$end.$spname.maf" >> qu/ucsc2maf/$ref_name.$counter.$spname.sh

    let "counter++";
    if (( $counter % $job_grouping == 0 ));then
	sbatch --account="palab" -t 4:00:00 -c 1 -J $ref_name.$counter.$spname -o out/ucsc2maf/$ref_name.$chrom.$counter.$spname --wrap="bash qu/ucsc2maf/$ref_name.$counter.$spname.sh";
    fi;
done < $region_file;

# Also last commands if necessary
if (( $counter % $job_grouping != 0 ));then
    sbatch --account="palab" -t 4:00:00 -c 1 -J $ref_name.$counter.$spname -o out/ucsc2maf/$ref_name.$chrom.$counter.$spname --wrap="bash qu/ucsc2maf/$ref_name.$counter.$spname.sh";
fi;

