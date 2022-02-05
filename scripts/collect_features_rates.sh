#!/usr/bin/env bash

# Arguments for script
while getopts w:r:c:mtsfn flag
do
    case "${flag}" in
        w) region_file=${OPTARG};;  # Bed files with coordinates of windows to extract
        r) ref_name=${OPTARG};;  # Genome to be used as reference sequence
        c) spname=${OPTARG};;  # Alias/code name for experiment, whatever string one likes/prefers, Must have been used in filter step
        s) sub_rates='true' ;; # Collect sub rates
        n) counts='true' ;; # Collect sub counts 
        m) methylation='true' ;; # Collect methylation
        t) rep_timing='true' ;; # Collect replication timing
        f) gc_sp_specific='true' ;; # Collect species-specific GC from maf file
    esac
done

maf_output_path="/moto/palab/projects/male_mutation_bias_XA/MAFs_2021";
features_path="/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes/features/$ref_name";

# Make sure to not repeatedly write on top of already existing command files
find qu/ -name collect.*.$ref_name.$spname.sh -delete;

# Iterate over windows and send prepare jobs for each chromosome
while read chunk;do 

    read chrom start end <<<$(echo $chunk);

    # Important files
    bed_file="$maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$spname.filtered.bed";
    bed_file=$(echo $bed_file | sed 's/Birds_g[1-9]/Birds/g');
    tr_tab="tr '\n' '\t'"
    cmd="";

    # Size of unfiltered sequence
    fasta_file="$features_path/$ref_name.fa.gz";
    if [ -s "$bed_file" ];then
	filt_size=$(awk '{size+=$3-$2+1}END{print size}' $bed_file);
    else
	filt_size=0;
    fi;
    unfilt_size=$(echo "$end-$start" | bc)
    cmd=$(echo "printf $chrom'\t'$start'\t'$end'\t'$filt_size'\t'$unfilt_size'\t'")

    # GC content
    #if [ "$gc_flag" == true ]; then
    #if [ -s "$bed_file" ];then
    #    cmd=$(echo $cmd"awk '{print \$1\":\"\$2\"-\"\$3}' $bed_file | sed \"s/$ref_name.//g\" | xargs samtools faidx $fasta_file | python gc_content_from_fasta.py -c $chrom -s $start -e $end | $tr_tab");
    #    cmd=$(echo $cmd";samtools faidx $fasta_file $chrom:${start}-${end} | python gc_content_from_fasta.py -c $chrom -s $start -e $end | cut -f4,5");
    #else
    #    filt_size=0;
    #    fasta_file="$features_path/$ref_name.fa.gz";
    #    cmd=$(echo $cmd"echo '' | python gc_content_from_fasta.py -c $chrom -s $start -e $end | $tr_tab");
    #    cmd=$(echo $cmd";samtools faidx $fasta_file $chrom:${start}-${end} | python gc_content_from_fasta.py -c $chrom -s $start -e $end | cut -f4,5");
    #fi;
    #fi;

    # Substitution rates in 1Mb windows
    if [ "$sub_rates" == true ]; then
    if [ -z "$cmd" ];then trail=""; else trail=" | $tr_tab;";fi
    maf="$maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$spname.filtered.maf";
    if compgen -G "$maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$spname.*.filtered*mod" > /dev/null;then
        high_mod=$(ls $maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$spname.*.filtered*mod | while read rep;do grep LNL $rep | awk -v r="$rep" '{print $NF"\t"r}';done | sort -k1 -V | head -1 | cut -f2);
        cmd=$(echo $cmd"$trail ./tree_from_mod.sh $high_mod")
    else
        cmd=$(echo $cmd"$trail printf 'nan\tnan\n'")
    fi;
    fi;

    # Substitution counts (numerator and denominator) in 1Mb windows
    if [ "$counts" == true ]; then
    if [ -z "$cmd" ];then trail=""; else trail=" | $tr_tab;";fi
    mut_types=( "C>A" "C>G" "C>T" "T>A" "T>C" "T>G" )
    mut_types_str=$( echo ${mut_types[@]} | tr ' ' ',' )
    ntypes=${#mut_types[@]}
    if compgen -G "$maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$spname.*.filtered*mod" > /dev/null;then
        high_mod=$(ls $maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$spname.*.filtered*mod | while read rep;do grep LNL $rep | awk -v r="$rep" '{print $NF"\t"r}';done | sort -k1 -V | head -1 | cut -f2);
        high_exptotsub=${high_mod%.mod}.exptotsub
        cmd=$(echo $cmd"$trail ./subcounts_from_exptotsub.R $high_mod $high_exptotsub '$mut_types_str'")
    else
        n_nans=$( echo "$ntypes + 2" | bc ) # How many nans to write
        nanstr=""
        for i in $( seq $n_nans ); do
            nanstr="${nanstr}nan\t"
        done
        nanstr="${nanstr%\\t}\n"
        cmd=$(echo $cmd"$trail printf '$nanstr'")
    fi;
    fi;

    # Replication timing
    if [ "$rep_timing" == true ]; then
    if [ -z "$cmd" ];then trail=""; else trail=" | $tr_tab;";fi
        reptiming_bed="$features_path/$ref_name.repTiming.bed";
        #reptiming_bed="$features_path/Koren_et_al_Table_S2.lift-hg38.bed";
        #reptiming_bed="$features_path/Mus_musculus.repTiming.Spermatogonia.bed";
        #reptiming_bed="$features_path/Mus_musculus.repTiming.PGC_female.bed";
        #reptiming_bed="$features_path/Mus_musculus.repTiming.PGC_male.bed";
        cmd=$(echo $cmd"$trail awk '{print \$1\"\t\"\$2-1\"\t\"\$3}' $bed_file | bedtools merge -i - | sed \"s/$ref_name.//g\" | bedtools intersect -a $reptiming_bed -b - | python weighted_average_replication_timing.py | $tr_tab");
        cmd=$(echo $cmd";printf '${chrom}\t${start}\t${end}\n' | bedtools intersect -a $reptiming_bed -b - | python weighted_average_replication_timing.py");
    fi;

    # CpG methylation
    if [ "$methylation" == true ]; then
    if [ -z "$cmd" ];then trail=""; else trail=" | $tr_tab;";fi
    if [ "$filt_size" -gt 0 ];then
        methylation_bed="$features_path/$ref_name.methylation.bed.gz";
        cmd=$(echo $cmd"$trail tabix $methylation_bed $chrom:${start}-${end} | bedtools intersect -a - -b <(sed \"s/$ref_name.//g\" $bed_file | bedtools merge -i - ) | awk '{if ((\$NF>0.5)) print}' | wc -l | awk '{print (\$1/2)/$filt_size}' | $tr_tab");
        cmd=$(echo $cmd";tabix $methylation_bed $chrom:${start}-${end} | awk '{if ((\$NF>0.5)) print}' | wc -l | awk '{print (\$1/2)/$unfilt_size}'");
    else
        cmd=$(echo $cmd"$trail printf 'nan\tnan\n'")
    fi;
    fi;

    # GC specific-specific
    if [ "$gc_sp_specific" == true ]; then
    if [ -z "$cmd" ];then trail=""; else trail=" | $tr_tab;";fi
        maf="$maf_output_path/$ref_name/$chrom/$chrom.$start.$end.$spname.filtered.maf";
	maf=$(echo $maf | sed 's/Birds_g[1-9]/Birds/g');
	maf=$(echo $maf | sed 's/Serpentes5/Serpentes/g');
        cmd=$(echo $cmd" $trail cat $maf | python gc_content_from_maf.py -s ../data/${spname}.txt");
    fi;

    # Write command to chromosome file
    echo $cmd >> qu/collect.$chrom.$ref_name.$spname.sh;

done < $region_file;

# Send jobs
for file in $(ls qu/collect.*.$ref_name.$spname.sh | sort -V);do
    chrom=$(echo $file | rev | cut -d"/" -f1 | rev | cut -d"." -f2);
    
    # Headers:
    printf "chrom\tstart\tend\tsize\tsize_unfiltered" | gzip > merged_features/$chrom.$spname.$ref_name.features.txt.gz;

    # If GC
    if [ "$gc_flag" == true ]; then
    printf "gc\tsize\tgc_unfiltered\tsize_unfiltered" | gzip > merged_features/$chrom.$spname.$ref_name.features.txt.gz;
    fi;

    # Substitution rates in 1Mb windows
    if [ "$sub_rates" == true ]; then
    printf "\ttree\ttree_var" | gzip >> merged_features/$chrom.$spname.$ref_name.features.txt.gz;
    fi;

    # Substitution counts in 1Mb windows
    if [ "$counts" == true ]; then
    outstr="\t"
    for mut in ${mut_types[@]}; do # Numerators
        tmpstr="tree_${mut}\t"
        outstr="${outstr}${tmpstr}"
    done
    for b in C T; do # Denominators
        tmpstr="tree_${b}\t"
        outstr="${outstr}${tmpstr}"
    done
    printf "${outstr%\\t}" | gzip >> merged_features/$chrom.$spname.$ref_name.features.txt.gz;
    fi;

    # If replication timing
    if [ "$rep_timing" == true ]; then
    printf "\treptiming\treptiming_unfiltered" | gzip >> merged_features/$chrom.$spname.$ref_name.features.txt.gz;
    fi;

    # If methylation
    if [ "$methylation" == true ]; then
        printf "\tmethylated_cpgs\tmethylated_cpgs_unfiltered" | gzip >> merged_features/$chrom.$spname.features.txt.gz;
    fi;
    
    # GC specific-specific
    if [ "$gc_sp_specific" == true ]; then
        printf "\t" | gzip >> merged_features/$chrom.$spname.$ref_name.features.txt.gz;
	cat ../data/${spname}.txt | tr "\n" "," | sed 's/,/-gc\t/g' | sed 's/[[:space:]]*$//' | gzip >> merged_features/$chrom.$spname.$ref_name.features.txt.gz;
    fi;

    # Add newline in any case
    printf "\n" | gzip >> merged_features/$chrom.$spname.$ref_name.features.txt.gz;
    
    # Send job by chromosome
    sbatch --account="palab" -t 01:00:00 -J $ref_name.$chrom.$spname.collect -o out/$ref_name.$chrom.$spname.collect --wrap="bash $file | gzip >> merged_features/$chrom.$spname.$ref_name.features.txt.gz";
done


