#!/usr/bin/env bash

spname="Mammals"
ref_name="Homo_sapiens"

out_spname="Mammals"

exclude=true # Exclude CpG-islands

if [ $exclude == "true" ]; then
    cgi_path="/moto/palab/projects/male_mutation_bias_XA/features/Homo_sapiens/cpgIslandExtUnmasked.inverse.bed"
    suffix=nonCGI
else    
    cgi_path="/moto/palab/projects/male_mutation_bias_XA/features/Homo_sapiens/cpgIslandExtUnmasked.bed"
    suffix=CGI
fi
maf_path="/moto/palab/projects/male_mutation_bias_XA/MAFs_2021/Homo_sapiens"
region_file="/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes/scripts/region_files/${ref_name}.1Mb.bed"


# iarr=( `seq 1 22` X )
# for i in ${iarr[@]}; do
#     file="qu/filter_CGI.chr${i}.$ref_name.$spname.sh"
#     cur_maf="${maf_path}/chr${i}/chr${i}.concat.Mammals.filtered.maf"
#     out_maf="${maf_path}/chr${i}/chr${i}.concat.Mammals.filtered.${suffix}.maf"
#     echo "maf_parse -o MAF -g ${cgi_path} ${cur_maf} > ${out_maf}" > $file
#     sbatch --account="palab" -t 11:50:00 -J $ref_name.$i.$out_spname.filter_CGI -o out/$ref_name.$i.$out_spname.filter_CGI --wrap="bash $file";
# done


find qu/ -name filter_CGI.*.$ref_name.$spname.sh -delete;
while read chunk;do 
    read chrom start end <<<$(echo $chunk);
    # if [ "$chrom" == "chrY" ]; then
    #     continue
    # fi

    file="qu/filter_CGI.${chrom}.$ref_name.$spname.sh"
    cur_maf="${maf_path}/${chrom}/${chrom}.${start}.${end}.Mammals.filtered.maf"
    out_maf="${maf_path}/${chrom}/${chrom}.${start}.${end}.Mammals.filtered.${suffix}.maf"
    echo "maf_parse -o MAF -g ${cgi_path} ${cur_maf} > ${out_maf}" >> $file
done < $region_file;

for file in $(ls qu/filter_CGI.*.$ref_name.$spname.sh);do
    chrom=$(echo $file | rev | cut -d"/" -f1 | rev | cut -d"." -f2);
    sbatch --account="palab" -t 11:50:00 -J $ref_name.$chrom.$out_spname.filter_CGI -o out/$ref_name.$chrom.$out_spname.filter_CGI --wrap="bash $file";
done

