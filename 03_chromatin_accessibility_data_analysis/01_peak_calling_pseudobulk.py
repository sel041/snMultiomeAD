#peak calling for pseudobulk bed files

#!/bin/bash

dir="../bed_files_new/"

counter=0
for x in ../bed_files_new/*.bed; do
    if [ $counter -lt 104 ]; then
    type=$(echo $x | sed 's=../bed_files_new/=  =g' | cut -f2 | sed 's=.bed==g'); 
   echo "$type"
        ((counter++))
    else
        break
    fi

mkdir q_05/${type}

grep '^chr' ${dir}${type}.bed > ${dir}stdchr/${type}_stdchr.bed

macs2 callpeak --treatment ${dir}stdchr/${type}_stdchr.bed --shift -75 --ext 150 -g mm --outdir q_05/${type} --name ${type} -q 0.05 --call-summits --nomodel -f BED

done
