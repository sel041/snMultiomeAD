#Fragment files to bigwig

#!/bin/bash

for x in bed_files_new/stdchr/*_stdchr.bed ; do 
  type=$(echo $x | sed 's=bed_files_new/stdchr/==g' | sed 's=_stdchr.bed==g');
if [[ -n "$type" ]]; then
echo $type
fi

~/bedtools2/bin/bedToBam -i bed_files_new/stdchr/${type}_stdchr.bed -g ~/mm10.chrom.sizes > bam/${type}.bam

~/samtools-1.19.2/samtools sort bam/${type}.bam > bam/${type}_sorted.bam

rm bam/${type}.bam

~/samtools-1.19.2/samtools index bam/${type}_sorted.bam

bamCoverage --normalizeUsing RPKM -b bam/${type}_sorted.bam -o bigwig/${type}.bw

done
