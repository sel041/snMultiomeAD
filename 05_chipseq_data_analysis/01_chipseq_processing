s=$1
genome=$2

fastqc_dir=${proj_dir}/fastqc
fastq_dir=${proj_dir}/rawdata
map_dir=${proj_dir}/mapping
trim_dir=${proj_dir}/trimmed
bw_dir=${proj_dir}/bw
log_dir=${proj_dir}/log
peak_dir=${proj_dir}/macs2

mm10_bowtie2=~/bowtie2_indexes/mm10
hg38_bowtie2=~/bowtie2_indexes/hg38
mm10_bed=~/gencode.vM25.annotation.gtf
hg38_bed=~/hg38.genes.gtf

if [ -z ${genome+x} ]; then echo "genome is not specified. Default to mm10"; genome="mm10"; fi
if [ $genome == "mm10" ]; then ref=${mm10_bowtie2}; macs2_genome="mm"; TSS_bed=${mm10_bed}; fi
if [ $genome == "hg38" ]; then ref=${hg38_bowtie2}; macs2_genome="hs"; TSS_bed=${hg38_bed}; fi

fastqc -t 16 -o ${fastqc_dir} ${fastq_dir}/${s}_R1_001.fastq.gz
fastqc -t 16 -o ${fastqc_dir} ${fastq_dir}/${s}_R2_001.fastq.gz

trim_galore -q 20 --paired ${fastq_dir}/${s}_R1_001.fastq.gz ${fastq_dir}/${s}_R2_001.fastq.gz -o ${trim_dir}

# mapping
(~/bowtie2-2.5.4-linux-x86_64/bowtie2 -x ${ref} -1 ${trim_dir}/${s}_R1_001_val_1.fq.gz -2 ${trim_dir}/${s}_R2_001_val_2.fq.gz --no-unal -p 16 -S ${map_dir}/${s}_${genome}.sam) 2>${map_dir}/${s}.log
samtools sort -@ 16 -T ${map_dir} -o ${map_dir}/${s}_${genome}_sorted.bam ${map_dir}/${s}_${genome}.sam
if [[ -f "${map_dir}/${s}_${genome}_sorted.bam" ]]; then rm ${map_dir}/${s}_${genome}.sam; fi

# deduplicate
java -Xmx8G -XX:ParallelGCThreads=3 -jar ~/picard-2.25.7/picard.jar MarkDuplicates I=${map_dir}/${s}_${genome}_sorted.bam TMP_DIR=${map_dir} METRICS_FILE=${log_dir}/${s}_${genome}_dup.qc VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=TRUE O=${map_dir}/${s}_${genome}_sorted_rmdup.bam REMOVE_DUPLICATES=true

# genreate bw, TSS enrichment
samtools index ${map_dir}/${s}_${genome}_sorted_rmdup.bam
bamCoverage -b ${map_dir}/${s}_${genome}_sorted_rmdup.bam -o ${bw_dir}/${s}_${genome}_sorted_rmdup.bw -p max --normalizeUsing RPKM
computeMatrix reference-point --referencePoint TSS -b 2000 -a 2000 -R ${TSS_bed} -S ${bw_dir}/${s}_${genome}_sorted_rmdup.bw --skipZeros -o ${bw_dir}/${s}_${genome}_sorted.matrix.gz -p 16
plotProfile -m ${bw_dir}/${s}_${genome}_sorted.matrix.gz -o ${bw_dir}/${s}_${genome}_sorted.profile.png --refPointLabel "TSS"
