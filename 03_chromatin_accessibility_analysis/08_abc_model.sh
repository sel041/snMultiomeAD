#!/bin/bash

for x in mcool_to_bedpe/* ; do 
  type=$(echo $x | grep '^mcool' | sed 's=mcool_to_bedpe/==g');
  echo ${type}


python ~/software/ABC-Enhancer-Gene-Prediction/src/run.neighborhoods.py \
--candidate_enhancer_regions ~/5xfad/atac_mouse/macs2_new/q_05/final_peaks/${type}.bed \
--genes mm10_gene_bounds_protein_coding.txt \
--ATAC ~/5xfad/atac_mouse/bam_files/${type}_sorted.bam \
--expression_table ~/5xfad/gene_expression/cpm/${type}.cpm.tsv \
--chrom_sizes mm10.chrom.sizes \
--ubiquitously_expressed_genes ubiqiutously_expressed_empty.txt \
--cellType ${type} \
--outdir neighborhoods/${type}

done

for x in mcool_to_bedpe/*/chrX ; do 
  type=$(echo $x | grep '^mcool' | sed 's=mcool_to_bedpe/==g' | sed 's=/chrX==g');
#  echo ${type}


mkdir abc_predictions/${type}

python ~/software/ABC-Enhancer-Gene-Prediction/src/predict.py \
--enhancers neighborhoods/${type}/EnhancerList.txt \
--genes neighborhoods/${type}/GeneListfilt.txt \
--HiCdir mcool_to_bedpe/${type} \
--hic_type bedpe \
--chrom_sizes mm10.chrom.sizes \
--hic_resolution 10000 \
--threshold .02 \
--cellType ${type} \
--outdir abc_predictions/${type} \
--make_all_putative

done
