# chromBPNet
bedtools slop -i ~/chrombpnet/data/Downloads/mm10.blacklist.bed  -g ~/chrombpnet/data/Downloads/mm10.chrom.sizes -b 1057 > ~/chrombpnet/data/Downloads/temp.bed

bedtools intersect -v -a Microglia_p_0.05_peaks.narrowPeak.bed -b ~/chrombpnet/data/Downloads/temp.bed > peaks_Microglia_p_0.05_peaks_no_blacklist.bed

head -n 21  ~/chrombpnet/data/Downloads/mm10.chrom.sizes  >  ~/chrombpnet/data/Downloads/mm10.chrom.subset.sizes

chrombpnet prep splits -c ~/chrombpnet/data/Downloads/mm10.chrom.subset.sizes -tcr chr1 chr3 chr6 -vcr chr8 chr19 -op splits/fold_0

chrombpnet prep nonpeaks -g ~/chrombpnet/data/Downloads/mm10.fa -p peaks_Microglia_p_0.05_peaks_no_blacklist.bed -c  ~/chrombpnet/data/Downloads/mm10.chrom.sizes  -fl splits/fold_0.json -br ~/chrombpnet/data/Downloads/mm10.blacklist.bed -o chrombpnet_nonpeaks_dam_ifn

chrombpnet pipeline \
        -ifrag 5XFAD_H151_DAM-IFN_stdchr.bed \
        -d "ATAC" \
        -g ~/chrombpnet/data/Downloads/mm10.fa \
        -c ~/chrombpnet/data/Downloads/mm10.chrom.sizes \
        -p peaks_Microglia_p_0.05_peaks_no_blacklist.bed   \
        -n chrombpnet_nonpeaks_dam_ifn_negatives.bed  \
        -fl splits/fold_0.json \
        -b ~/chrombpnet/data/Downloads/bias_model/models/bias_model.h5 \
        -o chrombpnet_DAM_IFN_H151/ \

chrombpnet pipeline \
        -ifrag 5XFAD_DMSO_DAM-IFN_stdchr.bed  \
        -d "ATAC" \
        -g ~/chrombpnet/data/Downloads/mm10.fa \
        -c ~/chrombpnet/data/Downloads/mm10.chrom.sizes \
        -p peaks_Microglia_p_0.05_peaks_no_blacklist.bed   \
        -n chrombpnet_nonpeaks_dam_ifn_negatives.bed  \
        -fl splits/fold_0.json \
        -b ~/chrombpnet/data/Downloads/bias_model/models/bias_model.h5 \
        -o chrombpnet_DMSO_IFN_H151/ \

chrombpnet pred_bw  -bm chrombpnet_DAM_IFN_H151/models/bias_model_scaled.h5 -cm chrombpnet_DAM_IFN_H151/models/chrombpnet.h5 -cmb chrombpnet_DAM_IFN_H151/models/chrombpnet_nobias.h5 -r Microglia_peaks.narrowPeak.bed -g ~/chrombpnet/data/Downloads/mm10.fa -c ~/chrombpnet/data/Downloads/mm10.chrom.sizes -op pred_bw_Micro_chrombpnet_DAM_IFN_H151

chrombpnet contribs_bw  -m chrombpnet_DAM_IFN_H151/models/chrombpnet_nobias.h5 -r Microglia_peaks.narrowPeak.bed -g ~/chrombpnet/data/Downloads/mm10.fa -c ~/chrombpnet/data/Downloads/mm10.chrom.sizes -op  contrib_bw_Micro_chrombpnet_DAM_IFN_H151_with_stats -os DAM_IFN_H151_stats

chrombpnet contribs_bw  -m chrombpnet_DMSO_IFN_H151/models/chrombpnet_nobias.h5 -r Microglia_peaks.narrowPeak.bed -g ~/chrombpnet/data/Downloads/mm10.fa -c ~/chrombpnet/data/Downloads/mm10.chrom.sizes -op  contrib_bw_Micro_chrombpnet_DAM_IFN_DMSO

chrombpnet contribs_bw  -m chrombpnet_DAM_IFN_H151/models/chrombpnet_nobias.h5 -r Microglia_peaks.narrowPeak.bed -g ~/chrombpnet/data/Downloads/mm10.fa -c ~/chrombpnet/data/Downloads/mm10.chrom.sizes -op  chrombpnet_stats/contrib_bw_Micro_chrombpnet_DAM_IFN_H151_with_stats -os DAM_IFN_H151_stats

chrombpnet contribs_bw  -m chrombpnet_DMSO_IFN_H151/models/chrombpnet_nobias.h5 -r Microglia_peaks.narrowPeak.bed -g ~/chrombpnet/data/Downloads/mm10.fa -c ~/chrombpnet/data/Downloads/mm10.chrom.sizes -op  chrombpnet_stats/contrib_bw_Micro_chrombpnet_DMSO_IFN_H151_with_stats -os DMSO_IFN_H151_stats


srun -N 1 -n 2  --gpus 1   --mem 100G -t 15:00:00 -p rtx6000  -q condo-gpu  -A csd788 --pty bash
  
chrombpnet pipeline \
        -ifrag 5XFAD_DMSO_DAM-IFN_downsampled_stdchr.bed  \
        -d "ATAC" \
        -g ~/chrombpnet/data/Downloads/mm10.fa \
        -c ~/chrombpnet/data/Downloads/mm10.chrom.sizes \
        -p peaks_Microglia_p_0.05_peaks_no_blacklist.bed   \
        -n chrombpnet_nonpeaks_dam_ifn_negatives.bed  \
        -fl splits/fold_0.json \
        -b ~/chrombpnet/data/Downloads/bias_model/models/bias_model.h5 \
        -o chrombpnet_DMSO_IFN_H151_downsampled/ \

chrombpnet contribs_bw  -m chrombpnet_DMSO_IFN_H151_downsampled/models/chrombpnet_nobias.h5 -r Microglia_peaks.narrowPeak.bed -g ~/chrombpnet/data/Downloads/mm10.fa -c ~/chrombpnet/data/Downloads/mm10.chrom.sizes -op  chrombpnet_stats/contrib_bw_Micro_chrombpnet_DMSO_IFN_H151_downsampled_with_stats -os DMSO_IFN_H151_stats_downsampled
