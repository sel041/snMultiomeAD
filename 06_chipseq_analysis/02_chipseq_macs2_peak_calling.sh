source activate ~/conda_envs/ATACseqQC

cd ~/macs2

macs2 callpeak -t ../mapping/EH_77_S35_mm10_sorted_rmdup.bam \
               -c ../mapping/EH_74_S32_mm10_sorted_rmdup.bam \
               -f BAMPE -g mm -n SIMA9_CTRL_IRF3 -q 0.05 --nomodel --outdir SIMA9_CTRL_IRF3
