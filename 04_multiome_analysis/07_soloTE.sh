# WT
python SoloTE_pipeline.py --threads 32 --bam ~/rna/bamfiles/Microglia/Microglia_WT_3M.bam --teannotation mm10_rmsk.bed --outputprefix Microglia_WT_3M  --outputdir ~/SoloTE_5xFAD
python SoloTE_pipeline.py --threads 32 --bam ~/rna/bamfiles/Microglia/Microglia_WT_9M.bam --teannotation mm10_rmsk.bed --outputprefix Microglia_WT_9M  --outputdir ~/SoloTE_5xFAD
python SoloTE_pipeline.py --threads 32 --bam ~/rna/bamfiles/Microglia//Microglia_WT_18M.bam --teannotation mm10_rmsk.bed --outputprefix Microglia_WT_18M  --outputdir ~/SoloTE_5xFAD
# 5xFAD
python SoloTE_pipeline.py --threads 32 --bam ~/rna/bamfiles/Microglia/Microglia_AD_3M.bam --teannotation mm10_rmsk.bed --outputprefix Microglia_AD_3M  --outputdir ~/SoloTE_5xFAD
python SoloTE_pipeline.py --threads 32 --bam ~/rna/bamfiles/Microglia/Microglia_AD_9M.bam --teannotation mm10_rmsk.bed --outputprefix Microglia_AD_9M  --outputdir ~/SoloTE_5xFAD
python SoloTE_pipeline.py --threads 32 --bam ~/rna/bamfiles/Microglia/Microglia_AD_18M.bam --teannotation mm10_rmsk.bed --outputprefix Microglia_AD_18M  --outputdir ~/SoloTE_5xFAD

# H-151
# WT
python SoloTE_pipeline.py --threads 32 --bam ~/rna/bamfiles/Microglia_h151/Microglia_C57_DMSO.bam --teannotation mm10_rmsk.bed --outputprefix Microglia_C57_DMSO  --outputdir ~/SoloTE_5xFAD
python SoloTE_pipeline.py --threads 32 --bam ~/rna/bamfiles/Microglia_h151/Microglia_C57_H151.bam --teannotation mm10_rmsk.bed --outputprefix Microglia_C57_H151  --outputdir ~/SoloTE_5xFAD
# 5xFAD
python SoloTE_pipeline.py --threads 32 --bam ~/rna/bamfiles/Microglia_h151/Microglia_5XFAD_DMSO.bam --teannotation mm10_rmsk.bed --outputprefix Microglia_5XFAD_DMSO  --outputdir ~/SoloTE_5xFAD
python SoloTE_pipeline.py --threads 32 --bam ~/rna/bamfiles/Microglia_h151/Microglia_5XFAD_H151.bam --teannotation mm10_rmsk.bed --outputprefix Microglia_5XFAD_H151  --outputdir ~/SoloTE_5xFAD
