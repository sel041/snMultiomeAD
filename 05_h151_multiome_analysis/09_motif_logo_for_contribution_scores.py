#!/bin/python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

%matplotlib inline

import logomaker
import pyBigWig

bw = pyBigWig.open("contrib_bw_Micro_chrombpnet_DAM_IFN_DMSO.counts_scores.bw")
for line in open("../../cgas_sting_genes/microglia_peaks_targeting_h151_microglia_up_p0.05.bed"):
    cols = line.strip().split()
    vals = bw.values(cols[0], int(cols[1]), int(cols[2]))
    

# Define the sequence
sequence = "TCTCTACACTTCCTTCTGCCTGCGGGATAAAAATAAAAGTCTTTCAAAACTCCCAGGTCCCGTGGCCTGC"

# Open the bigWig file and fetch values
bw = pyBigWig.open("contrib_bw_Micro_chrombpnet_DAM_IFN_DMSO.counts_scores.bw")
vals = bw.values("chr7", 128128110, 128128180)

# Close the bigWig file
bw.close()

# Create a DataFrame for the output
data = []

# Map the values to the appropriate column based on the sequence
for pos, (val, base) in enumerate(zip(vals, sequence)):
    row = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    if val is not None:  # Check if the value is not None (in case of gaps)
        row[base] = val
    row['pos'] = pos
    data.append(row)

# Convert the list of dictionaries into a DataFrame
df = pd.DataFrame(data)

# Set 'pos' as the index
df.set_index('pos', inplace=True)

# create Logo object
crp_logo = logomaker.Logo(df)#,
                          #shade_below=.5,
                          #fade_below=.5)#,
                          #font_name='Arial Rounded MT Bold')

# style using Logo methods
crp_logo.style_spines(visible=False)
crp_logo.style_spines(spines=['left', 'bottom'], visible=True)
crp_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

# style using Axes methods
crp_logo.ax.set_ylabel("Contribution score", labelpad=-1)
crp_logo.ax.xaxis.set_ticks_position('none')
crp_logo.ax.xaxis.set_tick_params(pad=-1)
crp_logo.ax.set_ylim(-0.04, 0.15)

plt.savefig("Itgax_enhancer_DMSO.pdf", format="pdf", bbox_inches="tight")
