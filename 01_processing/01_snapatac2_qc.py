import snapatac2 as snap
import h5py
import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
import polars as pl
import itertools
import pathlib
import re
import os
from matplotlib import pyplot as plt

# load metadata containing library name 
metadata=pd.read_excel('~/metadata.xlsx')

name = metadata['library_id'].tolist()
h5ad_files=[]
for name in name:
    h5ad_files.append((name,'~/h5_files/'+name+'.h5ad'))

data = snap.AnnDataSet(h5ad_files, filename='mouseAD_snapatac.h5ads')
data

data.obs['tsse'] = data.adatas.obs['tsse']
data.obs['frac_dup'] = data.adatas.obs['frac_dup']
data.obs['n_fragment'] = data.adatas.obs['n_fragment']
data.obs['frac_mito'] = data.adatas.obs['frac_mito']
data.obs['doublet_score'] = data.adatas.obs['doublet_score']

df=pd.DataFrame(data.obs_names)
df['sample']=data.obs['sample'] 
df['tsse']=data.obs['tsse'] 
df['frac_dup']=data.obs['frac_dup'] 
df['n_fragment']=data.obs['n_fragment']
df['frac_mito']=data.obs['frac_mito'] 
df['doublet_score']=data.obs['doublet_score']

df.to_csv('~/atac_metadata.csv.gz',index=False)

data.close()
data
