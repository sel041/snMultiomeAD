import mudata
import scenicplus
from scenicplus.plotting.dotplot import heatmap_dotplot
import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import anndata

scplus_mdata = mudata.read("scplusmdata.h5mu")
scplus_mdata.uns["direct_e_regulon_metadata"]

cluster_order = ["HMG", "DAM", "IFN"]
group_order = ["WT_3M", "WT_9M", "WT_18M", "AD_3M", "AD_9M", "AD_18M"]

fig = heatmap_dotplot(
  scplus_mudata=scplus_mdata,
  color_modality="direct_gene_based_AUC",
  size_modality="direct_region_based_AUC",
  group_variable="scRNA_counts:final_cluster",
  group_variable_order=cluster_order[::-1],
  eRegulon_metadata_key="direct_e_regulon_metadata",
  color_feature_key="Gene_signature_name",
  size_feature_key="Region_signature_name",
  feature_name_key="eRegulon_name",
  sort_data_by="direct_gene_based_AUC",
  orientation="horizontal",
  figsize=(16, 3.5)
)
