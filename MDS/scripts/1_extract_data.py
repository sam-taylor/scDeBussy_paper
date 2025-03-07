import os
import pandas as pd
import scanpy as sc
import anndata as ad

# Read all the h5 files
base_dir = "../data/MDS/raw"
out_dir = "../data/MDS/results"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
data_files = os.listdir(base_dir)
h5_files = [f for f in data_files if f.endswith('.h5')]
sample_names = [x.split('_')[1] for x in h5_files]

# Read the h5 files
adata_list = []
for i, sample_name in enumerate(sample_names):
    adata = sc.read_10x_h5(base_dir + h5_files[i])
    adata.obs['sample'] = sample_name
    adata.obs_names = adata.obs_names.str.replace('-1', '') + '_' + sample_name
    adata.var_names_make_unique()
    adata_list.append(adata)
adata_combined = ad.concat(adata_list)
print(adata_combined.shape) #(131709, 22164)

# Load the metadata
metadata_files = [f for f in data_files if f.endswith('.txt.gz')]
metadata_list = []
for metadata_file in metadata_files:
    metadata = pd.read_csv(base_dir + metadata_file, sep='\t', index_col=0)
    metadata_list.append(metadata)
metadata_combined = pd.concat(metadata_list)
print(metadata_combined.shape) #(116955, 3)

# Subset the dataset to only include the cells that have metadata
adata_combined = adata_combined[adata_combined.obs.index.isin(metadata_combined.index)].copy()
adata_combined.obs = adata_combined.obs.merge(metadata_combined, left_index=True, right_index=True)

# Save the data
adata_combined.write_h5ad(out_dir + "MDS_raw.h5ad")
