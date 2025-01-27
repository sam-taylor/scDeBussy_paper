import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import os

# load the tumor and normal data
print("Loading the tumor and normal data")
path = "../data/Nadeu2022_NatMed_scRNAseq_data/filtered"
adata = sc.read_10x_mtx(path)
cell_metadata = pd.read_table(os.path.join(path, "metadata.tsv"))
adata.obs = cell_metadata
# drop problematic columns
adata.layers['counts'] = adata.X.copy()
adata.obs.donor_id = adata.obs.donor_id.astype('category')
#adata.obs.is_richter = adata.obs.annotation_final.apply(lambda x: "RT" in x)


# load the annotated tumor data
print("Loading the annotated tumor data")
path = "../data/Nadeu2022_NatMed_scRNAseq_data/combined"
adata_annotated = sc.read_10x_mtx(path)
cell_metadata = pd.read_table(os.path.join(path, "metadata.tsv"))
adata_annotated.obs = cell_metadata
adata_annotated.obs.is_richter = adata_annotated.obs.annotation_final.apply(lambda x: "RT" in x)
annotation_mapping = dict(zip(adata_annotated.obs_names, adata_annotated.obs.annotation_final))
adata.obs.loc[:,'annotation'] = adata.obs_names.map(annotation_mapping)
adata.obs.loc[:,'annotation_coarse'] = adata.obs.loc[:,'annotation'].apply(lambda x: "normal" if x is np.nan else "RT" if "RT" in x else "CLL")
adata.var["mito"] = adata.var_names.str.upper().str.startswith(("MT-"))
adata.obs["mito_frac"] = np.sum(adata[:, adata.var["mito"]].X, axis=1) / np.sum(
    adata.X, axis=1
)

print("Filtering highly variable genes...")
sc.pp.highly_variable_genes(adata, n_top_genes=10000, batch_key='donor_id', flavor='seurat_v3')
print("Normalize the data...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print("Running PCA...")
sc.tl.pca(adata)
print("Computing neighbors...")
sc.pp.neighbors(adata)
print("Running UMAP...")
sc.tl.umap(adata)
print("Drawing graph...")
sc.tl.draw_graph(adata)
adata.obs = adata.obs.drop(columns = ["keep_cells", "scrublet_predicted_doublet", "is_hashed"])

print("Saving the annotated tumor and normal data")
outpath = "../data/Nadeu2022_NatMed_scRNAseq_data/results/"
if not os.path.exists(outpath):
    os.makedirs(outpath)
adata.write_h5ad(os.path.join(outpath, 'combined.h5ad'))
