import scanpy as sc
import pandas as pd
import numpy as np
import os

path = "../data/HCA-BM/20492a4b-0def-457b-9574-60dfdde2a0f2/BM_standard_design.h5ad"
print("Loading adata from path:", path)
adata = sc.read_h5ad(path)
is_b = adata.obs['anno'].isin(['Pro-B cells', 'Pre-B cells', 'Naive B cells'])
b_adata = adata[is_b].copy()
b_adata.X = b_adata.raw.X.copy()
del adata

print("Filtering highly variable genes...")
sc.pp.highly_variable_genes(b_adata, n_top_genes=10000, flavor='seurat_v3', batch_key='Donor')
print("Normalize the data...")
sc.pp.normalize_total(b_adata, target_sum=1e4)
sc.pp.log1p(b_adata)
print("Running PCA...")
sc.tl.pca(b_adata)
print("Computing neighbors...")
sc.pp.neighbors(b_adata)
print("Running UMAP...")
sc.tl.umap(b_adata)
print("Drawing graph...")
sc.tl.draw_graph(b_adata)

print("Saving the data")
outpath = "../data/HCA-BM/results/"
if not os.path.exists(outpath):
    os.makedirs(outpath)
b_adata.write_h5ad(os.path.join(outpath, 'combined.h5ad'))
