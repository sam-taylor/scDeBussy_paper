import os
import scanpy as sc

# Load the data
base_dir = "../data/MDS/results/"
out_dir = base_dir
adata = sc.read_h5ad(base_dir + "MDS_raw.h5ad")
adata = adata[adata.obs['CellType'] != 'not assigned'].copy()
adata.layers['counts'] = adata.X.copy()

print("Filtering highly variable genes...")
n_top_genes = 2000
sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key='sample', flavor='seurat_v3')
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


print("Saving the processed data")
adata.write_h5ad(os.path.join(out_dir, 'MDS_combined.h5ad'))
