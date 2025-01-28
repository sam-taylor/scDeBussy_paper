import os
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import palantir
import scvi
import torch


def determine_terminal_cells(adata_patient, cell_type, dm_mat, annotation_column):
    init_pdt = adata_patient.obs[annotation_column].map(lambda x: 1 if x == cell_type else 0)
    rho = [pearsonr(dm_mat.iloc[:,i], init_pdt)[0] for i in np.arange(dm_mat.shape[1])]
    dc_select = np.argmax(np.abs(rho))
    adata_patient.uns['DC_select'] = dc_select 
    start_cell_idx = (np.sign(rho[dc_select]) * dm_mat.iloc[:, dc_select]).argmax()
    start_cell_bc = adata_patient.obs.index[start_cell_idx]
    return start_cell_bc

def run_palantir_analysis(adata, start_cell_type, end_cell_type, cell_type_col, save_path):
    print("Preparing DM matrix...")
    dm_key = 'DM_EigenVectors'
    dm_mat = pd.DataFrame(adata.obsm[dm_key], index=adata.obs.index)
    
    print("Determining start cell...")
    start_cell_bc = determine_terminal_cells(adata, start_cell_type, dm_mat, cell_type_col)
    
    print("Determining end cell...")
    end_cell_bc = determine_terminal_cells(adata, end_cell_type, dm_mat, cell_type_col)
    
    print("Running Palantir...")
    terminal_states = pd.Series([end_cell_type], index=[end_cell_bc])
    pr_res = palantir.core.run_palantir(
        dm_mat, start_cell_bc, terminal_states=terminal_states, use_early_cell_as_start=True
    )
    adata.obs['palantir_pseudotime'] = pr_res.pseudotime
    adata.write_h5ad(save_path)

def run_cell_aligndtw(path, batch_key, start_cell_type, end_cell_type, cell_type_col):
    # Load data
    adata = sc.read_h5ad(path)
    adata.X = adata.layers['counts'].copy()
    patients = adata.obs.loc[:,batch_key].unique()
    n_components = [10, 10, 10, 10, 10]

    for i, patient in enumerate(patients):
        print(patient)
        adata_patient = adata[adata.obs.loc[:,batch_key] == patient].copy()
        print("Filtering highly variable genes...")
        sc.pp.highly_variable_genes(adata_patient, n_top_genes=10000, flavor='seurat_v3')
        sc.pp.normalize_total(adata_patient, target_sum=1e4)
        sc.pp.log1p(adata_patient)
        print("Running PCA...")
        sc.pp.pca(adata_patient)
        print("Computing neighbors...")
        sc.pp.neighbors(adata_patient)
        print("Running UMAP...")
        sc.tl.umap(adata_patient)
        print("Running PAGA...")
        sc.tl.paga(adata_patient, groups=cell_type_col)
        sc.pl.paga(adata_patient, plot=True)
        print("Drawing graph...")
        sc.tl.draw_graph(adata_patient, init_pos='paga')
        print("Running Magic Imputation...")
        dm_res = palantir.utils.run_diffusion_maps(adata_patient, n_components=n_components[i])
        imputed_X = palantir.utils.run_magic_imputation(adata_patient)
        print("Saving tumor cells adata for ", patient)
        adata_patient.write_h5ad(os.path.join(os.path.dirname(path), f'{patient}.h5ad'))
        outpath = os.path.join(os.path.dirname(path), f'{patient}_palantir.h5ad')
        adata_patient = run_palantir_analysis(adata_patient, start_cell_type, end_cell_type, cell_type_col, outpath)


def run_harmony(path, batch_key, start_cell_type, end_cell_type, cell_type_col):
    adata = sc.read_h5ad(path)
    sc.external.pp.harmony_integrate(adata, batch_key, basis="X_pca")
    sc.pp.neighbors(adata, use_rep="X_pca_harmony")
    print("Running UMAP...")
    sc.tl.umap(adata)
    print("Running PAGA...")
    sc.tl.paga(adata, groups=cell_type_col)
    sc.pl.paga(adata, plot=True)
    print("Drawing graph with Harmony...")
    sc.tl.draw_graph(adata, use_rep='X_pca_harmony', init_pos='paga')
    print("Running Magic Imputation...")
    dm_res = palantir.utils.run_diffusion_maps(adata, n_components=10, pca_key='X_pca_harmony')
    imputed_X = palantir.utils.run_magic_imputation(adata)
    print("Saving tumor cells adata to path:", path)
    adata.write_h5ad(os.path.join(os.path.dirname(path), 'harmony_combined.h5ad'))
    outpath = os.path.join(os.path.dirname(path), f'harmony_palantir.h5ad')
    adata = run_palantir_analysis(adata, start_cell_type, end_cell_type, cell_type_col, outpath)


def run_scvi(path, batch_key, start_cell_type, end_cell_type, cell_type_col):
    adata = sc.read_h5ad(path)
    torch.set_float32_matmul_precision("high")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    scvi.model.SCVI.setup_anndata(adata_hvg, batch_key=batch_key, layer="counts")
    model = scvi.model.SCVI(adata_hvg, n_latent=50)
    model.train(early_stopping=True, early_stopping_monitor="reconstruction_loss_validation")
    latent = model.get_latent_representation()
    denoised = model.get_normalized_expression(library_size=1e4)
    integrated = model.get_normalized_expression(
    transform_batch=adata.obs[:,batch_key].unique().tolist(), library_size=1e4
    )
    adata.obsm["X_scVI"] = latent
    adata.obsm["scvi_normalized"] = denoised
    adata.obsm["scvi_normalized_integrated"] = integrated
    
    print("Computing neighbors with scVI...")
    sc.pp.neighbors(adata, use_rep="X_scVI")
    print("Running UMAP...")
    sc.tl.umap(adata)
    print("Running PAGA...")
    sc.tl.paga(adata, groups=cell_type_col)
    sc.pl.paga(adata, plot=True)
    print("Drawing graph with Harmony...")
    sc.tl.draw_graph(adata, use_rep='X_scVI', init_pos='paga')
    print("Running Magic Imputation...")
    dm_res = palantir.utils.run_diffusion_maps(adata, n_components=10, pca_key='X_scVI')
    imputed_X = palantir.utils.run_magic_imputation(adata)
    print("Saving results...")
    adata.write_h5ad(os.path.join(os.path.dirname(path), 'scvi_combined.h5ad'))
    outpath = os.path.join(os.path.dirname(path), f'scvi_palantir.h5ad')
    adata = run_palantir_analysis(adata, start_cell_type, end_cell_type, cell_type_col, outpath)


def main():
    parser = argparse.ArgumentParser(description='Preprocess data with different correction methods.')
    parser.add_argument('path', type=str, help='Path to the data directory.')
    parser.add_argument('batch_correction', choices=['no', 'harmony', 'scvi'], 
                        help='Correction method to use: no, harmony, or scvi.')
    parser.add_argument('batch_key', type=str, 
                        help='Whether to separate process individual samples. Only available when batch_correction is no.')
    parser.add_argument('start_cell_type', type=str, 
                        help='Cell type to start the trajectory inference.')
    parser.add_argument('end_cell_type', type=str, 
                        help='Cell type to end the trajectory inference.')
    parser.add_argument('cell_type_col', type=str, 
                        help='Column in the adata.obs dataframe that contains the cell type information.')
    args = parser.parse_args()

    if args.batch_correction == 'no':
        run_cell_aligndtw(args.path, args.batch_key, args.start_cell_type, args.end_cell_type, args.cell_type_col)
    elif args.batch_correction == 'harmony':
        run_harmony(args.path, args.batch_key, args.start_cell_type, args.end_cell_type, args.cell_type_col)
    elif args.batch_correction == 'scvi':
        run_scvi(args.path, args.batch_key, args.start_cell_type, args.end_cell_type, args.cell_type_col)

if __name__ == "__main__":
    main()
