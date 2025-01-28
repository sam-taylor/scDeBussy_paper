import scanpy as sc
import pandas as pd
import os
import numpy as np
from CellAlignDTW.pp import create_cellrank_probability_df
from CellAlignDTW import CellAlignDTW, gam_smooth_expression
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

path = "../data/Nadeu2022_NatMed_scRNAseq_data/results/"
patients = ['63', '12', '19', '365', '3299']
cluster_ordering = ['CLL', 'RT']
df_combined = pd.DataFrame()

h5ad_files = [os.path.join(path, f"{patient}_palantir.h5ad") for patient in patients]
logger.info("Starting to create cellrank probability dataframe...")
combined_adata, df = create_cellrank_probability_df(h5ad_files, 'annotation_coarse', patients, cluster_ordering, 'palantir_pseudotime')
logger.info("Cellrank probability dataframe created successfully.")

aligner = CellAlignDTW(df=df, 
                       cluster_ordering=cluster_ordering,
                       subject_col='subject',
                       score_col='score',
                       cell_id_col='cell_id',
                       cell_type_col='cell_type',
                       verbose=True)
logger.info("Starting alignment process...")
aligner.align()
logger.info("Alignment process completed.")

combined_adata.write_h5ad(os.path.join(path, 'combined_palantir.h5ad'))
logger.info("Combined adata written to file.")

aligner.df.to_csv(os.path.join(path, 'aligned_pseudotime.csv'), index=False)
logger.info("Aligned pseudotime dataframe written to file.")

logger.info("Starting GAM smoothing expression...")
df = aligner.df
combined_adata = sc.read_h5ad(os.path.join(path, 'combined_palantir.h5ad'))
genes = combined_adata.var_names.values.tolist()
summary_df, gene_curves, scores_df = gam_smooth_expression(df, genes, n_splines=6, lam=4)
logger.info("GAM smoothing expression completed.")

summary_df.to_csv(os.path.join(path, 'gam_summary.csv'), index=False)
gene_curves.to_csv(os.path.join(path, 'gam_gene_curves.csv'), index=False)
scores_df.to_csv(os.path.join(path, 'gam_scores.csv'), index=False)
