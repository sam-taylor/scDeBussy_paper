import os
import pandas as pd
import scanpy as sc
import anndata as ad
import seaborn as sns
from matplotlib.cm import get_cmap
from matplotlib import cycler
from scipy.stats import spearmanr
from scipy.stats import wilcoxon
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec  # Import GridSpec


tab10_colors = get_cmap('tab10').colors  # Extract colors from tab10
color_cycler = cycler(color=tab10_colors)

# Load un-batch corrected data
print('Loading data...')
path = "../data/Nadeu2022_NatMed_scRNAseq_data/results/"
adata = sc.read_h5ad(os.path.join(path, 'combined.h5ad'))
adata = adata[adata.obs['annotation_coarse'] != 'normal',:]
adata_harmony = sc.read_h5ad(os.path.join(path, 'harmony_palantir.h5ad'))
adata_harmony = adata_harmony[adata_harmony.obs['annotation_coarse'] != 'normal']

df_aligned = pd.read_csv(os.path.join(path, 'aligned_pseudotime.csv'))
pseudotime_map = dict(zip(df_aligned['cell_id'], df_aligned['aligned_score']))
adata.obs['aligned_score'] = adata.obs_names.map(lambda x: pseudotime_map[x])

# Spearman correlation of cell type labels
print('Calculating Spearman correlation of cell type labels...')
correlations = []
for i, donor in enumerate(adata.obs['donor_id'].unique()):
    donor_adata = adata[adata.obs['donor_id'] == donor,:]
    corr, _ = spearmanr(donor_adata.obs['aligned_score'], donor_adata.obs['annotation_coarse'].map(lambda x: 0 if x == 'CLL' else 1))
    correlations.append({
        'method': 'CellAlignDTW',
        'donor': donor,
        'correlation': corr
    })
aligned_correlations = pd.DataFrame(correlations)

correlations = []
for i, donor in enumerate(adata.obs['donor_id'].unique()):
    donor_adata = adata_harmony[adata.obs['donor_id'] == donor,:]
    corr, _ = spearmanr(donor_adata.obs['palantir_pseudotime'], donor_adata.obs['annotation_coarse'].map(lambda x: 0 if x == 'CLL' else 1))
    correlations.append({
        'method': 'Harmony',
        'donor': donor,
        'correlation': corr
    })
harmony_correlations = pd.DataFrame(correlations)
correlations = pd.concat([harmony_correlations, aligned_correlations])

print(wilcoxon(correlations[correlations['method'] == 'Harmony']['correlation'],
         correlations[correlations['method'] == 'CellAlignDTW']['correlation'],
         ))

# Create the main figure
fig = plt.figure(figsize=(7.1, 3))  # Adjust figure size
legend_fontsize=8
label_size=9
# Define GridSpec (2 rows, 4 columns)
gs = gridspec.GridSpec(2, 4, width_ratios=[1, 1, 1, 1.5])  # Right panel is slightly wider

# ---- Left Panel (2×3 Grid) ----
colors = ['annotation_coarse', 'donor_id', 'palantir_pseudotime']
color_names = ['Cell type', 'Donor', 'Pseudotime']

left_axes = []
for i in range(2):  # Rows
    for j in range(3):  # Columns
        ax = fig.add_subplot(gs[i, j])  # Assign each grid space
        left_axes.append(ax)

# First dataset (Top row)
for idx, color in enumerate(colors):
    #legend_loc = 'lower right'
    color = 'aligned_score' if color == 'palantir_pseudotime' else color
    sc.pl.draw_graph(
        adata, 
        color=color, 
        ax=left_axes[idx], 
        show=False, 
        frameon=False, 
        alpha=0.6,
        #legend_loc=legend_loc,
        palette=color_cycler, cmap='Spectral_r', vmin='p0.5', vmax='p99.5',
        legend_fontsize=legend_fontsize
    )
    if color == 'aligned_score':
        cb_ax = left_axes[idx].collections[-1].colorbar.ax  # Get colorbar axis
        cb_ax.tick_params(labelsize=legend_fontsize)
    left_axes[idx].set_title(color_names[idx], fontsize=label_size)

# Second dataset (Bottom row)
for idx, color in enumerate(colors):
    sc.pl.draw_graph(
        adata_harmony, 
        color=color, 
        ax=left_axes[idx + 3], 
        show=False, 
        frameon=False, 
        sort_order=True,
        legend_loc=None,
        palette=color_cycler, cmap='Spectral_r', vmin='p1', vmax='p99',
        legend_fontsize=legend_fontsize
    )
    if color == 'palantir_pseudotime':
        cb_ax = left_axes[idx + 3].collections[-1].colorbar.ax  # Get colorbar axis
        cb_ax.tick_params(labelsize=legend_fontsize)
    left_axes[idx + 3].set_title('')  # Remove second row titles

# Add left-side labels
fig.text(0.0, 0.7, 'No Batch Correction', ha='center', va='center', rotation=90, fontsize=label_size)
fig.text(0.0, 0.25, 'Harmony', ha='center', va='center', rotation=90, fontsize=label_size)

# ---- Right Panel (Spearman Correlation Scatter) ----
right_ax = fig.add_subplot(gs[:, 3])  # Span both rows in last column

sns.scatterplot(data=correlations, x='method', y='correlation', hue='donor', palette='tab10', ax=right_ax)
sns.lineplot(data=correlations, x='method', y='correlation', hue='donor', palette='tab10', marker='o', legend=False, ax=right_ax)

right_ax.set_ylim(0, 1)
right_ax.set_ylabel('Spearman correlation', fontsize=label_size)
right_ax.legend(loc='lower right', title='Donor', fontsize=legend_fontsize, title_fontsize=legend_fontsize)
right_ax.set_xlabel('')
right_ax.tick_params(axis='both', which='major', labelsize=legend_fontsize)
# ---- Labels "A" and "B" ----
fig.text(-0.01, 0.95, "A", fontsize=label_size, fontweight="bold")  # Label for left side
fig.text(0.68, 0.95, "B", fontsize=label_size, fontweight="bold")  # Label for right side

plt.tight_layout()
plt.show()

# Save figure
fig.savefig('../figures/figure_4.png', dpi=300)
