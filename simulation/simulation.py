# The simulation is based on Genes2Genes Gaussian Process simulation for pseudotime data (matching case)
# https://github.com/Teichlab/G2G_notebooks/blob/main/Simulation_experiments/DataGeneration_TrajectorySimulatorUsingGPs.ipynb
import torch
import os
import anndata
import gpytorch
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
import warnings
import pickle
from tqdm import tqdm
import pandas as pd
import scipy

warnings.simplefilter(action='ignore', category=FutureWarning)
torch.cuda.is_available() 
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

def plot_simulated_dataset(simulated_dataset, T_f, figsize = (15,3)):
    l = len(simulated_dataset)
    ncol = 4
    nrow = int(np.ceil(l/ncol)) 
    plt.subplots(nrow,ncol,figsize=figsize) 
    k=1
    for pair in simulated_dataset:
        plt.subplot(nrow,ncol,k)
        plt.scatter(T_f,pair[0], color='midnightblue')
        plt.scatter(T_f,pair[1], color='forestgreen')
        plt.scatter(T_f,pair[2], color='orange')
        k=k+1


def generate_matched_functions(f, σ = 0.1):
    sample1 = gpytorch.distributions.MultivariateNormal(f, (σ**2)*torch.eye(len(f))).rsample()
    sample2 = gpytorch.distributions.MultivariateNormal(f, (σ**2)*torch.eye(len(f))).rsample()
    sample3 = gpytorch.distributions.MultivariateNormal(f, (σ**2)*torch.eye(len(f))).rsample()
    return sample1,sample2,sample3


def run_all_match_simulator(T_f, n_simulations, base_mean_factor=0.0):
    matched_datasets = []
    for i in tqdm(range(n_simulations)):
        σ = np.random.uniform(0.05, 1.0)  # Random variation
        base_mean_factor = np.random.uniform(0.5, 9.0)
        base_mean = torch.zeros(len(T_f)) + base_mean_factor
        # Define the kernel and compute covariance matrix
        base_kernel = gpytorch.kernels.RBFKernel()
        cov_matrix = base_kernel(T_f).evaluate()
        # Add a small jitter term to the diagonal to ensure positive definiteness
        jitter = 1e-4 * torch.eye(cov_matrix.shape[0], device=cov_matrix.device)
        cov_matrix += jitter
        # Sample from the multivariate normal distribution
        f = gpytorch.distributions.MultivariateNormal(base_mean, cov_matrix).rsample()
        # Generate matched functions
        sample1, sample2, sample3 = generate_matched_functions(f, σ)
        dataset = [sample1.detach().numpy(), sample2.detach().numpy(), sample3.detach().numpy()]
        matched_datasets.append(dataset)
    return matched_datasets


def convert2adata(simulated_datasets,T_f, write=False):
    # if const_std mode, G2G uses a constant of std=0.1 for all time points 
    n_sim_genes = len(simulated_datasets) 
    sample_1_df = []
    sample_2_df = []
    sample_3_df = []
    time = T_f.detach().numpy()
    
    for dataset in simulated_datasets:
        sample_1_df.append( dataset[0])
        sample_2_df.append( dataset[1])
        sample_3_df.append( dataset[2])

    sample_1_df = pd.DataFrame(sample_1_df).transpose() 
    sample_2_df = pd.DataFrame(sample_2_df).transpose() 
    sample_3_df = pd.DataFrame(sample_3_df).transpose() 
    sample_1_df.columns = np.asarray(['Gene' + str(x) for x in np.arange(n_sim_genes)]) 
    sample_2_df.columns = sample_1_df.columns
    sample_3_df.columns = sample_1_df.columns
    gene_list = list(sample_1_df.columns)

    adata_sample_1 = anndata.AnnData(X=scipy.sparse.csr_matrix((sample_1_df)))
    adata_sample_1.var_names = sample_1_df.columns
    adata_sample_1.obs_names = sample_1_df.index
    adata_sample_1.obs['time'] = time
    adata_sample_2 = anndata.AnnData(X=scipy.sparse.csr_matrix(sample_2_df))
    adata_sample_2.var_names = sample_2_df.columns
    adata_sample_2.obs_names = sample_2_df.index
    adata_sample_2.obs['time'] = time
    adata_sample_3 = anndata.AnnData(X=scipy.sparse.csr_matrix(sample_3_df))
    adata_sample_3.var_names = sample_3_df.columns
    adata_sample_3.obs_names = sample_3_df.index
    adata_sample_3.obs['time'] = time
    
    if(write):
        dataset_dir = 'data/simulation/'
        adata_sample_1.write_h5ad(  dataset_dir   + 'sample_1' +'.h5ad')
        adata_sample_2.write_h5ad(dataset_dir   + 'sample_2' +'.h5ad')
        adata_sample_3.write_h5ad(dataset_dir   + 'sample_3' +'.h5ad')

    return adata_sample_1, adata_sample_2, adata_sample_3


torch.manual_seed(1) 
np.random.seed(1)
T_f = torch.tensor(np.linspace(0,1,300,dtype=np.float64))
matched_datasets = run_all_match_simulator(T_f, n_simulations=500) # genes x sample x time
plot_simulated_dataset(matched_datasets[10:26], T_f, figsize=(15,10))
if not os.path.exists('data/simulation'):
    os.makedirs('data/simulation')
pickle.dump(matched_datasets, open('data/simulation/matched_datasets.pkl', 'wb'))

adata_sample_1, adata_sample_2, adata_sample_3 = convert2adata(matched_datasets, T_f, write=False)
