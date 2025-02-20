import torch
import os
import gpytorch
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
import warnings
import pickle
from tqdm import tqdm

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


torch.manual_seed(1) 
np.random.seed(1)
T_f = torch.tensor(np.linspace(0,1,300,dtype=np.float64))
matched_datasets = run_all_match_simulator(T_f, n_simulations=500) # genes x sample x time
plot_simulated_dataset(matched_datasets[10:26], T_f, figsize=(15,10))
if not os.path.exists('data/simulation'):
    os.makedirs('data/simulation')
pickle.dump(matched_datasets, open('data/simulation/matched_datasets.pkl', 'wb'))
