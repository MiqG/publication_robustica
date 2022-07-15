#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Evaluate how robust is the computation of robust independent components
# applying random noise based on the centroid's mean and standard deviation.

import argparse
import os
import pandas as pd
import numpy as np
from sklearn.metrics import pairwise_distances
from tqdm import tqdm
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
fdrtool = importr("fdrtool")

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}

"""
Development
-----------
ROOT = '/home/miquel/projects/publication_robustica'
RESULTS_DIR = os.path.join(ROOT,'results','benchmark_sign_inference')
algorithm = 'robustica_pca'
S_mean_file = os.path.join(RESULTS_DIR,'files','Sastry2019','robustica_pca-S.tsv.gz')
S_std_file = os.path.join(RESULTS_DIR,'files','Sastry2019','robustica_pca-S_std.tsv.gz')
n_samples=10
index_col=0 #'log-TPM'
"""


##### FUNCTIONS #####
def load_data(S_mean_file, S_std_file):
    """
    We consider the first iteration as our ground truth.
    """
    # ICA runs
    S_mean = pd.read_table(S_mean_file)
    S_std = pd.read_table(S_std_file)

    return S_mean, S_std


def define_modules(x, cutoff=0.01):
    fdr = fdrtool.fdrtool(
        ro.FloatVector(x), cutoff_method="fndr", plot=False, verbose=False
    )
    fdr = np.array(fdr.rx2("qval"))
    return fdr < cutoff


def get_modules(X):
    modules = np.apply_along_axis(define_modules, axis=0, arr=X)
    return modules


def compare_modules(A, B):
    jaccard = 1 - pairwise_distances(A.T, B.T, metric="jaccard")
    return jaccard.max(axis=1)


def map_sampled_modules(S_mean, S_std, index_col, modules_ref):
    # prep inputs
    mat_mean = S_mean.drop(columns=S_mean.columns[index_col]).values
    mat_std = S_std.drop(columns=S_mean.columns[index_col]).values
    
    # sample each weight
    mat_sampled = np.full(mat_mean.shape, np.nan)
    for i in range(mat_mean.shape[0]):
        for j in range(mat_mean.shape[1]):
            mu, sigma = mat_mean[i,j], mat_std[i,j]
            mat_sampled[i,j] = np.random.normal(mu, sigma, 1)
    
    # obtain modules
    modules_sampled = get_modules(mat_sampled)

    # compare modules
    jaccard = compare_modules(modules_sampled, modules_ref)
    
    result = pd.Series({
        'mean': np.mean(jaccard),
        'median': np.median(jaccard),
        'std': np.std(jaccard),
        'q25': np.quantile(jaccard, 0.25),
        'q75': np.quantile(jaccard, 0.75),
        'min': np.min(jaccard),
        'max': np.max(jaccard)
    })
    
    return result


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--S_mean_file", type=str)
    parser.add_argument("--S_std_file", type=str)
    parser.add_argument("--output_file", type=str)
    parser.add_argument("--algorithm", type=str)
    parser.add_argument("--n_samples", type=int)
    parser.add_argument("--index_col", type=int, default=0)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    S_mean_file = args.S_mean_file
    S_std_file = args.S_std_file
    output_file = args.output_file
    algorithm = args.algorithm
    n_samples = args.n_samples
    index_col = args.index_col

    # load data
    S_mean, S_std = load_data(S_mean_file, S_std_file)
    
    # sample modules jaccards
    modules_ref = get_modules(S_mean.drop(columns=S_mean.columns[index_col]))
    results = [
        map_sampled_modules(
            S_mean, S_std, index_col, modules_ref
        )
        for it in tqdm(range(n_samples))
    ]
    results = pd.DataFrame(results)
    results['algorithm'] = algorithm
    results['n_samples'] = n_samples

    # save outputs
    results.to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
