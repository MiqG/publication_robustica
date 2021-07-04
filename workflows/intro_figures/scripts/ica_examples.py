#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Evaluate differenti clustering metrics on AgglomerativeClustering and different
# clustering methods using absolute pearson correlation dissimilarity.

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import FastICA

# variables
N_FEATURES = 20
N_SAMPLES = 15
N_COMPONENTS = 8

"""
Development
-----------
n_features = N_FEATURES
n_samples = N_SAMPLES
n_components = N_COMPONENTS
output_dir='.'
"""

##### FUNCTIONS #####
def make_data(n_features, n_samples):
    return np.random.randn(n_features, n_samples)
    
    
def make_ica_output(data, n_components):
    ica = FastICA(n_components)
    S = ica.fit_transform(data)
    A = ica.mixing_
    
    # make S extreme
    idx = np.abs(S) > S.std(axis=0) * 1
    S_extreme = np.zeros(S.shape)
    S_extreme[idx] = np.sign(S[idx])
    
    # make A extreme
    idx = np.abs(A) > A.std(axis=0) * 1
    A_extreme = np.zeros(A.shape)
    A_extreme[idx] = np.sign(A[idx])
    
    # data
    data_extreme = S.dot(A.T)
    
    return S, A, data_extreme, S_extreme, A_extreme

 
def make_plots(data, S, A, output_dir):
    # overview data
    f = plt.figure(figsize=(4,7))
    sns.heatmap(np.round(data,2), cmap='coolwarm', linewidths=.1, linecolor='black', cbar=False, square=True, center=0, annot_kws={"fontsize":7})
    f.savefig(os.path.join(output_dir,"data.pdf"))
    
    # overview of A
    f = plt.figure(figsize=(4,2.5))
    sns.heatmap(A.T, cmap='coolwarm', linewidths=.1, linecolor='black', square=True, cbar=False, center=0, annot_kws={"fontsize":7})
    f.savefig(os.path.join(output_dir,"A.pdf"))
    
    # overview of S
    f = plt.figure(figsize=(2.5,7))
    sns.heatmap(S, cmap='coolwarm', linewidths=.1, linecolor='black', square=True, cbar=False, center=0, annot_kws={"fontsize":7})
    f.savefig(os.path.join(output_dir,"S.pdf"))
    
    # pearson correlation
    f = plt.figure(figsize=(5,5))
    sns.heatmap(np.round(np.corrcoef(S.T),2), cmap='coolwarm', linewidths=.1, linecolor='black', square=True, cbar=True, center=0, vmin=-1, vmax=1, annot_kws={"fontsize":7})
    f.savefig(os.path.join(output_dir,"S_pearson.pdf"))
    
    # absolute pearson correlation
    f = plt.figure(figsize=(5,5))
    sns.heatmap(np.abs(np.round(np.corrcoef(S.T),2)), cmap='coolwarm', linewidths=.1, linecolor='black', square=True, cbar=False, center=0, vmin=-1, vmax=1, annot_kws={"fontsize":7})
    f.savefig(os.path.join(output_dir,"S_abs_pearson.pdf"))
        
    
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_samples", type=int, default=N_SAMPLES)
    parser.add_argument("--n_features", type=int, default=N_FEATURES)
    parser.add_argument("--n_components", type=int, default=N_COMPONENTS)
    parser.add_argument("--output_dir", type=str)
    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    n_features = args.n_features
    n_samples = args.n_samples
    n_components = args.n_components
    output_dir = args.output_dir
    
    os.makedirs(output_dir)
    data = make_data(n_features, n_samples)
    
    S_extremes = []
    A_extremes = []
    for i in range(2):
        S, A, data_extreme, S_extreme, A_extreme = make_ica_output(data, n_components)
        S_extremes.append(S_extreme)
        A_extremes.append(A_extreme)
        
    S_extreme = np.concatenate(S_extremes, axis=1)
    A_extreme = np.concatenate(A_extremes, axis=1)
    make_plots(data_extreme, S_extreme, A_extreme, output_dir)
    
    
    
##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")