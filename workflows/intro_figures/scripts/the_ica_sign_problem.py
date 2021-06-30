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
N_FEATURES = 15
N_SAMPLES = 20
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
    
    return S, A, S_extreme

 
def make_plots(X, output_dir):
    # overview of S
    f = plt.figure(figsize=(2,4))
    sns.heatmap(X, cmap='coolwarm', linewidths=.1, linecolor='black', annot=True, cbar=False, center=0, annot_kws={"fontsize":7})
    f.savefig(os.path.join(output_dir,"S.pdf"))
    
    # pearson correlation
    f = plt.figure(figsize=(3,3))
    sns.heatmap(np.round(np.corrcoef(X.T),2), cmap='coolwarm', linewidths=.1, linecolor='black', annot=True, cbar=False, center=0, annot_kws={"fontsize":7})
    f.savefig(os.path.join(output_dir,"pearson.pdf"))
    
    # absolute pearson correlation
    f = plt.figure(figsize=(3,3))
    sns.heatmap(np.abs(np.round(np.corrcoef(X.T),2)), cmap='coolwarm', linewidths=.1, linecolor='black', annot=True, cbar=False, center=0, annot_kws={"fontsize":7})
    f.savefig(os.path.join(output_dir,"abs_pearson.pdf"))
        
    
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
    S, A, S_extreme = make_ica_output(data, n_components)
    make_plots(S_extreme, output_dir)
    
    
    
##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")