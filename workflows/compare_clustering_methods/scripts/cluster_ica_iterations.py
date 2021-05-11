#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Load ICA iterations, cluster them with every algorithm.

import argparse
import pandas as pd
import os
import numpy as np
from robustica import RobustICA
from joblib import parallel_backend
from sklearn.metrics import silhouette_samples
import time

# variables
SAVE_PARAMS = {"sep": "\t", "index": False}


"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_robustica'
PREP_DIR = os.path.join(ROOT,'data','prep')
RESULTS_DIR = os.path.join(ROOT,'results','compare_clustering_methods')
input_file = os.path.join(PREP_DIR,'TCGA','genexpr','LGG.tsv.gz')
S_all_file = os.path.join(RESULTS_DIR,'files','ica_iterations-LGG-S.pickle')
A_all_file = os.path.join(RESULTS_DIR,'files','ica_iterations-LGG-A.pickle')
method = 'DBSCAN'
iterations = 100
n_jobs = 15
"""


##### FUNCTIONS #####
def process_inputs(input_file):
    # load
    X = pd.read_table(input_file, index_col=0)

    # drop genes with no variation
    X = X.loc[X.std(1) > 0]

    # normalize
    X = (X - X.mean(1).values.reshape(-1, 1)) / X.std(1).values.reshape(-1, 1)
    return X


def load_data(input_file, S_file, A_file):
    X = process_inputs(input_file)
    S_all = pd.read_pickle(S_file)
    A_all = pd.read_pickle(A_file)

    return X, S_all, A_all


def evaluate_clustering(S_all, A_all, labels, signs):
    """
    - std within cluster by row
    - silhouette
    """
    print("Evaluating the clustering...")

    # consider signs of components
    S = np.multiply(S_all, signs)
    A = np.multiply(A_all, signs)

    # silhouette of components
    silhouette = silhouette_samples(S.T.values, labels)
    silhouette = pd.Series(silhouette, index=S_all.columns)

    # standard deviation of features within cluster
    feature_stats = {"S": {}, "A": {}}
    for cluster in np.unique(labels):
        idx = labels == cluster
        feature_stats["S"][cluster] = S.iloc[:, idx].std(1)
        feature_stats["A"][cluster] = A.iloc[:, idx].std(1)
    S_std = pd.concat(feature_stats["S"], axis=1)
    A_std = pd.concat(feature_stats["A"], axis=1)

    return silhouette, S_std, A_std


def cluster_ica_iterations(X, S_all, A_all, iterations, method, n_jobs):
    # prepare
    n_components = int(S_all.shape[1] / iterations)

    # initialize class
    rica = RobustICA(
        n_components=n_components,
        robust_method=method,
        robust_iter=iterations,
        robust_dimreduce=True,
    )
    rica.S_all = S_all.values
    rica.A_all = A_all.values
    
    start_time = time.time()
    with parallel_backend("loky", n_jobs=n_jobs):
        rica._cluster_components(X.values)
    seconds = time.time() - start_time
    labels = rica.clustering.labels_
    signs = rica.iteration_signs_

    silhouette, S_std, A_std = evaluate_clustering(S_all, A_all, labels, signs)

    # prepare output
    clustering_summary = pd.DataFrame(
        {
            "component": S_all.columns,
            "label": labels,
            "sign": signs,
            "silhouette": silhouette,
            "time": seconds,
            "clustering_method": method
        }
    )
    S_mean = pd.DataFrame(rica.S, index=S_all.index, columns=np.unique(labels))
    A_mean = pd.DataFrame(rica.A, index=A_all.index, columns=np.unique(labels))

    return clustering_summary, S_mean, A_mean, S_std, A_std


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str)
    parser.add_argument("--S_all_file", type=str)
    parser.add_argument("--A_all_file", type=str)
    parser.add_argument("--method", type=str)
    parser.add_argument("--iterations", type=int)
    parser.add_argument("--output_dir", type=str)
    parser.add_argument("--n_jobs", type=int)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    input_file = args.input_file
    S_all_file = args.S_all_file
    A_all_file = args.A_all_file
    method = args.method
    iterations = args.iterations
    output_dir = args.output_dir
    n_jobs = args.n_jobs

    os.makedirs(output_dir)

    X, S_all, A_all = load_data(input_file, S_all_file, A_all_file)
    clustering_summary, S_mean, A_mean, S_std, A_std = cluster_ica_iterations(
        X, S_all, A_all, iterations, method, n_jobs
    )

    # save
    clustering_summary.to_csv(
        os.path.join(output_dir, "clustering_summary.tsv"), **SAVE_PARAMS
    )
    S_mean.reset_index().to_csv(
        os.path.join(output_dir, "S_mean.tsv.gz"), compression="gzip", **SAVE_PARAMS
    )
    A_mean.reset_index().to_csv(
        os.path.join(output_dir, "A_mean.tsv.gz"), compression="gzip", **SAVE_PARAMS
    )
    S_std.reset_index().to_csv(
        os.path.join(output_dir, "S_std.tsv.gz"), compression="gzip", **SAVE_PARAMS
    )
    A_std.reset_index().to_csv(
        os.path.join(output_dir, "A_std.tsv.gz"), compression="gzip", **SAVE_PARAMS
    )


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
