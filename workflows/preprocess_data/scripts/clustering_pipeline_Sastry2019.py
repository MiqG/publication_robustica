#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Reproduce clustering step according to Sastry (2019).

import argparse
import numpy as np
import pandas as pd
import sys
from sklearn.cluster import DBSCAN
from scipy import stats
import time
import os

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_robustica'
PREP_DIR = os.path.join(ROOT,'data','prep')
S_all_file = os.path.join(PREP_DIR,'ica_runs','Sastry2019','S.pickle.gz')
A_all_file = os.path.join(PREP_DIR,'ica_runs','Sastry2019','A.pickle.gz')
distance_func = "original"
eps = 0.5
min_frac = 0.5
n_iters = 100
"""

##### FUNCTIONS ####
def timeit(start):
    end = time.time()
    t = end - start
    if t < 60:
        print("{:.2f} seconds elapsed".format(t))
    elif t < 3600:
        print("{:.2f} minutes elapsed".format(t / 60))
    else:
        print("{:.2f} hours elapsed".format(t / 3600))
    return end


def load_data(
    S_all_file, A_all_file,
):
    # ICA runs
    S_all = pd.read_pickle(S_all_file)
    A_all = pd.read_pickle(A_all_file)

    return S_all, A_all


def compute_distances(S_all, distance_func):
    """
    Compute distances as in Sastry 2019 or as we do.
    """
    if distance_func == "original":
        D = 1 - abs(np.dot(S_all.T, S_all))
        D[D < 0.5] = 0
    elif distance_func == "full":
        D = 1 - abs(np.dot(S_all.T, S_all))
    elif distance_func == "robustica":
        D = 1 - abs(np.corrcoef(S_all.T))
    else:
        print("Error: Set a valid 'distance_func'.")
        sys.exit(1)

    D = np.clip(D, 0, 1)

    return D


def cluster_components(D, eps, min_frac, n_iters, S_all, A_all):

    # Run DBSCAN
    min_samples = int(round(min_frac * n_iters)) + 1
    dbscan = DBSCAN(eps=eps, min_samples=min_samples, metric="precomputed")
    labels = dbscan.fit_predict(D)
    n_clusters = max(labels) + 1
    print("Number of clusters:", n_clusters)

    # Gather final S and A matrices
    print("\nGathering final S and A matrices")

    S_final = pd.DataFrame(columns=range(n_clusters), index=S_all.index)
    A_final = pd.DataFrame(columns=range(n_clusters), index=A_all.index)
    df_stats = pd.DataFrame(
        columns=["S_mean_std", "A_mean_std", "count"], index=range(n_clusters)
    )

    t = time.time()

    for lab in range(n_clusters):

        idx = labels == lab

        S_clust = S_all.loc[:, idx].values.T
        A_clust = A_all.loc[:, idx].values.T

        # First item is base component
        Svec0 = S_clust[0]
        Avec0 = A_clust[0]

        # Make sure base component is facing positive
        if abs(min(Svec0)) > max(Svec0):
            Svec0 = -Svec0
            Avec0 = -Avec0

        S_single = [Svec0]
        A_single = [Avec0]

        # Add in rest of components
        for j in range(1, len(S_clust)):
            Svec = S_clust[j]
            Avec = A_clust[j]
            if stats.pearsonr(Svec, Svec0)[0] > 0:
                S_single.append(Svec)
                A_single.append(Avec)
            else:
                S_single.append(-Svec)
                A_single.append(-Avec)

        # Add centroid of cluster to final S matrix
        S_final[lab] = np.array(S_single).T.mean(axis=1)
        A_final[lab] = np.array(A_single).T.mean(axis=1)

        # Get component stats
        df_stats.loc[lab, "S_mean_std"] = np.array(S_single).std(axis=1).mean()
        df_stats.loc[lab, "A_mean_std"] = np.array(A_single).std(axis=1).mean()
        df_stats.loc[lab, "count"] = len(S_single)

    print("\nFinal components created")
    t = timeit(t)

    return S_final, A_final, df_stats


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--S_all_file", type=str)
    parser.add_argument("--A_all_file", type=str)
    parser.add_argument("--distance_func", type=str)
    parser.add_argument("--eps", type=float)
    parser.add_argument("--min_frac", type=float)
    parser.add_argument("--n_iters", type=int)
    parser.add_argument("--output_dir", type=str)
    args = parser.parse_args()

    return args


def main():
    # unpack arguments
    args = parse_args()
    S_all_file = args.S_all_file
    A_all_file = args.A_all_file
    distance_func = args.distance_func
    eps = args.eps
    min_frac = args.min_frac
    n_iters = args.n_iters
    output_dir = args.output_dir

    # load data
    print("Loading data..")
    S_all, A_all = load_data(S_all_file, A_all_file)

    # cluster
    print("Clustering ICs...")
    D = compute_distances(S_all, distance_func)
    S, A, df_stats = cluster_components(D, eps, min_frac, n_iters, S_all, A_all)

    # save
    print("Saving...")
    os.makedirs(output_dir, exist_ok=True)
    S.to_csv(os.path.join(output_dir, "S.csv"))
    A.to_csv(os.path.join(output_dir, "A.csv"))
    df_stats.to_csv(os.path.join(output_dir, "component_stats.csv"))


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
