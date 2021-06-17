#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Run different algorithms to perform robust ICA.

import argparse
import pandas as pd
import numpy as np
from sklearn.decomposition import FastICA, PCA
from sklearn.cluster import AgglomerativeClustering
from scipy.stats import pearsonr
from memory_profiler import memory_usage
import time

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}
MEM_KWS = {
    "interval": 0.1,
    "timestamps": True,
    "include_children": True,
    "retval": True,
}

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_robustica'
RESULTS_DIR = os.path.join(ROOT,'results','benchmark_performance')
input_file = os.path.join(RESULTS_DIR,'files','sample_data.tsv.gz')
iterations = 10
algorithm = 'icasso'
"""

##### FUNCTIONS #####
def load_data(input_file):
    return pd.read_table(input_file, index_col=0)


def iterate_ica(data, n_components, iterations):
    S_all = []
    A_all = []
    ica = FastICA(n_components=n_components, tol=0.1)
    for it in range(iterations):
        S = ica.fit_transform(data)
        A = ica.mixing_
        S_all.append(S)
        A_all.append(A)

    S_all = np.hstack(S_all)
    A_all = np.hstack(A_all)
    return S_all, A_all


def compute_distance(X):
    D = 1 - np.abs(np.corrcoef(X.T))
    return D


def cluster_components(X, n_components, **kws):
    model = AgglomerativeClustering(n_clusters=n_components, linkage="average", **kws)
    labels = model.fit_predict(X.T)  # features are in rows
    return labels


def get_iteration_signs(data, S_all, A_all, n_components, iterations):
    """
        Correct direction of every iteration of ICA with respect to the first one.
        """
    # get component that explains the most variance of X from first iteration
    S0 = S_all[:, 0:n_components]
    A0 = A_all[:, 0:n_components]
    tss = []
    for i in range(n_components):
        pred = np.dot(S0[:, i].reshape(-1, 1), A0[:, i].reshape(1, -1))
        tss.append(np.sum((data - pred) ** 2))  # total sum of squares
    best_comp = S0[:, np.argmax(tss)]

    # correlate best component with the rest of iterations to decide signs
    signs = np.full(S_all.shape[1], np.nan)
    signs[0:n_components] = 1
    for it in range(1, iterations):
        start = n_components * it
        end = start + n_components
        S_it = S_all[:, start:end]

        correl = np.apply_along_axis(
            lambda x: pearsonr(x, best_comp)[0], axis=0, arr=S_it
        )
        best_correl = correl[np.argmax(np.abs(correl))]
        signs[start:end] = np.sign(best_correl)

    return signs


def compute_centroids_w_signs(S_all, A_all, labels, signs):
    # correct signs
    S_all = (signs * S_all).T
    A_all = (signs * A_all).T

    # put clusters together
    S = []
    A = []
    for label in np.unique(labels):
        # subset
        idx = np.where(labels == label)[0]
        S_clust = S_all[idx, :]
        A_clust = A_all[idx, :]

        # save centroids
        S.append(np.array(S_clust).T.mean(axis=1))
        A.append(np.array(A_clust).T.mean(axis=1))

    # prepare output
    S = np.vstack(S).T
    A = np.vstack(A).T

    return S, A


def compute_centroids_wo_signs(S_all, A_all, labels):
    """
    Based on https://github.com/SBRG/precise-db/blob/master/scripts/cluster_components.py
    """
    S = []
    A = []
    for label in np.unique(labels):
        # subset
        idx = labels == label
        S_clust = S_all[:, idx]
        A_clust = A_all[:, idx]

        # first item is base component
        Svec0 = S_clust[:, 0]
        Avec0 = A_clust[:, 0]

        # Make sure base component is facing positive
        if abs(min(Svec0)) > max(Svec0):
            Svec0 = -Svec0
            Avec0 = -Avec0

        S_single = [Svec0]
        A_single = [Avec0]
        # Add in rest of components
        for j in range(1, S_clust.shape[1]):
            Svec = S_clust[:, j]
            Avec = A_clust[:, j]
            if pearsonr(Svec, Svec0)[0] > 0:
                S_single.append(Svec)
                A_single.append(Avec)
            else:
                S_single.append(-Svec)
                A_single.append(-Avec)

        # save centroids
        S.append(np.array(S_single).T.mean(axis=1))
        A.append(np.array(A_single).T.mean(axis=1))

    # prepare output
    S = np.vstack(S).T
    A = np.vstack(A).T

    return S, A


def get_centroids(S_all, A_all, labels, signs=None):
    if signs is None:
        S_robust, A_robust = compute_centroids_wo_signs(S_all, A_all, labels)
    else:
        S_robust, A_robust = compute_centroids_w_signs(S_all, A_all, labels, signs)
    return S_robust, A_robust


def run_pca(n_components, S_all):
    return PCA(n_components).fit_transform(S_all.T).T


def icasso(data, iterations):
    n_components = int(0.5 * data.shape[1])
    mem = {}
    t = {}
    start = time.time()
    mem["iterate_ica"], (S_all, A_all) = memory_usage(
        (iterate_ica, (data, n_components, iterations)), **MEM_KWS
    )
    t["iterate_ica"] = time.time() - start
    start = time.time()
    mem["compute_distance"], D = memory_usage((compute_distance, (S_all,)), **MEM_KWS)
    t["compute_distance"] = time.time() - start
    start = time.time()
    mem["cluster_components"], labels = memory_usage(
        (cluster_components, (D, n_components), {"affinity": "precomputed"}), **MEM_KWS
    )
    t["cluster_components"] = time.time() - start
    start = time.time()
    mem["get_centroids"], (S_robust, A_robust) = memory_usage(
        (get_centroids, (S_all, A_all, labels)), **MEM_KWS
    )
    t["get_centroids"] = time.time() - start

    return S_robust, A_robust, mem, t


def robustica(data, iterations, do_pca=False):
    n_components = int(0.5 * data.shape[1])
    mem = {}
    t = {}

    start = time.time()
    mem["iterate_ica"], (S_all, A_all) = memory_usage(
        (iterate_ica, (data, n_components, iterations)), **MEM_KWS
    )
    t["iterate_ica"] = time.time() - start

    if do_pca:
        start = time.time()
        mem["run_pca"], pcs = memory_usage((run_pca, (n_components, S_all)), **MEM_KWS)
        t["run_pca"] = time.time() - start

        start = time.time()
        mem["cluster_components"], labels = memory_usage(
            (cluster_components, (pcs, n_components)), **MEM_KWS
        )
        t["cluster_components"] = time.time() - start
    else:
        start = time.time()
        mem["cluster_components"], labels = memory_usage(
            (cluster_components, (S_all, n_components)), **MEM_KWS
        )
        t["cluster_components"] = time.time() - start

    start = time.time()
    mem["get_iteration_signs"], signs = memory_usage(
        (get_iteration_signs, (data.values, S_all, A_all, n_components, iterations)),
        **MEM_KWS
    )
    t["get_iteration_signs"] = time.time() - start

    start = time.time()
    mem["get_centroids"], (S_robust, A_robust) = memory_usage(
        (get_centroids, (S_all, A_all, labels, signs)), **MEM_KWS
    )
    t["get_centroids"] = time.time() - start

    return S_robust, A_robust, mem, t


def robustica_pca(data, iterations):
    return robustica(data, iterations, do_pca=True)


def get_performance(mem, t):
    performance = pd.concat(
        [
            pd.DataFrame(mem[key], columns=["memory", "timestamp"]).assign(function=key)
            for key in mem.keys()
        ]
    )
    performance = pd.merge(
        performance,
        pd.Series(t).reset_index().rename(columns={"index": "function", 0: "time"}),
        on="function",
    )

    return performance


def evaluate_performance(data, algorithm, iterations):
    func = eval(algorithm)
    S_robust, A_robust, mem, t = func(data, iterations)
    performance = get_performance(mem, t)
    performance['algorithm'] = algorithm
    return S_robust, A_robust, performance


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str)
    parser.add_argument("--output_file", type=str)
    parser.add_argument("--iterations", type=int)
    parser.add_argument("--algorithm", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    input_file = args.input_file
    output_file = args.output_file
    iterations = args.iterations
    algorithm = args.algorithm

    data = load_data(input_file)
    _, _, result = evaluate_performance(data, algorithm, iterations)
    result.to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
