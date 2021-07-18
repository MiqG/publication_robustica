#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
#

import argparse
import os
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_samples, pairwise_distances
from memory_profiler import memory_usage
import time
from robustica import RobustICA, compute_iq, corrmats
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

fdrtool = importr("fdrtool")

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}
MEM_KWS = {
    "interval": 0.1,
    "timestamps": True,
    "include_children": True,
    "retval": True,
}

ALGORITHMS = {
    "icasso": {
        "robust_infer_signs": False,
        "robust_dimreduce": False,
        "robust_kws": {"n_clusters": 100},
    },
    "robustica_nosign": {
        "robust_infer_signs": False,
        "robust_dimreduce": False,
        "robust_kws": {"affinity": "euclidean", "linkage": "average"},
    },
    "robustica": {
        "robust_infer_signs": True,
        "robust_dimreduce": False,
        "robust_kws": {"affinity": "euclidean", "linkage": "average"},
    },
    "robustica_pca": {
        "robust_infer_signs": True,
        "robust_dimreduce": True,
        "robust_kws": {"affinity": "euclidean", "linkage": "average"},
    },
}

"""
Development
-----------
ROOT = '/home/miquel/projects/publication_robustica'
PREP_DIR = os.path.join(ROOT,'data','prep','ica_runs')
S_all_file = os.path.join(PREP_DIR,'Sastry2019','S.pickle.gz')
A_all_file = os.path.join(PREP_DIR,'Sastry2019','A.pickle.gz')
iterations = 100
n_components = 100
algorithms = 'icasso,robustica_pca'.split(',')
"""


##### FUNCTIONS #####
def load_data(S_all_file, A_all_file):
    """
    We consider the first iteration as our ground truth.
    """
    # ICA runs
    S_all = pd.read_pickle(S_all_file)
    A_all = pd.read_pickle(A_all_file)

    return S_all, A_all


class icasso:
    def __init__(self, **kws):
        self.clustering = AgglomerativeClustering(
            affinity="precomputed", linkage="average", **kws
        )

    def compute_distance(self, X):
        return 1 - np.abs(np.corrcoef(X))

    def fit(self, X):
        self.clustering.fit(X)
        self.labels_ = self.clustering.labels_


ALGORITHMS["icasso"]["robust_method"] = icasso


def compute_robust_components(rica):
    mem = {}
    t = {}

    # get objects
    S_all = rica.S_all
    A_all = rica.A_all

    # infer component signs across runs
    if rica.robust_infer_signs:
        start = time.time()
        mem["infer_components_signs"], signs = memory_usage(
            (
                rica._infer_components_signs,
                (S_all, rica.n_components, rica.robust_runs),
            ),
            **MEM_KWS
        )
        t["infer_components_signs"] = time.time() - start
    else:
        signs = np.ones(S_all.shape[1])

    Y = S_all * signs

    # reduce dimensions
    if rica.robust_dimreduce != False:
        start = time.time()
        mem["run_pca"], Y = memory_usage(
            (rica.dimreduce.fit_transform, (Y.T,)), **MEM_KWS
        )
        t["run_pca"] = time.time() - start
        Y = Y.T

    # compute dissimilarity
    if "icasso" in str(rica.clustering):
        start = time.time()
        mem["compute_distance"], Y = memory_usage(
            (rica.clustering.compute_distance, (Y.T,)), **MEM_KWS
        )
        t["compute_distance"] = time.time() - start

    # cluster
    start = time.time()
    mem["cluster_components"], _ = memory_usage(
        (rica.clustering.fit, (Y.T,)), **MEM_KWS
    )
    labels = rica.clustering.labels_
    t["cluster_components"] = time.time() - start

    # compute robust components
    start = time.time()
    (
        mem["compute_centroids"],
        (S, A, S_std, A_std, clustering_stats, orientation),
    ) = memory_usage(
        (rica._compute_centroids, (S_all * signs, A_all * signs, labels)), **MEM_KWS
    )
    t["compute_centroids"] = time.time() - start

    # add clustering info
    clustering_info = pd.DataFrame(
        {
            "component": np.arange(S_all.shape[1]),
            "cluster_id": labels,
            "sign": signs,
            "orientation": orientation,
        }
    )

    return S, A, S_std, A_std, mem, t, clustering_info


def prep_performance(mem, t):
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


def evaluate_performance(rica):
    start = time.time()
    tmp_mem, (S, A, S_std, A_std, mem, t, clustering_info) = memory_usage(
        (compute_robust_components, (rica,)), **MEM_KWS
    )
    tmp_t = time.time() - start

    # prepare performance evaluation
    mem["evaluate_performance"] = tmp_mem
    t["evaluate_performance"] = tmp_t
    performance = prep_performance(mem, t)

    # prepare clustering info
    labels = clustering_info["cluster_id"].values
    X = (
        rica.S_all
        * clustering_info["sign"].values
        * clustering_info["orientation"].values
    ).T
    clustering_info = pd.merge(clustering_info, compute_iq(X, labels), on="cluster_id")
    clustering_info["silhouette_euclidean"] = silhouette_samples(X, labels)
    D = 1 - np.abs(np.corrcoef(rica.S_all.T))
    clustering_info["silhouette_pearson"] = silhouette_samples(
        D, labels, metric="precomputed"
    )

    return S, A, S_std, A_std, performance, clustering_info


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


def get_relative_timestamp(x):
    x0 = x[0]
    return x - x0


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--S_all_file", type=str)
    parser.add_argument("--A_all_file", type=str)
    parser.add_argument("--output_file", type=str)
    parser.add_argument("--n_components", type=int)
    parser.add_argument("--iterations", type=int)
    parser.add_argument("--algorithms", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    S_all_file = args.S_all_file
    A_all_file = args.A_all_file
    output_file = args.output_file
    n_components = args.n_components
    iterations = args.iterations
    algorithms = args.algorithms.split(",")

    # load data
    S_all, A_all = load_data(S_all_file, A_all_file)

    # evaluate performance
    mappings = []
    performances = []
    clustering_infos = []
    for algorithm in algorithms:
        print(algorithm)

        modules = {}
        for it in np.arange(5,iterations+5,5)[::-1]:
            print(it)

            rica_kws = ALGORITHMS[algorithm].copy()
            rica_kws["n_components"] = n_components
            rica_kws["robust_runs"] = it
            rica = RobustICA(**rica_kws)
            rica.S_all = S_all.values[:, : it * n_components]
            rica.A_all = A_all.values[:, : it * n_components]

            S, A, S_std, A_std, performance, clustering_info = evaluate_performance(
                rica
            )

            # summarize
            ## mean clustering stats
            clustering_info = (
                clustering_info.groupby(["cluster_id"])[
                    ["silhouette_euclidean", "silhouette_pearson", "iq"]
                ]
                .mean()
                .reset_index()
            )
            ## max time and memory
            performance = performance.loc[
                performance["function"] != "evaluate_performance"
            ].copy()
            performance["rel_timestamp"] = get_relative_timestamp(
                performance["timestamp"]
            )
            performance = pd.DataFrame(performance[["rel_timestamp", "memory"]].max()).T
            performance.columns = ["time", "max_memory"]

            # obtain modules
            modules[it] = get_modules(S)
            mapping = pd.DataFrame(
                {
                    "cluster_id": np.arange(S.shape[1]),
                    "algorithm": algorithm,
                    "runs": it,
                    "jaccard": compare_modules(modules[it], modules[iterations]),
                }
            )

            # add info
            performance["algorithm"] = algorithm
            clustering_info["algorithm"] = algorithm
            performance["runs"] = it
            clustering_info["runs"] = it

            # save
            mappings.append(mapping)
            performances.append(performance)
            clustering_infos.append(clustering_info)

    mappings = pd.concat(mappings)
    performances = pd.concat(performances)
    clustering_infos = pd.concat(clustering_infos)

    # prepare output
    output = pd.merge(mappings, clustering_infos, on=["algorithm", "runs", "cluster_id"])
    output = pd.merge(output, performances, on=["algorithm", "runs"])

    # save outputs
    output.to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
