#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Run different algorithms to perform robust ICA and prove that our
# sign correction allows for different clustering methods.
#

import argparse
import os
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering, DBSCAN
from sklearn.metrics import silhouette_samples, r2_score, pairwise_distances
from memory_profiler import memory_usage
import time
from robustica import RobustICA, compute_iq, corrmats

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}
MEM_KWS = {
    "interval": 0.1,
    "timestamps": True,
    "include_children": True,
    "retval": True,
}

METHOD_OI = "DBSCAN"
ALGORITHMS = {
    "DBSCAN": {
        "icasso": {
            "robust_infer_signs": False,
            "robust_dimreduce": False,
            "robust_kws": {"min_samples": int(100 * 0.5)},
        },
        "icasso_infer_signs": {
            "robust_infer_signs": True,
            "robust_dimreduce": False,
            "robust_kws": {"min_samples": int(100 * 0.5)},            
        },
        "icasso_pca": {
            "robust_infer_signs": True,
            "robust_dimreduce": True,
            "robust_kws": {"min_samples": int(100 * 0.5)},            
        },
        "icasso_pca_nosign": {
            "robust_infer_signs": False,
            "robust_dimreduce": True,
            "robust_kws": {"min_samples": int(100 * 0.5)},            
        },
        "robustica_nosign": {
            "robust_infer_signs": False,
            "robust_dimreduce": False,
            "robust_method": "DBSCAN",
            "robust_kws": {"metric": "euclidean"},
        },
        "robustica": {
            "robust_infer_signs": True,
            "robust_dimreduce": False,
            "robust_method": "DBSCAN",
            "robust_kws": {"metric": "euclidean"},
        },
        "robustica_pca": {
            "robust_infer_signs": True,
            "robust_dimreduce": True,
            "robust_method": "DBSCAN",
            "robust_kws": {"metric": "euclidean"},
        }
    },
    "AgglomerativeClustering": {
        "icasso": {
            "robust_infer_signs": False,
            "robust_dimreduce": False,
            "robust_kws": {"n_clusters": 100},
        },
        "icasso_infer_signs": {
            "robust_infer_signs": True,
            "robust_dimreduce": False,
            "robust_kws": {"n_clusters": 100},
        },
        "icasso_pca": {
            "robust_infer_signs": True,
            "robust_dimreduce": True,
            "robust_kws": {"n_clusters": 100},
        },
        "icasso_pca_nosign": {
            "robust_infer_signs": False,
            "robust_dimreduce": True,
            "robust_kws": {"n_clusters": 100},
        },
        "robustica_nosign": {
            "robust_infer_signs": False,
            "robust_dimreduce": False,
            "robust_method": "AgglomerativeClustering",
            "robust_kws": {"affinity": "euclidean", "linkage": "average"},
        },
        "robustica": {
            "robust_infer_signs": True,
            "robust_dimreduce": False,
            "robust_method": "AgglomerativeClustering",
            "robust_kws": {"affinity": "euclidean", "linkage": "average"},
        },
        "robustica_pca": {
            "robust_infer_signs": True,
            "robust_dimreduce": True,
            "robust_method": "AgglomerativeClustering",
            "robust_kws": {"affinity": "euclidean", "linkage": "average"},
        },
    },
}

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_robustica'
PREP_DIR = os.path.join(ROOT,'data','prep','ica_runs')
S_all_file = os.path.join(PREP_DIR,'Sastry2019','S.pickle.gz')
A_all_file = os.path.join(PREP_DIR,'Sastry2019','A.pickle.gz')
iterations = 100
algorithms = 'icasso,icasso_infer_signs,robustica_pca'.split(',')
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
        if METHOD_OI == "DBSCAN":
            self.clustering = DBSCAN(metric="precomputed", **kws)
        elif METHOD_OI == "AgglomerativeClustering":
            self.clustering = AgglomerativeClustering(
                affinity="precomputed", linkage="average", **kws
            )

    def compute_distance(self, X):
        return 1 - np.abs(np.corrcoef(X))

    def fit(self, X):
        self.clustering.fit(X)
        self.labels_ = self.clustering.labels_


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
    clustering_info["silhouette_euclidean"] = silhouette_samples(X, labels)
    D = 1 - np.abs(np.corrcoef(rica.S_all.T))
    clustering_info["silhouette_pearson"] = silhouette_samples(
        D, labels, metric="precomputed"
    )
    clustering_info = pd.merge(clustering_info, compute_iq(X, labels), on="cluster_id")

    return S, A, S_std, A_std, performance, clustering_info


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--S_all_file", type=str)
    parser.add_argument("--A_all_file", type=str)
    parser.add_argument("--output_file", type=str)
    parser.add_argument("--iterations", type=int)
    parser.add_argument("--algorithms", type=str)

    args = parser.parse_args()

    return args

# add custom method class
ALGORITHMS[METHOD_OI]["icasso"]["robust_method"] = icasso
ALGORITHMS[METHOD_OI]["icasso_infer_signs"]["robust_method"] = icasso
ALGORITHMS[METHOD_OI]["icasso_pca"]["robust_method"] = icasso
ALGORITHMS[METHOD_OI]["icasso_pca_nosign"]["robust_method"] = icasso

def main():
    args = parse_args()
    S_all_file = args.S_all_file
    A_all_file = args.A_all_file
    output_file = args.output_file
    iterations = args.iterations
    algorithms = args.algorithms.split(",")

    # load data
    S_all, A_all = load_data(S_all_file, A_all_file)

    # evaluate performance
    performances = []
    clustering_infos = []
    S_robusts = {}
    S_stds = {}
    A_robusts = {}
    A_stds = {}
    for algorithm in algorithms:
        print(algorithm, METHOD_OI)

        rica_kws = ALGORITHMS[METHOD_OI][algorithm].copy()
        rica_kws["n_components"] = int(S_all.shape[1] / iterations)
        rica_kws["robust_runs"] = iterations
        rica = RobustICA(**rica_kws)
        rica.S_all = S_all.values
        rica.A_all = A_all.values

        S, A, S_std, A_std, performance, clustering_info = evaluate_performance(rica)

        # add algorithm
        performance["algorithm"] = algorithm
        clustering_info["algorithm"] = algorithm

        # prepare outpus
        labels = np.unique(rica.clustering.labels_)
        S_robusts[algorithm] = pd.DataFrame(S, index=S_all.index, columns=labels)
        S_stds[algorithm] = pd.DataFrame(S_std, index=S_all.index, columns=labels)
        A_robusts[algorithm] = pd.DataFrame(A, index=A_all.index, columns=labels)
        A_stds[algorithm] = pd.DataFrame(A_std, index=A_all.index, columns=labels)
        performances.append(performance)
        clustering_infos.append(clustering_info)

    performances = pd.concat(performances)
    clustering_infos = pd.concat(clustering_infos)

    # save outputs
    performances.to_csv(output_file, **SAVE_PARAMS)
    output_dir = os.path.dirname(output_file)
    clustering_infos.to_csv(
        os.path.join(output_dir, "clustering_info.tsv.gz"), **SAVE_PARAMS
    )
    for algorithm in algorithms:
        S_robusts[algorithm].reset_index().to_csv(
            os.path.join(output_dir, "%s-S.tsv.gz") % algorithm, **SAVE_PARAMS
        )
        S_stds[algorithm].reset_index().to_csv(
            os.path.join(output_dir, "%s-S_std.tsv.gz") % algorithm, **SAVE_PARAMS
        )
        A_robusts[algorithm].reset_index().to_csv(
            os.path.join(output_dir, "%s-A.tsv.gz") % algorithm, **SAVE_PARAMS
        )
        A_stds[algorithm].reset_index().to_csv(
            os.path.join(output_dir, "%s-A_std.tsv.gz") % algorithm, **SAVE_PARAMS
        )


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
