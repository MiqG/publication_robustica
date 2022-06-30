#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Evaluate differenti clustering metrics on AgglomerativeClustering and different
# clustering methods using absolute pearson correlation dissimilarity.

import argparse
import os
import time
import pandas as pd
import numpy as np
from sklearn.metrics import silhouette_samples
from memory_profiler import memory_usage
from robustica import RobustICA, compute_iq, corrmats

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}
MEM_KWS = {
    "interval": 0.1,
    "timestamps": True,
    "include_children": True,
    "retval": True,
}

RICA_KWS = {
    "metrics": {
        "euclidean": {"robust_infer_signs": True, "robust_dimreduce": True},
        "l1": {"robust_infer_signs": True, "robust_dimreduce": True},
        "l2": {"robust_infer_signs": True, "robust_dimreduce": True},
        "manhattan": {"robust_infer_signs": True, "robust_dimreduce": True},
        "cosine": {"robust_infer_signs": True, "robust_dimreduce": True},
        "precomputed": {"robust_infer_signs": False, "robust_dimreduce": False},
    },
    "methods": {
        "AgglomerativeClustering": {
            "robust_kws": {"affinity": "precomputed", "linkage": "average"},
        },
        "AffinityPropagation": {"robust_kws": {"affinity": "precomputed"}}, # the higher the more similar (not distant)
        "DBSCAN": {"robust_kws": {"metric": "precomputed"}},
        "FeatureAgglomeration": {"robust_kws": {"affinity": "precomputed", "linkage":"average"}},
        "OPTICS": {"robust_kws": {"metric": "precomputed"}},
        "KMedoids": {"robust_kws": {"metric": "precomputed"}},
        "CommonNNClustering": {"robust_kws": {"metric": "precomputed"}},
    },
    "methods_defaults": {}
}

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_robustica'
PREP_DIR = os.path.join(ROOT,'data','prep')
S_all_file = os.path.join(PREP_DIR,'ica_runs','Sastry2019','S.pickle.gz')
A_all_file = os.path.join(PREP_DIR,'ica_runs','Sastry2019','A.pickle.gz')
iterations = 100
property_type = 'methods'
properties_oi = 'OPTICS'.split(',')
"""


##### FUNCTIONS #####
def load_data(
    S_all_file, A_all_file,
):
    # ICA runs
    S_all = pd.read_pickle(S_all_file)
    A_all = pd.read_pickle(A_all_file)

    return S_all, A_all
    

def compute_pearson_affinity(X):
    print("using custom affinity")
    # similarly to negative euclidean distance
    D = np.clip((1 - np.abs(np.corrcoef(X.T))), 0, 1)
    D = (-1) * D
    return D
RICA_KWS["methods"]["AffinityPropagation"]["robust_precompdist_func"] = compute_pearson_affinity
    
    
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
    if np.isin(["precomputed"], list(rica.robust_kws.values()))[0]:
        start = time.time()
        mem["compute_distance"], Y = memory_usage(
            (rica.robust_precompdist_func, (Y,)), **MEM_KWS
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
    try:
        clustering_info["silhouette_euclidean"] = silhouette_samples(X, labels)
    except:
        clustering_info["silhouette_euclidean"] = np.nan
        
    D = 1 - np.abs(np.corrcoef(rica.S_all.T))
    try:
        clustering_info["silhouette_pearson"] = silhouette_samples(
            D, labels, metric="precomputed"
        )
    except:    
        clustering_info["silhouette_pearson"] = np.nan
        
    clustering_info = pd.merge(
        clustering_info, compute_iq(X, labels), on="cluster_id"
    )
    return S, A, S_std, A_std, performance, clustering_info


def evaluate_clustering(property_type, property_oi, S_all, A_all, iterations):

    # update rica_kws accordingly
    rica_kws = RICA_KWS[property_type][property_oi].copy()
    rica_kws["n_components"] = int(S_all.shape[1] / iterations)
    rica_kws["robust_runs"] = iterations
    if property_type == "metrics":
        rica_kws["robust_kws"] = {"affinity": property_oi, "linkage": "average"}
    elif property_type == "methods":
        rica_kws["robust_infer_signs"] = False
        rica_kws["robust_dimreduce"] = False
        rica_kws["robust_method"] = property_oi

    print(rica_kws)

    rica = RobustICA(**rica_kws)
    rica.S_all = S_all.values
    rica.A_all = A_all.values

    S, A, S_std, A_std, performance, clustering_info = evaluate_performance(rica)

    # add info current loop
    performance["property_type"] = property_type
    clustering_info["property_type"] = property_type
    performance["property_oi"] = property_oi
    clustering_info["property_oi"] = property_oi

    return S, A, S_std, A_std, performance, clustering_info


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--S_all_file", type=str)
    parser.add_argument("--A_all_file", type=str)
    parser.add_argument("--output_file", type=str)
    parser.add_argument("--iterations", type=int)
    parser.add_argument("--property_type", type=str)
    parser.add_argument("--properties_oi", type=str)
    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    S_all_file = args.S_all_file
    A_all_file = args.A_all_file
    output_file = args.output_file
    iterations = args.iterations
    property_type = args.property_type
    properties_oi = args.properties_oi.split(",")

    # load data
    S_all, A_all = load_data(S_all_file, A_all_file)

    # evaluate performance
    performances = []
    clustering_infos = []
    S_robusts = {}
    S_stds = {}
    A_robusts = {}
    A_stds = {}
    for property_oi in properties_oi:
        # run
        S, A, S_std, A_std, p, c = evaluate_clustering(property_type, property_oi, S_all, A_all, iterations)
        
        # save
        performances.append(p)
        clustering_infos.append(c)
        
        labels = np.unique(c["cluster_id"])
        S_robusts[property_oi] = pd.DataFrame(S, index=S_all.index, columns=labels)
        S_stds[property_oi] = pd.DataFrame(S_std, index=S_all.index, columns=labels)
        A_robusts[property_oi] = pd.DataFrame(A, index=A_all.index, columns=labels)
        A_stds[property_oi] = pd.DataFrame(A_std, index=A_all.index, columns=labels)

    performances = pd.concat(performances)
    clustering_infos = pd.concat(clustering_infos)

    # save outputs
    performances.to_csv(output_file, **SAVE_PARAMS)
    output_dir = os.path.dirname(output_file)
    clustering_infos.to_csv(
        os.path.join(output_dir, "clustering_info.tsv.gz"), **SAVE_PARAMS
    )
    for algorithm in properties_oi:
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
