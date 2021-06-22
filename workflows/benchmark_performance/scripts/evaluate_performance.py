#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Run different algorithms to perform robust ICA and prove that our 
# sign correction allows for different clustering methods.
#
# Outline
# -------
# We will compare 4 different algorithms:
# - icasso: the standard
# - robustica: allows fully customizable clustering algorithm thanks to our
#   method to guess component signs.
# - robustica_nosign: to check what happens when we don't correct component
#   signs and use a different clustering metric.
# - robustica_pca: to check if compressing the feature space speeds up the 
#   algorithm.
# 
# We will record:
# - Performance
#     - memory usage of every step of the algorithm
#     - time spent in every step of the algorithm
# - Clustering
#     - labels
#     - signs
#     - final orientation
#     - Silhouette scores w/ absolute pearson and euclidean metrics
#     - Validation scores with the ground truth set
#        - maximum R2
#        - maximum absolute pearson
#        - minimum distance

import argparse
import os
import pandas as pd
import numpy as np
from sklearn.decomposition import FastICA, PCA
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_samples, r2_score, pairwise_distances
from scipy.stats import pearsonr
from memory_profiler import memory_usage
import time
import itertools

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
S_true_file = os.path.join(RESULTS_DIR,'files','validation_data','Sastry2019','S.tsv.gz')
A_true_file = os.path.join(RESULTS_DIR,'files','validation_data','Sastry2019','A.tsv.gz')
S_all_file = os.path.join(RESULTS_DIR,'files','ica_iterations','validation_data','Sastry2019','S.pickle')
A_all_file = os.path.join(RESULTS_DIR,'files','ica_iterations','validation_data','Sastry2019','A.pickle')
iterations = 100
algorithms = 'icasso,robustica_nosign,robustica_pca'.split(',')
"""


##### FUNCTIONS #####
def load_data(S_all_file, A_all_file, S_true_file, A_true_file):
    """
    We consider the first iteration as our ground truth.
    """
    # ICA runs
    S_all = pd.read_pickle(S_all_file).values
    A_all = pd.read_pickle(A_all_file).values
    
    # Ground truth
    S_true = pd.read_table(S_true_file, index_col=0).values
    A_true = pd.read_table(A_true_file, index_col=0).values

    return S_true, A_true, S_all, A_all


def compute_distance(X):
    D = 1 - np.abs(np.corrcoef(X.T))
    return D


def corrmats(X,Y):
    """
    Correlation matrix between rows in X and rows in Y.
    """
    X_cent = X - X.mean(axis=1).reshape(-1,1)
    Y_cent = Y - Y.mean(axis=1).reshape(-1,1)
    
    num = X_cent.dot(Y_cent.T)
    den = np.sqrt(
        (X_cent**2).sum(axis=1).reshape(-1,1) * (Y_cent**2).sum(axis=1).reshape(1,-1)
    )
    r = num / den
    return r
    

def cluster_components(X, n_components, **kws):
    model = AgglomerativeClustering(n_clusters=n_components, linkage="average", **kws)
    labels = model.fit_predict(X.T)  # features are in rows
    return labels


def compute_iq(D, labels):
    """
    Compute cluster quality index.
    """
    abs_r = 1 - D
    iqs = []
    for label in np.unique(labels):
        idx_cluster = labels == label
        
        avg_in = abs_r[idx_cluster,:][:,idx_cluster].mean()
        avg_out = abs_r[idx_cluster,:][:,~idx_cluster].mean()
        iq_cluster = avg_in - avg_out
        iqs.append(iq_cluster)
    
    df = pd.DataFrame({
        'label': np.unique(labels),
        'iq': iqs
    })
    
    return df
        

def infer_components_signs(S_all, n_components, iterations):
    """
    Correct direction of every run of ICA with respect to the first one.
    """
    # components from first run are the reference
    S0 = S_all[:,0:n_components]

    # correlate best component with the rest of iterations to decide signs
    signs = np.full(S_all.shape[1], np.nan)
    signs[0:n_components] = 1
    for it in range(1, iterations):
        start = n_components * it
        end = start + n_components
        S_it = S_all[:, start:end]

        correl = corrmats(S0.T,S_it.T)
        rows_oi = np.abs(correl).argmax(axis=0)
        cols_oi = np.arange(correl.shape[1])
        best_correls = correl[(rows_oi,cols_oi)]
        signs[start:end] = np.sign(best_correls)

    return signs


def get_centroids(S_all, A_all, labels):
    """
    Based on https://github.com/SBRG/precise-db/blob/master/scripts/cluster_components.py
    """
    S = []
    A = []
    orientation = np.full(S_all.shape[1], np.nan)
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
            ori = [-1]
        else:
            ori = [1]
        
        S_single = [Svec0]
        A_single = [Avec0]
        
        # Add in rest of components
        for j in range(1, S_clust.shape[1]):
            Svec = S_clust[:, j]
            Avec = A_clust[:, j]
            if pearsonr(Svec, Svec0)[0] > 0:
                S_single.append(Svec)
                A_single.append(Avec)
                ori.append(1)
            else:
                S_single.append(-Svec)
                A_single.append(-Avec)
                ori.append(-1)
                
        # save centroids
        S_single = np.array(S_single).T
        A_single = np.array(A_single).T
        S.append(S_single.mean(axis=1))
        A.append(A_single.mean(axis=1))
        
        # save orientation
        orientation[idx] = ori

    # prepare output
    S = np.vstack(S).T
    A = np.vstack(A).T

    return S, A, orientation


def run_pca(n_components, S_all):
    return PCA(n_components).fit_transform(S_all.T).T


def icasso(S_all, A_all, iterations):
    # init
    n_components = int(S_all.shape[1]/iterations)
    mem = {}
    t = {}
    
    # compute distance
    start = time.time()
    mem["compute_distance"], D = memory_usage((compute_distance, (S_all,)), **MEM_KWS)
    t["compute_distance"] = time.time() - start
    
    # cluster components
    start = time.time()
    mem["cluster_components"], labels = memory_usage(
        (cluster_components, (D, n_components), {"affinity": "precomputed"}), **MEM_KWS
    )
    t["cluster_components"] = time.time() - start
    
    # get centroids
    start = time.time()
    mem["get_centroids"], (S_robust, A_robust, orientation) = memory_usage(
        (get_centroids, (S_all, A_all, labels)), **MEM_KWS
    )
    t["get_centroids"] = time.time() - start
    
    # add clustering info
    clustering_info = pd.DataFrame({
        'component': np.arange(S_all.shape[1]),
        'label': labels,
        'sign': 1, # we don't deal with signs in icasso
        'orientation': orientation
    })

    
    return S_robust, A_robust, mem, t, clustering_info


def robustica(S_all, A_all, iterations, do_pca=False, correct_signs=True):
    n_components = int(S_all.shape[1]/iterations)
    mem = {}
    t = {}
    
    # correct component signs across runs
    if correct_signs:
        start = time.time()
        mem["infer_components_signs"], signs = memory_usage(
            (infer_components_signs, (S_all, n_components, iterations)),
            **MEM_KWS
        )
        t["infer_components_signs"] = time.time() - start
        S_all = np.multiply(S_all, signs)
        A_all = np.multiply(A_all, signs)
    else:
        signs = np.ones(S_all.shape[1])
    
    # compress feature dimensions with PCA
    if do_pca:
        start = time.time()
        mem["run_pca"], X = memory_usage((run_pca, (n_components, S_all)), **MEM_KWS)
        t["run_pca"] = time.time() - start
    else:
        X = S_all
    
    # cluster components
    start = time.time()
    mem["cluster_components"], labels = memory_usage(
        (cluster_components, (X, n_components)), **MEM_KWS
    )
    t["cluster_components"] = time.time() - start

    # get centroids
    start = time.time()
    mem["get_centroids"], (S_robust, A_robust, orientation) = memory_usage(
        (get_centroids, (S_all, A_all, labels)), **MEM_KWS
    )
    t["get_centroids"] = time.time() - start

    # add clustering info
    clustering_info = pd.DataFrame({
        'component': np.arange(S_all.shape[1]),
        'label': labels,
        'sign': signs,
        'orientation': orientation
    })
    
    return S_robust, A_robust, mem, t, clustering_info


def robustica_nosign(S_all, A_all, iterations):
    return robustica(S_all, A_all, iterations, do_pca=True, correct_signs=False)


def robustica_pca(S_all, A_all, iterations):
    return robustica(S_all, A_all, iterations, do_pca=True)


def compute_r2(Y_true, Y_pred):
    assert Y_true.shape == Y_pred.shape
    
    r2 = np.full((Y_true.shape[1], Y_pred.shape[1]), np.nan)
    for i in range(Y_true.shape[1]):
        for j in range(Y_pred.shape[1]):
            r2[i,j] = r2_score(Y_true[:,i], Y_pred[:,j])
    return r2.max(axis=0)
    
    
def score_prediction(Y_true, Y_pred):
    r2 = compute_r2(Y_true, Y_pred)
    p = np.abs(corrmats(Y_true.T, Y_pred.T)).max(axis=0)
    d = pairwise_distances(Y_true.T, Y_pred.T, metric='euclidean').min(axis=0)
    df = pd.DataFrame({
        'max_r2': r2, 
        'max_pearson': p,
        'min_dist': d,
        'label': np.arange(Y_pred.shape[1])
    })
    return df 
    
    
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


def evaluate_performance(S_all, A_all, algorithm, iterations):
    func = eval(algorithm)

    start = time.time()
    tmp_mem, (S_robust, A_robust, mem, t, clustering_info) = memory_usage(
        (func, (S_all, A_all, iterations)), **MEM_KWS
    )
    tmp_t = time.time() - start
    
    # prepare performance evaluation
    mem["evaluate_performance"] = tmp_mem
    t["evaluate_performance"] = tmp_t
    performance = prep_performance(mem, t)
    performance['algorithm'] = algorithm
    
    # prepare clustering info
    labels = clustering_info['label'].values
    X = (S_all*clustering_info['sign'].values*clustering_info['orientation'].values).T
    clustering_info['silhouette_euclidean'] = silhouette_samples(X, labels)
    D = 1 - np.abs(np.corrcoef(S_all.T))
    clustering_info['silhouette_pearson'] = silhouette_samples(D, labels, metric='precomputed')
    clustering_info = pd.merge(clustering_info, compute_iq(D, labels), on='label')
    clustering_info['algorithm'] = algorithm 
    
    return S_robust, A_robust, performance, clustering_info


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--S_all_file", type=str)
    parser.add_argument("--A_all_file", type=str)
    parser.add_argument("--S_true_file", type=str)
    parser.add_argument("--A_true_file", type=str)
    parser.add_argument("--output_file", type=str)
    parser.add_argument("--iterations", type=int)
    parser.add_argument("--algorithms", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    S_all_file = args.S_all_file
    A_all_file = args.A_all_file
    S_true_file = args.S_true_file
    A_true_file = args.A_true_file
    output_file = args.output_file
    iterations = args.iterations
    algorithms = args.algorithms.split(',')
    
    # load data
    S_true, A_true, S_all, A_all = load_data(S_all_file, A_all_file, S_true_file, A_true_file)
    
    # evaluate performance
    performances = []
    clustering_infos = []
    S_robusts = {}
    A_robusts = {}
    for algorithm in algorithms:
        print(algorithm)
        S_robust, A_robust, performance, clustering_info = evaluate_performance(
            S_all, A_all, algorithm, iterations
        )
        
        # add robust components validation scores
        clustering_info = pd.merge(
            clustering_info,
            score_prediction(S_true, S_robust).assign(algorithm=algorithm),
            on=['algorithm','label']
        )
        
        # prepare outpus
        S_robusts[algorithm] = pd.DataFrame(S_robust)
        A_robusts[algorithm] = pd.DataFrame(A_robust)
        performances.append(performance)
        clustering_infos.append(clustering_info)
        
    performances = pd.concat(performances)
    clustering_infos = pd.concat(clustering_infos)
    
    # save outputs
    performances.to_csv(output_file, **SAVE_PARAMS)
    output_dir = os.path.dirname(output_file)
    clustering_infos.to_csv(
        os.path.join(output_dir,'clustering_info.tsv.gz'),
        **SAVE_PARAMS
    )
    for algorithm in algorithms:
        S_robusts[algorithm].reset_index().to_csv(
            os.path.join(output_dir,'%s-S.tsv.gz') % algorithm,
            **SAVE_PARAMS
        )
        A_robusts[algorithm].reset_index().to_csv(
            os.path.join(output_dir,'%s-A.tsv.gz') % algorithm,
            **SAVE_PARAMS
        )
        

##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
