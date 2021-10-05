#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Cluster and evaluate ICA runs.

import argparse
import pandas as pd
import numpy as np
import robustica

# variables
SAVE_PARAMS = {'sep':'\t', 'compression':'gzip', 'index':False}
ALGORITHMS = {
    "icasso": {
        "robust_infer_signs": False,
        "robust_dimreduce": False,
        "robust_kws": {"n_clusters": np.nan},
    },
    "robustica_pca": {
        "robust_infer_signs": True,
        "robust_dimreduce": True,
        "robust_kws": {"affinity": "euclidean", "linkage": "average"},
    },
}
RANDOM_SEED = 1234

##### FUNCTIONS #####
def load_data(S_all_file, A_all_file, n_genes):
    S_all = pd.read_pickle(S_all_file)
    A_all = pd.read_pickle(A_all_file)
    
    if n_genes is not None:
        print('Picking %s genes at random...' % n_genes)
        np.random.seed(RANDOM_SEED)
        idx = np.random.randint(0, S_all.shape[0], size=n_genes)
        S_all = S_all.iloc[idx].copy()
    
    return S_all, A_all


def cluster_ica_runs(S_all, A_all, n_components, robust_runs, algorithm):
    # cluster
    ALGORITHMS['icasso']['robust_kws']['n_clusters'] = n_components
    
    rica = robustica.RobustICA(
        n_components = n_components, 
        robust_runs = robust_runs,
        **ALGORITHMS[algorithm]
    )
    S, A, S_std, A_std, clustering_stats, signs, orientation = rica._compute_robust_components(S_all.values, A_all.values) 

    # evaluate clusters
    ## silhouette euclidean
    evaluation = rica.evaluate_clustering(S_all.values, rica.clustering.labels_, signs, orientation)
    evaluation.columns = ['cluster_id','silhouette_euclidean','iq']
    ## silhouette pearson
    D = 1 - np.abs(np.corrcoef(S_all.values.T))
    eval_pearson = rica.evaluate_clustering(D, rica.clustering.labels_, signs, orientation, metric='precomputed')
    eval_pearson = eval_pearson[['cluster_id','mean_silhouette']]\
                        .rename(columns={'mean_silhouette':'silhouette_pearson'})
    evaluation = pd.merge(evaluation, eval_pearson, on='cluster_id')

    # prepare outputs
    S = pd.DataFrame(S, index=S_all.index).reset_index()
    A = pd.DataFrame(A, index=A_all.index).reset_index()
    stats = pd.merge(clustering_stats, evaluation, on='cluster_id')
    stats['n_genes'] = S.shape[0]
    stats['n_samples'] = A.shape[0]
    
    return S, A, stats
    
    
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--S_all_file", type=str)
    parser.add_argument("--A_all_file", type=str)
    parser.add_argument("--S_file", type=str)
    parser.add_argument("--A_file", type=str)
    parser.add_argument("--stats_file", type=str)
    parser.add_argument("--n_components", type=int)
    parser.add_argument("--robust_runs", type=int)
    parser.add_argument("--algorithm", type=str)
    parser.add_argument("--n_genes", type=int, default=None)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    S_all_file = args.S_all_file
    A_all_file = args.A_all_file
    S_file = args.S_file
    A_file = args.A_file
    stats_file = args.stats_file
    n_components = args.n_components
    robust_runs = args.robust_runs
    algorithm = args.algorithm
    n_genes = args.n_genes
    
    print('Loading data...')
    S_all, A_all = load_data(S_all_file, A_all_file, n_genes)
    
    print('Clustering ICA runs...')
    S, A, stats = cluster_ica_runs(
        S_all, A_all, n_components, robust_runs, algorithm
    )
    
    # save
    S.to_csv(S_file, **SAVE_PARAMS)
    A.to_csv(A_file, **SAVE_PARAMS)
    stats.to_csv(stats_file, **SAVE_PARAMS) 
    
    
##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")