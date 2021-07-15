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

##### FUNCTIONS #####
def cluster_ica_runs(S_all, A_all, n_components, robust_runs):
    # cluster
    rica = robustica.RobustICA(
        n_components = n_components, 
        robust_runs = robust_runs,
        robust_kws = {"affinity": "euclidean", "linkage": "average"},
        robust_infer_signs = True,
        robust_dimreduce = True
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
        
    # load data
    S_all = pd.read_pickle(S_all_file)
    A_all = pd.read_pickle(A_all_file)

    S, A, stats = cluster_ica_runs(S_all, A_all, n_components, robust_runs)
    
    # save
    S.to_csv(S_file, **SAVE_PARAMS)
    A.to_csv(A_file, **SAVE_PARAMS)
    stats.to_csv(stats_file, **SAVE_PARAMS) 
    
    
##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")