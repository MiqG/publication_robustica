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
"""

##### FUNCTIONS ####
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
    if distance_func=="original":
        D = 1 - abs(np.dot(S_all.T,S_all))
        D[D < 0.5] = 0
    elif distance_func=="full":
        D = 1 - abs(np.dot(S_all.T,S_all))
    elif distance_func=="robustica":
        D = 1 - abs(np.corrcoef(S_all.T))
    else:
        print("Error: Set a valid 'distance_func'.")
        sys.exit(1)
    
    D = np.clip(D, 0, 1)

    return D


def cluster_components(D, eps, min_samples, S_all, A_all):
    
    # Run DBSCAN
    dbscan = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed')
    labels = dbscan.fit_predict(D)
    n_clusters = max(labels)+1
    
    # Place clustered components into correct bins
    if rank == 0:
        print('Loading individual S and A matrices')
        start = 0
        end = 0
        S_bins = {i:[] for i in range(n_clusters)}
        A_bins = {i:[] for i in range(n_clusters)}

        for i in range(nWorkers):
            # Get labels for each partial matrix
            end += block_size[i]
            proc_labels = labels[start:end]
            start = end

            # Add parts of matrix to full table
            S_part = pd.read_csv(os.path.join(tmp_dir,'proc_{}_S.csv'.format(i)),index_col=0)
            A_part = pd.read_csv(os.path.join(tmp_dir,'proc_{}_A.csv'.format(i)),index_col=0)
            S_part.columns = range(S_part.shape[1])
            A_part.columns = range(A_part.shape[1])
            for i,lab in enumerate(proc_labels):
                if lab != -1:
                    S_bins[lab].append(S_part[i])
                    A_bins[lab].append(A_part[i])
    
    # Gather final S and A matrices
    if rank == 0:

        t = timeit(t)
        print('\nGathering final S and A matrices')

        S_final = pd.DataFrame(columns=range(n_clusters),index=S_part.index)
        A_final = pd.DataFrame(columns=range(n_clusters),index=A_part.index)
        df_stats = pd.DataFrame(columns=['S_mean_std','A_mean_std','count'],index=range(n_clusters))

        for lab in range(n_clusters):
            S_clust = S_bins[lab]
            A_clust = A_bins[lab]

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
            for j in range(1,len(S_clust)):
                Svec = S_clust[j]
                Avec = A_clust[j]
                if stats.pearsonr(Svec,Svec0)[0] > 0:
                    S_single.append(Svec)
                    A_single.append(Avec)
                else:
                    S_single.append(-Svec)
                    A_single.append(-Avec)

            # Add centroid of cluster to final S matrix        
            S_final[lab] = np.array(S_single).T.mean(axis=1)
            A_final[lab] = np.array(A_single).T.mean(axis=1)

            # Get component stats
            df_stats.loc[lab,'S_mean_std'] = np.array(S_single).std(axis=1).mean()
            df_stats.loc[lab,'A_mean_std'] = np.array(A_single).std(axis=1).mean()
            df_stats.loc[lab,'count'] = len(S_single)

        print('\nFinal components created')
        t = timeit(t)    
        
        
    return S_final, A_final, df_stats


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--S_all_file", type=str)
    parser.add_argument("--A_all_file", type=str)
    parser.add_argument("--output_file", type=str)
    parser.add_argument("--distance_func", type=str)
    args = parser.parse_args()

    return args

def main():
    # unpack arguments
    args = parse_args()
    S_all_file = args.S_all_file
    A_all_file = args.A_all_file
    output_file = args.output_file
    n_iters = args.n_iters
    
    # load data
    S_all, A_all = load_data(S_all_file, A_all_file)
    
    # cluster
    D = compute_distances(S_all, distance_func)
    S, A, df_stats = cluster_components(D, S_all, A_all)
    

##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")