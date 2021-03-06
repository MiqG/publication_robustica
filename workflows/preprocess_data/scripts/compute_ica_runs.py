#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# compute ICA iterations.

import argparse
import pandas as pd
import numpy as np
from robustica import RobustICA

# variables
SEED = 1234

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_robustica'
PREP_DIR = os.path.join(ROOT,'data','prep')
input_file = os.path.join(PREP_DIR,'TCGA','genexpr','BRCA.tsv.gz')
iterations = 3
transpose = False
n_jobs = 10
"""

##### FUNCTIONS #####
def process_inputs(input_file, transpose):
    # load
    X = pd.read_table(input_file, index_col=0)

    # drop genes with missing values
    X = X.loc[~X.isnull().any(axis=1)]
    
    # drop genes with no variation
    X = X.loc[X.std(1) > 0]

    # transpose
    if transpose:
        X = X.T

    # normalize
    X = (X - X.mean(1).values.reshape(-1, 1)) / X.std(1).values.reshape(-1, 1)
    return X


def compute_ica_iterations(X, n_components, iterations, n_jobs):
    print("Selected %s components" % n_components)

    # iterate
    rica = RobustICA(n_components=n_components, robust_runs=iterations, n_jobs=n_jobs)
    S_all, A_all, time = rica._iterate_ica(X.values)
    S_all = pd.DataFrame(S_all, index=X.index)
    A_all = pd.DataFrame(A_all, index=X.columns)
    return S_all, A_all


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str)
    parser.add_argument("--S_file", type=str)
    parser.add_argument("--A_file", type=str)
    parser.add_argument("--n_components", type=int)
    parser.add_argument("--iterations", type=int)
    parser.add_argument("--transpose", type=bool, default=False)
    parser.add_argument("--n_jobs", type=int)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    input_file = args.input_file
    S_file = args.S_file
    A_file = args.A_file
    n_components = args.n_components
    iterations = args.iterations
    transpose = args.transpose
    n_jobs = args.n_jobs

    print("Loading data...")
    X = process_inputs(input_file, transpose)

    print("Computing ICA...")
    S_all, A_all = compute_ica_iterations(X, n_components, iterations, n_jobs)

    print("Saving...")
    S_all.to_pickle(S_file)
    A_all.to_pickle(A_file)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
