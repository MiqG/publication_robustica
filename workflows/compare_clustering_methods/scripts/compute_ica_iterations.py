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
from robustica import RobustICA, InferComponents

# variables
MAX_VARIANCE_EXPLAINED = 0.9
MAX_ITER = int(1e5)
TOL = 1e-4 # default

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
    
    # drop genes with no variation
    X = X.loc[X.std(1) > 0]
    
    # transpose
    if transpose: X = X.T
    
    # normalize
    X = (X - X.mean(1).values.reshape(-1,1)) / X.std(1).values.reshape(-1,1)
    return X
    
    
def compute_ica_iterations(X, iterations, n_jobs):
    # infer components
    ncomp = InferComponents(MAX_VARIANCE_EXPLAINED).fit_predict(X)
    print('Selected %s components' % ncomp)
    
    # iterate
    rica = RobustICA(n_components=ncomp, 
                     robust_iter=iterations, 
                     n_jobs=n_jobs,
                     max_iter=MAX_ITER, 
                     tol=TOL)
    rica._iterate_ica(X.values)
    S_all = pd.DataFrame(rica.S_all, index=X.index)
    A_all = pd.DataFrame(rica.A_all, index=X.columns)
    return S_all, A_all
    

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str)
    parser.add_argument("--S_file", type=str)
    parser.add_argument("--A_file", type=str)
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
    iterations = args.iterations
    transpose = args.transpose
    n_jobs = args.n_jobs
    
    print('Loading data...')    
    X = process_inputs(input_file, transpose)
    
    print('Computing ICA...')
    S_all, A_all = compute_ica_iterations(X, iterations, n_jobs)
    
    print('Saving...')
    S_all.to_pickle(S_file)
    A_all.to_pickle(A_file)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")