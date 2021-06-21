#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Create a sample dataset to evaluate the performance of different algorithms
# to carry out robust ICA.

import os
import argparse
import pandas as pd
import numpy as np
from sklearn.decomposition import FastICA

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}
SEED = 2021

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_robustica'
PREP_DIR = os.path.join(ROOT,'data','prep')
input_file = os.path.join(PREP_DIR,'genexpr','Sastry2019.tsv.gz')
n_components = 200
"""

##### FUNCTIONS #####
def load_data(input_file):
    return pd.read_table(input_file, index_col=0)


def create_validation_data(data, n_components):
    """
    Compute one ground truth ICA and re-calculate the data so we know 
    these components are the ground truth of the new dataset.
    """
    # compute ICA
    ica = FastICA(n_components=n_components, random_state=SEED)
    
    # recover ground truth data set
    S = ica.fit_transform(data.values)
    A = ica.mixing_
    X = S.dot(A.T) + ica.mean_
    
    # prepare outputs
    X = pd.DataFrame(
        X, 
        index=['row%s'%r for r in range(X.shape[0])],
        columns=['col%s'%c for c in range(X.shape[1])]
    )
    
    S = pd.DataFrame(
        S, 
        index=['row%s'%r for r in range(S.shape[0])],
        columns=['col%s'%c for c in range(S.shape[1])]
    )
    
    A = pd.DataFrame(
        A, 
        index=['row%s'%r for r in range(A.shape[0])],
        columns=['col%s'%c for c in range(A.shape[1])]
    )
    
    
    return X, S, A


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str)
    parser.add_argument("--n_components", type=int)
    parser.add_argument("--output_dir", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    input_file = args.input_file
    n_components = args.n_components
    output_dir = args.output_dir
    
    data = load_data(input_file)
    X, S, A = create_validation_data(data, n_components)
    
    X.reset_index().to_csv(os.path.join(output_dir,'X.tsv.gz'), **SAVE_PARAMS)
    S.reset_index().to_csv(os.path.join(output_dir,'S.tsv.gz'), **SAVE_PARAMS)
    A.reset_index().to_csv(os.path.join(output_dir,'A.tsv.gz'), **SAVE_PARAMS)

    
##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")