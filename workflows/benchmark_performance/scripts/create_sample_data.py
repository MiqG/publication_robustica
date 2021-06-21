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

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}
SEED = 2021

"""
Development
-----------
n_samples = 500
n_features = 20000
"""

##### FUNCTIONS #####
def create_sample_data(n_samples, n_features):
    """
    Based on https://scikit-learn.org/stable/auto_examples/decomposition/plot_ica_blind_source_separation.html#sphx-glr-auto-examples-decomposition-plot-ica-blind-source-separation-py
    """
    
    np.random.seed(SEED)
    
    time = np.linspace(0, 10, n_features)
    
    s1 = np.sin(2 * time)  # Signal 1 : sinusoidal signal
    s2 = np.sign(np.sin(3 * time))  # Signal 2 : square signal
    s3 = np.cos(2 * time)  # Signal 3: saw tooth signal
    s4 = np.random.beta(0.5, 0.25, n_features)
    s5 = np.random.standard_t(1.5, n_features)
    s6 = np.random.uniform(-1, 1, n_features)
    s7 = np.random.chisquare(2, n_features)
    
    S = np.vstack([s1, s2, s3, s4, s5, s6, s7]).T # features x samples
    n_components = S.shape[1]
    
    # add noise
    S = S + np.random.normal(size=S.shape)
    
    # standardize every feature
    S = S / S.std(axis=1).reshape(-1,1)
    
    # mix
    A = np.random.choice(np.linspace(0,10,20), size=(n_samples, n_components))
    X = S.dot(A.T)
    
    # prepare outputs
    X = pd.DataFrame(
        X, 
        index=['row%s'%r for r in range(n_features)],
        columns=['col%s'%c for c in range(n_samples)]
    )
    
    S = pd.DataFrame(
        S, 
        index=['row%s'%r for r in range(n_features)],
        columns=['col%s'%c for c in range(n_components)]
    )
    
    A = pd.DataFrame(
        A, 
        index=['row%s'%r for r in range(n_samples)],
        columns=['col%s'%c for c in range(n_components)]
    )
    
    
    return X, S, A


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_samples", type=int)
    parser.add_argument("--n_features", type=int)
    parser.add_argument("--output_dir", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    n_samples = args.n_samples
    n_features = args.n_features
    output_dir = args.output_dir
    
    X, S, A = create_sample_data(n_samples, n_features)
    
    os.makedirs(output_dir)
    X.reset_index().to_csv(os.path.join(output_dir,'X.tsv.gz'), **SAVE_PARAMS)
    S.reset_index().to_csv(os.path.join(output_dir,'S.tsv.gz'), **SAVE_PARAMS)
    A.reset_index().to_csv(os.path.join(output_dir,'A.tsv.gz'), **SAVE_PARAMS)

    
##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")