#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# compute ICA iterations and keep:
#  - clustering summary + params: combined parameters label
#  - run summary:
#     - ICA
#       - score at convergence
#       - iteration at convergence
#       - time at convergence
#       - max_iter
#       - tol
#       - params: combined parameters label
#     - RobustICA
#       - time to run all iterations
#       - n_jobs
#       - iteration

import argparse
import pandas as pd
import numpy as np
from robustica import RobustICA, InferComponents
from joblib import parallel_backend
import os
import time

# variables
MAX_VARIANCE_EXPLAINED = 0.9
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_robustica'
PREP_DIR = os.path.join(ROOT,'data','prep')
input_file = os.path.join(PREP_DIR,'TCGA','genexpr','LGG.tsv.gz')
iterations = 3
transpose = False
n_jobs = 10
n_components = 10
max_iter = int(1e3)
tol = 1e-2
"""

##### FUNCTIONS #####
def process_inputs(input_file):
    # load
    X = pd.read_table(input_file, index_col=0)

    # drop genes with no variation
    X = X.loc[X.std(1) > 0]

    # normalize
    X = (X - X.mean(1).values.reshape(-1, 1)) / X.std(1).values.reshape(-1, 1)
    return X


def evaluate_ica_iterations(X, n_components, iterations, n_jobs, max_iter, tol):
    # infer components
    if n_components is None:
        n_components = InferComponents(MAX_VARIANCE_EXPLAINED).fit_predict(X)
    print("Selected %s components" % n_components)

    # run robust ICA
    rica = RobustICA(
        n_components=n_components, max_iter=max_iter, tol=tol, robust_iter=iterations
    )
    start_time = time.time()
    with parallel_backend("loky", n_jobs=n_jobs):
        S, A = rica.fit_transform(X.values)
    seconds = time.time() - start_time

    # prepare output
    S = pd.DataFrame(S, index=X.index)
    A = pd.DataFrame(A, index=X.columns)
    S_std = pd.DataFrame(rica.S_std, index=X.index)
    A_std = pd.DataFrame(rica.A_std, index=X.columns)
    clustering_stats = rica.clustering_stats_
    summary = rica.prepare_summary()
    summary["time_robustica"] = seconds
    summary["n_jobs"] = n_jobs
    summary["n_components"] = n_components

    # add params label
    summary["params"] = "max_iter=%s & tol=%s" % (max_iter, tol)
    clustering_stats["params"] = "max_iter=%s & tol=%s" % (max_iter, tol)

    # keep only last row of each RobustICA iteration
    idx = (
        summary.groupby(["iteration_robustica"])["iteration_ica"].transform(max)
        == summary["iteration_ica"]
    )
    summary = summary.loc[idx].copy()

    return S, A, S_std, A_std, summary, clustering_stats


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str)
    parser.add_argument("--output_dir", type=str)
    parser.add_argument("--iterations", type=int)
    parser.add_argument("--n_jobs", type=int)
    parser.add_argument("--max_iter", type=float)
    parser.add_argument("--tol", type=float)
    parser.add_argument("--n_components", type=int, default=None)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    input_file = args.input_file
    output_dir = args.output_dir
    iterations = args.iterations
    n_jobs = args.n_jobs
    max_iter = int(args.max_iter)
    tol = args.tol
    n_components = args.n_components

    os.makedirs(output_dir)

    print("Loading data...")
    X = process_inputs(input_file)

    print("Computing ICA...")
    S, A, S_std, A_std, summary, clustering_stats = evaluate_ica_iterations(
        X, n_components, iterations, n_jobs, max_iter, tol
    )

    print("Saving...")
    S.reset_index().to_csv(os.path.join(output_dir, "S.tsv.gz"), **SAVE_PARAMS)
    A.reset_index().to_csv(os.path.join(output_dir, "A.tsv.gz"), **SAVE_PARAMS)
    S_std.reset_index().to_csv(os.path.join(output_dir, "S_std.tsv.gz"), **SAVE_PARAMS)
    A_std.reset_index().to_csv(os.path.join(output_dir, "A_std.tsv.gz"), **SAVE_PARAMS)
    summary.to_csv(os.path.join(output_dir, "summary.tsv.gz"), **SAVE_PARAMS)
    clustering_stats.to_csv(
        os.path.join(output_dir, "clustering_stats.tsv.gz"), **SAVE_PARAMS
    )


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
