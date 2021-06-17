#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Create a sample dataset to evaluate the performance of different algorithms
# to carry out robust ICA.

import argparse
import pandas as pd
import numpy as np

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}

seed = 1234
##### FUNCTIONS #####
def create_sample_data(n_cols, n_rows):
    data = pd.DataFrame(
        np.vstack([np.random.normal(size=n_cols) for i in range(n_rows)]), 
        index=['row%s'%r for r in range(n_rows)],
        columns=['col%s'%c for c in range(n_cols)]
    )
    
    return data


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_cols", type=int)
    parser.add_argument("--n_rows", type=int)
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    n_cols = args.n_cols
    n_rows = args.n_rows
    output_file = args.output_file
    
    data = create_sample_data(n_cols, n_rows)
    
    data.reset_index().to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")