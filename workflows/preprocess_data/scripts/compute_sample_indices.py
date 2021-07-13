#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Compute sample indices

import argparse
import pandas as pd

# variables
MITOTIC_INDEX_SYMBOLS = [
    "CDKN3",
    "ILF2",
    "KDELR2",
    "RFC4",
    "TOP2A",
    "MCM3",
    "KPNA2",
    "CKS2",
    "CDK1",
]  # CDC2 was the previous name for CDK1
SAVE_PARAMS = {"sep": "\t", "index": False}

"""
Development
-----------
import os
ROOT = '~/projects/publication_robustica'
PREP_DIR = os.path.join(ROOT,'data','prep')
genexpr_file = os.path.join(PREP_DIR,'genexpr','LGG.tsv.gz')
"""

##### FUNCTION #####
def compute_mitotic_index(genexpr):
    mitotic_index = genexpr.loc[MITOTIC_INDEX_SYMBOLS].mean()
    mitotic_index.name = "mitotic_index"
    return mitotic_index


def compute_sample_indices(genexpr):
    # indices
    mitotic_index = compute_mitotic_index(genexpr)

    # prepare output
    sample_indices = pd.concat([mitotic_index], axis=1)

    sample_indices.index.name = "sampleID"

    return sample_indices


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genexpr_file", type=str)
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    genexpr_file = args.genexpr_file
    output_file = args.output_file

    genexpr = pd.read_table(genexpr_file, index_col=0)
    sample_indices = compute_sample_indices(genexpr)

    # save
    sample_indices.reset_index().to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
