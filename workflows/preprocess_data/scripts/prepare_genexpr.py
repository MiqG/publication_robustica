#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Preprocess raw metadata.
#

import os
import argparse
import pandas as pd

# variables
SAVE_PARAMS = {"sep": "\t", "index": True, "compression": "gzip"}

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_robustica'
DATA_DIR = os.path.join(ROOT,'data')
RAW_DIR = os.path.join(DATA_DIR,'raw')
PREP_DIR= os.path.join(DATA_DIR,'prep')

genexpr_file = os.path.join(RAW_DIR,'UCSCXena','TCGA','rnaseq','AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz')
metadata_file = os.path.join(DATA_DIR,'prep','TCGA','metadata','PANCAN.tsv')
sample_col = 'sampleID'
subset_col = 'cancer_type'
subset_values = 'BRCA'.split(',')

genexpr_file = os.path.join(RAW_DIR,'GTEx','v8','rnaseq','GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz')
metadata_file = os.path.join(PREP_DIR,'metadata','GTEx','muscle.tsv')
sample_col = 'sampleID'
subset_col = 'tissue_type'
subset_values = 'muscle'.split(",")
gene_col = "Description"
skiprows = 2
"""

##### FUNCTIONS #####
def get_file_header(file, skiprows=0, nrows=1):
    """returns the first line of a file"""
    return list(pd.read_table(file, nrows=nrows, skiprows=skiprows).columns)


def get_samples_oi(metadata_file, sample_col, subset_col, subset_values):
    metadata = pd.read_table(metadata_file)
    idx = metadata[subset_col].isin(subset_values)
    samples_oi = list(metadata.loc[idx, sample_col])

    return samples_oi


def load_data(genexpr_file, metadata_file, sample_col, subset_col, subset_values, gene_col, skiprows):
    # get samples
    samples_oi = get_samples_oi(metadata_file, sample_col, subset_col, subset_values)
    samples_avail = get_file_header(genexpr_file, skiprows)
    common_samples = list(set(samples_oi).intersection(samples_avail))

    # load
    genexpr = pd.read_table(
        genexpr_file, usecols=[gene_col] + common_samples, skiprows=skiprows
    ).set_index(gene_col)
    return genexpr


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genexpr_file", type=str)
    parser.add_argument("--gene_col", type=str, default="sample")
    parser.add_argument("--metadata_file", type=str)
    parser.add_argument("--sample_col", type=str)
    parser.add_argument("--subset_col", type=str)
    parser.add_argument("--skiprows", type=int, default=0)
    parser.add_argument("--subset_values", type=str)
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    genexpr_file = args.genexpr_file
    metadata_file = args.metadata_file
    sample_col = args.sample_col
    subset_col = args.subset_col
    subset_values = args.subset_values.split(",")
    output_file = args.output_file
    gene_col = args.gene_col
    skiprows = args.skiprows

    genexpr = load_data(
        genexpr_file, metadata_file, sample_col, subset_col, subset_values, gene_col, skiprows
    )

    genexpr.to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
