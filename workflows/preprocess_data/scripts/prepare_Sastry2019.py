#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Preprocess data from Sastry et al. 2019 (DOI: https://doi.org/10.1038/s41467-019-13483-w) as in the paper.
#
# Outline
# -------
# 1. load metadata and gene expression
# 2. for the gene expression, substract the mean from the control conditions
# 3. Save final S and A matrices

import os
import argparse
import pandas as pd

# variables
SAVE_PARAMS = {"sep": "\t", "index": True}
CONTROLS = ['control__wt_glc__1', 'control__wt_glc__2']

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_robustica'
DATA_DIR = os.path.join(ROOT,'data')
RAW_DIR = os.path.join(DATA_DIR,'raw')
suptab_file = os.path.join(RAW_DIR,'articles','Sastry2019','ica_analysis.xlsx')
"""

##### FUNCTIONS #####
def load_data(suptab_file):
    data = pd.ExcelFile(suptab_file)
    genexpr = data.parse('Expression Data')
    metadata = data.parse('Metadata')
    S = data.parse('S matrix')
    A = data.parse('A matrix')
    
    return genexpr, metadata, S, A


def preprocess_genexpr(genexpr):
    # index
    genexpr = genexpr.set_index(genexpr.columns[0])
    
    # compute logFC w.r.t. controls
    ctrl = genexpr[CONTROLS].mean(1).values.reshape(-1,1)
    genexpr = genexpr - ctrl
    
    return genexpr.reset_index()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--suptab_file", type=str)
    parser.add_argument("--genexpr_file", type=str)
    parser.add_argument("--metadata_file", type=str)
    parser.add_argument("--S_file", type=str)
    parser.add_argument("--A_file", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    suptab_file = args.suptab_file
    genexpr_file = args.genexpr_file
    metadata_file = args.metadata_file
    S_file = args.S_file
    A_file = args.A_file
    
    genexpr, metadata, S, A = load_data(suptab_file)
    
    genexpr = preprocess_genexpr(genexpr)
    S = S.loc[S.iloc[:,0] != 'Threshold'].rename(columns={'Unnamed: 0':'index'})
    A = A.rename(columns={'Unnamed: 0':'index'})
    
    genexpr.to_csv(genexpr_file, compression='gzip', **SAVE_PARAMS)
    metadata.to_csv(metadata_file, **SAVE_PARAMS)
    S.to_csv(S_file, compression='gzip', **SAVE_PARAMS)
    A.to_csv(A_file, compression='gzip', **SAVE_PARAMS)
    

##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")