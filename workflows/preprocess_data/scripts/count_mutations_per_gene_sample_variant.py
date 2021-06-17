#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# Count mutations per sample by effect. 
#

import pandas as pd
import argparse

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_robustica'
RAW_DIR = os.path.join(ROOT,'data','raw')
snv_file = os.path.join(RAW_DIR,'UCSCXena','TCGA','snv','mc3.v0.2.8.PUBLIC.xena.gz')
metadata_file = os.path.join(ROOT,'data','prep','metadata','tcga.tsv')
"""

##### FUNCTIONS #####
def load_data(snv_file, metadata_file):
    snv = pd.read_table(snv_file, low_memory=False)
    snv = snv.rename(columns={'sample':'sampleID'})
    metadata = pd.read_table(metadata_file)
    
    # add cancer types
    snv = pd.merge(snv, metadata[['sampleID','cancer_type']], how='left', on='sampleID')
    
    return snv


def count_mutations_per_gene_and_effect(snv):
    # mutations per gene and effect
    # if a sample has multiple mutations in 
    # the same gene it is counted only once
    counts = snv[['cancer_type','sampleID','gene','effect']]\
                .drop_duplicates()\
                .groupby(['cancer_type','sampleID','gene','effect'])\
                .size()\
                .reset_index()\
                .rename(columns={0:'n'})
    
    return counts


def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--snv_file',type=str)
    parser.add_argument('--metadata_file',type=str)
    parser.add_argument('--output_file',type=str)
    
    args = parser.parse_args()
    return args


def main():
    # parse arguments
    args = parse_args()
    
    snv_file = args.snv_file
    metadata_file = args.metadata_file
    output_file = args.output_file
    
    # load
    snv = load_data(snv_file, metadata_file)
    result = count_mutations_per_gene_and_effect(snv)
    
    # save
    result.to_csv(output_file, sep='\t', index=False, compression='gzip')


##### SCRIPT #####
if __name__ == '__main__':
    main()
    print('Done!')

