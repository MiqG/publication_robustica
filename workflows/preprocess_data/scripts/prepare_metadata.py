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
FEATURES_OI = {
    'tcga':{
        'sample': 'sampleID',
        '_PATIENT': 'subjectID',
        'sample_type': 'sample_type',
        '_primary_disease': 'primary_disease',
        'cancer type abbreviation': 'cancer_type',
        'age_at_initial_pathologic_diagnosis': 'subject_age',
        'gender': 'subject_sex'
    } 
}

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_robustica'
DATA_DIR = os.path.join(ROOT,'data')
RAW_DIR = os.path.join(DATA_DIR,'raw')
raw_tcga_clinical = os.path.join(RAW_DIR,'UCSCXena','TCGA','phenotype','TCGA_phenotype_denseDataOnlyDownload.tsv.gz')
raw_tcga_survival = os.path.join(RAW_DIR,'UCSCXena','TCGA','phenotype','Survival_SupplementalTable_S1_20171025_xena_sp.gz')
"""

##### FUNCTIONS #####
def preprocess_tcga(raw_tcga_clinical, raw_tcga_survival):
    # load and preprocess
    clinical = pd.read_table(raw_tcga_clinical)
    survival = pd.read_table(raw_tcga_survival)
    
    # combine
    df = pd.merge(clinical, survival, on='sample', how='outer')
    
    # rename
    df = df.rename(columns=FEATURES_OI['tcga'])
    
    # fill up cancer types
    cancer_types = df[['primary_disease','cancer_type']]\
                        .drop_duplicates()\
                        .sort_values('primary_disease')\
                        .dropna()
    cancer_types = {row['primary_disease']: row['cancer_type'] 
                    for idx, row in cancer_types.iterrows()}
    df['cancer_type'] = [cancer_types[primary_disease] 
                          for primary_disease in df['primary_disease']]
    
    df['sample_type'] = df['sample_type'].str.replace('Primary Blood Derived Cancer - Peripheral Blood', 'PBDC-PB')
    return df


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw_tcga_clinical", type=str)
    parser.add_argument("--raw_tcga_survival", type=str)
    parser.add_argument("--output_dir", type=str)

    args = parser.parse_args()
    
    return args


def main():
    args = parse_args()
    raw_tcga_clinical = args.raw_tcga_clinical
    raw_tcga_survival = args.raw_tcga_survival
    output_dir = args.output_dir
    
    # load
    tcga = preprocess_tcga(raw_tcga_clinical,raw_tcga_survival)

    # save
    for cancer_type in tcga['cancer_type'].unique():
        df = tcga.loc[tcga['cancer_type']==cancer_type]
        filename = os.path.join(output_dir,'%s.tsv') % cancer_type
        df.to_csv(filename, sep="\t", index=False)
    
    filename = os.path.join(output_dir,'%s.tsv') % 'PANCAN'
    tcga.to_csv(filename, sep="\t", index=False)
    

##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
