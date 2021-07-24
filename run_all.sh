# 
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Run full pipeline.

set -o errexit
set -o nounset

NCORES=$1

echo "Downloading data..."
snakemake -s workflows/download_data/snakefile --cores=$NCORES
echo "Finished downloading data."


echo "Preprocessing data..."
snakemake -s workflows/preprocess_data/snakefile --cores=$NCORES
echo "Finished preprocessing data."


echo "Comparing clustering methods for robust ICA..."
snakemake -s workflows/benchmark_clustering/snakefile --cores=$NCORES
echo "Finished comparing clustering methods for robust ICA."


echo "Evaluating sign inference for robust ICA..."
snakemake -s workflows/benchmark_sign_inference/snakefile --cores=$NCORES
echo "Finished evaluating sign inference for robust ICA."


echo "Analyzing case study..."
snakemake -s workflows/case_study/snakefile --cores=$NCORES
echo "Finished analyzing case study."

echo "Making intro figures..."
snakemake -s workflows/case_study/snakefile --cores=$NCORES
echo "Finished making intro figures."

echo "Done!"