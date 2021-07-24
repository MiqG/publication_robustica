# `robustica`: customizable robust independent component analysis

This is the repository to reproduce the analysis carried out in our manuscript, from the data download to the figures. 
The article describes how we developed [`robustica`](), a Python package to carry out robust independent component analysis.

## Overview
<p align="center">
  <img src="images/fig-main.png" width="">
</p>

**Development and implementation of `robustica` to carry out robust Independent Component Analysis (ICA).** (**A**)`robustica` enables fully customized robust ICA. We built `robustica` following `scikit-learn`’s programming conventions to enable full control of both the iterative and clustering steps of robust ICA facilitating customization and optimization. In particular, `robustica` includes a subroutine to infer and correct the signs of components across ICA runs that improves the precision and efficiency of the clustering step by enabling us to use Euclidean distance metrics and to compress the feature space. (**B**) Comparison of clustering algorithms for robust ICA using [Sastry (2019)](https://doi.org/10.1038/s41467-019-13483-w)‘s dataset. Maximum memory usage and time for each clustering algorithm to cluster 100 ICA runs with 100 components each. Dot sizes indicate silhouette scores (the larger the better). In boldface, the best clustering algorithm ([*AgglomerativeClustering*](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html)). (**C**) Development steps to improve the precision (Weight Std.) and efficiency while reducing the bias of robust ICA through our sign inference-and-correction subroutine combined with PCA and Euclidean distances, using [Sastry (2019)](https://doi.org/10.1038/s41467-019-13483-w)‘s dataset. (**D**) Case study workflow for robust ICA. We dissected >500 tumor samples from [LGG patients](https://xenabrowser.net/datapages/?dataset=EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) with `robustica` into 100 robust independent components. Component 67 was simultaneously associated with multiple sample features (*IDH1* mutation status, mitotic index, and overall survival) and contained genes known to be mechanistically associated with mutations in *IDH1* that modulate tumor aggressiveness.

## Structure
- `images`: images for the repository.
- `src`: contains scripts with helper functions used throughout the project.
- `workflows`: contains the `snakemake` workflows and scripts to reproduce the whole analysis
    - `download_data`: outputs the datasets from [Sastry (2019)](https://doi.org/10.1038/s41467-019-13483-w) and [The Cancer Genome Atlas](https://xenabrowser.net/datapages/?dataset=EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) into `data/raw`.
    - `preprocess_data`: outputs clean data into `data/prep`.
    - `benchmark_clustering`: outputs the results of comparing different clustering algorithms to perform robust ICA into `results/benchmark_clustering`
    - `benchmark_sign_inference`: outputs the results of introducing our subroutine to infer-and-correct the sign of components across ICA runs into `results/benchmark_sign_inference`
    - `case_study`: outputs the results of dissecting LGG tumor samples through robust ICA into `results/case_study`
    - `intro_figures`: generates pretty heatmaps used to make figures for the manuscript. Outputs into `results/intro_figures`
    
After running all workflows (indicated below), the figures used in the manuscript and their accompaining data will be in the subdirectory `figures` of each workflow outputting into `results`.

## Requirements
- Python
    - snakemake
    - rpy2
    - pandas
    - numpy
    - matplotlib
    - seaborn
    - scikit-learn
    - scikit-learn-extra
    - memory_profiler
- R
    - tidyverse
    - ggpubr
    - pheatmap
    - reshape2
    - ggrepel
    - writexl
    - proxy
    - frdtool
    - latex2exp
    - ggnewscale
    - enrichplot
    - clusterProfiler
    - org.Hs.eg.db
    - survival
    - extrafont
    - sruvminer

## Usage
This will generate `data` and `results`.

```shell
bash run_all.sh 4 # number of cores available
```
