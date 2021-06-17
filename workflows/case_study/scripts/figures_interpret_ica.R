#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Interpret the independent components.
#
# Outline
# -------
# - Which components differentiate IDH1 mutant from WT samples?
#     - barplot of mut vs wt sample counts
#     - wilcoxon test + FDR, select components of interest
#     - volcano plot of statistical test
#     - heatmap of A all components
#     - heatmap of A components of interest
# - Which biological processes are enriched in those components?
# - How certain are we of these components?



require(tidyverse)
require(ggpubr)
require(clusterProfiler)
require(pheatmap)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variable
GENE_OI = 'IDH1'
EFFECT_OI = c('Missense_Mutation','Nonsense_Mutation')

# Development
# -----------
PREP_DIR = file.path(ROOT,'data','prep')
RESULTS_DIR = file.path(ROOT,'results','compare_clustering_methods')
S_file = file.path(RESULTS_DIR,'files','cluster_iterations-LGG-DBSCAN','S_mean.tsv.gz')
A_file = file.path(RESULTS_DIR,'files','cluster_iterations-LGG-DBSCAN','A_mean.tsv.gz')
genexpr_file = file.path(PREP_DIR,'genexpr','LGG.tsv.gz')
snv_file = file.path(PREP_DIR,'snv','LGG.tsv.gz')
metadata_file = file.path(PREP_DIR,'metadata','LGG.tsv')
figs_dir = file.path(RESULTS_DIR,'figures','LGG-DBSCAN')


##### FUNCTIONS #####
plot_mut_components = function(S, A, metadata){
    
    # which components have different weights in samples with and without mutation?
    X = A %>% column_to_rownames('index')
    samples = metadata[['sampleID']]
    is_mut = metadata[['is_mut']]
    
    result = apply(X[samples,], 2, function(x){
        a = x[is_mut]
        b = x[!is_mut]

        result = data.frame(
            a_median = median(a), 
            b_median = median(b),
            med_diff = median(a)-median(b),
            pvalue = wilcox.test(a, b)$p.value
        )    
        
        return(result)
    })
    result = do.call(rbind,result)
    result[['fdr']] = p.adjust(result[['pvalue']], method='fdr')
    
    plts = list()
    
}


make_plots = function(clustering_stats){
    plts = list(
    )
    plts = do.call(c,plts)
    return(plts)
}


save_plot = function(plt, plt_name, extension='.pdf', 
                      directory='', dpi=350, 
                      width = par("din")[1], height = par("din")[2], units='cm'){
        filename = file.path(directory,paste0(plt_name,extension))
        ggsave(filename, 
               plt, 
               width=width, height=height, dpi=dpi, limitsize=FALSE, units=units)
}


save_plots = function(plts, figs_dir){
    lapply(names(plts), function(plt_name){
        save_plot(plts[[plt_name]] + 
                  theme_pubr(base_size=10, x.text.angle=70),
                  plt_name, '.png', figs_dir, width=12, height=12)
    })
}


main = function(){
    args = getParsedArgs()
    S_file = args$S_file
    A_file = args$A_file
    genexpr_file = args$genexpr_file
    snv_file = args$snv_file
    metadata_file = args$metadata_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load data
    S = read_tsv(S_file)
    A = read_tsv(A_file)
    genexpr = read_tsv(genexpr_file)
    snv = read_tsv(snv_file)
    metadata = read_tsv(metadata_file)
    
    # drop unwanted components: -1, silhouette
    S = S %>% dplyr::select(-one_of('-1'))
    A = A %>% dplyr::select(-one_of('-1'))
    
    # add IDH1 mutation in metadata
    mut_samples = snv %>% 
        filter(gene==GENE_OI & effect%in%EFFECT_OI) %>% 
        pull(sampleID)
    metadata = metadata %>% mutate(is_mut = sampleID %in% mut_samples)
    
    
    plts = make_plots(data[['clustering_stats']])
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
