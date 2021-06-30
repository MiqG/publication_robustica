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
require(latex2exp)
require(ggrepel)

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
clustering_summary_file = file.path(RESULTS_DIR,'files','cluster_iterations-LGG-DBSCAN','clustering_summary.tsv')
genexpr_file = file.path(PREP_DIR,'genexpr','LGG.tsv.gz')
snv_file = file.path(PREP_DIR,'snv','LGG.tsv.gz')
metadata_file = file.path(PREP_DIR,'metadata','LGG.tsv')
figs_dir = file.path(RESULTS_DIR,'figures','LGG-DBSCAN')


##### FUNCTIONS #####
make_heatmap = function(X, metadata){
    # get parameters
    metadata_sample_col = 'sampleID'
    metadata_cols_oi = c('is_mut')
    palettes_cols_oi = 'lancet'
    
    # prepare matrix
    mat = X
    
    # prepare annotation
    annot = metadata %>% 
        column_to_rownames(metadata_sample_col) %>%
        dplyr::select(all_of(metadata_cols_oi))
    annot = annot[(rownames(annot) %in% colnames(X)),,drop=FALSE]
    
    # prepare annotation colors
    colors = sapply(1:length(metadata_cols_oi), simplify = FALSE, 
                    function(i){
                        x = unique(annot[,metadata_cols_oi[i]])
                        palette = get_palette(palettes_cols_oi[i], length(x))
                        color = setNames(palette, x)
                        return(color)
                    })
    names(colors) = metadata_cols_oi
    
    # plot
    plt = pheatmap(t(mat), 
                     show_rownames = FALSE, 
                     annotation_row = annot, 
                     annotation_colors = colors, 
                     color = rev(get_palette('RdYlBu',25)),
                     silent = TRUE)
    return(plt)
}


plot_mut_components = function(S, A, metadata, clustering_summary){
    plts = list()
    
    X = A %>% column_to_rownames('index') %>% t()
    comp_sil = clustering_summary %>% 
        mutate(label=as.character(label)) %>%
        group_by(label) %>%
        summarise(mean_silhouette = mean(silhouette),
                  n=n()) %>%
        dplyr::rename(component=label) 
    
    # which components have different weights in samples with and without mutation?
    samples = metadata[['sampleID']]
    is_mut = metadata[['is_mut']]=='TRUE'
    
    result = apply(X[,samples], 1, function(x){
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
    result = do.call(rbind,result) %>% 
        rownames_to_column('component') %>%
        mutate(fdr = p.adjust(pvalue, method='fdr'),
               log10_fdr = -log10(fdr)) %>%
        left_join(comp_sil, by='component')
    
    components_oi = result %>% 
        filter(fdr<0.01 & abs(med_diff)>5) %>% 
        pull(component)
    
    
    # mut vs wt counts
    plts[['counts_mutated']] = metadata %>%
        count(is_mut) %>%
        ggbarplot(x='is_mut', y='n', label=TRUE, fill='is_mut', color=NA)
    
    
    # overview sample weights and mutation status
    plts[['A-overview_heatmap']] = make_heatmap(X, metadata)
    
    
    plts[['A-differential_weights_mutated-volcano']] = result %>%
        ggscatter(x='med_diff', y='log10_fdr', color='mean_silhouette') +
        geom_hline(yintercept = -log10(0.01), linetype='dashed') +
        gradient_color(c("grey", "red")) +
        #geom_vline(xintercept = c(-5,5), linetype='dashed') + 
        geom_text_repel(
            aes(label=component), 
            result %>% filter(component %in% components_oi)
        ) +
        labs(x=TeX('$\\Delta W_{mut-wt}$'), 
             y=TeX('$-log_{10}(FDR)$'),
             color='Mean Silhouette')
    
    # overview components of interest
    plts[['A-components_oi-heatmap']] = make_heatmap(X[components_oi,], metadata)
    
    corr = cor(t(X[components_oi,]))
    plts[['A-components_oi_corr-heatmap']] = pheatmap(corr, cutree_cols = 2, cutree_rows = 2, silent=TRUE)
    
    plts[['S-components_oi_corr-heatmap']] = pheatmap(cor(S[,components_oi]), cutree_cols = 2, cutree_rows = 2, silent=TRUE)
    
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
    clustering_summary_file = args$clustering_summary_file
    genexpr_file = args$genexpr_file
    snv_file = args$snv_file
    metadata_file = args$metadata_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load data
    S = read_tsv(S_file)
    A = read_tsv(A_file)
    clustering_summary = read_tsv(clustering_summary_file)
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
    metadata = metadata %>% mutate(is_mut = as.character(sampleID %in% mut_samples))
    
    
    plts = make_plots(data[['clustering_stats']])
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
