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
require(reshape2)
require(survival)
require(survminer)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variable
GENE_OI = 'IDH1'
EFFECT_OI = c('Missense_Mutation','Nonsense_Mutation')

# Development
# -----------
PREP_DIR = file.path(ROOT,'data','prep')
RESULTS_DIR = file.path(ROOT,'results','case_study')

S_file = file.path(RESULTS_DIR,'files','cluster_iterations','LGG','S.tsv.gz')
A_file = file.path(RESULTS_DIR,'files','cluster_iterations','LGG','A.tsv.gz')
stats_file = file.path(RESULTS_DIR,'files','cluster_iterations','LGG','stats.tsv.gz')

genexpr_file = file.path(PREP_DIR,'genexpr','LGG.tsv.gz')
snv_file = file.path(PREP_DIR,'snv','LGG.tsv.gz')
metadata_file = file.path(PREP_DIR,'metadata','LGG.tsv')
figs_dir = file.path(RESULTS_DIR,'figures','LGG')

sample_properties_file = '~/projects/sso_targets/results/exon_associations_ignoring_confounders/LGG/files/sample_properties.tsv'

##### FUNCTIONS #####
make_heatmap = function(X, metadata, cluster_rows=TRUE, cluster_cols=TRUE){
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
                     silent = TRUE, 
                     cluster_rows = cluster_rows,
                     cluster_cols = cluster_cols)
    return(plt)
}


plot_silhouettes = function(stats){
    X = stats

    plts=list()
    
    plts[['silhouettes-scatter']] = X %>% 
        ggscatter(x='silhouette_euclidean', y='silhouette_pearson') + 
        geom_text_repel(aes(label=cluster_id))
    
    return(plts)
}


plot_sample_properties = function(genexpr, S, A, sample_properties){
    
    # do any weights correlate with any sample property?
    X = A %>% column_to_rownames('index')
    properties = sample_properties %>% column_to_rownames('sampleID')
    properties_oi = setdiff(colnames(properties),c('OS','OS.time','sample_type','ABSOLUTE'))
    samples = rownames(X)

    # spearman correlations
    correls = sapply(properties_oi, function(property_oi){
        y = properties[samples,property_oi]
        correls = apply(X[samples,], 2, function(x){
            cor(x,y, method='spearman', use='pairwise.complete.obs')
        })
    }, simplify=FALSE)
    correls = do.call(rbind, correls)
    
    # survival associationss
    survs = apply(X[samples,], 2, function(x){
        dat = cbind(weight=x, properties[samples,c('OS','OS.time')])
        result = summary(coxph(Surv(OS.time, OS) ~ weight, dat))
        assoc = result$coefficients['weight','z']
        return(assoc)
    })
    survs = data.frame(coxph = survs)
    
    plts = list()
    
    plts[['sample_properties-coxph']] = pheatmap(survs, cluster_cols = FALSE, silent = TRUE)
    plts[['sample_properties-correls']] = pheatmap(t(correls), silent = TRUE)
    
    # focus on interesting components
    ## ESTIMATE
    component_oi = '83'
    ### scatter
    plts[['components_oi-estimate-scatter']] = data.frame(
        component = X[samples,component_oi], 
        index = properties[samples,'ESTIMATE']) %>% 
    ggscatter(x='component', y='index', alpha=0.5) + 
    labs(x=sprintf('Weights Component %s',component_oi), y='ESTIMATE') + 
    stat_cor(method='spearman', label.y = 0.6)
    ### heatmap of important genes
    genes_oi = S[,c('sample',component_oi)] %>% 
        filter(abs(get(component_oi)) > (4*sd(get(component_oi)) + mean(get(component_oi)))) %>% 
        pull(sample)
    plts[['components_oi-estimate-heatmap_genexpr']] = pheatmap(
        genexpr %>% filter(sample %in% genes_oi) %>% column_to_rownames('sample'), 
        main = paste0(
            sprintf('Gene expression from %s genes in component %s',length(genes_oi),component_oi)
        ), show_colnames = FALSE, silent=TRUE)
    
    ## Mitotic Index
    component_oi = '95'
    ### scatter
    plts[['components_oi-mitotic_index-scatter']] = data.frame(
            component = X[samples,component_oi], 
            index = properties[samples,'mitotic_index']
        ) %>% 
        ggscatter(x='component', y='index', alpha=0.5) + 
        labs(x=sprintf('Weights Component %s',component_oi), y='Mitotic Index') + 
        stat_cor(method='spearman')
    ### heatmap of important genes
    genes_oi = S[,c('sample',component_oi)] %>% 
        filter(abs(get(component_oi)) > (4*sd(get(component_oi)) + mean(get(component_oi)))) %>% 
        pull(sample)
    plts[['components_oi-mitotic_index-heatmap_genexpr']] = pheatmap(
        genexpr %>% filter(sample %in% genes_oi) %>% column_to_rownames('sample'), 
        main = paste0(
            sprintf('Gene expression from %s genes in component %s',length(genes_oi),component_oi)
        ), show_colnames = FALSE, silent=TRUE)
    
    ## survival
    component_oi = '70'
    dat = cbind(component = X[samples,component_oi], properties[samples,c('OS.time','OS')])
    ### component weights vs OS.time
    plts[['components_oi-coxph-scatter']] =  dat %>% 
        ggscatter(x='OS.time', y='component', alpha=0.5) + 
        ggtitle(sprintf('Cox PH Z-score = %s', round(survs[component_oi,'coxph'],2))) +
        labs(y=sprintf('Weights Component %s',component_oi), x='Time (Days)')
    
    ### Kaplan Meier plot
    dat = dat %>% mutate(weight = ifelse(component < median(component), 'Low', 'High'))
    fit = survfit(Surv(OS.time, OS) ~ weight, data = dat)
    plts[['components_oi-coxph-kaplan_meier']] = fit %>% ggsurvplot(
            xlab='Time (Days)', ylab="Overall Survival Probability",
            legend.title = "Component Weight",
            pval = TRUE,
            conf.int = TRUE,
            palette = 'Dark2'
        ) + ggtitle(sprintf('KM from Weights Component %s',component_oi))
    
    ### heatmap of important genes
   genes_oi = S[,c('sample',component_oi)] %>% 
        filter(abs(get(component_oi)) > (4*sd(get(component_oi)) + mean(get(component_oi)))) %>% 
        pull(sample)
   plts[['components_oi-coxph-heatmap_genexpr']] = pheatmap(
        genexpr %>% filter(sample %in% genes_oi) %>% column_to_rownames('sample'), 
        main = paste0(
            sprintf('Gene expression from %s genes in component %s',length(genes_oi),component_oi)
        ), show_colnames = FALSE, silent=TRUE)
    
    
    return(plts)
}


plot_mut_components = function(genexpr, S, A, metadata){
    
    X = A %>% column_to_rownames('index') %>% t()
    
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
       
    plts = list()

    # mut vs wt counts
    plts[['mutation-counts_mutated']] = metadata %>%
        count(is_mut) %>%
        ggbarplot(x='is_mut', y='n', label=TRUE, fill='is_mut', color=NA)
    
    
    # overview sample weights and mutation status
    plts[['mutation-A-overview_heatmap']] = make_heatmap(X, metadata)
    
    
    plts[['mutation-A-volcano']] = result %>%
        ggscatter(x='med_diff', y='log10_fdr', color='silhouette_euclidean') +
        geom_hline(yintercept = -log10(0.01), linetype='dashed') +
        gradient_color(c("grey", "red")) +
        #geom_vline(xintercept = c(-5,5), linetype='dashed') + 
        geom_text_repel(
            aes(label=component), 
            result %>% filter(fdr<0.05 & abs(med_diff)>10)
        ) +
        labs(x=TeX('$\\Delta W_{mut-wt}$'), 
             y=TeX('$-log_{10}(FDR)$'),
             color='Mean Silhouette')
    
    component_oi = '62'
    
    # overview components of interest
    plts[['mutation-A-sel_component-heatmap']] = make_heatmap(X[component_oi,,drop=FALSE], metadata, cluster_cols=FALSE)
    
    # violins of components oi
    plts[['mutation-A-sel_component-violin']] = X %>% 
        melt(varnames = c('component','sampleID'), value.name='weight') %>% 
        left_join(metadata, by='sampleID') %>% 
        filter(component==component_oi) %>% 
        mutate(is_mut = ifelse(is_mut, 'Mut.', 'WT')) %>%
        ggviolin(x='component', y='weight', fill='is_mut', 
                 palette='lancet', color=NA) + 
        geom_boxplot(aes(group=is_mut), width=0.1, position = position_dodge(0.8)) +
        labs(x='Component', y='Weight', fill=TeX('\\textit{IDH1}'))
    
    ### heatmap of important genes
    genes_oi = S[,c('sample',component_oi)] %>% 
        filter(abs(get(component_oi)) > (4*sd(get(component_oi)) + mean(get(component_oi)))) %>% 
        pull(sample)
    plts[['mutation-sel_component-heatmap_genexpr']] = make_heatmap(genexpr %>% filter(sample %in% genes_oi) %>% column_to_rownames('sample'), metadata)
    
    return(plts)
}


make_plots = function(genexpr, S, A, sample_properties, metadata){
    plts = list(
        plot_sample_properties(genexpr, S, A, sample_properties),
        plot_silhouettes(stats),
        plot_mut_components(genexpr, S, A, metadata)
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
    save_plot(plts[['sample_properties-coxph']], 'sample_properties-coxph', '.pdf', figs_dir, width=5, height=25)
    save_plot(plts[['sample_properties-correls']], 'sample_properties-correls', '.pdf', figs_dir, width=15, height=25)
    
    save_plot(plts[['components_oi-estimate-scatter']], 'components_oi-estimate-scatter', '.png', figs_dir, width=12, height=12)
    save_plot(plts[['components_oi-mitotic_index-scatter']], 'components_oi-mitotic_index-scatter', '.png', figs_dir, width=12, height=12)
    save_plot(plts[['components_oi-coxph-scatter']], 'components_oi-coxph-scatter', '.png', figs_dir, width=12, height=12)
    ggsave(file.path(figs_dir,'components_oi-coxph-kaplan_meier.pdf'), print(plts[['components_oi-coxph-kaplan_meier']]), units='cm', width=12, height=12)
    
    save_plot(plts[['components_oi-estimate-heatmap_genexpr']], 'components_oi-estimate-heatmap_genexpr', '.pdf', figs_dir, width=12, height=12)
    save_plot(plts[['components_oi-mitotic_index-heatmap_genexpr']], 'components_oi-mitotic_index-heatmap_genexpr', '.pdf', figs_dir, width=12, height=12)
    save_plot(plts[['components_oi-coxph-heatmap_genexpr']], 'components_oi-coxph-heatmap_genexpr', '.pdf', figs_dir, width=12, height=12)
    
    
    save_plot(plts[['silhouettes-scatter']], 'silhouettes-scatter', '.pdf', figs_dir, width=12, height=12)
    save_plot(plts[['silhouettes-scatter-unfiltered']], 'silhouettes-scatter-unfiltered', '.pdf', figs_dir, width=12, height=12)
    
   save_plot(plts[['mutation-counts_mutated']], 'mutation-counts_mutated', '.pdf', figs_dir, width=12, height=12) 
    save_plot(plts[['mutation-A-volcano']], 'mutation-A-volcano', '.pdf', figs_dir, width=12, height=12)
    save_plot(plts[['mutation-A-sel_component-violin']], 'mutation-A-sel_component-violin', '.pdf', figs_dir, width=12, height=12)
    
    save_plot(plts[['mutation-A-overview_heatmap']], 'mutation-A-overview_heatmap', '.pdf', figs_dir,width=15, height=17)
    save_plot(plts[['mutation-A-sel_component-heatmap']], 'mutation-A-sel_component-heatmap', '.pdf', figs_dir, width=15, height=17)
    save_plot(plts[['mutation-A-sel_component-heatmap_genexpr']], 'mutation-A-sel_component-heatmap_genexpr', '.pdf', figs_dir, width=15, height=17)
    
}


main = function(){
    args = getParsedArgs()
    S_file = args$S_file
    A_file = args$A_file
    stats_file = args$stats_file
    genexpr_file = args$genexpr_file
    snv_file = args$snv_file
    metadata_file = args$metadata_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load data
    S = read_tsv(S_file)
    A = read_tsv(A_file)
    stats = read_tsv(stats_file)
    genexpr = read_tsv(genexpr_file)
    snv = read_tsv(snv_file)
    metadata = read_tsv(metadata_file)
    sample_properties = read_tsv(sample_properties_file)
    
    # checkout silhouettes
    plt = plot_silhouettes(stats)[[1]]
    
    # drop unwanted components: -1, silhouette
    clean_components = stats %>% filter(silhouette_euclidean>0.8) %>% pull(cluster_id)
    S = S[,c(1,1+clean_components)]
    A = A[,c(1,1+clean_components)]
    stats = stats %>% filter(cluster_id %in% clean_components)
    
    # add IDH1 mutation in metadata
    mut_samples = snv %>% 
        filter(gene==GENE_OI & effect%in%EFFECT_OI) %>% 
        pull(sampleID)
    metadata = metadata %>% mutate(is_mut = as.character(sampleID %in% mut_samples))
    
    # make plots
    plts = make_plots(genexpr, S, A, sample_properties, metadata)
    plts[['silhouettes-scatter-unfiltered']] = plt
    
    # save
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
