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
require(grid)

require(clusterProfiler)
require(org.Hs.eg.db)
require(enrichplot)

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


analyze_components = function(A, sample_properties){
    
    # do any weights correlate with any sample property?
    X = A %>% column_to_rownames('index')
    properties = sample_properties %>% column_to_rownames('sampleID')
    properties_oi = setdiff(colnames(properties),c('OS','OS.time','sample_type','ABSOLUTE'))
    samples = rownames(X)

    correls = sapply(properties_oi, function(property_oi){
        y = properties[samples,property_oi]
        correls = apply(X[samples,], 2, function(x){
            cor(x,y, method='spearman', use='pairwise.complete.obs')
        })
    }, simplify=FALSE)
    correls = do.call(rbind, correls)
    
    # Are there any survival associations?
    survs = apply(X[samples,], 2, function(x){
        dat = cbind(weight=x, properties[samples,c('OS','OS.time')])
        result = summary(coxph(Surv(OS.time, OS) ~ weight, dat))
        assoc = result$coefficients['weight','z']
        return(assoc)
    })
    survs = data.frame(coxph = survs)
    
    # which components have different weights in samples with and without mutation?
    samples = metadata[['sampleID']]
    is_mut = metadata[['is_mut']]=='TRUE'
    
    result = apply(X[samples,], 2, function(x){
        a = x[is_mut]
        b = x[!is_mut]

        result = data.frame(
            a_median = median(a, na.rm=TRUE), 
            b_median = median(b, na.rm=TRUE),
            med_diff = median(a, na.rm=TRUE) - median(b, na.rm=TRUE),
            pvalue = wilcox.test(a, b)$p.value
        )    
        
        return(result)
    })
    result = do.call(rbind,result) %>% 
        rownames_to_column('component') %>%
        mutate(fdr = p.adjust(pvalue, method='fdr'),
               log10_fdr = -log10(fdr))
    
    
    results = list(
        sample_indices = correls,
        survival = survs,
        mutation = result
    )
    
    return(results)
}


plot_silhouettes = function(stats){
    X = stats

    plts=list()
    
    plts[['silhouettes-scatter']] = X %>% 
        ggplot(aes(x=silhouette_euclidean, y=silhouette_pearson)) +
        geom_text(aes(label=cluster_id)) + 
        theme_pubr() + 
        geom_abline(intercept = 0, slope = 1, linetype='dashed') +
        labs(x='Euclidean', y='1 - abs(Pearson Corr.)')
    
    return(plts)
}


run_enrichment = function(genes_oi){
    result = enrichGO(genes_oi, OrgDb = org.Hs.eg.db, keyType='SYMBOL', ont='BP')
    return(result)
}


plot_component_oi = function(genexpr, S, A, metadata, sample_properties, component_oi){    
    # weights, mitotic index, survival, mutation status
    X = A[,c('index',component_oi)] 
    colnames(X) = c('sampleID', 'component')
    X = X %>%
        left_join(metadata[,c('sampleID','OS','OS.time','is_mut')], 
                  by='sampleID') %>%
        left_join(sample_properties[,c('sampleID','mitotic_index')],
                  by='sampleID') %>%
        mutate(OS = OS == 1,
               is_mut = ifelse(is_mut, 'Mutated', 'WT'),
               weight = ifelse(component < median(component), 'Low', 'High')) %>%
        drop_na(OS)
        
    
    plts = list()
    plts[['weights_vs_mitotic_index']] = X %>%
        ggscatter(x='mitotic_index', y='component', alpha=0.5) + 
        labs(x='Mitotic Index', y=sprintf('Weights Component %s',component_oi)) + 
        stat_cor(method='spearman')
    
    plts[['weights_vs_surv_time']] = X %>%
        ggscatter(x='OS.time', y='component', palette='Dark2',
                  color='OS', shape='OS', alpha=0.5) + 
        labs(x='Time (Days)', y=sprintf('Weights Component %s',component_oi), 
             color='Censored', shape='Censored')
    
    fit = survfit(Surv(OS.time, OS) ~ weight, data = X)
    plts[['weights_vs_surv_km']] = fit %>% 
        ggsurvplot(
            xlab='Time (Days)', ylab="Overall Survival Probability",
            legend.title = "Component Weight",
            pval = TRUE,
            conf.int = TRUE,
            palette = 'Dark2'
        ) + ggtitle(sprintf('KM from Weights Component %s',component_oi))
    
    plts[['weights_vs_mutation']] = X %>%
        ggviolin(x='is_mut', y='component', 
                 fill='is_mut', color=NA, palette='lancet') +
        geom_boxplot(width=0.1) +
        stat_compare_means(method='wilcox.test') +
        guides(fill=FALSE) +
        labs(x=TeX('\\textit{IDH1}'), y=sprintf('Weights Component %s',component_oi))
    
    # genes info
    df = S[,c('sample',component_oi)]
    colnames(df) = c('gene', 'weights')
    df[['component']] = component_oi
    plts[['weights_genes']] = df %>% 
        ggviolin(x='component', y='weights', fill='orange', color=NA) +
        geom_boxplot(width=0.1, outlier.size = 0.1) +
        labs(x='Component', y= sprintf('Weights Component %s',component_oi)) +
        geom_text_repel(
            aes(label=gene), 
            df %>% slice_max(order_by = abs(weights), n = 20)
        )
    
    genes_oi = df %>% filter(abs(weights) > 4*sd(weights)) %>% pull(gene)
    plts[['weights_vs_genexpr']] = make_heatmap(
        genexpr %>% filter(sample%in%genes_oi) %>% column_to_rownames('sample'), 
        metadata
    )
    
    enrichment = run_enrichment(genes_oi)
    plts[['weights_vs_enrichment_dotplot']] = dotplot(enrichment) +
        ggtitle(sprintf('n=%s',length(genes_oi)))
    
    plts[['weights_vs_enrichment_cnetplot']] = cnetplot(enrichment) +
        ggtitle(sprintf('n=%s',length(genes_oi)))
    
    return(plts)  
}


make_plots = function(genexpr, S, A, sample_properties, metadata){
    plts = list(
        plot_component_oi(genexpr, S, A, metadata, sample_properties, component_oi),
        plot_silhouettes(stats)
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
    save_plot(plts[['weights_vs_mitotic_index']], 'weights_vs_mitotic_index', '.png', figs_dir, width=12, height=12)
    
    save_plot(plts[['weights_vs_surv_time']], 'weights_vs_surv_time', '.png', figs_dir, width=12, height=12)
    save_plot(print(plts[['weights_vs_surv_km']]), 'weights_vs_surv_km', '.pdf', figs_dir, width=12, height=12)
    
    save_plot(plts[['weights_vs_mutation']], 'weights_vs_mutation', '.pdf', figs_dir, width=12, height=12)
    
    save_plot(plts[['weights_genes']], 'weights_genes', '.pdf', figs_dir, width=12, height=12)
    save_plot(plts[['weights_vs_genexpr']], 'weights_vs_genexpr', '.pdf', figs_dir, width=24, height=12)
    
    save_plot(plts[['weights_vs_enrichment_dotplot']], 'weights_vs_enrichment_dotplot', '.pdf', figs_dir, width=20, height=20)
    save_plot(plts[['weights_vs_enrichment_cnetplot']], 'weights_vs_enrichment_cnetplot', '.pdf', figs_dir, width=20, height=20)
    
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
    
    # analysis - 95 stands out
    results = analyze_components(A, sample_properties)
    comps_correls = sort(abs(results[['sample_indices']]['mitotic_index',]))
    comps_survs = results[['survival']] %>% arrange(abs(coxph))
    comps_mitotic_index = results[['mutation']] %>% filter(fdr < 0.05) %>% arrange(abs(med_diff))
    
    # make plots for this component
    plts = plot_component_oi(genexpr, S, A, metadata, sample_properties, '95')
    plts[['silhouettes-scatter-unfiltered']] = plt
    
    # save
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
