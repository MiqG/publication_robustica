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
require(fdrtool)
require(writexl)

require(clusterProfiler)
require(org.Hs.eg.db)
require(enrichplot)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variable
GENE_OI = 'IDH1'
EFFECT_OI = c('Missense_Mutation','Nonsense_Mutation')

MSIGDB_HALLMARKS = file.path(ROOT,'data','raw','MSigDB','msigdb_v7.4','msigdb_v7.4_files_to_download_locally','msigdb_v7.4_GMTs','h.all.v7.4.symbols.gmt')
MSIGDB_IMMUNE = file.path(ROOT,'data','raw','MSigDB','msigdb_v7.4','msigdb_v7.4_files_to_download_locally','msigdb_v7.4_GMTs','c7.immunesigdb.v7.4.symbols.gmt')

THRESH_FDR = 0.01

# Development
# -----------
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','case_study')

# S_file = file.path(RESULTS_DIR,'files','cluster_iterations','LGG','S.tsv.gz')
# A_file = file.path(RESULTS_DIR,'files','cluster_iterations','LGG','A.tsv.gz')
# stats_file = file.path(RESULTS_DIR,'files','cluster_iterations','LGG','stats.tsv.gz')

# genexpr_file = file.path(PREP_DIR,'genexpr','LGG.tsv.gz')
# snv_file = file.path(PREP_DIR,'snv','LGG.tsv.gz')
# metadata_file = file.path(PREP_DIR,'metadata','LGG.tsv')
# figs_dir = file.path(RESULTS_DIR,'figures','LGG')

# sample_indices_file = file.path(PREP_DIR,'sample_indices','LGG.tsv')

##### FUNCTIONS #####
define_module = function(x, cutoff=0.01){
    fdr = fdrtool(x, plot=FALSE, cutoff.method="fndr", verbose=FALSE)[['qval']]
    x = fdr < cutoff
    return(x)
}


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


analyze_components = function(A, metadata){
    
    # do any weights correlate with any sample property?
    X = A %>% column_to_rownames('index')
    metadata = metadata %>% column_to_rownames('sampleID')
    properties_oi = 'mitotic_index'
    samples = rownames(X)

    correls = sapply(properties_oi, function(property_oi){
        y = metadata[samples,property_oi]
        correls = apply(X[samples,,drop=FALSE], 2, function(x){
            cor(x,y, method='spearman', use='pairwise.complete.obs')
        })
    }, simplify=FALSE)
    correls = do.call(rbind, correls)
    
    # Are there any survival associations?
    survs = apply(X[samples,,drop=FALSE], 2, function(x){
        dat = cbind(weight=x, metadata[samples,c('OS','OS.time'),drop=FALSE])
        result = summary(coxph(Surv(OS.time, OS) ~ weight, dat))
        assoc = result$coefficients['weight','z']
        return(assoc)
    })
    survs = data.frame(coxph = survs)
    
    # which components have different weights in samples with and without mutation?
    samples = rownames(metadata)
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


run_enrichments = function(genes_oi, universe, msigdb){
    print('Running enrichments...')
    result = list(
        'GO_BP' = enrichGO(genes_oi, OrgDb = org.Hs.eg.db, keyType='SYMBOL', ont='BP', universe = universe),
        'MSigDB_Hallmarks' = enricher(genes_oi, TERM2GENE=msigdb[['h']], universe = universe)
    )
    return(result)
}


plot_component_oi = function(genexpr, S, A, metadata, sample_indices, 
                             component_oi, enrichments){    
    # weights, mitotic index, survival, mutation status
    X = A[,c('index',component_oi)] 
    colnames(X) = c('sampleID', 'component')
    X = X %>%
        left_join(metadata[,c('sampleID','OS','OS.time','is_mut')], 
                  by='sampleID') %>%
        left_join(sample_indices[,c('sampleID','mitotic_index')],
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
    plts[['weights_vs_surv_km']] = ggsurvplot(
            fit = fit, data = X,
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
    df$fdr = fdrtool(df$weights, plot = FALSE, verbose=FALSE)[['qval']]
    genes_oi = df %>% filter(fdr < THRESH_FDR) %>% pull(gene)
    cutoffs = df %>% 
        filter(fdr < THRESH_FDR) %>% 
        mutate(is_pos = weights > 0) 
    cutoffs = c(
        cutoffs %>% filter(is_pos) %>% pull(weights) %>% min(),
        cutoffs %>% filter(!is_pos) %>% pull(weights) %>% max()
    )
    plts[['weights_genes']] = df %>% 
        ggviolin(x='component', y='weights', fill='orange', color=NA) +
        geom_boxplot(width=0.1, outlier.size = 0.1) +
        labs(x='Component', y= sprintf('Weights Component %s',component_oi)) +
        geom_hline(yintercept = cutoffs, linetype='dashed')

    plts[['weights_vs_genexpr']] = make_heatmap(
        genexpr %>% filter(sample%in%genes_oi) %>% column_to_rownames('sample'), 
        metadata
    )

    plts[['weights_vs_enrichment_dotplot-GO_BP']] = dotplot(enrichments[['GO_BP']]) +
        ggtitle(sprintf('GO Biological Processes | n=%s',length(genes_oi)))
    plts[['weights_vs_enrichment_cnetplot-GO_BP']] = cnetplot(enrichments[['GO_BP']]) +
        ggtitle(sprintf('GO Biological Processes | n=%s',length(genes_oi)))
    
    plts[['weights_vs_enrichment_dotplot-MSigDB_Hallmarks']] = dotplot(enrichments[['MSigDB_Hallmarks']]) +
        ggtitle(sprintf('MSigDB Hallmarks | n=%s',length(genes_oi)))
    plts[['weights_vs_enrichment_cnetplot-MSigDB_Hallmarks']] = cnetplot(enrichments[['MSigDB_Hallmarks']]) +
        ggtitle(sprintf('MSigDB Hallmarks | n=%s',length(genes_oi)))
    
        
    return(plts) 
    
    # component 95 correlates well with E2F transcription factors
    # genexpr %>% filter(grepl('^E2F', sample)) %>% column_to_rownames('sample') %>% t() %>% as.data.frame() %>% rownames_to_column('index') %>% left_join(A %>% dplyr::select(one_of(c('index','95'))), by='index') %>% column_to_rownames('index') %>% cor(method='spearman') %>% pheatmap()
    
    # upstream regulator of E2F differentially expressed in IDH mut
    # genexpr %>% filter(grepl('CDKN1A',sample)) %>% column_to_rownames('sample') %>% t() %>% as.data.frame() %>% rownames_to_column('sampleID') %>% left_join(metadata[,c('sampleID','is_mut')], by='sampleID') %>% ggviolin(x='is_mut', y='CDKN1A') + stat_compare_means(method = 'wilcox.test')
    
    # clustering the genes by their gene expression
    # guilt-by-association
    # plt = genexpr %>% filter(sample %in% genes_oi) %>% column_to_rownames('sample') %>% t() %>% cor(method='spearman') %>% pheatmap(cutree_cols = 4, cutree_rows = 4)
 
}


make_figdata = function(S, A, metadata, sample_indices, enrichments){
    # make gene modules
    is_selected = S %>% mutate_at(vars(-('sample')), define_module)
    modules = S %>% dplyr::rename(gene=sample)
    modules[,-1][!is_selected[,-1]] = NA
    
    figdata = list(
        'case_study-analysis' = list(
            'source_matrix_S' = S %>% dplyr::rename(gene=sample),
            'mixing_matrix_A' = A,
            'sample_metadata' = metadata,
            'sample_indices' = sample_indices,
            'gene_modules' = modules,
            'enrichment-GO_BP' = as.data.frame(enrichments[['GO_BP']]),
            'enrichment-MSigDB_Hallmarks' = as.data.frame(enrichments[['MSigDB_Hallmarks']])
        )
    )
    return(figdata)
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
    save_plot(plts[['weights_vs_genexpr']], 'weights_vs_genexpr', '.pdf', figs_dir, width=70, height=15)
    
    save_plot(plts[['weights_vs_enrichment_dotplot-GO_BP']], 'weights_vs_enrichment_dotplot-GO_BP', '.pdf', figs_dir, width=20, height=20)
    save_plot(plts[['weights_vs_enrichment_cnetplot-GO_BP']], 'weights_vs_enrichment_cnetplot-GO_BP', '.pdf', figs_dir, width=20, height=20)
    
    save_plot(plts[['weights_vs_enrichment_dotplot-MSigDB_Hallmarks']], 'weights_vs_enrichment_dotplot-MSigDB_Hallmarks', '.pdf', figs_dir, width=20, height=20)
    save_plot(plts[['weights_vs_enrichment_cnetplot-MSigDB_Hallmarks']], 'weights_vs_enrichment_cnetplot-MSigDB_Hallmarks', '.pdf', figs_dir, width=20, height=20)
    
    save_plot(plts[['silhouettes-scatter-unfiltered']], 'silhouettes-scatter-unfiltered', '.pdf', figs_dir, width=12, height=12)
}


save_figdata = function(figdata, dir){
    lapply(names(figdata), function(x){
        filename = file.path(dir,'figdata',paste0(x,'.xlsx'))
        dir.create(dirname(filename), recursive=TRUE)
        write_xlsx(figdata[[x]], filename)
    })
}


main = function(){
    args = getParsedArgs()
    S_file = args$S_file
    A_file = args$A_file
    stats_file = args$stats_file
    genexpr_file = args$genexpr_file
    snv_file = args$snv_file
    metadata_file = args$metadata_file
    stats_file = args$stats_file
    sample_indices_file = args$sample_indices_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load data
    msigdb = list(
        h = read.gmt(MSIGDB_HALLMARKS),
        c7 = read.gmt(MSIGDB_IMMUNE)
    )    
    
    S = read_tsv(S_file)
    A = read_tsv(A_file)
    stats = read_tsv(stats_file)
    genexpr = read_tsv(genexpr_file)
    snv = read_tsv(snv_file)
    metadata = read_tsv(metadata_file)
    sample_indices = read_tsv(sample_indices_file)
    
    # checkout silhouettes
    plt = plot_silhouettes(stats)[[1]]
    
    # drop unwanted components: -1, silhouette
    clean_components = stats %>% 
        filter(silhouette_pearson>0.9) %>% 
        pull(cluster_id)
    S = S[,c(1,1+clean_components)]
    A = A[,c(1,1+clean_components)]
    stats = stats %>% filter(cluster_id %in% clean_components)
    
    # add IDH1 mutation in metadata
    mut_samples = snv %>% 
        filter(gene==GENE_OI & effect%in%EFFECT_OI) %>% 
        pull(sampleID)
    metadata = metadata %>% 
        mutate(is_mut = as.character(sampleID %in% mut_samples)) %>%
        left_join(sample_indices, by='sampleID')
    
    # analysis - 95 stands out
    results = analyze_components(A, metadata)
    comps_correls = sort(abs(results[['sample_indices']]['mitotic_index',]))
    comps_survs = results[['survival']] %>% arrange(abs(coxph))
    comps_mitotic_index = results[['mutation']] %>% filter(fdr < 0.05) %>% arrange(abs(med_diff))
    
    # enrichment for component of interest
    component_oi = '95'
    df = S[,c('sample',component_oi)]
    colnames(df) = c('gene', 'weights')
    df[['component']] = component_oi
    df$fdr = fdrtool(df$weights, plot = FALSE, verbose=FALSE)[['qval']]
    genes_oi = df %>% filter(fdr < THRESH_FDR) %>% pull(gene)
    enrichments = run_enrichments(genes_oi, df %>% pull(gene), msigdb)
    
    # make plots for this component
    plts = plot_component_oi(genexpr, S, A, metadata, sample_indices, '95', enrichments)
    plts[['silhouettes-scatter-unfiltered']] = plt
    
    # make figdata
    figdata = make_figdata(S, A, metadata, sample_indices, enrichments)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
