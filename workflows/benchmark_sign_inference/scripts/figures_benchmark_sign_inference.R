#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Evaluate the performance in terms of memory and execution time of the
# different algorithms to compute robust independent components.
#
# Outline
# -------
# - memory usage and time spent for every algorithm
# - silhouette scores of every algorithm


require(tidyverse)
require(ggpubr)
require(pheatmap)
require(reshape2)
require(ggrepel)
require(proxy)
require(fdrtool)
require(latex2exp)
require(writexl)
require(ggnewscale)

require(org.EcK12.eg.db)
require(clusterProfiler)
require(enrichplot)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
ALGORITHMS = c('icasso','robustica_nosign','robustica','robustica_pca')
ORGDB = org.EcK12.eg.db
THRESH_JACCARD = 0.9

# Development
# -----------
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','benchmark_sign_inference','files')
# algorithms = c('icasso','robustica_nosign','robustica','robustica_pca')

# performance_evaluation_file = file.path(RESULTS_DIR,'Sastry2019','performance_evaluation.tsv.gz')
# clustering_info_file = file.path(RESULTS_DIR,'Sastry2019','clustering_info.tsv.gz')

# S_means_files = unlist(strsplit(file.path(RESULTS_DIR,'Sastry2019',paste0(algorithms,'-','S.tsv.gz')),','))
# A_means_files = unlist(strsplit(file.path(RESULTS_DIR,'Sastry2019',paste0(algorithms,'-','A.tsv.gz')),','))

# S_stds_files = unlist(strsplit(file.path(RESULTS_DIR,'Sastry2019',paste0(algorithms,'-','S_std.tsv.gz')),','))
# A_stds_files = unlist(strsplit(file.path(RESULTS_DIR,'Sastry2019',paste0(algorithms,'-','A_std.tsv.gz')),','))

# gene_annotation_file=file.path(PREP_DIR,'gene_annotation','Sastry2019.tsv')
# genexpr_file = file.path(PREP_DIR,'genexpr','Sastry2019.tsv.gz')

# figs_dir = file.path(ROOT,'results','benchmark_sign_inference','figures','Sastry2019')
# files_dir = file.path(ROOT,'results','benchmark_sign_inference','files','Sastry2019')


##### FUNCTIONS #####
get_relative_timestamp = function(x){
    x0=x[1]
    return(x-x0)
}


load_mats = function(files, rwnms, idvar, valuename){
    mats = lapply(files, function(file){
        read_tsv(file) %>% 
        column_to_rownames(rwnms) %>%
        rownames_to_column(idvar) %>%
        melt(id.vars = idvar, 
             variable.name = 'component', 
             value.name = valuename) %>% 
        mutate(algorithm = gsub('-.*','',basename(file)))
    })
    mats = do.call(rbind, mats)
    mats = mats %>% 
        mutate(
            component = paste0('comp',component),
            algorithm = factor(algorithm, levels = ALGORITHMS)
        )
    return(mats)
}


plot_corr_vs_euc = function(len=10){
    set.seed(123)
    compA = rnorm(len)
    compB = rnorm(len)
    compC = compA - 6
    df = data.frame(
        feature = 1:length(compA),
        compA,
        compB,
        compC
    ) %>% melt(id.vars = 'feature')
    
    metrics = data.frame(
        cor_AB = cor(compA, compB, method = 'pearson'),
        cor_AC = cor(compA, compC, method = 'pearson'),
        cor_BC = cor(compB, compC, method = 'pearson'),
        dist_AB = dist(rbind(compA, compB), method = 'euclidean')[1],
        dist_AC = dist(rbind(compA, compC), method = 'euclidean')[1],
        dist_BC = dist(rbind(compB, compC), method = 'euclidean')[1]
    ) %>% melt()
    
    palette = get_palette('simpsons', 3)
    plts = list()
    
    plts[['corr_vs_euc-scatter']] = df %>% 
        ggline(x='feature', y='value', color='variable', 
               linetype='dashed', palette=palette) + 
        labs(x="Components' Feature", y='Weight', color='Component') + 
        theme(aspect.ratio=1)
    
    plts[['corr_vs_euc-loess_compAB']] = df %>% 
        mutate(pair = variable %in% c('compA','compB')) %>% 
        ggscatter(x='feature', y='value', color='pair', conf.int = TRUE, 
                  add='loess', add.params = list(linetype='dashed'))
    
    plts[['corr_vs_euc-loess_compAC']] = df %>% 
        mutate(pair = variable %in% c('compA','compC')) %>% 
        ggscatter(x='feature', y='value', color='pair', conf.int = TRUE, 
                  add='loess', add.params = list(linetype='dashed'), 
                  palette='Dark2')
    
    plts[['corr_vs_euc-table']] = ggtexttable(metrics, rows = NULL)
    
    return(plts)
}


plot_performance_profile = function(df){
    X = df
    
    plts = list()
    plts[['mem_time-scatter']] = X %>% 
        ggline(x='rel_timestamp', y='memory',numeric.x.axis = TRUE, 
               color='function', linetype='dashed', 
               shape=20, point.size=0.0001, 
               palette='uchicago',alpha=0.5) + 
    facet_wrap(~algorithm, scales='free_x', ncol=1) + 
    guides(color=guide_legend(ncol=1)) + 
    labs(x='Time (s)', y='Memory (MiB)', color='Function') + 
    theme_pubr(border = TRUE) 
    
    return(plts)
}


plot_clustering_evaluation = function(df){
    # get mean silhouette per cluster
    X = df %>%
        group_by(algorithm,component,silhouette_euclidean) %>%
        summarize(mean_std_cluster = mean(weight_std))
    
    plts = list()
    plts[['clustering-silhouettes-violins']] = X %>%
        ggviolin(x='algorithm', y='silhouette_euclidean', trim = TRUE,
                 fill='algorithm', color=NA, palette='Paired') + 
        guides(fill=FALSE) +
        geom_boxplot(width=0.1, outlier.size = 0.1) +
        labs(x='Algorithm', y='Silhouette Score') +
        theme_pubr(x.text.angle = 70)
     
    plts[['clustering-S_stds-violins']] = X %>% 
        ggviolin(x='algorithm', y='mean_std_cluster', trim = TRUE,
                 fill='algorithm', color=NA, palette='Paired') + 
        geom_boxplot(width=0.1, outlier.size = 0.1) +
        guides(fill=FALSE) +
        labs(x='Algorithm', y='Cluster Mean Std.') +
        theme_pubr(x.text.angle = 70)
    
     plts[['clustering-silhouettes_vs_stds-scatter']] = X %>%
        ggplot(aes(x=silhouette_euclidean, y=mean_std_cluster)) +
        geom_text(aes(label=component)) +
        facet_wrap(~algorithm, scales='free') + 
        theme_pubr(border=TRUE) +
        labs(x='Silhouette Score', y='Cluster Mean Std.')
    
    return(plts)
}


plot_weights_precision = function(S_info){
    
    correls = S_info %>% 
        group_by(algorithm, component, silhouette_euclidean) %>% 
        summarize(correlation = cor(abs(weight_mean), weight_std))
    X = correls %>% 
        group_by(algorithm) %>%
        arrange(correlation) %>% 
        filter(row_number()==1 | row_number()==n()) %>%
        left_join(S_info, by=c('algorithm','component'))
    

    plts = list()
    
    plts[['clustering-weightS_means_vs_std-violins']] = correls %>%
        ggviolin(x='algorithm', y='correlation', trim = TRUE,
                 fill='algorithm', color=NA, palette='Paired') + 
        geom_boxplot(width=0.05, outlier.size = 0.1) +
        guides(fill=FALSE) +
        labs(x='Algorithm', y='Correlation Weights Mean vs Std.')
    
    plts[['clustering-weightS_means_vs_std-scatters']] = X %>%
        ggscatter(x='weight_mean', y='weight_std', 
                  size=0.5, alpha=0.5) +
        facet_wrap(~algorithm+component, ncol=2) +
        geom_text(aes(x=-0.05,y=0.1,label=sprintf('R=%s',round(correlation,2))), 
                  X %>% distinct(algorithm, component, correlation)) +
        theme_pubr(border = TRUE)
    
    return(plts)
}


define_module = function(x, cutoff=0.01){
    fdr = fdrtool(x, plot=FALSE, cutoff.method="fndr", verbose=FALSE)[['qval']]
    x = fdr < cutoff
    return(x)
}


run_enrichments = function(genes_oi){
    print('Running enrichments...')
    symbols = gene_annotation %>% 
        filter(b_number %in% genes_oi) %>% 
        pull(gene_name)
    result = enrichGO(symbols, OrgDb = ORGDB, keyType='SYMBOL', ont='BP')
    return(result)
}


make_module_comparisons = function(S_info){
    # get modules from source matrix
    modules = S_info %>%
        group_by(algorithm, component) %>%
        mutate(in_module=define_module(weight_mean))
    
    # module sizes
    module_size = modules %>% 
        filter(in_module) %>% 
        count(algorithm, component)
    
    # use jaccard similarity to map modules
    icasso = modules %>% 
        filter(algorithm == 'icasso') %>% 
        pivot_wider(id_cols = gene, names_from = component, values_from = in_module)
    robustica_pca = modules %>% 
        filter(algorithm == 'robustica_pca') %>% 
        pivot_wider(id_cols = gene, names_from = component, values_from = in_module)
    sim = simil(icasso[,-1], robustica_pca[,-1], 
                method='Jaccard', by_rows=FALSE) %>% as.matrix()
    
    
    # map modules
    rows2cols = paste0(rownames(sim),'_',colnames(sim)[apply(sim,1,which.max)])
    cols2rows = paste0(rownames(sim)[apply(sim,2,which.max)],'_',colnames(sim))
    mapped = intersect(rows2cols, cols2rows) %>% str_split_fixed('_', 2) %>% as.data.frame()
    unmapped_icasso = setdiff(rows2cols, cols2rows) %>% str_split_fixed('_',2) %>% as.data.frame() %>% mutate(V2=NA)
    unmapped_robustica_pca = setdiff(cols2rows, rows2cols) %>% str_split_fixed('_',2) %>% as.data.frame() %>% mutate(V1=NA)
    mapping = rbind(
        mapped,
        unmapped_icasso,
        unmapped_robustica_pca
    )
    colnames(mapping) = c('icasso','robustica_pca')
    mapping = mapping %>% mutate(jaccard = sim[cbind(icasso, robustica_pca)])
        
    # add module sizes
    module_sizes = module_size %>% 
        filter(algorithm %in% c('icasso','robustica_pca')) %>%
        pivot_wider(id_cols = component, names_from = algorithm, values_from = n) %>%
        column_to_rownames('component')
    
    mapping = mapping %>% 
        mutate(
            icasso_module_size = module_sizes[icasso,'icasso'],
            robustica_pca_module_size = module_sizes[robustica_pca,'robustica_pca'],
            module_size_diff = icasso_module_size - robustica_pca_module_size,
            abs_module_size_diff = abs(module_size_diff)
        ) %>%
        mutate(
            diff_summary = ifelse(module_size_diff == 0, 'equal', 'diff'),
            diff_summary = ifelse(module_size_diff > 0, 'icasso', diff_summary),
            diff_summary = ifelse(module_size_diff < 0, 'robustica_pca', diff_summary)
        )
    
    
    # add module silhouettes and S_stds
    mapping = mapping %>% 
        left_join(
            S_info %>% filter(algorithm=='icasso') %>% 
            group_by(component, silhouette_euclidean) %>% 
            summarize(avg_weight_std=mean(weight_std)) %>% 
            rename_all(~ paste0('icasso_',.x)) %>%
            dplyr::rename(icasso=icasso_component),
            by = 'icasso'
        ) %>% 
        left_join(
            S_info %>% filter(algorithm=='robustica_pca') %>% 
            group_by(component, silhouette_euclidean) %>% 
            summarize(avg_weight_std=mean(weight_std)) %>% 
            rename_all(~ paste0('robustica_pca_',.x)) %>%
            dplyr::rename(robustica_pca=robustica_pca_component),
            by = 'robustica_pca'
        )   

#     # enrichments in modules with low jaccard
#     enrichments = list()
#     ## icasso
#     components_oi = mapping %>% 
#         filter(!is.na(icasso)) %>%
#         filter(jaccard < THRESH_JACCARD | is.na(jaccard)) %>%
#         pull(icasso)
#     genes_oi = modules %>% 
#         filter(algorithm=='icasso' & component%in%components_oi & in_module)
#     genes_oi = split(genes_oi$gene, genes_oi$component)
    
#     enrichments[['icasso']] = sapply(genes_oi, simplify=FALSE, function(x){run_enrichments(x)})
    
#     ## robustica
#     components_oi = mapping %>% 
#         filter(!is.na(robustica_pca)) %>%
#         filter(jaccard < THRESH_JACCARD | is.na(jaccard)) %>%
#         pull(robustica_pca)
#     genes_oi = modules %>% 
#         filter(algorithm=='robustica_pca' & component%in%components_oi & in_module)
#     genes_oi = split(genes_oi$gene, genes_oi$component)
    
#     enrichments[['robustica_pca']] = sapply(genes_oi, simplify=FALSE, function(x){run_enrichments(x)})
    
#     mapping = mapping %>% 
#         left_join(
#             sapply(enrichments[['icasso']], nrow) %>% 
#             enframe(name = 'icasso', value = 'icasso_enriched'),
#             by = 'icasso'
#         )%>% 
#         left_join(
#             sapply(enrichments[['robustica_pca']], nrow) %>% 
#             enframe(name = 'robustica_pca', value = 'robustica_pca_enriched'),
#             by = 'robustica_pca'
#         )
    
    module_comparisons = list(
        'sim' = sim,
        'mapping' = mapping
    )
    return(module_comparisons)
}


plot_module_comparisons = function(module_comparisons){
    # unlist data
    sim = module_comparisons[['sim']]
    mapping = module_comparisons[['mapping']]
    
    # make plots
    plts = list()
    
    plts[['module_comparisons-jaccard']] = sim %>% pheatmap(silent=TRUE, cluster_cols = FALSE, cluster_rows=FALSE, fontsize = 2)
    
    plts[['module_comparisons-module_sizes-scatter']] = mapping %>% 
        arrange(jaccard) %>%
        ggscatter(x='icasso_module_size',
                  y='robustica_pca_module_size',
                  color='jaccard', alpha=0.5) +
        scale_color_gradient2(low = 'blue', mid='grey', high = 'red', midpoint = 0.5) + 
        guides(size=FALSE) + 
        geom_abline(intercept = 0, slope = 1, linetype='dashed') +
        new_scale("color") +
        geom_text_repel(
            aes(label=value, color=name), 
            mapping %>% 
                filter(abs_module_size_diff > 8) %>% 
                pivot_longer(cols = c(icasso,robustica_pca)),
            box.padding = 0.7
        ) +
        scale_color_manual("Algorithm", values=c('#A6CEE3','#33A02C')) +
        labs(x='Module Size Icasso', 
             y='Module Size robustica', 
             color='Jaccard Sim.',
             title=TeX(sprintf('Mapped gene modules between \\textit{Icasso} and \\textit{robustica} | n=%s',sum(mapping[['mapped_icasso']]))))
    
    
    plts[['module_comparisons-module_sizes_diffs-hist']] = mapping %>%
        gghistogram(x='module_size_diff', fill='darkblue', color=NA, bins=25) +
        labs(x='Module Size Icasso - Module Size Robustica', y='Count')
    
    
    plts[['module_comparisons-module_sizes_diffs-barplot']] = mapping %>% 
        count(diff_summary) %>% drop_na() %>%
        ggbarplot(x='diff_summary', y='n', palette='Dark2', 
                  fill='diff_summary', color=NA, label = TRUE) +
        guides(fill=FALSE) +
        labs(x='Summary of Module Size Difference')
    
    plts[['module_comparisons-module_sizes_diffs-strip']] = mapping %>% 
        drop_na(diff_summary) %>%
        ggstripchart(x='diff_summary', y='module_size_diff', alpha=0.5,
                     palette='Dark2', color='diff_summary') +
        labs(x='Summary of Module Size Difference') +
        guides(color=FALSE)
    
    plts[['module_comparisons-module_sizes_vs_jaccard']] = mapping %>% 
        ggscatter(x='robustica_pca_module_size', y='jaccard', add='reg.line', 
                  add.params=list(linetype='dashed'), conf.int = TRUE, 
                  size=1, alpha=0.5) + 
        stat_cor(method='spearman', label.x = 20)
    
    plts[['module_comparisons-low_jaccard']] = mapping %>% 
        mutate(is_low=jaccard<THRESH_JACCARD) %>% 
        count(is_low) %>%
        drop_na() %>%
        ggbarplot(x='is_low', y='n', fill='is_low', 
                  color=NA, label=TRUE, palette='jco') +
        labs(x=sprintf('Jaccard < %s',THRESH_JACCARD), y='Count') +
        guides(fill=FALSE)
        
    return(plts)
}


make_plots = function(performance_evaluation, S_info, module_comparisons){
    plts = list(
        plot_corr_vs_euc(),
        plot_performance_profile(performance_evaluation),
        plot_clustering_evaluation(S_info),
        plot_weights_precision(S_info),
        plot_module_comparisons(module_comparisons)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(performance_evaluation, clustering_evaluation, S_info, A_info, module_comparisons){
    
    # prep
    ## gene module evaluation
    modules = S_info %>%
        group_by(algorithm, component) %>%
        mutate(in_module=define_module(weight_mean))
    sim = module_comparisons[['sim']]
    sim = matrix(sim, nrow = nrow(sim), ncol=ncol(sim), dimnames = dimnames(sim)) %>% 
        as.data.frame()
    mapping = module_comparisons[['mapping']]
    
    # prepare figdata object
    figdata = list()
    figdata[['benchmark_sign_inference-ICA-source_matrices-weight_means']] = S_info
    figdata[['benchmark_sign_inference-ICA-source_matrices-weight_stds']] = S_info
    igdata[['benchmark_sign_inference-ICA-mixing_matrices-weight_means']] = A_info
    figdata[['benchmark_sign_inference-ICA-mixing_matrices-weight_stds']] = A_info
    
    figdata[['benchmark_sign_inference-evaluation']] = list(
        performance_evaluation = performance_evaluation,
        clustering_evaluation = clustering_evaluation,
        gene_modules = modules,
        gene_modules_jaccard_similarity = sim %>% rownames_to_column('component'),
        gene_modules_mapping = mapping
    )
    
    return(figdata)
}


save_plot = function(plt, plt_name, extension='.pdf', 
                      directory='', dpi=350, 
                      width = par("din")[1], height = par("din")[2], units='cm'){
        filename = file.path(directory,paste0(plt_name,extension))
        ggsave(filename, plt, width=width, height=height, 
               dpi=dpi, limitsize=FALSE, units=units)
}


save_plots = function(plts, figs_dir){
    # example
    save_plot(plts[['corr_vs_euc-scatter']],'corr_vs_euc-scatter','.pdf',figs_dir, width=12,height=12)
    save_plot(plts[['corr_vs_euc-loess_compAB']],'corr_vs_euc-loess_compAB','.pdf',figs_dir, width=12,height=12)
    save_plot(plts[['corr_vs_euc-loess_compAC']],'corr_vs_euc-loess_compAC','.pdf',figs_dir, width=12,height=12)
    save_plot(plts[['corr_vs_euc-table']],'corr_vs_euc-table','.pdf',figs_dir, width=12,height=12)
    
    # memory time profiles
    save_plot(plts[['mem_time-scatter']],'mem_time-scatter','.png',figs_dir, width=12,height=20)
    
    # clustering
    save_plot(plts[['clustering-silhouettes-violins']],'clustering-silhouettes-violins','.pdf',figs_dir, width=12,height=12)
    save_plot(plts[['clustering-S_stds-violins']],'clustering-S_stds-violins','.pdf',figs_dir, width=12,height=12)
    save_plot(plts[['clustering-silhouettes_vs_stds-scatter']],'clustering-silhouettes_vs_stds-scatter','.png',figs_dir, width=20,height=20)
    save_plot(plts[['clustering-weightS_means_vs_std-violins']],'clustering-weightS_means_vs_std-violins','.pdf',figs_dir, width=12,height=12)
    save_plot(plts[['clustering-weightS_means_vs_std-scatters']],'clustering-weightS_means_vs_std-scatters','.png',figs_dir, width=20,height=20)
    
    # module comparisons
    save_plot(plts[['module_comparisons-jaccard']],'module_comparisons-jaccard','.png',figs_dir, width=30, height=30)
    
    save_plot(plts[['module_comparisons-module_sizes-scatter']],'module_comparisons-module_sizes-scatter','.pdf',figs_dir, width=12, height=12)
    save_plot(plts[['module_comparisons-module_sizes_diffs-hist']],'module_comparisons-module_sizes_diffs-hist','.pdf',figs_dir, width=12, height=12)
    save_plot(plts[['module_comparisons-module_sizes_diffs-barplot']],'module_comparisons-module_sizes_diffs-barplot','.pdf',figs_dir, width=12, height=12)
    save_plot(plts[['module_comparisons-module_sizes_diffs-strip']],'module_comparisons-module_sizes_diffs-strip','.pdf',figs_dir, width=12, height=12)
    save_plot(plts[['module_comparisons-module_sizes_vs_jaccard']],'module_comparisons-module_sizes_vs_jaccard','.pdf',figs_dir, width=12, height=12)
    save_plot(plts[['module_comparisons-low_jaccard']],'module_comparisons-low_jaccard','.pdf',figs_dir, width=12, height=12)
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

    performance_evaluation_file = args$performance_evaluation_file
    clustering_info_file = args$clustering_info_file
    S_stds_files = unlist(strsplit(args$S_stds_files,','))
    S_means_files = unlist(strsplit(args$S_means_files,','))
    A_means_files = unlist(strsplit(args$A_means_files,','))
    A_stds_files = unlist(strsplit(args$A_stds_files,','))
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load data
    # gene_annotation = read_tsv(gene_annotation_file) %>% dplyr::rename(b_number=`b-number`)
    performance_evaluation = read_tsv(performance_evaluation_file) %>% 
        filter(`function` != 'evaluate_performance') %>%
        group_by(algorithm) %>% 
        mutate(rel_timestamp = get_relative_timestamp(timestamp),
               `function` = paste0(`function`,'()'),
               algorithm = factor(algorithm, levels = ALGORITHMS))
    
    clustering_info = read_tsv(clustering_info_file) %>%
        mutate(
            cluster_id = paste0('comp',cluster_id),
            algorithm = factor(algorithm, levels = ALGORITHMS)
        )
    
    S_means = load_mats(S_means_files, 'log-TPM', 'gene', 'weight_mean')
    S_stds = load_mats(S_stds_files, 'log-TPM', 'gene', 'weight_std')
    
    A_means = load_mats(A_means_files, 'index', 'sample', 'weight_mean')
    A_stds = load_mats(A_stds_files, 'index', 'sample', 'weight_std')
    
    ## process clustering info
    clustering_evaluation = clustering_info %>%
        group_by(algorithm, cluster_id) %>%
        summarise(
            silhouette_pearson = mean(silhouette_pearson),
            silhouette_euclidean = mean(silhouette_euclidean)
        ) %>% 
        dplyr::rename(component=cluster_id) %>% ungroup() 
    
    ## combine clustering_info, means and stds
    S_info = S_means %>%
        left_join(S_stds, by = c('algorithm','component','gene')) %>%
        left_join(clustering_evaluation, by = c('algorithm','component'))
    
    A_info = A_means %>%
        left_join(A_stds, by = c('algorithm','component','sample'))
    
    # analysis
    module_comparisons = make_module_comparisons(S_info)
    figdata = make_figdata(performance_evaluation, clustering_evaluation, 
                           S_info, A_info, module_comparisons)
    
    # visualize
    plts = make_plots(performance_evaluation, S_info, module_comparisons)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata[1], figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}