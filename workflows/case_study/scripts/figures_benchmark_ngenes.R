#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
#
# Outline
# -------


require(tidyverse)
require(ggpubr)
require(ComplexHeatmap)
require(circlize)
require(ggplotify)
require(writexl)
require(extrafont)
require(fdrtool)
require(proxy)

loadfonts()

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
PARAMS = list(
    'LGG' = list(
        'n_genes_oi' = seq(500,30000,1000),
        'thresh_silh' = 0.9,
        'dodge_val' = 500
    ),
    'Sastry2019' = list(
        'n_genes_oi' = seq(200,3900,200),
        'thresh_silh' = 0.9,
        'dodge_val' = 100
    )
)


# Development
# -----------
RESULTS_DIR = file.path(ROOT,'results','case_study')
dataset = 'LGG'
stats_file = file.path(RESULTS_DIR,'files','benchmark_ica_ngenes',dataset,'stats_merged.tsv.gz')
S_icasso_file = file.path(RESULTS_DIR,'files','benchmark_ica_ngenes',dataset,'20224','icasso','S.tsv.gz')
S_robustica_file = file.path(RESULTS_DIR,'files','benchmark_ica_ngenes',dataset,'20224','robustica_pca','S.tsv.gz')
figs_dir = file.path(RESULTS_DIR,'figures','benchmark_ngenes',dataset)


##### FUNCTIONS #####
plot_silhouettes = function(stats, params){
    # unpack params
    n_genes_oi = params[['n_genes_oi']]
    thresh_silh = params[['thresh_silh']]
    dodge_val = params[['dodge_val']]
    
    X = stats %>% 
        filter(n_genes %in% n_genes_oi | n_genes==max(n_genes))
    
    plts = list()
    
    # how do silhouette scores distributions change with different numbers of genes
    meds = X %>% 
        group_by(algorithm,n_genes) %>%
        summarize(
            silh_max=quantile(silhouette_euclidean,0.75,na.rm=TRUE),
            silh_min=quantile(silhouette_euclidean,0.25,na.rm=TRUE),
            silhouette_euclidean=median(silhouette_euclidean,na.rm=TRUE))
    
    palette = get_palette('Paired', 4)[c(1,4)]
    dodge = position_dodge(dodge_val)
    plts[['silhouette-violins']] = X %>% 
        ggplot(aes(x=n_genes, y=silhouette_euclidean,
                   group=interaction(n_genes,algorithm))) + 
        geom_violin(aes(fill=algorithm), position=dodge, color=NA, trim=TRUE) +
        fill_palette(palette) +
        geom_line(data=meds, linetype='dashed', position=dodge) +
        geom_point(data=meds, position=dodge, size=0.5) +
        geom_errorbar(data=meds, 
                      mapping=aes(ymin=silh_min, ymax=silh_max), 
                      position=dodge, width=dodge_val/2) +
        geom_hline(yintercept=thresh_silh, linetype='dashed', size=0.1) +
        #facet_wrap(~algorithm, ncol=1) +
        theme_pubr(border=TRUE) +
        labs(x='No. Genes Randomly Picked', y='Silhouette Score') +
        ylim(NA,1)
    
    # how many confident components?
    plts[['silhouette-confident_components']] = X %>%
        mutate(is_confident=silhouette_euclidean>thresh_silh) %>%
        group_by(algorithm, n_genes, is_confident) %>%
        summarize(n=n()) %>%
        filter(is_confident) %>%
        ggplot(aes(x=n_genes, y=n, fill=algorithm)) +
        geom_bar(stat='identity', position=position_dodge(dodge_val)) +
        fill_palette(palette) +
        labs(x='No. Genes Randomly Picked', 
             y=sprintf('Count Silh. Score > %s',thresh_silh)) +
        theme_pubr()
    
    return(plts)
}


make_module_comparisons = function(S_info){
    # get modules from source matrix
    modules = S_info %>%
        group_by(algorithm, component) %>%
        mutate(in_module=define_module(weight_mean))
    
    
    # get modules
    S = S_robustica
    is_selected = S %>% mutate_at(vars(-('sample')), define_module)
    modules = S %>% dplyr::rename(gene=`sample`)
    modules[,-1][!is_selected[,-1]] = NA
    robustica = modules[,-1] %>% is.na() 
    
    S = S_icasso
    is_selected = S %>% mutate_at(vars(-('sample')), define_module)
    modules = S %>% dplyr::rename(gene=`sample`)
    modules[,-1][!is_selected[,-1]] = NA
    icasso = modules[,-1] %>% is.na() 
    
    # module sizes
    module_size = rbind(
        robustica %>% colSums() %>% 
            enframe('component','n') %>% mutate(algorithm='robustica_pca'),
        icasso %>% colSums() %>% 
            enframe('component','n') %>% mutate(algorithm='icasso')
    )
    
    # use jaccard similarity to map modules
#     icasso = modules %>% 
#         filter(algorithm == 'icasso') %>% 
#         pivot_wider(id_cols = gene, names_from = component, values_from = in_module)
#     robustica_pca = modules %>% 
#         filter(algorithm == 'robustica_pca') %>% 
#         pivot_wider(id_cols = gene, names_from = component, values_from = in_module)
    sim = simil(icasso[,-1], robustica[,-1], 
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
    
    module_comparisons = list(
        'sim' = sim,
        'mapping' = mapping
    )
    return(module_comparisons)
}


plot_module_comparisons = function(module_comparisons){
    # unlist data
    sim = module_comparisons[['sim']]
    mapping = module_comparisons[['mapping']] %>%
        mutate(is_low=jaccard<THRESH_JACCARD)

    
    # make plots
    plts = list()
    
    plts[['module_comparisons-jaccard']] = sim %>% pheatmap(silent=TRUE, cluster_cols = FALSE, cluster_rows=FALSE, fontsize = 5)
    
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
                filter(abs_module_size_diff > 15) %>% 
                pivot_longer(cols = c(icasso,robustica_pca)),
            box.padding = 0.7,
            family="Arial",
            size=2
        ) +
        scale_color_manual("Algorithm", values=c('#A6CEE3','#33A02C')) +
        labs(x='Module Size Icasso', 
             y='Module Size robustica', 
             color='Jaccard Sim.',
             title=TeX(sprintf('Mapped gene modules between \\textit{Icasso} and \\textit{robustica} | n=%s',sum(!is.na(mapping[['icasso']])))))
    
    
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
        drop_na(jaccard, is_low) %>%
        ggscatter(x='robustica_pca_module_size', y='jaccard', 
                  color='is_low', palette='jco', add='reg.line', 
                  add.params=list(linetype='dashed',color='black'), 
                  conf.int = TRUE, size=1, alpha=0.5) + 
        stat_cor(method='spearman', label.x = 20) +
        ylim(NA,1)+
        labs(y='Jaccard', color=sprintf('Jaccard < %s',THRESH_JACCARD))
    
    plts[['module_comparisons-low_jaccard']] = mapping %>% 
        count(is_low) %>%
        drop_na() %>%
        ggbarplot(x='is_low', y='n', fill='is_low', 
                  color=NA, label=TRUE, palette='jco') +
        labs(x=sprintf('Jaccard < %s',THRESH_JACCARD), y='Count') +
        guides(fill=FALSE)
        
    return(plts)
}


plot_compare_components = function(S_robustica, S_icasso){
    
    col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    
    plts = list()
    # icasso does not find robust independent components
    corr = cor(S_icasso[,-1])
    plts[['comparison-corr_icasso']] = corr %>% 
        Heatmap(name = 'Correlation', col = col_fun, 
                column_names_gp = grid::gpar(fontsize = 4, fontfamily='Arial'),
                row_names_gp = grid::gpar(fontsize = 4, fontfamily='Arial'))
    # robustica does
    corr = cor(S_robustica[,-1])
    plts[['comparison-corr_robustica_pca']] = corr %>% 
        Heatmap(name = 'Correlation', col = col_fun,
                column_names_gp = grid::gpar(fontsize = 4, fontfamily='Arial'),
                row_names_gp = grid::gpar(fontsize = 4, fontfamily='Arial'))
    
    # 
    corr = cor(S_icasso[,-1],S_robustica[,-1])
    corr %>% 
        Heatmap(name='Correlation', col = col_fun,
                column_names_gp = grid::gpar(fontsize = 4, fontfamily='Arial'),
                row_names_gp = grid::gpar(fontsize = 4, fontfamily='Arial'))
    
    
    
    
    
    
    # map modules with each other
    
    
    
    
    return(plts)
}


define_module = function(x, cutoff=0.01){
    fdr = fdrtool(x, plot=FALSE, cutoff.method="fndr", verbose=FALSE)[['qval']]
    x = fdr < cutoff
    return(x)
}


make_plots = function(stats, S_robustica, S_icasso, params){
    plts = list(
        plot_silhouettes(stats, params),
        plot_compare_components(S_robustica, S_icasso)
    )
    plts = do.call(c,plts)
    return(plts)
}


save_plot = function(plt, plt_name, extension='.pdf', 
                      directory='', dpi=350, change_params=TRUE,
                      width = par("din")[1], height = par("din")[2], units='cm'){
    if (change_params){
        plt = ggpar(plt, font.title=11, font.subtitle=10, font.caption=10, 
                    font.x=10, font.y=10, font.legend=10,
                    font.tickslab=8, font.family='Arial')
    }
    
    filename = file.path(directory,paste0(plt_name,extension))
    ggsave(filename, 
           plt, 
           width=width, height=height, dpi=dpi, limitsize=FALSE, units=units)
}


save_plots = function(plts, figs_dir){
    save_plot(plts[['silhouette-violins']], 'silhouette-violins', '.pdf', figs_dir, width=8, height=5)
    save_plot(plts[['silhouette-confident_components']], 'silhouette-confident_components', '.pdf', figs_dir, width=8, height=5)
    
    save_plot(as.ggplot(grid.grabExpr(draw(plts[['comparison-corr_icasso']]))), 'comparison-corr_icasso', '.pdf', figs_dir, width=10, height=8, change_params=FALSE)
    save_plot(as.ggplot(grid.grabExpr(draw(plts[['comparison-corr_robustica_pca']]))), 'comparison-corr_robustica_pca', '.pdf', figs_dir, width=10, height=8, change_params=FALSE)
}


make_figdata = function(stats, S_robustica, S_icasso){
    figdata = list(
        'benchmark_ngenes' = list(
            'clustering_summary'= stats,
            'S_robustica' = S_robustica,
            'S_icasso' = S_icasso
        )
    )
    return(figdata)
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
    stats_file = args$stats_file
    S_icasso_file = args$S_icasso_file
    S_robustica_file = args$S_robustica_file
    dataset = args$dataset
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)

    # load data
    stats = read_tsv(stats_file)
    S_icasso = read_tsv(S_icasso_file)
    S_robustica = read_tsv(S_robustica_file)
    
    # make plots
    plts = make_plots(stats, S_robustica, S_icasso, PARAMS[[dataset]])
    
    # make figdata
    figdata = make_figdata(stats, S_robustica, S_icasso)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}