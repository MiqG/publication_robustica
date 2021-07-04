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

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables

# Development
# -----------
RESULTS_DIR = file.path(ROOT,'results','benchmark_clustering','files')
metrics_performance_file = file.path(RESULTS_DIR,'metrics','Sastry2019','performance_evaluation_merged.tsv.gz')
metrics_clustering_file = file.path(RESULTS_DIR,'metrics','Sastry2019','clustering_info_merged.tsv.gz')
methods_performance_file = file.path(RESULTS_DIR,'methods','Sastry2019','performance_evaluation_merged.tsv.gz')
methods_clustering_file = file.path(RESULTS_DIR,'methods','Sastry2019','clustering_info_merged.tsv.gz')
figs_dir = file.path(ROOT,'results','benchmark_clustering','figures','Sastry2019')

##### FUNCTIONS #####
get_relative_timestamp = function(x){
    x0=x[1]
    return(x-x0)
}


plot_performance_profile = function(df){
    X = df
    
    plt = X %>% 
        ggline(x='rel_timestamp', y='memory',numeric.x.axis = TRUE, 
               color='function', linetype='dashed', 
               shape=20, point.size=0.0001, 
               palette='uchicago',alpha=0.5) + 
    facet_wrap(~property_oi, scales='free_x', ncol=1) + 
    guides(color=guide_legend(ncol=1)) + 
    labs(x='Time (s)', y='Memory (MiB)', color='Function') + 
    theme_pubr(border = TRUE) 
    
    return(plt)
}


plot_silhouettes = function(df, lab){
    X = df %>% 
        dplyr::select(cluster_id, property_oi, time, 
                      max_memory, silhouette_pearson) %>%
        melt(id.vars = c('cluster_id','property_oi','time','max_memory')) %>%
        mutate(time=as.numeric(time))
    
    palette = get_palette('Paired',length(unique(X[['property_oi']])))
    
    plts = list()    
    
    # silhouettes vs time
    med = X %>% 
        group_by(property_oi, time, max_memory) %>%
        summarize(median_value=median(value))
    plts[['silhouettes_vs_time']] = ggplot(
        X, aes(x=time, y=value, fill=NULL, color=property_oi)
        ) + 
        geom_boxplot(outlier.alpha = 0, width=20) + 
        geom_point(aes(x=time, y=median_value, label=property_oi), med, size=0.5) +
        theme_pubr() + 
        geom_text_repel(
            aes(x=time, y=median_value, label=property_oi), med) + 
        guides(color=FALSE) +
        labs(x='Time (s)', y='Silhouette Score')
    
    # silhouettes vs memory
    med = X %>% 
        group_by(property_oi, time, max_memory) %>%
        summarize(median_value=median(value))
    plts[['silhouettes_vs_max_memory']] = ggplot(
        X, aes(x=max_memory, y=value, fill=NULL, color=property_oi)
        ) + 
        geom_boxplot(outlier.alpha = 0, width=20) + 
        geom_point(aes(x=max_memory, y=median_value, label=property_oi), med, size=0.5) +
        theme_pubr() + 
        geom_text_repel(
            aes(x=max_memory, y=median_value, label=property_oi), med) +
        guides(color=FALSE) +
        labs(x='Memory (MiB)', y='Silhouette Score')
    
    plts = sapply(plts, function(plt){ set_palette(plt, palette) }, simplify=FALSE)
    
    names(plts) = paste0(lab, names(plts))
    
    return(plts)
}


make_plots = function(
    metrics_performance, metrics_clustering, 
    methods_performance, methods_clustering
){
    plts = list(
        'metrics-mem_time-scatter' = plot_performance_profile(metrics_performance),
        'methods-mem_time-scatter' = plot_performance_profile(methods_performance)
    )
    plts = list(
        plts,
        plot_silhouettes(methods_clustering,'methods-'),
        plot_silhouettes(metrics_clustering,'metrics-')
    )
    plts = do.call(c,plts)
    return(plts)
}


save_plot = function(plt, plt_name, extension='.pdf', 
                     directory='', 
                     dpi=350, 
                     width = par("din")[1], 
                     height = par("din")[2], 
                     device = NULL,
                     units='cm'){
        filename = file.path(directory,paste0(plt_name,extension))
        ggsave(filename, 
               plt, 
               width=width, height=height, dpi=dpi, 
               limitsize=FALSE, useDingbats = FALSE,
               units=units, device=device)
}


save_plots = function(plts, figs_dir){
    
    save_plot(plts[['metrics-mem_time-scatter']],'metrics-mem_time-scatter','.png',figs_dir, width=12,height=20)
    save_plot(plts[['methods-mem_time-scatter']] + facet_wrap(~property_oi, ncol=2, scales='free'),'methods-mem_time-scatter','.png',figs_dir, width=20,height=20)
    
    save_plot(plts[['methods-silhouettes_vs_time']], 'methods-silhouettes_vs_time','.pdf',figs_dir, width=12, height=12, device=cairo_pdf)
    save_plot(plts[['methods-silhouettes_vs_max_memory']], 'methods-silhouettes_vs_max_memory','.pdf',figs_dir, width=12, height=12, device=cairo_pdf)
    
    save_plot(plts[['metrics-silhouettes_vs_time']], 'metrics-silhouettes_vs_time','.pdf',figs_dir, width=12, height=12, device=cairo_pdf)
    save_plot(plts[['metrics-silhouettes_vs_max_memory']], 'metrics-silhouettes_vs_max_memory','.pdf',figs_dir, width=12, height=12, device=cairo_pdf)
}


main = function(){
    args = getParsedArgs()

    metrics_performance_file = args$metrics_performance_file
    metrics_clustering_file = args$metrics_clustering_file
    methods_performance_file = args$methods_performance_file
    methods_clustering_file = args$methods_clustering_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load data
    ## performances
    metrics_performance = read_tsv(metrics_performance_file) %>% 
        filter(`function` != 'evaluate_performance') %>%
        group_by(property_oi) %>% 
        mutate(rel_timestamp = get_relative_timestamp(timestamp),
               `function` = paste0(`function`,'()')
              )
    methods_performance = read_tsv(methods_performance_file) %>% 
        filter(`function` != 'evaluate_performance') %>%
        group_by(property_oi) %>% 
        mutate(rel_timestamp = get_relative_timestamp(timestamp),
               `function` = paste0(`function`,'()')
              )
    ## clustering time and evaluation
    ### metrics
    clustering_time = metrics_performance %>% 
        filter(`function`=='cluster_components()') %>% 
        dplyr::select(time,property_oi) %>% 
        distinct()
    max_memory = metrics_performance %>% 
        group_by(property_oi) %>% 
        summarize(max_memory=max(memory))
    metrics_clustering = read_tsv(metrics_clustering_file) %>%
        left_join(clustering_time, by='property_oi') %>%
        left_join(max_memory, by='property_oi')
    
    ### methods
    clustering_time = methods_performance %>% 
        filter(`function`=='cluster_components()') %>% 
        dplyr::select(time,property_oi) %>% 
        distinct()
    max_memory = methods_performance %>% 
        group_by(property_oi) %>% 
        summarize(max_memory=max(memory))
    methods_clustering = read_tsv(methods_clustering_file) %>%
        left_join(clustering_time, by='property_oi') %>%
        left_join(max_memory, by='property_oi')
    
    plts = make_plots(metrics_performance, metrics_clustering, 
                      methods_performance, methods_clustering)
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
