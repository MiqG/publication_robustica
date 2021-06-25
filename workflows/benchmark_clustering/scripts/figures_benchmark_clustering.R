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


plot_silhouettes = function(df){
    # get mean silhouette per cluster
    X = df %>%
        group_by(property_oi, cluster_id) %>%
        summarise(
            `1 - abs(Pearson)` = mean(silhouette_pearson),
            `Euclidean` = mean(silhouette_euclidean)
        ) %>%
        melt(id.vars = c('property_oi','cluster_id'))
    
    palette = get_palette('Paired',length(unique(X[['property_oi']])))
    
    plt = X %>%
        ggviolin(x='variable', y='value', fill='property_oi', 
                 color=NA, palette=palette) + 
        geom_boxplot(aes(fill=property_oi), width=.1, alpha=1, 
                     position = position_dodge(0.8)) + 
        labs(x='Metric', y='Silhouette')
    
    return(plt)
}


plot_silhouette_time = function(df){
    X = df %>%
        group_by(property_oi, time, cluster_id) %>%
        summarise(
            silhouette_pearson = mean(silhouette_pearson),
            silhouette_euclidean = mean(silhouette_euclidean)
        ) %>%
        group_by(property_oi, time) %>%
        summarise(
            silhouette_euclidean = mean(silhouette_euclidean),
            silhouette_pearson = mean(silhouette_pearson)
        ) %>%
        melt(id.vars = c('property_oi','time'))
    
    palette = get_palette('Paired',length(unique(X[['property_oi']])))
    
    plt = X %>%
        ggscatter(x='time', y='value', color='property_oi', 
                  facet.by='variable', palette=palette) + 
        geom_text_repel(aes(label=property_oi)) + 
        guides(color=FALSE) +
        labs(x = 'Time (s)', y = 'Silhouette')
}


make_plots = function(
    metrics_performance, metrics_clustering, 
    methods_performance, methods_clustering
){
    plts = list(
        'metrics-mem_time-scatter' = plot_performance_profile(metrics_performance),
        'methods-mem_time-scatter' = plot_performance_profile(methods_performance),
        'metrics-clustering_silhouettes' = plot_silhouettes(metrics_clustering) + labs(fill='Clustering Metric'),
        'methods-clustering_silhouettes' = plot_silhouettes(methods_clustering) + labs(fill='Clustering Method'),
        'methods-clustering_silhouettes_time' = plot_silhouette_time(methods_clustering)
    )
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
    
    save_plot(plts[['metrics-mem_time-scatter']],'metrics-mem_time-scatter','.png',figs_dir, width=12,height=20)
    save_plot(plts[['methods-mem_time-scatter']] + facet_wrap(~property_oi, ncol=2, scales='free'),'methods-mem_time-scatter','.png',figs_dir, width=20,height=20)
    
    save_plot(plts[['metrics-clustering_silhouettes']],'metrics-clustering_silhouettes','.pdf',figs_dir, width=15,height=12)
    save_plot(plts[['methods-clustering_silhouettes']]+guides(fill=guide_legend(ncol=3)),'methods-clustering_silhouettes','.pdf',figs_dir, width=20,height=12)
    
    save_plot(plts[['methods-clustering_silhouettes_time']], 'methods-clustering_silhouettes_time','.pdf',figs_dir, width=16, height=9)
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
    
    clustering_time = metrics_performance %>% 
        filter(`function`=='cluster_components()') %>% 
        dplyr::select(time,property_oi) %>% 
        distinct()
    metrics_clustering = read_tsv(metrics_clustering_file) %>%
        left_join(clustering_time, by='property_oi')
    clustering_time = methods_performance %>% 
        filter(`function`=='cluster_components()') %>% 
        dplyr::select(time,property_oi) %>% 
        distinct()
    methods_clustering = read_tsv(methods_clustering_file) %>%
        left_join(clustering_time, by='property_oi')
    
    plts = make_plots(metrics_performance, metrics_clustering, 
                      methods_performance, methods_clustering)
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
