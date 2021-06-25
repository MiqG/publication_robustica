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

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables

# Development
# -----------
# RESULTS_DIR = file.path(ROOT,'results','benchmark_performance','files')
# performance_evaluation_file = file.path(RESULTS_DIR,'Sastry2019','performance_evaluation.tsv.gz')
# clustering_info_file = file.path(RESULTS_DIR,'Sastry2019','clustering_info.tsv.gz')
# figs_dir = file.path(ROOT,'results','benchmark_performance','figures','Sastry2019')

##### FUNCTIONS #####
get_relative_timestamp = function(x){
    x0=x[1]
    return(x-x0)
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


plot_S_correlations = function(Ss){
    plts = list()
    
    # are all components found independent?
    for (algorithm in names(Ss)){
        corr = cor(Ss[[algorithm]])
        plt_name = paste0('self_corr-S-',algorithm)
        plts[[plt_name]] = pheatmap(
            corr, show_colnames = FALSE, show_rownames = FALSE, 
            silent=TRUE, main=plt_name
        )
    }
    
    # Do we get the same output?
    corr = cor(Ss[['icasso']],Ss[['robustica_pca']])
    plts[['corr-S-icasso_vs_robustica_pca']] = pheatmap(
        corr, show_colnames = FALSE, show_rownames = FALSE, 
        main='icasso vs robustica_pca', silent=TRUE
    )
    
    corr = cor(Ss[['robustica']],Ss[['robustica_pca']])
    plts[['corr-S-robustica_vs_robustica_pca']] = pheatmap(
        corr, show_colnames = FALSE, show_rownames = FALSE, 
        main='robustica vs robutica_pca', silent=TRUE
    )
    
    return(plts)
}


plot_silhouettes = function(df){
    # get mean silhouette per cluster
    X = df %>%
        group_by(algorithm, cluster_id) %>%
        summarise(
            `1 - abs(Pearson)` = mean(silhouette_pearson),
            `Euclidean` = mean(silhouette_euclidean)
        ) %>%
        melt(id.vars = c('algorithm','cluster_id'))
    
    palette = get_palette('Paired',length(unique(X[['algorithm']])))
    
    plts = list()
    plts[['clustering-silhouettes']] = X %>%
        ggviolin(x='variable', y='value', fill='algorithm', 
                 color=NA, palette=palette) + 
        geom_boxplot(aes(fill=algorithm), width=.1, alpha=1, 
                     position = position_dodge(0.8)) + 
        labs(x='Metric', y='Silhouette', fill='Algorithm')
    
    return(plts)
}


make_plots = function(performance_evaluation, clustering_info){
    plts = list(
        plot_performance_profile(performance_evaluation),
        plot_silhouettes(clustering_info)
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
    save_plot(plts[['mem_time-scatter']],'mem_time-scatter','.png',figs_dir, width=12,height=20)
    
    save_plot(plts[['clustering-silhouettes']],'clustering-silhouettes','.pdf',figs_dir, width=15,height=12)
}


main = function(){
    args = getParsedArgs()

    performance_evaluation_file = args$performance_evaluation_file
    clustering_info_file = args$clustering_info_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load data
    performance_evaluation = read_tsv(performance_evaluation_file) %>% 
        filter(`function` != 'evaluate_performance') %>%
        group_by(algorithm) %>% 
        mutate(rel_timestamp = get_relative_timestamp(timestamp),
               `function` = paste0(`function`,'()')
              )
    clustering_info = read_tsv(clustering_info_file)
    
    plts = make_plots(performance_evaluation, clustering_info)
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}