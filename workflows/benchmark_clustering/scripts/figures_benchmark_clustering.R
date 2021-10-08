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
require(writexl)
require(extrafont)
require(ggbreak)

loadfonts()

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables

# Development
# -----------
# RESULTS_DIR = file.path(ROOT,'results','benchmark_clustering','files')
# performance_file = file.path(RESULTS_DIR,'Sastry2019','methods','performance_evaluation_merged.tsv.gz')
# clustering_file = file.path(RESULTS_DIR,'Sastry2019','methods','clustering_info_merged.tsv.gz')
# figs_dir = file.path(ROOT,'results','benchmark_clustering','figures','Sastry2019','methods')


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


plot_silhouettes = function(df, lab='', labsize=0.1){
    X = df %>% 
        dplyr::select(cluster_id, property_oi, time, 
                      max_memory, silhouette_euclidean) %>%
        melt(id.vars = c('cluster_id','property_oi','time','max_memory')) %>%
        mutate(time=as.numeric(time)) %>%
        group_by(property_oi, time, max_memory, cluster_id) %>%
        summarize(value=mean(value)) # keep the mean silhouette per cluster
    
    palette = get_palette('Paired',length(unique(X[['property_oi']])))
    
    plts = list()    
    
    # silhouettes vs time
    plts[['silhouettes_vs_time']] = X %>%
        ggplot(aes(x=time, y=value, fill=property_oi)) + 
        geom_boxplot(outlier.size = 0.1, width=100) + 
        theme_pubr() + 
        labs(x='Time (s)', y='Silhouette Score') +
        theme(axis.text = element_text(size=labsize, family='Arial')) +
        #scale_x_break(c(200,500)) +
        guides(fill='none') +
        fill_palette(palette)

    
    # silhouettes vs memory
    plts[['silhouettes_vs_max_memory']] = X %>%
        ggplot(aes(x=max_memory, y=value, fill=property_oi)) + 
        geom_boxplot(outlier.size = 0.1, width=100) + 
        theme_pubr() + 
        labs(x='Max. Memory (MiB)', y='Silhouette Score') +
        theme(text = element_text(size=labsize,family='Arial')) +
        #scale_x_break(c(4500,5500)) +
        guides(fill='none') +
        fill_palette(palette)
    
    # memory vs time vs silhouettes
    plts[['silhouettes_vs_max_memory_vs_time']] = X %>% 
        group_by(property_oi, time, max_memory) %>% 
        summarize(silhouette = median(value)) %>% 
        ggscatter(x='time', y='max_memory', size='silhouette', 
                  repel = TRUE, color='property_oi', label='property_oi', 
                  palette=palette) +
        guides(color='none') +
        labs(x='Time (s)', y='Max. Memory (MiB)', size='Silhouette Score') +
        theme(text = element_text(size=labsize,family='Arial'))

    names(plts) = paste0(lab, names(plts))
    return(plts)
}


make_plots = function(performance, clustering){
    plts = list(
        'mem_time-scatter' = plot_performance_profile(performance)
    )
    plts = list(
        plts,
        plot_silhouettes(clustering)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(performance, clustering){
    summary_clustering = clustering %>% 
        group_by(property_oi, time, max_memory, cluster_id) %>% 
        summarize(mean_silhouette = mean(silhouette_euclidean)) %>%
        group_by(property_oi, time, max_memory) %>%
        summarize(
            mean=mean(mean_silhouette), 
            median=median(mean_silhouette), 
            std=sd(mean_silhouette), 
            range=max(mean_silhouette) - min(mean_silhouette)
        )
    
    figdata = list(
        'benchmark_clustering-evaluation' = list(
            'performance_evaluation' = performance,
            'clustering_evaluation' = clustering,
            'clustering_evaluation_summary' = summary_clustering
        )
    )
    return(figdata)
}


save_plot = function(plt, plt_name, extension='.pdf', 
                     directory='', 
                     dpi=350, 
                     width = par("din")[1], 
                     height = par("din")[2], 
                     device = NULL,
                     units='cm'){

    plt = ggpar(plt, font.title=11, font.subtitle=10, font.caption=10, 
                font.x=10, font.y=10, font.legend=10,
                font.tickslab=8, font.family='Arial')
    filename = file.path(directory,paste0(plt_name,extension))
    ggsave(filename, 
           plt, 
           width=width, height=height, dpi=dpi, 
           limitsize=FALSE,
           units=units, device=device)
}


save_plots = function(plts, figs_dir){
    save_plot(plts[['mem_time-scatter']] + facet_wrap(~property_oi, ncol=2, scales='free'),'mem_time-scatter','.png',figs_dir, width=12, height=12)
    save_plot(plts[['silhouettes_vs_time']], 'silhouettes_vs_time','.pdf',figs_dir, width=6, height=6, device=cairo_pdf)
    save_plot(plts[['silhouettes_vs_max_memory']], 'silhouettes_vs_max_memory','.pdf',figs_dir, width=6, height=6, device=cairo_pdf)
    save_plot(plts[['silhouettes_vs_max_memory_vs_time']], 'silhouettes_vs_max_memory_vs_time','.pdf',figs_dir, width=6, height=6, device=cairo_pdf)
    
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

    performance_file = args$performance_file
    clustering_file = args$clustering_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load data
    ## performances
    performance = read_tsv(performance_file) %>% 
        filter(`function` != 'evaluate_performance') %>%
        group_by(property_oi) %>% 
        mutate(rel_timestamp = get_relative_timestamp(timestamp),
               `function` = paste0(`function`,'()'))
    ## clustering time and evaluation    
    clustering_time = performance %>% 
        filter(`function`=='cluster_components()') %>% 
        dplyr::select(time,property_oi) %>% 
        distinct()
    max_memory = performance %>% 
        group_by(property_oi) %>% 
        summarize(max_memory=max(memory))
    clustering = read_tsv(clustering_file) %>%
        left_join(clustering_time, by='property_oi') %>%
        left_join(max_memory, by='property_oi')
    
    plts = make_plots(performance, clustering)
    figdata = make_figdata(performance, clustering)
    
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
