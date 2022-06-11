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
require(cowplot)
require(reshape2)
require(ggrepel)
require(writexl)
require(extrafont)
require(scattermore)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
FONT_FAMILY = "Arial"

# Development
# -----------
# RESULTS_DIR = file.path(ROOT,'results','benchmark_clustering','files')
# performance_file = file.path(RESULTS_DIR,'Sastry2019','methods','performance_evaluation_merged.tsv.gz')
# clustering_file = file.path(RESULTS_DIR,'Sastry2019','methods','clustering_info_merged.tsv.gz')
# pca_components_file = file.path(RESULTS_DIR,'Sastry2019','pca','components.tsv.gz')
# figs_dir = file.path(ROOT,'results','benchmark_clustering','figures','Sastry2019','methods')


##### FUNCTIONS #####
get_relative_timestamp = function(x){
    x0=x[1]
    return(x-x0)
}


plot_performance_profile = function(performance){
    X = performance
    
    plts = list()
    
    plts[["mem_time-scatter"]] = X %>% 
        ggplot(aes(x=rel_timestamp, y=memory)) +
        geom_scattermore(aes(color=`function`), pointsize=8, 
                         pixels=c(1000,1000), shape=20) +
        geom_line(aes(color=`function`), linetype="dashed")+
        color_palette("uchicago") +
        facet_wrap(~property_oi, ncol=2, scales='free') +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color=guide_legend(ncol=1)) + 
        labs(x='Time (s)', y='Memory (MiB)', color='Function') + 
        theme_pubr(border = TRUE)
    
    return(plts)
}


plot_silhouettes = function(df, lab='', labsize=0.1){
    X = df %>% 
        dplyr::select(cluster_id, property_oi, time, 
                      max_memory, silhouette_euclidean) %>%
        melt(id.vars = c('cluster_id','property_oi','time','max_memory')) %>%
        mutate(time=as.numeric(time)) %>%
        group_by(property_oi, time, max_memory, cluster_id) %>%
        summarize(value=mean(value)) %>% # keep the mean silhouette per cluster
        ungroup()
    
    palette = get_palette('Paired',length(unique(X[['property_oi']])))
    
    plts = list()    
    
    # silhouettes vs time
    plts[['silhouettes_vs_time']] = X %>%
        ggplot(aes(x=time, y=value, fill=property_oi)) + 
        geom_boxplot(outlier.size = 0.1, width=100) + 
        fill_palette(palette) +
        theme_pubr() +
        labs(x='Time (s)', y='Silhouette Score')

    # silhouettes vs memory
    plts[['silhouettes_vs_max_memory']] = X %>%
        ggplot(aes(x=max_memory, y=value, fill=property_oi)) + 
        geom_boxplot(outlier.size = 0.1, width=100) + 
        theme_pubr() + 
        labs(x='Max. Memory (MiB)', y='Silhouette Score') +
        fill_palette(palette)
    
    # memory vs time vs silhouettes
    plts[['silhouettes_vs_max_memory_vs_time']] = X %>% 
        group_by(property_oi, time, max_memory) %>% 
        summarize(silhouette = median(value)) %>% 
        ggscatter(x='time', y='max_memory', size='silhouette', 
                  repel = TRUE, color='property_oi', label='property_oi', 
                  palette=palette) +
        labs(x='Time (s)', y='Max. Memory (MiB)', size='Silhouette Score') +
        guides(color="none")

    names(plts) = paste0(lab, names(plts))
    return(plts)
}


plot_pca = function(pca_components, clustering){
    plts = list()
    
    X = clustering %>% 
        distinct(component, cluster_id, property_oi) %>% 
        mutate(component = paste0("comp",component),
               cluster_id = paste0("clust",cluster_id),
               is_noise = cluster_id == "clust-1") %>%
        left_join(pca_components, 
                  by="component") %>%
        arrange(component)
    
    plts[["clustering_pca-overview"]] = X %>% 
        distinct(PC1, PC2) %>%
        ggplot(aes(x=PC1, y=PC2)) + 
        geom_scattermore(pointsize=3, alpha=0.5, pixels=c(1000,1000)) + 
        theme_pubr()
    
    # highlight clusters
    for (property in unique(X[["property_oi"]])) {
        x = X %>%
            filter(property_oi == property)
        
        palette = get_palette("Paired", length(unique(x[["cluster_id"]])))
        is_noise = x[["cluster_id"]] == "clust-1"
        if (any(is_noise)) {palette[1] = "black"}
        
        clean_clusters = setdiff(unique(x[["cluster_id"]]),"clust-1")
        
        plts[[sprintf("clustering_pca-%s",property)]] = x %>%
            ggplot(aes(x=PC1, y=PC2)) + 
            geom_scattermore(aes(color=cluster_id), pointsize=3, 
                             alpha=0.8, pixels=c(1000,1000)) + 
            color_palette(palette) + 
            guides(color="none") +
            theme(aspect.ratio = 1) +
            theme_pubr() +
            labs(title=sprintf("%s | n. clusters=%s | n. noisy=%s", 
                               property, length(clean_clusters), sum(is_noise)))
        
    }
    
    return(plts)
}

make_plots = function(performance, clustering, pca_components){
    plts = list(
        plot_performance_profile(performance),
        plot_silhouettes(clustering),
        plot_pca(pca_components, clustering)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(performance, clustering, pca_components){
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
            'clustering_evaluation_summary' = summary_clustering,
            'pca_components' = pca_components
        )
    )
    return(figdata)
}


save_plt = function(plts, plt_name, extension=".pdf", 
                    directory="", dpi=350, format=TRUE,
                    width = par("din")[1], height = par("din")[2]){
    print(plt_name)
    plt = plts[[plt_name]]
    if (format){
        plt = ggpar(plt, font.title=8, font.subtitle=8, font.caption=8, 
                    font.x=8, font.y=8, font.legend=6,
                    font.tickslab=6, font.family=FONT_FAMILY, device=cairo_pdf)
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units="cm")
}


save_plots = function(plts, figs_dir){
    save_plt(plts, 'mem_time-scatter','.pdf',figs_dir, width=12, height=12)
    save_plt(plts, 'silhouettes_vs_time','.pdf',figs_dir, width=6, height=6)
    save_plt(plts, 'silhouettes_vs_max_memory','.pdf',figs_dir, width=6, height=6)
    save_plt(plts, 'silhouettes_vs_max_memory_vs_time','.pdf',figs_dir, width=6, height=6)
    
    save_plt(plts, 'clustering_pca-overview','.pdf',figs_dir, width=5, height=5)
    save_plt(plts, 'clustering_pca-AgglomerativeClustering','.pdf',figs_dir, width=5, height=5)
    save_plt(plts, 'clustering_pca-DBSCAN','.pdf',figs_dir, width=5, height=5)
    save_plt(plts, 'clustering_pca-KMedoids','.pdf',figs_dir, width=5, height=5)
    save_plt(plts, 'clustering_pca-CommonNNClustering','.pdf',figs_dir, width=5, height=5)
    save_plt(plts, 'clustering_pca-OPTICS','.pdf',figs_dir, width=5, height=5)
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
    performance = read_tsv(performance_file)
    clustering = read_tsv(clustering_file)
    pca_components = read_tsv(pca_components_file)
    
    ## keep those algorithms that we could evaluate
    avail_algorithms = intersect(
        performance[["property_oi"]], clustering[["property_oi"]]
    )
    
    performance = performance %>% filter(property_oi %in% avail_algorithms)
    clustering = clustering %>% filter(property_oi %in% avail_algorithms)
    
    ## performances
    performance = performance %>% 
        filter(`function` != 'evaluate_performance') %>%
        group_by(property_oi) %>% 
        mutate(rel_timestamp = timestamp - min(timestamp),
               `function` = paste0(`function`,'()')) %>%
        ungroup()
    
    ## clustering time and evaluation    
    clustering_time = performance %>% 
        filter(`function`=='cluster_components()') %>% 
        group_by(property_oi) %>%
        summarize(time=max(time)) %>%
        distinct()
    
    max_memory = performance %>% 
        group_by(property_oi) %>% 
        summarize(max_memory=max(memory))
    
    clustering = clustering %>%
        left_join(clustering_time, by='property_oi') %>%
        left_join(max_memory, by='property_oi')
    
    ## prepare pca
    pca_components = pca_components %>% 
      mutate(component = row_number() - 1,
             component = paste0("comp",component))
    
    plts = make_plots(performance, clustering, pca_components)
    figdata = make_figdata(performance, clustering, pca_components)
    
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
