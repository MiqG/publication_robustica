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
require(ComplexHeatmap)
require(ggplotify)
require(gridExtra)
require(circlize)
require(proxy)
require(fdrtool)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
FONT_FAMILY = "Arial"
PAL_ALGOS = setNames(
    get_palette("Paired",6),
    c("AffinityPropagation", "AgglomerativeClustering",
      "CommonNNClustering", "DBSCAN", "KMedoids", "OPTICS"))

# Development
# -----------
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','benchmark_clustering','files')
# performance_file = file.path(RESULTS_DIR,'Sastry2019','methods','performance_evaluation_merged.tsv.gz')
# clustering_file = file.path(RESULTS_DIR,'Sastry2019','methods','clustering_info_merged.tsv.gz')
# pca_components_file = file.path(RESULTS_DIR,'Sastry2019','pca','components.tsv.gz')
# components_robustica_file = file.path(RESULTS_DIR,'Sastry2019','methods','DBSCAN','DBSCAN-S.tsv.gz')
# components_Sastry2019_file = file.path(PREP_DIR,'original_pipeline_Sastry2019','original','S.csv')
# figs_dir = file.path(ROOT,'results','benchmark_clustering','figures','Sastry2019','methods')


##### FUNCTIONS #####
define_module = function(x, cutoff=0.01){
    fdr = fdrtool(x, plot=FALSE, cutoff.method="fndr", verbose=FALSE)[['qval']]
    x = fdr < cutoff
    return(x)
}

plot_performance_profile = function(performance){
    X = performance
    
    plts = list()
    
    plts[["mem_time-scatter"]] = X %>% 
        ggplot(aes(x=rel_timestamp, y=memory)) +
        geom_scattermore(aes(color=`function`), pointsize=10, 
                         pixels=c(1000,500)) +
        geom_line(aes(color=`function`), linetype="dashed", size=0.25)+
        color_palette("uchicago") +
        facet_wrap(~property_oi, ncol=2, scales='free') +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color=guide_legend(ncol=1)) + 
        labs(x='Time (s)', y='Memory (MiB)', color='Function') + 
        theme_pubr(border = TRUE)
    
    return(plts)
}


plot_silhouettes = function(clustering, lab='', labsize=0.1){
    X = clustering %>% 
        dplyr::select(cluster_id, property_oi, time, 
                      max_memory, silhouette_euclidean) %>%
        melt(id.vars = c('cluster_id','property_oi','time','max_memory')) %>%
        mutate(time=as.numeric(time)) %>%
        group_by(property_oi, time, max_memory, cluster_id) %>%
        summarize(value=mean(value), # keep the mean silhouette per cluster
                  n = n()) %>% 
        ungroup()
    
    plts = list()    
    
    # silhouettes vs time
    plts[['silhouettes_vs_time']] = X %>%
        ggplot(aes(x=time, y=value, color=property_oi)) + 
        geom_boxplot(outlier.size = 0.1, width=100) + 
        color_palette(PAL_ALGOS) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x='Time (s)', y='Silhouette Score')

    # silhouettes vs memory
    plts[['silhouettes_vs_max_memory']] = X %>%
        ggplot(aes(x=max_memory, y=value, color=property_oi)) + 
        geom_boxplot(outlier.size = 0.1, width=100) + 
        theme_pubr() + 
        theme(aspect.ratio=1) +
        labs(x='Max. Memory (MiB)', y='Silhouette Score') +
        color_palette(PAL_ALGOS)
    
    # memory vs time vs silhouettes
    plts[['silhouettes_vs_max_memory_vs_time']] = X %>% 
        group_by(property_oi, time, max_memory) %>% 
        summarize(silhouette = median(value)) %>% 
        ggscatter(x='time', y='max_memory', size='silhouette', 
                  color='property_oi', palette=PAL_ALGOS) +
        geom_text_repel(aes(label=property_oi, color=property_oi), size=2) +
        labs(x='Time (s)', y='Max. Memory (MiB)', size='Silhouette Score') +
        theme(aspect.ratio=1) +
        guides(color="none")
    
    # silhouettes vs n components in cluster
    x = X %>%
        group_by(property_oi) %>%
        mutate(property_lab = paste0(property_oi," | n=",n())) %>%
        ungroup()
    plts[["silhouettes_vs_n_components"]] = x %>%
        ggplot(aes(x=n, y=value)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=10, alpha=0.5) +
        geom_text_repel(aes(label=cluster_id), x %>% filter(n>300)) +
        facet_wrap(~property_lab) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Cluster Size", y="Silhouette Score") +
        theme_pubr()
    
    
    names(plts) = paste0(lab, names(plts))
    return(plts)
}


plot_pca = function(pca_components, clustering){
    plts = list()
    
    X = clustering %>% 
        distinct(component, cluster_id, property_oi, silhouette_euclidean) %>% 
        group_by(property_oi, cluster_id) %>%
        mutate(
            cluster_size = n(),
            component = paste0("comp",component),
            cluster_id = paste0("clust",cluster_id),
            is_noise = cluster_id=="clust-1" | mean(silhouette_euclidean) < 0.5
        ) %>%
        ungroup() %>%
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
        
        cluster_ids = unique(x[["cluster_id"]])
        palette = setNames(get_palette("Paired", length(cluster_ids)), cluster_ids)
        is_noise = x %>% filter(is_noise) %>% pull(cluster_id) %>% unique()
        palette[is_noise] = "black"
        
        clean_clusters = setdiff(cluster_ids, is_noise)
        
        plts[[sprintf("clustering_pca-%s",property)]] = x %>%
            ggplot(aes(x=PC1, y=PC2)) + 
            geom_scattermore(aes(color=cluster_id), pointsize=3, 
                             alpha=0.8, pixels=c(1000,1000)) + 
            color_palette(palette) + 
            guides(color="none") +
            theme_pubr() +
            theme(aspect.ratio = 1) +
            labs(title=sprintf("%s | n. clusters=%s | n. low qual.=%s", 
                               property, length(clean_clusters), sum(x[["is_noise"]])))
        
    }
    
    return(plts)
}


plot_comparison_methods = function(comparison_pearson, comparison_jaccard){
    # Icasso algorithm run using DBSCAN with either robustica or 
    # Sastry2019's original implementation using same parameters
    
    plts = list()
    
    # do we find the same number of components of Sastry2019?
    X = data.frame(
        dataset = c("Sastry2019","robustica"),
        n_components = c(ncol(comparison_pearson), nrow(comparison_pearson))
    )
    plts[["comparison_methods-n_components"]] = X %>%
        ggbarplot(x="dataset", y="n_components", fill="orange", color=NA,
                  label=TRUE, lab.size=2, lab.family=FONT_FAMILY) +
        labs(x="Icasso Variant", y="N. Components")
    
    # create common palette
    col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    
    # check components found in either implementation
    X = comparison_pearson
    plts[["comparison_methods-pearson-heatmap"]] = X %>% 
        Heatmap(
            name="Pearson\nCorrelation", cluster_columns=TRUE, cluster_rows=TRUE,
            col=col_fun, rect_gp = gpar(col = "white", lwd = 0.5),
            row_names_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY),
            column_names_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY),
            heatmap_legend_param = list(legend_gp = gpar(fontsize = 6, fontfamily=FONT_FAMILY))
        ) %>% 
        draw() %>%
        grid.grabExpr() %>%
        as.ggplot()
    ## robustica found 49; Sastry2019 found 45+1 (considering noisy component -1)
    ## robustica creates components not found by Sastry2019
    
    # make gene sets using 3 stds as threhsold
    # compare gene sets through jaccard distances
    X = comparison_jaccard
    plts[["comparison_methods-jaccard-heatmap"]] = X %>% 
        Heatmap(
            name="Jaccard\nDistance", cluster_columns=TRUE, cluster_rows=TRUE,
            col=col_fun, rect_gp = gpar(col = "white", lwd = 0.5),
            row_names_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY),
            column_names_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY),
            heatmap_legend_param = list(legend_gp = gpar(fontsize = 6, fontfamily=FONT_FAMILY))
        ) %>% 
        draw() %>%
        grid.grabExpr() %>%
        as.ggplot()
    
    return(plts)
}


make_plots = function(performance, clustering, pca_components,comparison_pearson, comparison_jaccard){
    plts = list(
        plot_performance_profile(performance),
        plot_silhouettes(clustering),
        plot_pca(pca_components, clustering),
        plot_comparison_methods(comparison_pearson, comparison_jaccard)
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
    
    save_plt(plts, 'silhouettes_vs_time','.pdf',figs_dir, width=7, height=7)
    save_plt(plts, 'silhouettes_vs_max_memory','.pdf',figs_dir, width=7, height=7)
    save_plt(plts, 'silhouettes_vs_max_memory_vs_time','.pdf',figs_dir, width=7, height=7)
    
    save_plt(plts, 'clustering_pca-overview','.pdf',figs_dir, width=4, height=4)
    save_plt(plts, 'clustering_pca-AgglomerativeClustering','.pdf',figs_dir, width=4, height=4)
    save_plt(plts, 'clustering_pca-AffinityPropagation','.pdf',figs_dir, width=4, height=4)
    save_plt(plts, 'clustering_pca-DBSCAN','.pdf',figs_dir, width=4, height=4)
    save_plt(plts, 'clustering_pca-KMedoids','.pdf',figs_dir, width=4, height=4)
    save_plt(plts, 'clustering_pca-CommonNNClustering','.pdf',figs_dir, width=4, height=4)
    save_plt(plts, 'clustering_pca-OPTICS','.pdf',figs_dir, width=4, height=4)
    
    save_plt(plts, 'comparison_methods-n_components','.pdf',figs_dir, width=4, height=4)
    save_plt(plts, 'comparison_methods-pearson-heatmap','.pdf',figs_dir, width=13, height=11)
    save_plt(plts, 'comparison_methods-jaccard-heatmap','.pdf',figs_dir, width=13, height=11)
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
    components_robustica = read_tsv(components_robustica_file)
    components_Sastry2019 = read_csv(components_Sastry2019_file)
    
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
    
    ## compare components Sastry vs. robustica
    components_robustica = components_robustica %>% 
        dplyr::select(-(`-1`)) %>%
        column_to_rownames("log-TPM")
    colnames(components_robustica) = paste0("rcomp",colnames(components_robustica))
    
    components_Sastry2019 = components_Sastry2019 %>% 
        column_to_rownames("log-TPM")
    colnames(components_Sastry2019) = paste0("scomp",colnames(components_Sastry2019))
    
    ### pearson correlation between components
    comparison_pearson = cor(components_robustica, components_Sastry2019, method="pearson")
    
    ### jaccard distances between gene sets
    gene_sets_robustica = apply(components_robustica, 2, define_module)
    gene_sets_Sastry2019 = apply(components_Sastry2019, 2, define_module)
    comparison_jaccard = simil(gene_sets_robustica, gene_sets_Sastry2019, 
                               method='Jaccard', by_rows=FALSE)
    
    
    # visualize
    plts = make_plots(performance, clustering, pca_components,comparison_pearson, comparison_jaccard)
    figdata = make_figdata(performance, clustering, pca_components)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
