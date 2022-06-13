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
require(ggrepel)
require(writexl)
require(extrafont)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
FONT_FAMILY = "Arial"
PAL_ALGOS = setNames(get_palette("Paired",6),
                     c("AffinityPropagation", "AgglomerativeClustering","CommonNNClustering",
                       "DBSCAN","KMedoids", "OPTICS"))

DATASETS_REF = rbind(
    data.frame(
        main = "",
        dataset = "Sastry2019"
    ),
    data.frame(
        main = "TCGA",
        dataset = c(
        'BRCA',
        'KIRC',
        'LUAD',
        'LUSC',
        'OV',
        'HNSC',
        'GBM',
        'UCEC',
        'THCA',
        'PRAD',
        'COAD',
        'LGG',
        'STAD',
        'SKCM',
        'LIHC',
        'BLCA',
        'KIRP',
        'CESC',
        'SARC',
        'ESCA',
        'LAML'
        )
    ),
    data.frame(
        main = "GTEx",
        dataset = c(
            'blood',
            'brain', 
            'skin',
            'esophagus',
            'blood_vessel',
            'adipose_tissue', 
            'heart',
            'muscle',
            'lung',
            'colon', 
            'thyroid',
            'nerve',
            'breast',
            'testis',
            'stomach', 
            'pancreas', 
            'pituitary',
            'adrenal_gland', 
            'prostate', 
            'spleen',
            'liver'
        )
    )
)

# Development
# -----------
# RESULTS_DIR = file.path(ROOT,'results','benchmark_clustering','files')
# performance_file = file.path(RESULTS_DIR,'clustering_performance_evaluation_merged.tsv.gz')
# clustering_file = file.path(RESULTS_DIR,'clustering_info_merged.tsv.gz')
# figs_dir = file.path(ROOT,'results','benchmark_clustering','figures','all')


##### FUNCTIONS #####
plot_clustering_eda = function(clustering){
    X = clustering
    
    plts = list()
    
    # how many components per cluster does each algorithm find?
    plts[["clustering_eda-components_per_cluster"]] = X %>%
        count(dataset, property_oi, cluster_id) %>%
        ggboxplot(x="dataset", y="n", fill="property_oi", 
                  palette=PAL_ALGOS, outlier.size=0.1) +
        facet_wrap(~property_oi, scales="free_x") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(fill="none") +
        labs(x="Dataset", y="N. Independent Components per Cluster") +
        coord_flip()
    
    # how many clusters per dataset and algorithm were found?
    plts[["clustering_eda-cluster_count"]] = X %>%
        distinct(dataset, property_oi, cluster_id) %>%
        count(dataset, property_oi) %>%
        ggbarplot(x="dataset", y="n", fill="property_oi", 
                  color=NA, palette=PAL_ALGOS) +
        facet_wrap(~property_oi) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(fill="none") +
        labs(x="Dataset", y="N. Robust Independent Components") +
        coord_flip()
    
    # how many components were assigned as noisy?
    plts[["clustering_eda-noisy_count"]] = X %>%
        distinct(dataset, property_oi, component, cluster_id) %>%
        mutate(is_noisy = as.factor(cluster_id == -1)) %>%
        count(dataset, property_oi, is_noisy, .drop=FALSE) %>%
        ggbarplot(x="dataset", y="n", fill="is_noisy", 
                  color=NA, palette=c("orange","black")) + 
        facet_wrap(~property_oi) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Dataset", y="N. Independent Components", fill="Is Noisy") +
        coord_flip()
        
    # time and memory to cluster
    plts[["clustering_eda-time"]] = X %>%
        distinct(property_oi, time, dataset) %>%
        ggplot(aes(x=property_oi, y=time)) +
        geom_boxplot(aes(fill=property_oi), outlier.size=0.1) +
        labs(x="Clustering Algorithm", y="Time (s)") +
        guides(fill="none") +
        fill_palette(PAL_ALGOS) +
        theme_pubr() +
        coord_flip()
    
    plts[["clustering_eda-max_memory"]] = X %>%
        distinct(property_oi, max_memory, dataset) %>%
        ggplot(aes(x=property_oi, y=max_memory)) +
        geom_boxplot(aes(fill=property_oi), outlier.size=0.1) +
        labs(x="Clustering Algorithm", y="Max. Memory (MiB)") +
        guides(fill="none") +
        fill_palette(PAL_ALGOS) +
        theme_pubr() +
        coord_flip()
    
    return(plts)
}

plot_silhouettes = function(clustering, lab='', labsize=0.1){
    # get median silhouette score for each dataset and clustering algorithm with
    # their corresponing runtime and max_memory usage
    
    X = clustering %>% 
        # mean silhouette per cluster
        group_by(dataset, property_oi, time, max_memory, cluster_id) %>%
        summarize(mean_silhouette = mean(silhouette_euclidean)) %>%
        ungroup() %>%
        # median silhouette accross clusters
        group_by(dataset, property_oi, time, max_memory) %>%
        summarize(median_silhouette = median(mean_silhouette)) %>%
        # summarize across datasets
        group_by(property_oi) %>%
        mutate(
            time=as.numeric(time),
            med_silhouette = median(median_silhouette, na.rm=TRUE),
            silh25 = quantile(median_silhouette, 0.25, na.rm=TRUE),
            silh75 = quantile(median_silhouette, 0.75, na.rm=TRUE),
            
            time25 = quantile(time, 0.25),
            time75 = quantile(time, 0.75),
            med_time = median(time),

            mem25 = quantile(max_memory, 0.25),
            mem75 = quantile(max_memory, 0.75),
            med_max_memory = median(max_memory)
        ) %>%
        ungroup()
        
    plts = list()    
    
    # silhouettes vs time
    plts[['silhouettes_vs_time']] = X %>%
        ggplot(aes(x=med_time, y=median_silhouette, fill=property_oi)) + 
        geom_boxplot(outlier.size = 0.1, width=100) +
        fill_palette(PAL_ALGOS) +
        theme(aspect.ratio = 1) +
        theme_pubr() +
        labs(x='Median Time (s)', y='Silhouette Score')

    # silhouettes vs memory
    plts[['silhouettes_vs_max_memory']] = X %>%
        ggplot(aes(x=med_max_memory, y=median_silhouette, fill=property_oi)) + 
        geom_boxplot(outlier.size = 0.1, width=100) + 
        theme(aspect.ratio = 1) +
        theme_pubr() + 
        labs(x='Median Max. Memory (MiB)', y='Silhouette Score') +
        fill_palette(PAL_ALGOS)
    
    # memory vs time vs silhouettes
    plts[['silhouettes_vs_max_memory_vs_time']] = X %>% 
        distinct(med_time, med_max_memory, med_silhouette, property_oi) %>%
        ggscatter(x='med_time', y='med_max_memory', size='med_silhouette', 
                  repel = TRUE, color='property_oi', label='property_oi', 
                  palette=PAL_ALGOS) +
        labs(x='Median Time (s)', y='Median Max. Memory (MiB)', size='Median Silhouette Score') +
        guides(color="none")

    names(plts) = paste0(lab, names(plts))
    return(plts)
}


make_plots = function(performance, clustering, pca_components){
    plts = list(
        plot_clustering_eda(clustering),
        plot_silhouettes(clustering)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(performance, clustering){
    summary_clustering = clustering %>% 
        group_by(dataset, property_oi, time, max_memory, cluster_id) %>% 
        summarize(mean_silhouette = mean(silhouette_euclidean)) %>%
        ungroup() %>%
        group_by(dataset, property_oi, time, max_memory) %>%
        summarize(
            mean=mean(mean_silhouette), 
            median=median(mean_silhouette), 
            std=sd(mean_silhouette), 
            range=max(mean_silhouette) - min(mean_silhouette)
        ) %>%
        ungroup()
    
    figdata = list(
        'benchmark_clustering-evaluation' = list(
            'performance_evaluation' = performance,
            'clustering_evaluation' = clustering,
            'clustering_evaluation_summary' = summary_clustering
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
    save_plt(plts, 'clustering_eda-components_per_cluster','.pdf',figs_dir, width=10, height=8)
    save_plt(plts, 'clustering_eda-cluster_count','.pdf',figs_dir, width=10, height=8)
    save_plt(plts, 'clustering_eda-noisy_count','.pdf',figs_dir, width=10, height=8)
    save_plt(plts, 'clustering_eda-time','.pdf',figs_dir, width=10, height=8)
    save_plt(plts, 'clustering_eda-max_memory','.pdf',figs_dir, width=10, height=8)
    
    save_plt(plts, 'silhouettes_vs_time','.pdf',figs_dir, width=6, height=6)
    save_plt(plts, 'silhouettes_vs_max_memory','.pdf',figs_dir, width=6, height=6)
    save_plt(plts, 'silhouettes_vs_max_memory_vs_time','.pdf',figs_dir, width=6, height=6)
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
    
    ## keep those algorithms that we could evaluate
    avail_algorithms = intersect(
        performance[["property_oi"]], clustering[["property_oi"]]
    )
    
    performance = performance %>% filter(property_oi %in% avail_algorithms)
    clustering = clustering %>% filter(property_oi %in% avail_algorithms)
    
    ## performances
    performance = performance %>% 
        filter(`function` != 'evaluate_performance') %>%
        group_by(dataset, property_oi) %>% 
        mutate(rel_timestamp = timestamp - min(timestamp),
               `function` = paste0(`function`,'()')) %>%
        ungroup()
    
    ## clustering time and evaluation    
    clustering_time = performance %>% 
        filter(`function`=='cluster_components()') %>% 
        group_by(dataset, property_oi) %>%
        summarize(time=max(time)) %>%
        distinct() %>%
        ungroup()
    
    max_memory = performance %>% 
        group_by(dataset, property_oi) %>% 
        summarize(max_memory=max(memory)) %>%
        ungroup()
    
    clustering = clustering %>%
        left_join(clustering_time, by=c('dataset','property_oi')) %>%
        left_join(max_memory, by=c('dataset','property_oi'))
    
    ## beautify names
    clustering = clustering %>%
        left_join(DATASETS_REF, by="dataset") %>%
        mutate(dataset = sprintf("%s | %s", main, dataset))
    
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
