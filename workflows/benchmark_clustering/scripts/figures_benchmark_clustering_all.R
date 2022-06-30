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
PAL_ALGOS = setNames(
    get_palette("Paired",6),
    c("AffinityPropagation", "AgglomerativeClustering",
      "CommonNNClustering", "DBSCAN", "KMedoids", "OPTICS"))

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
# dataset_info_file = file.path(ROOT,'results','preprocess_data','files','dataset_info.tsv')
# figs_dir = file.path(ROOT,'results','benchmark_clustering','figures','all')

##### FUNCTIONS #####
plot_clustering_eda = function(clustering){
    X = clustering 
    
    plts = list()
    
    # how many components per cluster does each algorithm find?
    plts[["clustering_eda-components_per_cluster"]] = X %>%
        count(label, property_oi, cluster_id) %>%
        ggboxplot(x="label", y="n", color="property_oi", 
                  palette=PAL_ALGOS, outlier.size=0.1) +
        facet_wrap(~property_oi, scales="free_x") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="Dataset", y="N. Independent Components per Cluster") +
        coord_flip()
    
    # how many clusters per dataset and algorithm were found?
    plts[["clustering_eda-cluster_count"]] = X %>%
        distinct(label, property_oi, cluster_id) %>%
        count(label, property_oi) %>%
        ggstripchart(x="label", y="n", color="property_oi", palette=PAL_ALGOS, size=1) +
        labs(x="Dataset", y="N. Robust Independent Components", color="Clustering Algorithm") +
        theme_pubr(x.text.angle = 70)
    
    # time and memory to cluster
    plts[["clustering_eda-time"]] = X %>%
        distinct(property_oi, time, label) %>%
        ggplot(aes(x=property_oi, y=time)) +
        geom_boxplot(aes(color=property_oi), outlier.size=0.1) +
        labs(x="Clustering Algorithm", y="log10(Dataset Time (s))") +
        yscale("log10", .format=TRUE) +
        guides(color="none") +
        color_palette(PAL_ALGOS) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        coord_flip()
    
    plts[["clustering_eda-max_memory"]] = X %>%
        distinct(property_oi, max_memory, dataset) %>%
        ggplot(aes(x=property_oi, y=max_memory)) +
        geom_boxplot(aes(color=property_oi), outlier.size=0.1) +
        labs(x="Clustering Algorithm", y="Dataset Max. Memory (MiB)") +
        guides(color="none") +
        color_palette(PAL_ALGOS) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        coord_flip()
    
    # time and memory to cluster with different number of samples and genes
    plts[["clustering_eda-n_samples_vs_time"]] = X %>%
        distinct(n_samples, time, property_oi) %>%
        ggscatter(x="n_samples", y="time", color="property_oi", palette=PAL_ALGOS, size=1) +
        labs(x="N. Samples", y="log10(Dataset Time (s))", color="Clustering Algorithm") +
        yscale("log10", .format=TRUE)
    
    plts[["clustering_eda-n_samples_vs_max_memory"]] = X %>%
        distinct(n_samples, max_memory, property_oi) %>%
        ggscatter(x="n_samples", y="max_memory", color="property_oi", palette=PAL_ALGOS, size=1) +
        labs(x="N. Samples", y="Dataset Max. Memory (MiB)", color="Clustering Algorithm") 
    
    plts[["clustering_eda-n_genes_vs_time"]] = X %>%
        distinct(n_genes, time, property_oi) %>%
        ggboxplot(x="n_genes", y="time", color="property_oi", palette=PAL_ALGOS, outlier.size=0.1) +
        labs(x="N. Genes", y="log10(Dataset Time (s))", color="Clustering Algorithm") +
        yscale("log10", .format=TRUE)
    
    plts[["clustering_eda-n_genes_vs_max_memory"]] = X %>%
        distinct(n_genes, max_memory, property_oi) %>%
        ggboxplot(x="n_genes", y="max_memory", color="property_oi", palette=PAL_ALGOS, outlier.size=0.1) +
        labs(x="N. Genes", y="Dataset Max. Memory (MiB)", color="Clustering Algorithm")
    
    # how many components were noisy? (noise = silhouette score < 0.5)
    plts[["clustering_eda-noisy_count"]] = X %>%
        distinct(label, property_oi, is_noisy, cluster_id, component) %>%
        count(label, property_oi, is_noisy, .drop=FALSE) %>%
        ggbarplot(x="label", y="n", fill="is_noisy", 
                  color=NA, palette=c("orange","black")) + 
        facet_wrap(~property_oi) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Dataset", y="N. Independent Components", fill="Is Noisy") +
        coord_flip()
    
    plts[["clustering_eda-low_silhouette_count"]] = X %>%
        distinct(label, property_oi, is_low_silhouette, cluster_id, component) %>%
        count(label, property_oi, is_low_silhouette, .drop=FALSE) %>%
        ggbarplot(x="label", y="n", fill="is_low_silhouette", 
                  color=NA, palette=c("orange","black")) + 
        facet_wrap(~property_oi) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Dataset", y="N. Independent Components", fill="Is Low Silhouette") +
        coord_flip()
    
    return(plts)
}


plot_silhouettes = function(clustering, lab='', labsize=0.1){
    # get median silhouette score for each dataset and clustering algorithm with
    # their corresponing runtime and max_memory usage
    
    X = clustering %>% 
        # mean silhouette per cluster
        group_by(dataset, label, property_oi, time, max_memory, cluster_id, n_genes, n_samples) %>%
        summarize(mean_silhouette = mean(silhouette_euclidean)) %>%
        ungroup() %>%
        # median silhouette accross clusters
        group_by(dataset, label, property_oi, time, max_memory, n_genes, n_samples) %>%
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
        ungroup() %>%
        mutate(
            med_silhouette = replace_na(med_silhouette, 0) # placeholder
        )
        
        
    plts = list()    
    
    # silhouettes vs time
    plts[['silhouettes_vs_time']] = X %>%
        ggplot(aes(x=med_time, y=median_silhouette, color=property_oi)) + 
        geom_boxplot(outlier.size = 0.1, width=1) +
        color_palette(PAL_ALGOS) +
        theme_pubr() +
        theme(aspect.ratio = 1) +
        labs(x='log10(Median Time (s))', y='Dataset Silhouette Score', color="Clustering Algorithm") +
        xscale("log10", .format=TRUE)

    # silhouettes vs memory
    plts[['silhouettes_vs_max_memory']] = X %>%
        ggplot(aes(x=med_max_memory, y=median_silhouette, color=property_oi)) + 
        geom_boxplot(outlier.size = 0.1, width=0.025) + 
        theme_pubr() + 
        theme(aspect.ratio = 1) +
        labs(x='log10(Median Max. Memory (MiB))', y='Dataset Silhouette Score', color="Clustering Algorithm") +
        color_palette(PAL_ALGOS)  +
        xscale("log10", .format=TRUE)
    
    # memory vs time vs silhouettes
    plts[['silhouettes_vs_max_memory_vs_time']] = X %>% 
        distinct(med_time, med_max_memory, med_silhouette, property_oi) %>%
        ggscatter(x='med_time', y='med_max_memory', size='med_silhouette', 
                  color='property_oi', palette=PAL_ALGOS) +
        geom_text_repel(aes(label=property_oi, color=property_oi),  
                        size=2, family=FONT_FAMILY, segment.size=0.1) +
        labs(x='log10(Median Time (s))', y='log10(Median Max. Memory (MiB))', 
             size='Median Silhouette Score') +
        theme(aspect.ratio = 1) +
        guides(color="none") +
        xscale("log10", .format=TRUE) +
        yscale("log10", .format=TRUE)
    
    plts[['silhouettes_vs_max_memory_vs_time_detailed']] = X %>% 
        ggscatter(x="time", y="max_memory", color="property_oi", 
                  palette=PAL_ALGOS, size=1, alpha=0.5) +
        geom_point(
            aes(x=med_time, y=med_max_memory, size=med_silhouette, color=property_oi),
            X %>% distinct(med_time, med_max_memory, med_silhouette, property_oi),
            shape=21
        ) +
        geom_text_repel(
            aes(x=med_time, y=med_max_memory, label=property_oi, color=property_oi),  
            X %>% distinct(med_time, med_max_memory, med_silhouette, property_oi),
            size=2, family=FONT_FAMILY, segment.size=0.1
        ) +
        labs(x='log10(Median Time (s))', y='log10(Median Max. Memory (MiB))', 
             size='Median Silhouette Score') +
        theme(aspect.ratio = 1) +
        guides(color="none") +
        xscale("log10", .format=TRUE) +
        yscale("log10", .format=TRUE)
    
    # silhouettes to cluster with different number of samples and genes
    plts[["clustering_eda-n_samples_vs_median_silhouette"]] = X %>%
        distinct(n_samples, median_silhouette, property_oi) %>%
        ggscatter(x="n_samples", y="median_silhouette", color="property_oi", palette=PAL_ALGOS, size=1) +
        labs(x="N. Samples", y="Dataset Silhouette Score", color="Clustering Algorithm")
    
    plts[["clustering_eda-n_genes_vs_median_silhouette"]] = X %>%
        distinct(n_genes, median_silhouette, property_oi) %>%
        ggboxplot(x="n_genes", y="median_silhouette", color="property_oi", palette=PAL_ALGOS, outlier.size=0.1) +
        labs(x="N. Genes", y="Dataset Silhouette Score", color="Clustering Algorithm")
    
    names(plts) = paste0(lab, names(plts))
    return(plts)
}


make_plots = function(performance, clustering){
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
    save_plt(plts, 'clustering_eda-components_per_cluster','.pdf',figs_dir, width=8, height=19)
    save_plt(plts, 'clustering_eda-noisy_count','.pdf',figs_dir, width=8, height=20.5)
    save_plt(plts, 'clustering_eda-low_silhouette_count','.pdf',figs_dir, width=8, height=20.5)
    save_plt(plts, 'clustering_eda-cluster_count','.pdf',figs_dir, width=11, height=8)
    save_plt(plts, 'clustering_eda-time','.pdf',figs_dir, width=6, height=6)
    save_plt(plts, 'clustering_eda-max_memory','.pdf',figs_dir, width=6, height=6)
    save_plt(plts, 'clustering_eda-n_samples_vs_time','.pdf',figs_dir, width=6, height=6)
    save_plt(plts, 'clustering_eda-n_samples_vs_max_memory','.pdf',figs_dir, width=6, height=6)
    save_plt(plts, 'clustering_eda-n_samples_vs_median_silhouette','.pdf',figs_dir, width=6, height=6)
    save_plt(plts, 'clustering_eda-n_genes_vs_time','.pdf',figs_dir, width=6, height=6)
    save_plt(plts, 'clustering_eda-n_genes_vs_max_memory','.pdf',figs_dir, width=6, height=6)
    save_plt(plts, 'clustering_eda-n_genes_vs_median_silhouette','.pdf',figs_dir, width=6, height=6)

    save_plt(plts, 'silhouettes_vs_time','.pdf',figs_dir, width=6, height=6)
    save_plt(plts, 'silhouettes_vs_max_memory','.pdf',figs_dir, width=6, height=6)
    save_plt(plts, 'silhouettes_vs_max_memory_vs_time','.pdf',figs_dir, width=6, height=6)
    save_plt(plts, 'silhouettes_vs_max_memory_vs_time_detailed','.pdf',figs_dir, width=6, height=6)
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
    dataset_info = read_tsv(dataset_info_file) %>%
        mutate(label = paste0(main, " | ", dataset),
               label = gsub("Sastry2019 \\|","",label))
    
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
        left_join(max_memory, by=c('dataset','property_oi')) %>%
        left_join(dataset_info, by="dataset") %>%
        # label noisy components/clusters
        group_by(property_oi, cluster_id) %>%
        mutate(is_low_silhouette = mean(silhouette_euclidean) < 0.5 | cluster_id==-1,
               is_noisy = cluster_id==-1) %>%
        ungroup()
    
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
