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
require(proxy)
require(fdrtool)
require(extrafont)
require(ComplexHeatmap)
require(circlize)
require(ggplotify)
require(optparse)

# variables
FONT_FAMILY = "Arial"
METRICS_OI = c("icasso","robustica_nosign","robustica","robustica_pca","icasso_pca")
PAL_METRICS = setNames(
    get_palette("npg",5),
    c('icasso','robustica_nosign','robustica','robustica_pca','icasso_pca')
)

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
# ROOT = here::here()
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','benchmark_sign_inference')
# performance_file = file.path(RESULTS_DIR,'files','merged_clustering_performance_evaluation.tsv.gz')
# clustering_file = file.path(RESULTS_DIR,'files','merged_clustering_info.tsv.gz')
# mapping_eval_file = file.path(RESULTS_DIR,'files','merged_module_mapping_evaluation.tsv.gz')
# mapping_robust_file = file.path(RESULTS_DIR,'files','merged_module_mapping_robustness.tsv.gz')
# dataset_info_file = file.path(ROOT,'results','preprocess_data','files','dataset_info.tsv')
# S_stds_file = file.path(RESULTS_DIR,'files','merged_summary_S_std.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,'figures','all')


##### FUNCTIONS #####
plot_clustering_eval = function(clustering_eval){
    X = clustering_eval %>%
        # compute average silhouette per cluster
        group_by(dataset, algorithm, cluster_id, time, max_memory, cluster_avg_std) %>%
        summarize(
            cluster_silhouette = mean(silhouette_euclidean),
            cluster_size = n()
        ) %>%
        ungroup() %>%
        # compute median of each feature across the clusters
        # of each dataset
        group_by(dataset, algorithm, time, max_memory)%>%
        summarize(
            dataset_cluster_silhouette = median(cluster_silhouette),
            dataset_cluster_size = median(cluster_size),
            dataset_cluster_std = median(cluster_avg_std),
            dataset_n_clusters = n()
        ) %>%
        ungroup() %>%
        distinct(dataset, algorithm, time, max_memory, dataset_cluster_silhouette, 
                 dataset_cluster_size, dataset_cluster_std, dataset_n_clusters)
    
    plts = list()
    # for each dataset and metric calculate the median of each feature
    # metric vs. time
    plts[["clustering_eval-metric_vs_time"]] = X %>%
        ggboxplot(x="algorithm", y="time", color="algorithm", 
                  palette=PAL_METRICS, fill=NA, outlier.size=0.1) +
        yscale("log10", .format=TRUE) +
        guides(color="none") +
        theme_pubr(x.text.angle=70) +
        labs(x="Metric", y="log10(Dataset Time (s))")
    
    # metric vs. max. memory
    plts[["clustering_eval-metric_vs_max_memory"]] = X %>% 
        ggboxplot(x="algorithm", y="max_memory", color="algorithm", 
                  palette=PAL_METRICS, fill=NA, outlier.size=0.1) +
        guides(color="none") +
        theme_pubr(x.text.angle=70) +
        labs(x="Metric", y="Dataset Max. Memory (MiB)")
    
    # metric vs. silhouette
    plts[["clustering_eval-metric_vs_silhouette"]] = X %>% 
        ggboxplot(x="algorithm", y="dataset_cluster_silhouette", color="algorithm",
                  palette=PAL_METRICS, fill=NA, outlier.size=0.1) +
        guides(color="none") +
        theme_pubr(x.text.angle=70) +
        labs(x="Metric", y="Dataset Silhouette Score")
    
    # metric vs. n. clusters
    plts[["clustering_eval-metric_vs_n_clusters"]] = X %>% 
        ggboxplot(x="algorithm", y="dataset_n_clusters", color="algorithm",
                  palette=PAL_METRICS, fill=NA, outlier.size=0.1) +
        guides(color="none") +
        theme_pubr(x.text.angle=70) +
        labs(x="Metric", y="N. Robust Independent Components")
    
    # metric vs. n. clusters vs. silhouette?
    plts[["clustering_eval-metric_vs_n_clusters_vs_silhouette"]] = X %>% 
        ggscatter(x="dataset_n_clusters", y="dataset_cluster_silhouette", color="algorithm",
                  palette=PAL_METRICS, fill=NA, size=1) +
        theme(aspect.ratio=1) +
        labs(x="N. Robust Independent Components", y="Dataset Silhouette Score", color="Metric")
    
    # metric vs. size clusters
    plts[["clustering_eval-metric_vs_cluster_size"]] = X %>% 
        ggboxplot(x="algorithm", y="dataset_cluster_size", color="algorithm",
                  palette=PAL_METRICS, fill=NA, outlier.size=0.1) +
        guides(color="none") +
        theme_pubr(x.text.angle=70) +
        labs(x="Metric", y="Dataset N. Components in Robust Comp.")

    # metric vs. cluster average std (add data!)
    plts[["clustering_eval-metric_vs_cluster_std"]] = X %>% 
        ggboxplot(x="algorithm", y="dataset_cluster_std", color="algorithm",
                  palette=PAL_METRICS, fill=NA, outlier.size=0.1) +
        guides(color="none") +
        theme_pubr(x.text.angle=70) +
        labs(x="Metric", y="Dataset Avg. Component Std.")
    
    return(plts)
}


plot_mapping_eval = function(mapping_eval){
    X = mapping_eval %>%
        group_by(dataset, algorithm, runs) %>%
        summarize(dataset_jaccard = mean(jaccard),
                  dataset_n_clusters = n()) %>%
        ungroup() 
    
    metrics_oi = c('icasso','robustica_pca')
    
    plts = list()
    plts[['mapping_eval-runs_vs_jaccard']] = X %>% 
        filter(algorithm%in%metrics_oi) %>%
        group_by(runs, algorithm) %>%
        summarize(
            median = median(dataset_jaccard),
            q25 = quantile(dataset_jaccard, 0.25),
            q75 = quantile(dataset_jaccard, 0.75)
        ) %>%
        ungroup() %>%
        ggplot(aes(x=runs, y=median, color=algorithm)) +
        geom_line(size=0.1, linetype="dashed") +
        geom_point(size=0.1) +
        geom_errorbar(aes(ymin=q25, ymax=q75), width=1, size=0.1) +
        color_palette(PAL_METRICS[metrics_oi]) +
        theme_pubr() +
        labs(x='Runs', y='Dataset Jaccard', color="Metric")
    
    plts[['mapping_eval-runs_vs_n_components']] = X %>% 
        filter(algorithm%in%metrics_oi) %>%
        group_by(runs, algorithm) %>%
        summarize(
            median = median(dataset_n_clusters),
            q25 = quantile(dataset_n_clusters, 0.25),
            q75 = quantile(dataset_n_clusters, 0.75)
        ) %>%
        ungroup() %>%
        ggplot(aes(x=runs, y=median, color=algorithm)) +
        geom_line(size=0.1, linetype="dashed") +
        geom_point(size=0.1) +
        geom_errorbar(aes(ymin=q25, ymax=q75), width=1, size=0.1) +
        color_palette(PAL_METRICS[metrics_oi]) +
        theme_pubr() + 
        labs(x='Runs', y='N. Independent Components', color="Metric")
    
    return(plts)
}


plot_mapping_robust = function(mapping_robust){
    X = mapping_robust %>%
        group_by(dataset, algorithm, n_samples) %>%
        summarize(dataset_jaccard = median(mean)) %>%
        ungroup()
    
    plts = list()
    
    plts[["mapping_robust-metric_vs_jaccard"]] = X %>%
        ggboxplot(x="algorithm", y="dataset_jaccard", color="algorithm", 
                  palette=PAL_METRICS, fill=NA, outlier.size=0.1) +
        stat_compare_means(method="wilcox.test", size=2, family=FONT_FAMILY) +
        theme_pubr(x.text.angle=70) +
        guides(color="none") +
        labs(x="Metric", y="Dataset Jaccard")
    
    return(plts)
}


make_plots = function(clustering_eval, mapping_eval, mapping_robust){
    plts = list(
        plot_clustering_eval(clustering_eval),
        plot_mapping_eval(mapping_eval),
        plot_mapping_robust(mapping_robust)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(clustering_eval, mapping_eval, mapping_robust){
    summary_clustering = clustering_eval %>% 
        group_by(dataset, algorithm, time, max_memory, cluster_id) %>% 
        summarize(mean_silhouette = mean(silhouette_euclidean)) %>%
        ungroup() %>%
        group_by(dataset, algorithm, time, max_memory) %>%
        summarize(
            mean=mean(mean_silhouette), 
            median=median(mean_silhouette), 
            std=sd(mean_silhouette), 
            range=max(mean_silhouette) - min(mean_silhouette)
        ) %>%
        ungroup()
    
    figdata = list(
        'benchmark_clustering-evaluation' = list(
            'clustering_evaluation' = clustering_eval,
            'clustering_evaluation_summary' = summary_clustering,
            'module_mapping_diff_runs' = mapping_eval,
            'module_mapping_robustness' = mapping_robust
        )
    )
    return(figdata)
}


save_figdata = function(figdata, dir){
    lapply(names(figdata), function(x){
        d = file.path(dir,'figdata',x)
        dir.create(d, recursive=TRUE)
        lapply(names(figdata[[x]]), function(nm){
            df = figdata[[x]][[nm]]
            filename = file.path(d, paste0(nm,'.tsv.gz'))
            write_tsv(df, filename)
            
            print(filename)
        })
    })
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
    save_plt(plts,'clustering_eval-metric_vs_time','.pdf',figs_dir, width=5, height=5)
    save_plt(plts,'clustering_eval-metric_vs_max_memory','.pdf',figs_dir, width=5.2, height=5)
    save_plt(plts,'clustering_eval-metric_vs_silhouette','.pdf',figs_dir, width=5, height=5)
    save_plt(plts,'clustering_eval-metric_vs_n_clusters','.pdf',figs_dir, width=5, height=5)
    save_plt(plts,'clustering_eval-metric_vs_n_clusters_vs_silhouette','.pdf',figs_dir, width=5, height=5)
    save_plt(plts,'clustering_eval-metric_vs_cluster_size','.pdf',figs_dir, width=5, height=5)
    save_plt(plts,'clustering_eval-metric_vs_cluster_std','.pdf',figs_dir, width=5.2, height=5)
    save_plt(plts,'mapping_eval-runs_vs_jaccard','.pdf',figs_dir, width=5, height=5)
    save_plt(plts,'mapping_eval-runs_vs_n_components','.pdf',figs_dir, width=5, height=5)
    save_plt(plts,'mapping_robust-metric_vs_jaccard','.pdf',figs_dir, width=5, height=5)
}


parseargs = function(){
    
    option_list = list( 
        make_option("--performance_file", type="character"),
        make_option("--clustering_file", type="character"),
        make_option("--S_stds_file", type="character"),
        make_option("--dataset_info_file", type="character"),
        make_option("--mapping_eval_file", type="character"),
        make_option("--mapping_robust_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    performance_file = args[["performance_file"]]
    clustering_file = args[["clustering_file"]]
    S_stds_file = args[["S_stds_file"]]
    dataset_info_file = args[["dataset_info_file"]]
    mapping_eval_file = args[["mapping_eval_file"]]
    mapping_robust_file = args[["mapping_robust_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load data
    performance = read_tsv(performance_file)
    clustering = read_tsv(clustering_file)
    S_stds = read_tsv(S_stds_file)
    dataset_info = read_tsv(dataset_info_file) %>%
        mutate(label = paste0(main, " | ", dataset),
               label = gsub("Sastry2019 \\|","",label))
        
    mapping_eval = read_tsv(mapping_eval_file)
    mapping_robust = read_tsv(mapping_robust_file)    
    
    # preprocess
    ## performance
    performance = performance %>% 
        filter(`function` != 'evaluate_performance') %>%
        group_by(dataset, algorithm) %>% 
        mutate(rel_timestamp = timestamp - min(timestamp),
               `function` = paste0(`function`,'()')) %>%
        ungroup()
    
    ## clustering time and evaluation    
    clustering_time = performance %>% 
        filter(`function`=='cluster_components()') %>% 
        group_by(dataset, algorithm) %>%
        summarize(time=max(time)) %>%
        distinct() %>%
        ungroup()
    
    max_memory = performance %>% 
        group_by(dataset, algorithm) %>% 
        summarize(max_memory=max(memory)) %>%
        ungroup()
    
    clustering_eval = clustering %>%
        left_join(clustering_time, by=c('dataset','algorithm')) %>%
        left_join(max_memory, by=c('dataset','algorithm')) %>%
        left_join(S_stds, by=c("cluster_id","dataset","algorithm")) %>%
        left_join(dataset_info, by="dataset") %>%
        # label noisy components/clusters
        group_by(dataset, algorithm, cluster_id) %>%
        mutate(is_low_silhouette = mean(silhouette_euclidean) < 0.5 | cluster_id==-1,
               is_noisy = cluster_id==-1) %>%
        ungroup() %>%
        # only evaluate Icasso, Sign Infer., PCA
        filter(algorithm %in% METRICS_OI) %>%
        mutate(algorithm = factor(algorithm, levels=METRICS_OI))
    
    # visualize
    plts = make_plots(clustering_eval, mapping_eval, mapping_robust)
    figdata = make_figdata(clustering_eval, mapping_eval, mapping_robust)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
