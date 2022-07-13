#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Interpret the independent components.
#
# Outline
# -------
# - Which components differentiate IDH1 mutant from WT samples?
#     - barplot of mut vs wt sample counts
#     - wilcoxon test + FDR, select components of interest
#     - volcano plot of statistical test
#     - heatmap of A all components
#     - heatmap of A components of interest
# - Which biological processes are enriched in those components?
# - How certain are we of these components?

require(tidyverse)
require(cowplot)
require(ggpubr)
require(extrafont)
require(scattermore)
require(fdrtool)
require(proxy)
require(ComplexHeatmap)
require(ggplotify)
require(circlize)
require(ggrepel)
require(survival)
require(survminer)
require(optparse)

# variables
DATASET_OI = "LGG"
METRICS_OI = c("icasso","robustica_nosign","robustica","robustica_pca","icasso_pca")

# formatting
FONT_FAMILY = "Arial"
PAL_METRICS = setNames(get_palette("npg",5), METRICS_OI)
PAL_FDR_DARK = "#005AB5"
PAL_FDR_LIGHT = "#DC3220"

# Development
# -----------
# ROOT = here::here()
# PREP_DIR = file.path(ROOT,'data','prep')

# RESULTS_SIGN_DIR = file.path(ROOT,'results','benchmark_sign_inference')
# performance_file = file.path(RESULTS_SIGN_DIR,'files','merged_clustering_performance_evaluation.tsv.gz')
# clustering_file = file.path(RESULTS_SIGN_DIR,'files','merged_clustering_info.tsv.gz')
# S_stds_file = file.path(RESULTS_SIGN_DIR,'files','merged_summary_S_std.tsv.gz')

# S_icasso_file = file.path(RESULTS_SIGN_DIR,'files',DATASET_OI,'icasso-S.tsv.gz')
# S_std_icasso_file = file.path(RESULTS_SIGN_DIR,'files',DATASET_OI,'icasso-S_std.tsv.gz')
# S_robustica_pca_file = file.path(RESULTS_SIGN_DIR,'files',DATASET_OI,'robustica_pca-S.tsv.gz')
# S_std_robustica_pca_file = file.path(RESULTS_SIGN_DIR,'files',DATASET_OI,'robustica_pca-S_std.tsv.gz')
# A_robustica_pca_file = file.path(RESULTS_SIGN_DIR,'files',DATASET_OI,'robustica_pca-A.tsv.gz')

# genexpr_file = file.path(PREP_DIR,'genexpr','TCGA','LGG.tsv.gz')
# snv_file = file.path(PREP_DIR,'snv','TCGA','LGG.tsv.gz')
# metadata_file = file.path(PREP_DIR,'metadata','TCGA','LGG.tsv')
# sample_indices_file = file.path(PREP_DIR,'sample_indices','LGG.tsv')

# RESULTS_DIR = file.path(ROOT,'results','case_study')
# enrichment_icasso_file = file.path(RESULTS_DIR,"files","gsea","LGG","icasso.tsv.gz")
# enrichment_robustica_pca_file = file.path(RESULTS_DIR,"files","gsea","LGG","robustica_pca.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,'figures',DATASET_OI)


###### FUNCTIONS #####
define_module = function(x, cutoff=0.01){
    fdr = fdrtool(x, plot=FALSE, cutoff.method="fndr", verbose=FALSE)[['qval']]
    x = fdr < cutoff
    return(x)
}


plot_performance = function(performance){
    X = performance
    
    plts = list()
    plts[["performance-time_vs_memory"]] = X %>%
        ggplot(aes(x=rel_timestamp, y=memory)) +
        geom_scattermore(aes(color=`function`), pointsize=10, 
                         pixels=c(1000,500)) +
        geom_line(aes(color=`function`), linetype="dashed", size=0.25)+
        color_palette("uchicago") +
        facet_wrap(~algorithm, ncol=1, scales='free') +
        theme_pubr(border = TRUE) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color=guide_legend(ncol=1)) + 
        labs(x='Time (s)', y='Memory (MiB)', color='Function')
    
    return(plts)
}


plot_clustering_eval = function(clustering_eval){
    X = clustering_eval
    
    plts = list()
    
    # metric vs cluster avg std.
    plts[["clustering_eval-metric_vs_cluster_std"]] = X %>%
        distinct(algorithm, cluster_id, cluster_avg_std) %>% 
        ggstripchart(x="algorithm", y="cluster_avg_std", size=0.1, 
                     color="algorithm", palette=PAL_METRICS) +
        geom_boxplot(aes(color=algorithm), fill=NA, outlier.shape=NA) +
        guides(color="none") +
        theme_pubr(x.text.angle = 70) +
        labs(x="Metric", y="Avg. Component Std.")
    
    # metric vs avg. silhouette score per cluster
    plts[["clustering_eval-metric_vs_silhouette"]] = X %>%
        distinct(algorithm, cluster_id, mean_silhouette) %>% 
        ggstripchart(x="algorithm", y="mean_silhouette", size=0.1, 
                     color="algorithm", palette=PAL_METRICS) +
        geom_boxplot(aes(color=algorithm), fill=NA, outlier.shape=NA) +
        guides(color="none") +
        theme_pubr(x.text.angle = 70) +
        labs(x="Metric", y="Silhouette Score")
    
    # metric vs n.components per cluster (+ label of total components per metric)
    x = X %>%
        distinct(algorithm, cluster_id) %>%
        count(algorithm) %>%
        mutate(label = paste0("n=",n))
    plts[["clustering_eval-metrics_vs_cluster_size"]] = X %>%
        distinct(algorithm, cluster_id, cluster_size) %>%
        ggstripchart(x="algorithm", y="cluster_size", size=0.1, 
                     color="algorithm", palette=PAL_METRICS) +
        geom_boxplot(aes(color=algorithm), fill=NA, outlier.shape=NA) +
        geom_text(aes(label=label, y=1e4), x, size=2, family=FONT_FAMILY) +
        guides(color="none") +
        theme_pubr(x.text.angle = 70) +
        yscale("log10", .format=TRUE) +
        labs(x="Metric", y="log10(N. Components in Robust Comp.)") 
    
    return(plts)
}


plot_comparison_methods = function(comparison_pearson, comparison_jaccard, mapping, clustering_eval, enrichments){
    # icasso vs robustica PCA
    
    plts = list()
    
    # Mapping of components
    ## create common palette
    col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    
    ## Pearson
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
    
    ## Jaccard
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
    
    ### there extra components found by robustica_pca contain new information
    
    # where do original components go?
    ## identify components in extra robust components identified by robustica_pca
    X = mapping
    ## mapping distribution
    plts[["comparison_methods-mapping-hist"]] = X %>%
        gghistogram(x="correlation", fill="grey", color=NA) +
        geom_vline(xintercept=0.5, linetype="dashed", size=0.1, color="black") +
        labs(x="Pearson Correlation", y="Count")
    
    ## silhouette of mapped vs unmapped in robustica_pca
    plts[["comparison_methods-mapping_vs_silhouette"]] = X %>%
        ggstripchart(x="mapped", y="mean_silhouette", size=0.1, color="mapped", palette="jco") +
        geom_boxplot(aes(color=mapped), fill=NA, width=0.4, outlier.shape=NA) +
        stat_compare_means(method="wilcox.test", size=2, family=FONT_FAMILY) +
        guides(color="none") +
        labs(x="Mapped Robust Component", y="Silhouette Score")
        
    ## cluster mean std of mapped vs unmapped in robustica_pca
    plts[["comparison_methods-mapping_vs_cluster_std"]] = X %>%
        ggstripchart(x="mapped", y="cluster_avg_std", size=0.1, color="mapped", palette="jco") +
        geom_boxplot(aes(color=mapped), fill=NA, width=0.4, outlier.shape=NA) +
        stat_compare_means(method="wilcox.test", size=2, family=FONT_FAMILY) +
        geom_text(aes(label=label, y=0.005), 
                  X %>%
                      count(mapped) %>%
                      mutate(label=paste0("n=",n)), 
                  size=2, family=FONT_FAMILY) +
        guides(color="none") +
        labs(x="Mapped Robust Component", y="Avg. Component Std.")

    ## find where those components are in the robust components in icasso
    unmapped_clusters = mapping %>%
        filter(!mapped) %>%
        pull(comp_robustica_pca)
    unmapped_components = clustering_eval %>%
        filter(algorithm=="robustica_pca" & cluster_id%in%unmapped_clusters) %>%
        pull(component)
    
    clustering_eval %>%
        filter(algorithm=="icasso" & component%in%unmapped_components) %>%
        count(cluster_id) %>%
        print() # all unmapped components are in cluster 0 in icasso
    
    ## TODO: number of components in icasso and robustica_pca clusters
    
    ## Are modules enriched in similar numbers of gene sets?
    plts[["comparison_methods-enrichments_vs_freq"]] = enrichments %>%
        count(dataset, ontology, algorithm) %>%
        ggstripchart(x="algorithm", y="n", color="algorithm", palette=PAL_METRICS, size=0.1) +
        geom_boxplot(aes(color=algorithm), fill=NA, outlier.shape=NA, width=0.4) +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~ontology, scales="free_y") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        stat_compare_means(method="wilcox.test", size=2, family=FONT_FAMILY) +
        guides(color="none") +
        labs(x="Metric", y="Counts")
    
    plts[["comparison_methods-mapped_vs_freq"]] = enrichments %>%
        count(dataset, ontology, algorithm) %>%
        filter(algorithm=="robustica_pca") %>%
        left_join(mapping, by=c("dataset"="cluster_id", "algorithm")) %>%
        ggstripchart(x="mapped", y="n", color="mapped", palette="jco", size=0.1) +
        geom_boxplot(aes(color=mapped), fill=NA, outlier.shape=NA, width=0.4) +
        facet_wrap(~ontology, scales="free_y") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        stat_compare_means(method="wilcox.test", size=2, family=FONT_FAMILY) +
        guides(color="none") +
        labs(x="Mapped Robust Component", y="Counts")
    
    return(plts)
}


analyze_components = function(A_robustica_pca, metadata){
    
    # prep
    X = A_robustica_pca %>% column_to_rownames("index")
    met = metadata %>% column_to_rownames("sampleID")
    samples = rownames(X)
    
    # correlation with continuous sample properties (mitotic index)
    properties_oi = "mitotic_index"
    correls = sapply(properties_oi, function(property_oi){
        y = met[samples,property_oi]
        correls = apply(X[samples,,drop=FALSE], 2, function(x){
            cor(x,y, method='spearman', use='pairwise.complete.obs')
        })
    }, simplify=FALSE)
    correls = do.call(rbind, correls) %>%
        t() %>%
        as.data.frame() %>% 
        rownames_to_column("component") %>%
        rename(spearman_mitotic_index=mitotic_index)
    
    # survival associations
    survs = apply(X[samples,,drop=FALSE], 2, function(x){
        dat = cbind(weight=x, met[samples,c('OS','OS.time'),drop=FALSE])
        result = summary(coxph(Surv(OS.time, OS) ~ weight, dat))
        assoc = result$coefficients['weight','z']
        return(assoc)
    })
    survs = data.frame(coxph_zscore = survs) %>% 
        rownames_to_column("component")
    
    # association with histological_type, histological_grade, mutations
    components_oi = colnames(X)
    properties_oi = c("histological_type","histological_grade","IDH1","TP53")
    result = lapply(components_oi, function(comp_oi){
        res = lapply(properties_oi, function(property_oi){
            x = X[samples, comp_oi]
            g = met[samples, property_oi]
            test = kruskal.test(x=x, g=g)
            out = data.frame(
                component = comp_oi,
                property = property_oi, 
                pvalue = test[["p.value"]]
            )
            return(out)
        })
        res = do.call(rbind, res)
        return(res)
    })
    assocs = do.call(rbind, result) %>%
        group_by(property) %>%
        mutate(fdr = p.adjust(pvalue, method="fdr")) %>%
        ungroup() %>%
        mutate(label=paste0("assoc_fdr_",property)) %>%
        pivot_wider(component, names_from="label", values_from="fdr")
    
    results = correls %>%
        left_join(survs, by="component") %>%
        left_join(assocs, by="component")
    
    return(results)
}


plot_components_oi = function(snv, results, metadata, A_robustica_pca, 
                              S_robustica_pca, enrichments, modules_robustica_pca){
    
    plts = list()
    
    # find interesting component in common and interesting component in euclidean
    # mutagenesis IDH1, correlation with mitotic index, 
    X = snv %>%
        count(gene, effect) %>%
        group_by(effect) %>%
        arrange(n) %>%
        mutate(index=row_number()) %>%
        ungroup()
    
    plts[["components_oi-eda-snv"]] = X %>%
        ggplot(aes(x=index, y=n)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=8) +
        theme_pubr() +
        facet_wrap(~effect, scales="free_x") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text_repel(aes(label=gene), X %>% slice_max(n, n=2), size=2, 
                        family=FONT_FAMILY, segment.size=0.1) +
        labs(x="Index", y="Mutation Frequency")
    
    ## find components of interest
    ### a component with strong mitotic index, differential IDH1, prognosis
    ### a component with differential TP53
    components_oi = c(
        results %>% 
            filter(component!="-1") %>% 
            slice_max(spearman_mitotic_index, n=1) %>% 
            pull(component),
        results %>% 
            filter(!mapped & component!="-1") %>% 
            slice_min(assoc_fdr_TP53, n=1) %>% 
            pull(component)
    )
    
    X = metadata %>%
        left_join(
            A_robustica_pca %>% 
            dplyr::select(all_of(c("index",components_oi))) %>%
            pivot_longer(cols=all_of(components_oi), names_to="component", values_to="weight"), 
            by=c("sampleID"="index")
        ) %>%
        mutate(component = factor(component, levels=components_oi)) %>%
        group_by(component) %>%
        mutate(weight_group= ifelse(weight < median(weight), "Low", "High")) %>%
        ungroup()
               

    # mutations
    x = X %>%
        pivot_longer(cols=c("IDH1","TP53"), names_to="gene", values_to="is_mut") %>%
        mutate(is_mut = ifelse(is_mut, "Mutated", "WT"))
    plts[["components_oi-mut_vs_weights"]] = x %>%
        ggplot(aes(x=gene, y=weight, group=interaction(gene, is_mut))) +
        geom_violin(aes(fill=is_mut), color=NA) +
        geom_boxplot(width=0.1, outlier.size = 0.1, position=position_dodge(0.9)) +
        stat_compare_means(method="wilcox.test", label="p.signif", size=2, family=FONT_FAMILY) +
        fill_palette("lancet") +
        theme_pubr() +
        facet_wrap(~component) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY), aspect.ratio = 1) +
        labs(x="Gene", y="Weights Component", fill="Status") +
        geom_text(aes(y=-45, label=label), 
                  x %>% count(component, gene, is_mut) %>% mutate(label=paste0("n=",n)),
                  size=2, family=FONT_FAMILY, position=position_dodge(0.9))
    
    # histological type
    x = X %>%
        drop_na(histological_type)
        
    plts[["components_oi-histological_type_vs_weights"]] = x %>%
        ggplot(aes(x=component, y=weight, group=interaction(component, histological_type))) +
        geom_violin(aes(fill=histological_type), color=NA) +
        geom_boxplot(width=0.1, outlier.size = 0.1, position=position_dodge(0.9)) +
        stat_compare_means(method="kruskal.test", label="p.signif", size=2, family=FONT_FAMILY) +
        fill_palette("Greens") +
        theme_pubr() +
        labs(x="Component", y="Weights Component", fill="Histological Type") +
        geom_text(aes(y=-45, label=label), 
                  x %>% count(component, histological_type) %>% mutate(label=paste0("n=",n)),
                  size=2, family=FONT_FAMILY, position=position_dodge(0.9))
    
    # survival
    ## component 12
    component_oi = "12"
    x = X %>% filter(component==component_oi)
    fit = x %>% survfit(Surv(OS.time, OS) ~ weight_group, data=.)
    plts[["components_oi-surv_vs_weights-comp12"]] = x %>%
        ggsurvplot(
            data = .,
            fit = fit, 
            xlab='Time (Days)', ylab="Overall Survival Probability",
            legend.title = "Component Weight",
            pval = TRUE,
            conf.int = TRUE,
            palette = 'Dark2',
            ggtheme = theme_classic2(base_size=10, base_family = "Arial"),
            font.family = "Arial"
        ) + 
        labs(title=sprintf('KM from Weights Component %s',component_oi))
    plts[["components_oi-surv_vs_weights-comp12"]] = plts[["components_oi-surv_vs_weights-comp12"]][["plot"]]
    
    ## component 72
    component_oi = "72"
    x = X %>% filter(component==component_oi)
    fit = x %>% survfit(Surv(OS.time, OS) ~ weight_group, data=.)
    plts[["components_oi-surv_vs_weights-comp72"]] = x %>%
        ggsurvplot(
            data = .,
            fit = fit, 
            xlab='Time (Days)', ylab="Overall Survival Probability",
            legend.title = "Component Weight",
            pval = TRUE,
            conf.int = TRUE,
            palette = 'Dark2',
            ggtheme = theme_classic2(base_size=10, base_family = "Arial"),
            font.family = "Arial"
        ) + 
        labs(title=sprintf('KM from Weights Component %s',component_oi))
    plts[["components_oi-surv_vs_weights-comp72"]] = plts[["components_oi-surv_vs_weights-comp72"]][["plot"]]

    # mitotic index
    plts[["components_oi-mitotic_index_vs_weights"]] = X %>%
        ggplot(aes(x=mitotic_index, y=weight)) +
        geom_scattermore(pixels=c(1000, 1000), pointsize=4, alpha=0.5) +
        theme_pubr() +
        facet_wrap(~component) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY), aspect.ratio = 1) +
        stat_cor(method="spearman", size=2, family=FONT_FAMILY) +
        labs(x="Mitotic Index", y="Weights Component")
    
    # weight distibutions in source matrices
    X = S_robustica_pca %>%
        dplyr::select(all_of(components_oi)) %>%
        pivot_longer(cols=all_of(components_oi), names_to="component", values_to="weight")
    plts[["components_oi-source_weights"]] = X %>%
        ggviolin(x="component", y="weight", fill="orange", color=NA) +
        geom_boxplot(width=0.1, outlier.size=0.1) +
        labs(x="Component", y="Weights Component")
    
    # enrichments
    X = enrichments %>%
        mutate(component = factor(as.character(dataset), levels=components_oi)) %>%
        filter(algorithm=="robustica_pca") %>%
        filter(component%in%components_oi) %>%
        rowwise() %>%
        mutate(gene_ratio = eval(parse(text=GeneRatio))) %>%
        ungroup() 
    
    gene_sets_oi = X %>%
        group_by(component, ontology) %>%
        slice_max(Count, n=5) %>%
        ungroup() %>%
        pull(Description) %>%
        unique()
    
    n_genes = modules_robustica_pca %>%
        as.data.frame() %>%
        dplyr::select(all_of(components_oi)) %>%
        pivot_longer(all_of(components_oi), names_to="component", values_to="in_module") %>%
        filter(in_module) %>%
        count(component)
    
    plts[["components_oi-enrichments"]] = X %>%
        filter(Description %in% gene_sets_oi) %>%
        left_join(n_genes, by="component") %>%
        mutate(component = sprintf("%s\n(n=%s)", component, n),
               component = factor(component, levels=sort(unique(component), decreasing=TRUE))) %>%
        ggplot(aes(x=component, y=Description)) +
        geom_point(aes(size=gene_ratio, color=p.adjust)) +
        scale_size(range=c(0.5,3)) + 
        scale_color_continuous(
            low=PAL_FDR_LIGHT, high=PAL_FDR_DARK, guide=guide_colorbar(reverse=TRUE)) +
        theme_pubr() +
        labs(x="Component", y="Enriched Term", size="Gene Ratio", color="FDR") +
        facet_wrap(~ontology, scales="free_y") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    return(plts)
}


make_plots = function(performance, clustering_eval, comparison_pearson, 
                      comparison_jaccard, mapping,
                      snv, results, metadata, A_robustica_pca, S_robustica_pca, 
                      enrichments, modules_robustica_pca){
    plts = list(
        plot_performance(performance),
        plot_clustering_eval(clustering_eval),
        plot_comparison_methods(comparison_pearson, comparison_jaccard, mapping, clustering_eval, enrichments),
        plot_components_oi(snv, results, metadata, A_robustica_pca, 
                           S_robustica_pca, enrichments, modules_robustica_pca)
    )
    plts = do.call(c,plts)
    return(plts)
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
    save_plt(plts,'performance-time_vs_memory','.pdf',figs_dir, width=5, height=15)
    
    save_plt(plts,'clustering_eval-metric_vs_cluster_std','.pdf',figs_dir, width=5, height=5)
    save_plt(plts,'clustering_eval-metric_vs_silhouette','.pdf',figs_dir, width=5, height=5)
    save_plt(plts,'clustering_eval-metrics_vs_cluster_size','.pdf',figs_dir, width=5, height=5)
    
    save_plt(plts,'comparison_methods-pearson-heatmap','.pdf',figs_dir, width=20, height=13)
    save_plt(plts,'comparison_methods-jaccard-heatmap','.pdf',figs_dir, width=20, height=13)
    save_plt(plts,'comparison_methods-mapping-hist','.pdf',figs_dir, width=5, height=5)
    save_plt(plts,'comparison_methods-mapping_vs_silhouette','.pdf',figs_dir, width=5, height=5)
    save_plt(plts,'comparison_methods-mapping_vs_cluster_std','.pdf',figs_dir, width=5, height=5)
    save_plt(plts,'comparison_methods-enrichments_vs_freq','.pdf',figs_dir, width=9, height=5)
    save_plt(plts,'comparison_methods-mapped_vs_freq','.pdf',figs_dir, width=9, height=5)
    
    save_plt(plts,'components_oi-eda-snv','.pdf',figs_dir, width=15, height=15)
    save_plt(plts,'components_oi-mut_vs_weights','.pdf',figs_dir, width=8, height=8)
    save_plt(plts,'components_oi-histological_type_vs_weights','.pdf',figs_dir, width=8, height=5)
    save_plt(plts,'components_oi-surv_vs_weights-comp12','.pdf',figs_dir, width=5, height=6.5)
    save_plt(plts,'components_oi-surv_vs_weights-comp72','.pdf',figs_dir, width=5, height=6.5)
    save_plt(plts,'components_oi-mitotic_index_vs_weights','.pdf',figs_dir, width=8, height=5)
    save_plt(plts,'components_oi-source_weights','.pdf',figs_dir, width=6, height=5)
    save_plt(plts,'components_oi-enrichments','.pdf',figs_dir, width=22, height=6)
}


parseargs = function(){
    
    option_list = list( 
        make_option("--performance_file", type="character"),
        make_option("--clustering_file", type="character"),
        make_option("--S_stds_file", type="character"),
        make_option("--S_icasso_file", type="character"),
        make_option("--S_std_icasso_file", type="character"),
        make_option("--S_robustica_pca_file", type="character"),
        make_option("--S_std_robustica_pca_file", type="character"),
        make_option("--A_robustica_pca_file", type="character"),
        make_option("--genexpr_file", type="character"),
        make_option("--snv_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--sample_indices_file", type="character"),
        make_option("--enrichment_icasso_file", type="character"),
        make_option("--enrichment_robustica_pca_file", type="character"),
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
    S_icasso_file = args[["S_icasso_file"]]
    S_std_icasso_file = args[["S_std_icasso_file"]]
    S_robustica_pca_file = args[["S_robustica_pca_file"]]
    S_std_robustica_pca_file = args[["S_std_robustica_pca_file"]]
    A_robustica_pca_file = args[["A_robustica_pca_file"]]
    genexpr_file = args[["genexpr_file"]]
    snv_file = args[["snv_file"]]
    metadata_file = args[["metadata_file"]]
    sample_indices_file = args[["sample_indices_file"]]
    enrichment_icasso_file = args[["enrichment_icasso_file"]]
    enrichment_robustica_pca_file = args[["enrichment_robustica_pca_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load data
    genexpr = read_tsv(genexpr_file)
    snv = read_tsv(snv_file)
    metadata = read_tsv(metadata_file)
    sample_indices = read_tsv(sample_indices_file)
    
    performance = read_tsv(performance_file) %>% 
        filter(dataset==DATASET_OI & algorithm%in%METRICS_OI) %>%
        mutate(algorithm = factor(algorithm, levels=METRICS_OI))
    clustering = read_tsv(clustering_file) %>% 
        filter(dataset==DATASET_OI & algorithm%in%METRICS_OI) %>%
        mutate(algorithm = factor(algorithm, levels=METRICS_OI))
    S_stds = read_tsv(S_stds_file) %>% 
        filter(dataset==DATASET_OI & algorithm%in%METRICS_OI) %>%
        mutate(algorithm = factor(algorithm, levels=METRICS_OI))
    
    S_icasso = read_tsv(S_icasso_file)
    S_robustica_pca = read_tsv(S_robustica_pca_file)

    S_std_icasso = read_tsv(S_std_icasso_file)
    S_std_robustica_pca = read_tsv(S_std_robustica_pca_file)
    
    A_robustica_pca = read_tsv(A_robustica_pca_file)
    
    enrichment_icasso = read_tsv(enrichment_icasso_file) %>% 
        mutate(algorithm="icasso")
    enrichment_robustica_pca = read_tsv(enrichment_robustica_pca_file) %>% 
        mutate(algorithm="robustica_pca")
    
    # preprocess
    ## relative time
    performance = performance %>% 
        filter(`function` != 'evaluate_performance') %>%
        group_by(algorithm) %>% 
        mutate(rel_timestamp = timestamp - min(timestamp),
               `function` = paste0(`function`,'()')) %>%
        ungroup()
    
    ## clustering time and evaluation    
    clustering_time = performance %>% 
        filter(`function`=='cluster_components()') %>% 
        group_by(algorithm) %>%
        summarize(time=max(time)) %>%
        distinct()
    
    max_memory = performance %>% 
        group_by(algorithm) %>% 
        summarize(max_memory=max(memory))
    
    clustering_eval = clustering %>%
        left_join(clustering_time, by='algorithm') %>%
        left_join(max_memory, by='algorithm') %>%
        left_join(S_stds, by=c("algorithm","cluster_id","dataset")) %>%
        group_by(algorithm, cluster_id) %>%
        mutate(mean_silhouette = mean(silhouette_euclidean),
               cluster_size = n()) %>%
        ungroup()
    
    ## icasso vs robustica_pca
    ### Pearson correlation between components
    comparison_pearson = cor(S_icasso[,-1], S_robustica_pca[,-1], method="pearson")

    ### Jaccard distances between modules
    modules_icasso = apply(S_icasso[,-1], 2, define_module)
    modules_robustica_pca = apply(S_robustica_pca[,-1], 2, define_module)
    comparison_jaccard = simil(modules_icasso, modules_robustica_pca, 
                               method='Jaccard', by_rows=FALSE)
    ### map robustica to icasso
    mapping = comparison_pearson %>% 
        as.data.frame() %>%
        rownames_to_column("comp_icasso") %>%
        pivot_longer(-comp_icasso, names_to="comp_robustica_pca", values_to="correlation") %>%
        group_by(comp_robustica_pca) %>%
        slice_max(abs(correlation), n=1) %>%
        ungroup() %>%
        mutate(mapped = abs(correlation) > 0.5,
               cluster_id = as.numeric(comp_robustica_pca)) %>%
        left_join(
            clustering_eval %>%
                distinct(algorithm, cluster_id, mean_silhouette, cluster_avg_std) %>%
                filter(algorithm=="robustica_pca"),
            by="cluster_id"
        )
    
    ## enrichments
    enrichments = enrichment_icasso %>%
        bind_rows(enrichment_robustica_pca)
    
    ## metadata
    ### add mutations
    genes_oi = c("TP53","IDH1")
    muts_oi = snv %>%
        filter(effect=="Missense_Mutation" & gene%in%genes_oi) %>%
        mutate(is_mut=TRUE) %>%
        pivot_wider(sampleID, names_from="gene", values_from="is_mut") 
    
    metadata = metadata %>%
        left_join(muts_oi, by="sampleID") %>%
        mutate(IDH1 = replace_na(IDH1, FALSE),
               TP53 = replace_na(TP53, FALSE))
    
    ### add sample indices
    metadata = metadata %>%
        left_join(sample_indices, by="sampleID")
    
    ### analyze weights mixing matrix
    results = analyze_components(A_robustica_pca, metadata)
    results = results %>%
        left_join(mapping %>% distinct(comp_robustica_pca, mapped), 
                  by=c("component"="comp_robustica_pca"))
    
    # visualize
    plts = make_plots(performance, clustering_eval, comparison_pearson, 
                      comparison_jaccard, mapping,
                      snv, results, metadata, A_robustica_pca, S_robustica_pca, 
                      enrichments, modules_robustica_pca)
    
    # save
    save_plots(plts, figs_dir)
    # save_figdata(figdata, figs_dir)
}

###### SCRIPT ######
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}