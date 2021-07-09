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
require(proxy)
require(fdrtool)
require(latex2exp)
require(writexl)


ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
ALGORITHMS = c('icasso','robustica_nosign','robustica','robustica_pca')

# Development
# -----------
# RESULTS_DIR = file.path(ROOT,'results','benchmark_sign_inference','files')
# performance_evaluation_file = file.path(RESULTS_DIR,'Sastry2019','performance_evaluation.tsv.gz')
# clustering_info_file = file.path(RESULTS_DIR,'Sastry2019','clustering_info.tsv.gz')

# algorithms = c('icasso','robustica_nosign','robustica','robustica_pca')
# S_stds_files = unlist(strsplit(file.path(RESULTS_DIR,'Sastry2019',paste0(algorithms,'-','S_std.tsv.gz')),','))

# Ss_files = unlist(strsplit(file.path(RESULTS_DIR,'Sastry2019',paste0(algorithms,'-','S.tsv.gz')),','))

# As_files = unlist(strsplit(file.path(RESULTS_DIR,'Sastry2019',paste0(algorithms,'-','A.tsv.gz')),','))

# figs_dir = file.path(ROOT,'results','benchmark_sign_inference','figures','Sastry2019')
# files_dir = file.path(ROOT,'results','benchmark_sign_inference','files','Sastry2019')


##### FUNCTIONS #####
get_relative_timestamp = function(x){
    x0=x[1]
    return(x-x0)
}


load_S_stds = function(files){
    Ss = lapply(files, function(file){
        read_tsv(file) %>% 
        column_to_rownames('log-TPM') %>%
        melt() %>% 
        mutate(algorithm = gsub('-.*','',basename(file)))
    })
    Ss = do.call(rbind, Ss)
    return(Ss)
}


load_mats = function(files,rownms){
    mats = lapply(files, function(file){
        df = read_tsv(file) %>% 
        column_to_rownames(rownms)
        colnames(df) = paste0('comp',colnames(df))
        return(df)
    })
    names(mats) = gsub('-.*','',basename(files))
    return(mats)
}


plot_corr_vs_euc = function(len=10){
    set.seed(123)
    compA = rnorm(len)
    compB = rnorm(len)
    compC = compA - 6
    df = data.frame(
        feature = 1:length(compA),
        compA,
        compB,
        compC
    ) %>% melt(id.vars = 'feature')
    
    metrics = data.frame(
        cor_AB = cor(compA, compB, method = 'pearson'),
        cor_AC = cor(compA, compC, method = 'pearson'),
        cor_BC = cor(compB, compC, method = 'pearson'),
        dist_AB = dist(rbind(compA, compB), method = 'euclidean')[1],
        dist_AC = dist(rbind(compA, compC), method = 'euclidean')[1],
        dist_BC = dist(rbind(compB, compC), method = 'euclidean')[1]
    ) %>% melt()
    
    plts = list()
    
    plts[['corr_vs_euc-scatter']] = df %>% 
        ggline(x='feature', y='value', color='variable', 
               linetype='dashed', palette='simpsons') + 
        labs(x="Components' Feature", y='Weight', color='Component') + 
        theme(aspect.ratio=1)
    
    plts[['corr_vs_euc-table']] = ggtexttable(metrics, rows = NULL)
    
    return(plts)
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


plot_clustering_evaluation = function(df){
    # get mean silhouette per cluster
    X = df
    
    plts = list()
    plts[['clustering-silhouettes-violins']] = X %>%
        ggviolin(x='algorithm', y='silhouette_pearson', trim = TRUE,
                 fill='algorithm', color=NA, palette='Paired') + 
        guides(fill=FALSE) +
        geom_boxplot(width=0.1) +
        labs(x='Algorithm', y='Silhouette Score') +
        theme_pubr(x.text.angle = 70)
     
    plts[['clustering-S_stds-violins']] = X %>% 
        ggviolin(x='algorithm', y='mean_std_cluster', trim = TRUE,
                 fill='algorithm', color=NA, palette='Paired') + 
        geom_boxplot(width=0.1) +
        guides(fill=FALSE) +
        labs(x='Algorithm', y='Cluster Mean Std.') +
        theme_pubr(x.text.angle = 70)
    
     plts[['clustering-silhouettes_vs_stds-scatter']] = X %>%
        ggplot(aes(x=silhouette_pearson, y=mean_std_cluster)) +
        geom_text(aes(label=cluster_id)) +
        facet_wrap(~algorithm, scales='free') + 
        theme_pubr(border=TRUE) +
        labs(x='Silhouette Score', y='Cluster Mean Std.')
    
    return(plts)
}


define_module = function(x, cutoff=0.01){
    fdr = fdrtool(x, plot=FALSE, cutoff.method="fndr", verbose=FALSE)[['qval']]
    x = fdr < cutoff
    return(x)
}


plot_module_comparisons = function(Ss, As){
    plts = list()
    
    # get modules from source matrix
    modules = sapply(Ss, simplify=FALSE, function(S){
        S %>% mutate_all(define_module)
    })
    
    # module sizes
    module_size = data.frame(
        module_icasso = 1:ncol(modules[['icasso']]),
        icasso = colSums(modules[['icasso']]),
        robustica = colSums(modules[['robustica_pca']])
    ) %>% melt() 
    
    # module sizes are similar
    plts[['module_comparisons-module_size-boxplot']] = module_size %>% 
        filter(variable %in% c('icasso','robustica')) %>%
        ggstripchart(x='variable', y='value') + 
        geom_boxplot(width=0.1, alpha=0) + 
        labs(x='Algorithm', y='Module Size') +
        stat_compare_means(method='wilcox.test')
    
    # use jaccard similarity to map modules
    sim = simil(modules[['icasso']], modules[['robustica_pca']], 
                method='Jaccard', by_rows=FALSE) %>% as.matrix()     
    plts[['module_comparisons-jaccard']] = sim %>% pheatmap(silent=TRUE)
    
    # map modules in icasso to robustica
    ## mapping
    modules_info = apply(sim,1,which.max) %>% 
        enframe(name = 'module_icasso', value = 'module_robustica') %>%
        mutate(module_robustica = paste0('comp',module_robustica-1)) %>%
        left_join(
            apply(sim,1,max) %>% 
            enframe(name='module_icasso', value='jaccard'),
            by='module_icasso'
        ) %>% left_join(
            data.frame(
                module_icasso = colnames(modules[['icasso']]),
                module_icasso_size = colSums(modules[['icasso']])
            ),
            by='module_icasso'
        ) %>% left_join(
            data.frame(
                module_robustica = colnames(modules[['robustica_pca']]),
                module_robustica_size = colSums(modules[['robustica_pca']])
            ),
            by='module_robustica'
        ) %>% arrange(-jaccard)
    
    ## modules in icasso not found by robustica
    notfound = modules_info %>% 
        count(module_robustica) %>% 
        filter(n>1) %>% 
        pull(module_robustica)
    modules_icasso = modules_info %>% 
        filter(module_robustica %in% notfound) %>% 
        group_by(module_robustica) %>%
        slice_min(order_by = jaccard, n = 1) %>%
        pull(module_icasso)
    
    plts[['module_comparisons-mapping_icasso-module_sizes']] = modules_info %>% 
        ggscatter(x='module_icasso_size','module_robustica_size',color='jaccard', 
                  size='jaccard', alpha=0.75) +
        geom_text_repel(
            aes(label=module_icasso), 
            modules_info %>% filter(module_icasso %in% modules_icasso),
            box.padding = 0.5) + 
        scale_color_gradient2(low = 'blue', mid='grey', high = 'red', midpoint = 0.5) + 
        guides(size=FALSE) + 
        labs(x='Module Size Icasso', y='Module Size robustica', color='Jaccard Sim.',
             title=TeX('Mapping gene modules in \\textit{Icasso} to \\textit{robustica}'))
    
    X = As[['icasso']][,modules_icasso] %>% 
        rownames_to_column('sample') %>% melt()
    plts[['module_comparisons-mapping_icasso-unmapped']] = X %>% 
        ggviolin(x='variable', y='value') + 
        geom_boxplot(width=0.1) + 
        geom_text_repel(aes(label=sample), 
                        X %>% group_by(variable) %>% 
                        slice_max(order_by = abs(value), n=10)) +
        labs(x='Component', y='Weight in Mixing Matrix (A)',
             title=TeX('Modules detected only by \\textit{Icasso}'))
    
    
    # map modules in robustica to icasso
    ## mapping
    modules_info = apply(sim,2,which.max) %>% 
        enframe(name = 'module_robustica', value = 'module_icasso') %>%
        mutate(module_icasso = paste0('comp',module_icasso-1)) %>%
        left_join(
            apply(sim,2,max) %>% 
            enframe(name='module_robustica', value='jaccard'),
            by='module_robustica'
        ) %>% left_join(
            data.frame(
                module_icasso = colnames(modules[['icasso']]),
                module_icasso_size = colSums(modules[['icasso']])
            ),
            by='module_icasso'
        ) %>% left_join(
            data.frame(
                module_robustica = colnames(modules[['robustica_pca']]),
                module_robustica_size = colSums(modules[['robustica_pca']])
            ),
            by='module_robustica'
        ) %>% arrange(-jaccard)

    ## modules in robustica not found by icasso
    notfound = modules_info %>% 
        count(module_icasso) %>% 
        filter(n>1) %>% 
        pull(module_icasso)
    modules_robustica = modules_info %>% 
        filter(module_icasso %in% notfound) %>% 
        group_by(module_icasso) %>%
        slice_min(order_by = jaccard, n = 1) %>%
        pull(module_robustica)
    
    ## icasso vs robustica: modules different when small.
    plts[['module_comparisons-mapping_robustica-module_sizes']] = modules_info %>% 
        ggscatter(x='module_icasso_size','module_robustica_size',color='jaccard', 
                  size='jaccard', alpha=0.75) + 
        geom_text_repel(
            aes(label=module_robustica), 
            modules_info %>% filter(module_robustica %in% modules_robustica),
            box.padding = 0.5) + 
        scale_color_gradient2(low = 'blue', mid='grey', high = 'red', midpoint = 0.5) + 
        guides(size=FALSE) + 
        labs(x='Module Size Icasso', y='Module Size robustica', color='Jaccard Sim.',
             title=TeX('Mapping gene modules in \\textit{robustica} to \\textit{Icasso}'))
    
    X = As[['robustica_pca']][,modules_robustica] %>% 
        rownames_to_column('sample') %>% melt()
    plts[['module_comparisons-mapping_robustica-unmapped']] = X %>% 
        ggviolin(x='variable', y='value') + 
        geom_boxplot(width=0.1) + 
        geom_text_repel(aes(label=sample), 
                        X %>% group_by(variable) %>% 
                        slice_max(order_by = abs(value), n=10)) +
        labs(x='Component', y='Weight in Mixing Matrix (A)',
             title=TeX('Modules detected only by \\textit{robustica}'))
    
    return(plts)
}


make_plots = function(performance_evaluation, clustering_evaluation, Ss, As){
    plts = list(
        plot_corr_vs_euc(),
        plot_performance_profile(performance_evaluation),
        plot_clustering_evaluation(clustering_evaluation),
        plot_module_comparisons(Ss, As)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(performance_evaluation, clustering_evaluation, Ss, As){
    
    # prep
    ## gene module evaluation
    modules = sapply(Ss, simplify=FALSE, function(S){
        S %>% 
        rownames_to_column("gene") %>%
        mutate_at(vars(-("gene")), define_module) %>%
        column_to_rownames("gene")
    })
    sim = simil(modules[['icasso']], modules[['robustica_pca']], 
                method='Jaccard', by_rows=FALSE) %>% as.matrix()  
    sim = matrix(sim, nrow = nrow(sim), ncol=ncol(sim), dimnames = dimnames(sim)) %>% 
        as.data.frame()
    
    # prepare figdata object
    figdata = list()
    figdata[['benchmark_sign_inference-ICA-source_matrices']] = sapply(
        Ss, simplify=FALSE, function(S){
            S %>% rownames_to_column('gene')
        })
    figdata[['benchmark_sign_inference-ICA-mixing_matrices']] = sapply(
        As, simplify=FALSE, function(A){
            A %>% rownames_to_column('sample')
        })
    figdata[['benchmark_sign_inference-evaluation']] = list(
        performance_evaluation = performance_evaluation,
        clustering_evaluation = clustering_evaluation,
        gene_modules_icasso = modules[['icasso']] %>% rownames_to_column('gene'),
        gene_modules_robustica_pca = modules[['robustica_pca']] %>% rownames_to_column('gene'),
        gene_modules_jaccard_similarity = sim %>% rownames_to_column('component')
    )
    
    return(figdata)
}


save_plot = function(plt, plt_name, extension='.pdf', 
                      directory='', dpi=350, 
                      width = par("din")[1], height = par("din")[2], units='cm'){
        filename = file.path(directory,paste0(plt_name,extension))
        ggsave(filename, plt, width=width, height=height, 
               dpi=dpi, limitsize=FALSE, units=units)
}


save_plots = function(plts, figs_dir){
    save_plot(plts[['corr_vs_euc-scatter']],'corr_vs_euc-scatter','.pdf',figs_dir, width=12,height=12)
    save_plot(plts[['corr_vs_euc-table']],'corr_vs_euc-table','.pdf',figs_dir, width=8,height=8)
    
    save_plot(plts[['mem_time-scatter']],'mem_time-scatter','.png',figs_dir, width=12,height=20)
    
    save_plot(plts[['clustering-silhouettes-violins']],'clustering-silhouettes-violins','.pdf',figs_dir, width=8,height=8)
    save_plot(plts[['clustering-S_stds-violins']],'clustering-S_stds-violins','.pdf',figs_dir, width=8,height=8)
    save_plot(plts[['clustering-silhouettes_vs_stds-scatter']],'clustering-silhouettes_vs_stds-scatter','.png',figs_dir, width=20,height=20)
    
    save_plot(plts[['module_comparisons-module_size-boxplot']],'module_comparisons-module_size-boxplot','.pdf',figs_dir, width=10,height=10)
    
    save_plot(plts[['module_comparisons-jaccard']],'module_comparisons-jaccard','.png',figs_dir, width=30, height=30)
    
    save_plot(plts[['module_comparisons-mapping_icasso-module_sizes']],'module_comparisons-mapping_icasso-module_sizes','.pdf',figs_dir, width=15, height=15)
    save_plot(plts[['module_comparisons-mapping_robustica-module_sizes']],'module_comparisons-mapping_robustica-module_sizes','.pdf',figs_dir, width=15, height=15)
    
    save_plot(plts[['module_comparisons-mapping_icasso-unmapped']],'module_comparisons-mapping_icasso-unmapped','.pdf',figs_dir, width=15, height=15)
    save_plot(plts[['module_comparisons-mapping_robustica-unmapped']],'module_comparisons-mapping_robustica-unmapped','.pdf',figs_dir, width=15, height=15)
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

    performance_evaluation_file = args$performance_evaluation_file
    clustering_info_file = args$clustering_info_file
    S_stds_files = unlist(strsplit(args$S_stds_files,','))
    Ss_files = unlist(strsplit(args$Ss_files,','))
    As_files = unlist(strsplit(args$As_files,','))
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load data
    performance_evaluation = read_tsv(performance_evaluation_file) %>% 
        filter(`function` != 'evaluate_performance') %>%
        group_by(algorithm) %>% 
        mutate(rel_timestamp = get_relative_timestamp(timestamp),
               `function` = paste0(`function`,'()'),
               algorithm = factor(algorithm, levels = ALGORITHMS))
    
    clustering_info = read_tsv(clustering_info_file) %>%
        mutate(algorithm = factor(
            algorithm, levels = ALGORITHMS)
        )
    
    S_stds = load_S_stds(S_stds_files) %>%
        mutate(algorithm = factor(algorithm, levels = ALGORITHMS))
    
    Ss = load_mats(Ss_files,'log-TPM')
    As = load_mats(As_files,'index')
    
    # combine clustering_info and S_stds
    clustering_evaluation = clustering_info %>%
        group_by(algorithm, cluster_id) %>%
        summarise(
            silhouette_pearson = mean(silhouette_pearson),
            silhouette_euclidean = mean(silhouette_euclidean)
        ) %>% 
        mutate(cluster_id = as.factor(cluster_id)) %>%
        left_join(
            S_stds %>%
            group_by(algorithm, variable) %>%
            summarise(mean_std_cluster = mean(value)) %>%
            dplyr::rename(cluster_id=variable),
            by = c('algorithm','cluster_id')
        )
    
    plts = make_plots(performance_evaluation, clustering_evaluation, Ss, As)
    figdata = make_figdata(performance_evaluation, clustering_evaluation, Ss, As)
    
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}