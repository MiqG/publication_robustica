#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# visualize benchmark of ICA iterations.
#

require(tidyverse)
require(reshape)
require(ggpubr)
require(proxy)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
REFERENCES = list(
    'LGG' = 'max_iter=10000 & tol=1e-1',
    'Sastry2019' = 'Sastry2019_ref'
)

# Development
# -----------
# RESULTS_DIR = file.path(ROOT,'results','benchmark_tol_vs_maxiter','files')
# dirnames = 'evaluate_iterations-Sastry2019_10000_1e-1,evaluate_iterations-Sastry2019_10000_1e-2,evaluate_iterations-Sastry2019_10000_1e-3,evaluate_iterations-Sastry2019_1000_1e-1,evaluate_iterations-Sastry2019_1000_1e-2,evaluate_iterations-Sastry2019_1000_1e-3,evaluate_iterations-Sastry2019_1000_1e-4,evaluate_iterations-Sastry2019_100_1e-1,evaluate_iterations-Sastry2019_100_1e-2,evaluate_iterations-Sastry2019_100_1e-3,evaluate_iterations-Sastry2019_100_1e-4,evaluate_iterations-Sastry2019_200_1e-2,evaluate_iterations-Sastry2019_200_1e-3,evaluate_iterations-Sastry2019_200_1e-4,evaluate_iterations-Sastry2019_500_1e-3,evaluate_iterations-Sastry2019_500_1e-4'
# dirnames = unlist(strsplit(dirnames, ','))
# evaluation_dirs = file.path(RESULTS_DIR,dirnames)
# figs_dir = file.path(ROOT,'results','benchmark_tol_vs_maxiter','figures')
# dataset = 'Sastry2019'
# S_ref_file = '~/projects/publication_robustica/data/prep/ica/Sastry2019-S.tsv.gz'


##### FUNCTIONS #####
get_label = function(dir_name, dataset){
    lab = gsub(paste0('evaluate_iterations-',dataset,'_'),'',basename(dir_name))
    lab = unlist(strsplit(lab,'_'))
    lab = apply(cbind(c('max_iter=','tol='),lab),1, paste, collapse='')
    lab = paste(lab, collapse=' & ')
    return(lab)
}


load_datatype = function(evaluation_dirs, filename, dataset, bind=TRUE){
    df = lapply(evaluation_dirs, function(dir_name){
        lab = get_label(dir_name, dataset)
        df = read_tsv(file.path(dir_name, filename))
        df$params = lab
        return(df)
    })
    if(bind){ df = do.call(rbind, df) }
    return(df)
}


load_components = function(evaluation_dirs, filename, dataset, bind=TRUE){
    idx_col = colnames(read_tsv(file.path(evaluation_dirs[1], filename), n_max = 1))[1]
    dfs = lapply(evaluation_dirs, function(dir_name){
        lab = get_label(dir_name, dataset)
        components = read_tsv(file.path(dir_name, filename)) %>% 
            rename_at(vars(-idx_col),function(x) paste0(x,'_',lab))
        metadata = data.frame(id=colnames(components), params=lab) %>% filter(id!=idx_col)
        gc()
        return(list(components=components, metadata=metadata))
    })
    if(bind){ 
        components = lapply(dfs, function(df){df[['components']]})
        components = do.call(cbind, components)
        to_drop = grep(idx_col,colnames(components))[-1]
        components = components[,-to_drop]
        metadata = lapply(dfs, function(df){df[['metadata']]})
        metadata = do.call(rbind, metadata)            
    }
    return(list(components=components, metadata=metadata))
}


load_data = function(evaluation_dirs, dataset, S_ref_file=NULL){
    # summary
    summary = load_datatype(evaluation_dirs, dataset, filename='summary.tsv.gz')
    summary = summary %>% mutate(has_converged = convergence_score < tol)
    
    # clustering stats
    clustering_stats = load_datatype(evaluation_dirs, dataset, filename='clustering_stats.tsv.gz')
    
    # S
    S = load_components(evaluation_dirs, dataset, filename='S.tsv.gz')
    
    data = list(
        summary = summary,
        clustering_stats = clustering_stats,
        components = S[['components']],
        metadata = S[['metadata']]
    )
    
     # S_ref
    if(dataset=='Sastry2019'){ 
        lab = 'Sastry2019_ref'
        idx_col = 'index'
        components = read_tsv(S_ref_file) %>% 
            rename_at(vars(-idx_col),function(x) paste0(x,'_',lab))
        metadata = data.frame(id=colnames(components), params=lab) %>% filter(id!=idx_col)
        data[['components']] = cbind(data[['components']], components[,setdiff(colnames(components),idx_col)])
        data[['metadata']] = rbind(data[['metadata']], metadata)
    }
    
                      
    # prepare output
    data[['summary']] = data[['summary']] %>% mutate(params=as.factor(params))
    data[['clustering_stats']] = data[['clustering_stats']] %>% mutate(params=as.factor(params))
    data[['metadata']] = data[['metadata']] %>% mutate(params=as.factor(params))
    return(data)
}


plot_summary = function(summary){
    X = summary
    plts = list()
    plts[['summary-has_converged']] = ggbarplot(
        X %>% count(params,has_converged), x='params', y='n', 
        fill='has_converged', label=TRUE, lab.size=3,
        color=NA, palette='npg', position=position_dodge(0.9)
    ) + 
    theme_pubr(x.text.angle=70) + 
    labs(x='Parameters', y='Count', fill='Has Converged')
    
    plts[['summary-iterations']] = ggstripchart(
        X, x='params', y='iteration_ica', 
        color='has_converged', palette='npg', alpha=0.5
    ) + 
    theme_pubr(x.text.angle=70) + 
    labs(x='Parameters', y='N. Iterations ICA', color='Has Converged')
    
    plts[['summary-convergence_score']] = ggstripchart(
        X, x='params', y='convergence_score', 
        color='has_converged', palette='npg', alpha=0.5
    ) + 
    yscale('log10') + 
    theme_pubr(x.text.angle=70) + 
    labs(x='Parameters', y='Convergence Score', color='Has Converged')
    
    plts[['summary-convergence_time']] = ggstripchart(
        X, x='params', y='time_ica', 
        color='has_converged', palette='npg', alpha=0.5
    ) + 
    theme_pubr(x.text.angle=70) + 
    labs(x='Parameters', y='Time ICA (s)', color='Has Converged')
    
    plts[['summary-robustica_time']] = ggstripchart(
        X %>% count(time_robustica,params,has_converged) %>% filter(has_converged), 
        x='params', y='time_robustica', 
        size='n'
    ) + 
    theme_pubr(x.text.angle=70) + 
    labs(x='Parameters', y='Time RobustICA (s)', size='N. Converged')
    
    return(plts)
}


plot_clustering_stats = function(clustering_stats){
    X = clustering_stats %>% mutate(is_noise=cluster_id==-1)
    plts = list()
    plts[['clustering-n_clusters']] = ggstripchart(
        X %>% count(params), x='params', y='n'
    ) + 
    theme_pubr(x.text.angle=70) +
    labs(x='Parameters', y='Number of Robust Components')
    
    plts[['clustering-cluster_size']] = ggstripchart(
        X, x='params', y='cluster_size',
        color='is_noise', palette='jco', alpha=0.5
    ) + 
    yscale('log10') +
    theme_pubr(x.text.angle=70) +
    labs(x='Parameters', y='Robust Component Size', color='Is Noise')
    
    plts[['clustering-cluster_mean_std']] = ggstripchart(
        X, x='params', y='S_mean_std',
        color='is_noise', palette='jco', alpha=0.5
    ) + 
    yscale('log10') +
    theme_pubr(x.text.angle=70) +
    labs(x='Parameters', y='Robust Component Mean Std.', color='Is Noise')
    
    
    return(plts)
}


define_module = function(x){
    # outliers
    mu = mean(x)
    std = sd(x)
    x = abs(x) > (mu+std*3)
    return(x)
}


plot_components_corr = function(components, metadata, dataset){
    
    # init
    reference = REFERENCES[[dataset]]
    idx_col = colnames(components)[1]
    queries = setdiff(unique(metadata$params),reference)
    plts = list()
    
    # correlation between robust components of best quality params and rest
    comparison = lapply(queries, function(query){
        cols_ref = metadata %>% filter(params %in% reference) %>% pull(id)
        cols_query = metadata %>% filter(params %in% query) %>% pull(id)
        sim = cor(components[,cols_ref], 
                  components[,cols_query], 
                  method = "pearson")
        
        # high similarity --> likely to be the "same" component
        # with respect to the reference
        max_correlations = apply(abs(sim), 1, max)
        
        df = data.frame(max_correlations=sort(max_correlations, decreasing=TRUE), 
                        reference=reference, 
                        query=query)
        
        # indicate only the possible maximum correlations
        max_detected = min(dim(sim))
        df$detected = TRUE
        if(max_detected < nrow(sim)){ df$detected[max_detected:nrow(df)] = FALSE }
        
        return(df)
    })
    comparison = do.call(rbind,comparison)
    
    plts[['components-detected']] = ggbarplot(
        comparison %>% count(query,detected), x='query', y='n', 
        fill='detected', label=TRUE, lab.size=3,
        color=NA, palette='uchicago', position=position_dodge(0.9)
    ) + 
    theme_pubr(x.text.angle=70) + 
    labs(x='Parameters', y='Count', fill='Detected')
    
    plts[['components-best_correlation']] = ggstripchart(
        comparison, x='query', y='max_correlations', 
        color='detected', palette='uchicago',
        alpha=0.5
    ) +
    theme_pubr(x.text.angle=70) + 
    labs(x='Parameters', y='Best Correlation', color='Detected')
    
    # module overlap
    ## define modules
    modules = components %>% 
        dplyr::select(-one_of(idx_col)) %>% 
        mutate_all(define_module)
    
    ## selected items (genes) per module per params
    X = enframe(apply(modules, 2, sum), 'id', 'n') %>% 
        left_join(metadata, by='id')
    
    plts[['components-module_size']] = ggstripchart(
        X, x='params', y='n', alpha=0.5
    ) +
    geom_text(data = X %>% count(params), 
              mapping = aes(x=params, y=max(X$n)+5, label=n)) +
    theme_pubr(x.text.angle=70) + 
    labs(x='Parameters', y='Module Size')
    
    ## what is the overlap between reference and query components?
    comparison = lapply(queries, function(query){
        cols_ref = metadata %>% filter(params %in% reference) %>% pull(id)
        cols_query = metadata %>% filter(params %in% query) %>% pull(id)
        sim = as.matrix(
            simil(modules[,cols_ref], modules[,cols_query], 
                  method = "Jaccard", by_rows=FALSE))
        
        # high similarity --> likely to be the "same" component
        # with respect to the reference
        max_similarities = apply(abs(sim), 1, max)
        df = data.frame(max_similarities=sort(max_similarities, decreasing=TRUE), 
                        reference=reference, 
                        query=query)
        # indicate only the possible maximum correlations
        max_detected = min(dim(sim))
        df$detected = TRUE
        if(max_detected < nrow(sim)){ df$detected[max_detected:nrow(df)] = FALSE }
        
        return(df)
    })
    comparison = do.call(rbind,comparison)
    
    plts[['components-best_jaccard']] = ggstripchart(
        comparison, x='query', y='max_similarities', 
        color='detected', palette='uchicago',
        alpha=0.5
    ) +
    theme_pubr(x.text.angle=70) + 
    labs(x='Parameters', y='Best Jaccard Similarity', color='Detected')
    
    return(plts)
}


make_plots = function(data, dataset){
    plts = list(
        plot_summary(data[['summary']]),
        plot_clustering_stats(data[['clustering_stats']]),
        plot_components_corr(data[['components']], data[['metadata']], dataset)
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
    lapply(names(plts), function(plt_name){
        save_plot(plts[[plt_name]] + 
                  theme_pubr(base_size=10, x.text.angle=70),
                  plt_name, '.png', figs_dir, width=12, height=12)
    })
}


main = function(){
    args = getParsedArgs()
            
    evaluation_dirs = unlist(strsplit(args$evaluation_dirs,' '))
    figs_dir = args$figs_dir
    dataset = args$dataset
    S_ref_file = args$S_ref_file
    
    dir.create(figs_dir, recursive = TRUE)
    
    data = load_data(evaluation_dirs, dataset, S_ref_file)
    plts = make_plots(data, dataset)
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}