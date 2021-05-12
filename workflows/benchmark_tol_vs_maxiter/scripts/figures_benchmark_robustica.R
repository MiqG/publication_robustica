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
require(pheatmap)
require(proxy)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# Development
# -----------
TOLS = paste0('1e-',c(1,2,3,4))
MAX_ITERS = as.character(c(100,200,500,1000,10000))
RESULTS_DIR = file.path(ROOT,'results','benchmark_tol_vs_maxiter','files')
dirnames = as.vector(sapply(TOLS, function(t){sapply(MAX_ITERS, function(m){ sprintf('evaluate_iterations-LGG_%s_%s',m,t) })}))

evaluation_dirs = file.path(RESULTS_DIR,dirnames)
figs_dir = file.path(ROOT,'results','benchmark_tol_vs_maxiter','figures')


##### FUNCTIONS #####
get_label = function(dir_name){
    lab = gsub('evaluate_iterations-LGG_','',basename(dir_name))
    lab = unlist(strsplit(lab,'_'))
    lab = apply(cbind(c('max_iter=','tol='),lab),1, paste, collapse='')
    lab = paste(lab, collapse=' & ')
    return(lab)
}


load_datatype = function(evaluation_dirs, filename, bind=TRUE){
    df = lapply(evaluation_dirs, function(dir_name){
        df = read_tsv(file.path(dir_name, filename))
        return(df)
    })
    if(bind){ df = do.call(rbind, df) }
    return(df)
}


load_components = function(evaluation_dirs, filename, bind=TRUE){
    dfs = lapply(evaluation_dirs, function(dir_name){
        lab = get_label(dir_name)
        components = read_tsv(file.path(dir_name, filename)) %>% 
            rename_at(vars(-sample),function(x) paste0(x,'_',lab))
        metadata = data.frame(id=colnames(components), params=lab) %>% filter(id!='sample')
        return(list(components=components, metadata=metadata))
    })
    if(bind){ 
        components = lapply(dfs, function(df) {df[['components']]})
        components = components %>% reduce(inner_join, by='sample')
        metadata = lapply(dfs, function(df) {df[['metadata']]})
        metadata = do.call(rbind, metadata)            
    }
    return(list(components=components, metadata=metadata))
}


load_data = function(evaluation_dirs){
    # summary
    summary = load_datatype(evaluation_dirs, filename='summary.tsv.gz')
    summary = summary %>% mutate(has_converged = convergence_score < tol)
    
    # clustering stats
    clustering_stats = load_datatype(evaluation_dirs, filename='clustering_stats.tsv.gz')
    
    # S
    S = load_components(evaluation_dirs, filename='S.tsv.gz')
    
    data = list(
        summary = summary,
        clustering_stats = clustering_stats,
        components = S[['components']],
        metadata = S[['metadata']]
    )
    
    return(data)
}


plot_summary = function(summary){
    X = summary
    plts = list()
    plts[['summary-has_converged']] = ggbarplot(
        X %>% count(params,has_converged), x='params', y='n', 
        fill='has_converged', label=TRUE, 
        color=NA, palette='npg', position=position_dodge(0.9)
    ) + 
    theme_pubr(x.text.angle=70)
    
    plts[['summary-iterations']] = ggstripchart(
        X, x='params', y='iteration_ica', 
        color='has_converged', palette='npg', alpha=0.5
    ) + 
    theme_pubr(x.text.angle=70)
    
    plts[['summary-convergence_score']] = ggstripchart(
        X, x='params', y='convergence_score', 
        color='has_converged', palette='npg', alpha=0.5
    ) + 
    yscale('log10') + 
    theme_pubr(x.text.angle=70)
    
    plts[['summary-convergence_time']] = ggstripchart(
        X, x='params', y='time_ica', 
        color='has_converged', palette='npg', alpha=0.5
    ) + 
    theme_pubr(x.text.angle=70)
    
    plts[['summary-robustica_time']] = ggstripchart(
        X %>% count(time_robustica,params,has_converged) %>% filter(has_converged), 
        x='params', y='time_robustica', 
        size='n'
    ) + 
    theme_pubr(x.text.angle=70)
    
    return(plts)
}


plot_clustering_stats = function(clustering_stats){
    X = clustering_stats %>% mutate(is_noise=cluster_id==-1)
    plts = list()
    plts[['clustering-n_clusters']] = ggstripchart(
        X %>% count(params), x='params', y='n'
    ) + 
    theme_pubr(x.text.angle=70)
    
    plts[['clustering-cluster_size']] = ggstripchart(
        X, x='params', y='cluster_size',
        color='is_noise', palette='jco', alpha=0.5
    ) + 
    yscale('log10') +
    theme_pubr(x.text.angle=70)
    
    plts[['clustering-cluster_mean_std']] = ggstripchart(
        X, x='params', y='S_mean_std',
        color='is_noise', palette='jco', alpha=0.5
    ) + 
    yscale('log10') +
    theme_pubr(x.text.angle=70)
    
    return(plts)
}


define_module = function(x){
    # outliers
    mu = mean(x)
    std = sd(x)
    x = abs(x) > (mu+std*3)
    return(x)
}


plot_components_corr = function(components, metadata){
    
    # init
    reference = 'max_iter=10000 & tol=1e-4'
    queries = setdiff(unique(metadata$params),reference)
    plts = list()
    
    # correlation between robust components of best quality params and rest
    comparison = lapply(queries, function(query){
        cols_ref = metadata %>% filter(params %in% reference) %>% pull(id)
        cols_query = metadata %>% filter(params %in% query) %>% pull(id)
        sim = cor(components[,cols_ref], components[,cols_query], 
                  method = "pearson")
        
        # high similarity --> likely to be the "same" component
        # with respect to the reference
        max_correlations = apply(abs(sim), 1, max)
        df = data.frame(max_correlations=max_correlations, 
                        reference=reference, 
                        query=query)
        return(df)
    })
    comparison = do.call(rbind,comparison)
    
    plts[['components-best_correlation']] = ggstripchart(
        comparison, x='query', y='max_correlations',
        alpha=0.5
    ) +
    theme_pubr(x.text.angle=70)
    
    # module overlap
    ## define modules
    modules = components %>% 
        dplyr::select(-one_of('sample')) %>% 
        mutate_all(define_module)
    
    ## selected items (genes) per module per params
    X = enframe(apply(modules, 2, sum), 'id', 'n') %>% 
        left_join(metadata, by='id')
    plts[['components-module_size']] = ggstripchart(
        X, x='params', y='n', alpha=0.5
    ) +
    geom_text(data = X %>% count(params), 
              mapping = aes(x=params, y=650, label=n)) +
    theme_pubr(x.text.angle=70)
    
    ## what is the overlap between reference and query components?
    comparison = lapply(queries, function(query){
        cols_ref = metadata %>% filter(params %in% reference) %>% pull(id)
        cols_query = metadata %>% filter(params %in% query) %>% pull(id)
        sim = as.matrix(
            simil(modules[,cols_ref], modules[,cols_query], 
                  method = "Jaccard", by_rows=FALSE))
        
        # high similarity --> likely to be the "same" component
        # with respect to the reference
        max_similarities = apply(sim, 1, max)
        df = data.frame(max_similarities=max_similarities, 
                        reference=reference, 
                        query=query)
    })
    comparison = do.call(rbind,comparison)
    
    plts[['components-best_jaccard']] = ggstripchart(
        comparison, x='query', y='max_similarities',
        alpha=0.5
    ) +
    theme_pubr(x.text.angle=70)
    
    return(plts)
}


make_plots = function(summary, clustering_stats, components, metadata){
    plts = list(
        plot_summary(summary),
        plot_clustering_stats(clustering_stats),
        plot_components_corr(components, metadata)
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


save_plots = function(plts, figs_dir, params){
    lapply(names(plts), function(plt_name){
        save_plot(plts[[plt_name]] + 
                  theme_pubr(base_size=10, x.text.angle=70),
                  plt_name, '.png', figs_dir, width=12, height=12)
    })
}


main = function(){
    args = getParsedArgs()
            
    evaluation_dirs = unlist(strsplit(args$evaluation_dirs,','))
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    data = load_data(evaluation_dirs)    
    plts = make_plots(data[['summary']], data[['clustering_stats']], 
                      data[['components']], data[['metadata']])
    save_plots(plts, figs_dir, params)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}