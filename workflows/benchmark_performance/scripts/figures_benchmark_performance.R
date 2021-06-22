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


require(tidyverse)
require(ggpubr)
require(pheatmap)


ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables

# Development
# -----------
#RESULTS_DIR = file.path(ROOT,'results','benchmark_performance','files')
#performance_evaluation_file = file.path(RESULTS_DIR,'performance_evaluation.tsv.gz')
#S_icasso_file = file.path(RESULTS_DIR,'icasso-S.tsv.gz')
#S_robustica_file = file.path(RESULTS_DIR,'robustica-S.tsv.gz')
#S_robustica_pca_file = file.path(RESULTS_DIR,'robustica_pca-S.tsv.gz')


figs_dir = file.path(ROOT,'results','benchmark_performance','figures')

##### FUNCTIONS #####
get_relative_timestamp = function(x){
    x0=x[1]
    return(x-x0)
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
    guides(color=guide_legend(ncol=2)) + 
    labs(x='Time (s)', y='Memory (MiB)', color='Function') + 
    theme_pubr(border = TRUE) 
    
    return(plts)
}


plot_S_correlations = function(Ss){
    plts = list()
    
    # are all components found independent?
    for (algorithm in names(Ss)){
        corr = cor(Ss[[algorithm]])
        plt_name = paste0('self_corr-S-',algorithm)
        plts[[plt_name]] = pheatmap(
            corr, show_colnames = FALSE, show_rownames = FALSE, 
            silent=TRUE, main=plt_name
        )
    }
    
    # Do we get the same output?
    corr = cor(Ss[['icasso']],Ss[['robustica_pca']])
    plts[['corr-S-icasso_vs_robustica_pca']] = pheatmap(
        corr, show_colnames = FALSE, show_rownames = FALSE, 
        main='icasso vs robustica_pca', silent=TRUE
    )
    
    corr = cor(Ss[['robustica']],Ss[['robustica_pca']])
    plts[['corr-S-robustica_vs_robustica_pca']] = pheatmap(
        corr, show_colnames = FALSE, show_rownames = FALSE, 
        main='robustica vs robutica_pca', silent=TRUE
    )
    
    return(plts)
}


make_plots = function(df, Ss){
    plts = list(
        plot_performance_profile(df),
        plot_S_correlations(Ss)
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
    save_plot(plts[['mem_time-scatter']],'mem_time-scatter','.png',figs_dir, width=12,height=15)
    
    lapply(grep('corr',names(plts),value=TRUE), function(plt_name){
        save_plot(plts[[plt_name]],
                  plt_name, '.pdf', figs_dir, width=12, height=12)
    })
}


main = function(){
    args = getParsedArgs()

    performance_evaluation_file = args$performance_evaluation_file
    S_icasso_file = args$S_icasso_file
    S_robustica_file = args$S_robustica_file
    S_robustica_pca_file = args$S_robustica_pca_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load data
    performance_evaluation = read_tsv(performance_evaluation_file) %>% 
        group_by(algorithm) %>% 
        mutate(rel_timestamp = get_relative_timestamp(timestamp))
    
    Ss = list(
        'icasso'=read_tsv(S_icasso_file),
        'robustica'=read_tsv(S_robustica_file),
        'robustica_pca'=read_tsv(S_robustica_pca_file)
    )
    
    plts = make_plots(performance_evaluation, Ss)
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}