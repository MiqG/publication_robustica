#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# visualize benchmark comparing multiple methods to cluster and compute
# robust ICA.
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
# RESULTS_DIR = file.path(ROOT,'results','compare_clustering_methods','files')
# dirnames = 'cluster_iterations-LGG-AgglomerativeClustering,cluster_iterations-LGG-CommonNNClustering,cluster_iterations-LGG-DBSCAN,cluster_iterations-LGG-KMeans,cluster_iterations-LGG-KMedoids'
# dirnames = unlist(strsplit(dirnames, ','))
# evaluation_dirs = file.path(RESULTS_DIR,dirnames)
# figs_dir = file.path(ROOT,'results','compare_clustering_methods','figures')


##### FUNCTIONS #####
get_label = function(dir_name){
    lab = gsub('cluster_iterations-LGG-','',basename(dir_name))
    return(lab)
}


load_datatype = function(evaluation_dirs, filename, bind=TRUE){
    df = lapply(evaluation_dirs, function(dir_name){
        lab = get_label(dir_name)
        df = read_tsv(file.path(dir_name, filename))
        if (!('time' %in% colnames(df))){ df$time=NA }
        df = df %>% mutate(params=lab)
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
        gc()
        return(list(components=components, metadata=metadata))
    })
    if(bind){ 
        components = lapply(dfs, function(df){df[['components']]})
        components = do.call(cbind, components)
        to_drop = grep('sample',colnames(components))[-1]
        components = components[,-to_drop]
        metadata = lapply(dfs, function(df){df[['metadata']]})
        metadata = do.call(rbind, metadata)            
    }
    return(list(components=components, metadata=metadata))
}


load_data = function(evaluation_dirs){
    # clustering stats
    clustering_stats = load_datatype(evaluation_dirs, filename='clustering_summary.tsv')
        
    data = list(
        clustering_stats = clustering_stats,
        components = S[['components']],
        metadata = S[['metadata']]
    )
    
    return(data)
}
                 

plot_clustering_stats = function(clustering_stats){
    X = clustering_stats %>% mutate(is_noise = label==-1)
    plts = list()
    plts[['clustering-n_clusters']] = ggstripchart(
        X %>% distinct(params, label) %>% count(params), 
        x='params', y='n'
    ) + 
    theme_pubr(x.text.angle=70)
    
    plts[['clustering-cluster_size']] = ggstripchart(
        X %>% count(params, label, is_noise), 
        x='params', y='n',
        color='is_noise', palette='jco', alpha=0.5
    ) + 
    yscale('log10') +
    theme_pubr(x.text.angle=70)
    
    plts[['clustering-silhouettes']] = ggviolin(
        X, x='params', y='silhouette', fill='params', color=NA
    ) + 
    theme_pubr(x.text.angle=70) +
    guides(fill=FALSE)
    
    plts[['clustering-silhouettes_mean_bycomponent']] = ggstripchart(
        X %>% group_by(params,label,is_noise) %>% summarize(sil=mean(silhouette)),
        x='params', y='sil',
        color='is_noise', palette='jco', alpha=0.5
    ) + 
    theme_pubr(x.text.angle=70)
    
    plts[['clustering-time']] = ggstripchart(
        X %>% distinct(params,time),
        x='params', y='time'
    ) + 
    theme_pubr(x.text.angle=70)
    
    return(plts)
}


make_plots = function(clustering_stats){
    plts = list(
        plot_clustering_stats(clustering_stats)
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
    
    dir.create(figs_dir, recursive = TRUE)
    
    data = load_data(evaluation_dirs)    
    plts = make_plots(data[['clustering_stats']])
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
