#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
# Last update: 2022-07-13
#
# Script purpose
# --------------
# EDA of datasets used in this project.
#

require(tidyverse)
require(ggpubr)
require(cowplot)
require(writexl)
require(extrafont)
require(optparse)

# formatting
FONT_FAMILY = "Arial"

# Development
# -----------
# ROOT = here::here()
# RESULTS_DIR = file.path(ROOT,'results','preprocess_data','files')
# dataset_info_file = file.path(RESULTS_DIR,'dataset_info.tsv')
# figs_dir = file.path(ROOT,'results','preprocess_data','figures','eda_datasets')


##### FUNCTIONS #####
plot_datasets_eda = function(dataset_info){
    X = dataset_info
    
    plts = list()
    
    plts[["datasets_eda-counts_samples"]] = X %>%
        ggbarplot(x="label", y="n_samples", fill="darkred", color=NA) +
        labs(x="Dataset", y="N. Samples") +
        coord_flip()
    
    plts[["datasets_eda-counts_genes"]] = X %>%
        ggbarplot(x="label", y="n_genes", fill="darkred", color=NA) +
        labs(x="Dataset", y="N. Genes") +
        coord_flip()
    
    return(plts)
}


make_plots = function(dataset_info){
    plts = list(
        plot_datasets_eda(dataset_info)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(dataset_info){
    figdata = list(
        'eda_datasets' = list(
            'dataset_info' = dataset_info
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
    save_plt(plts, 'datasets_eda-counts_samples','.pdf',figs_dir, width=10, height=10)
    save_plt(plts, 'datasets_eda-counts_genes','.pdf',figs_dir, width=10, height=10)
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


parseargs = function(){
    
    option_list = list( 
        make_option("--dataset_info_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    dataset_info_file = args[["dataset_info_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load data
    dataset_info = read_tsv(dataset_info_file) %>%
        mutate(label = paste0(main, " | ", dataset),
               label = gsub("Sastry2019 \\|","",label))
    
    plts = make_plots(dataset_info)
    figdata = make_figdata(dataset_info)
    
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}