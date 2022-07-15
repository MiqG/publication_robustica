#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# EDA of gene dependencies regressed on event PSI and gene TPMs.
#

require(tidyverse)
require(fdrtool)
require(clusterProfiler)
require(optparse)

# variables
THRESH_FDR = 0.05
THRESH_MEDIAN_DIFF = 5

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_SIGN_DIR = file.path(ROOT,'results','benchmark_sign_inference')
# msigdb_dir = file.path(RAW_DIR,'MSigDB','msigdb_v7.4','msigdb_v7.4_files_to_download_locally','msigdb_v7.4_GMTs')
# S_file = file.path(RESULTS_SIGN_DIR,'files',"LGG",'icasso-S.tsv.gz')
# features_col = "sample"

##### FUNCTIONS #####
define_module = function(x, cutoff=0.01){
    fdr = fdrtool(x, plot=FALSE, cutoff.method="fndr", verbose=FALSE)[['qval']]
    x = fdr < cutoff
    return(x)
}


get_genes_lists = function(S, features_col){
    components = setdiff(colnames(S), features_col)
    df = S %>% mutate_at(components, define_module)
    
    genes_lists = df %>% 
        pivot_longer(all_of(components), names_to='component', values_to='in_module') %>%
        filter(in_module) %>%
        drop_na() %>%
        dplyr::select(all_of(c(features_col, "component"))) %>%
        distinct() %>%
        with(., split(get(features_col), component))
        
    return(genes_lists)
}


get_universe = function(S, features_col){
    universe = S[[features_col]] %>% unique()
    return(universe)
}


run_enrichment = function(genes, universe, ontologies){
    enrichments = list()
    if(length(genes)>0){
        enrichments[['hallmarks']] = enricher(genes, TERM2GENE=ontologies[['hallmarks']], universe=universe)
        enrichments[['oncogenic_signatures']] = enricher(genes, TERM2GENE=ontologies[['oncogenic_signatures']], universe=universe)
        enrichments[['GO_BP']] = enricher(genes, TERM2GENE=ontologies[['GO_BP']], universe=universe)
    }
    
    return(enrichments)
}


run_enrichments = function(genes_lists, universe, ontologies){
    enrichments = sapply(names(genes_lists), function(cond){
            genes = genes_lists[[cond]]
            run_enrichment(genes, universe, ontologies)
        }, simplify=FALSE)
    return(enrichments)
}


get_enrichment_result = function(enrich_list, thresh){
    ## datasets are extracted from names
    if (length(enrich_list)>0){
        ontos = names(enrich_list[[1]])
        datasets = names(enrich_list)
        results = sapply(
            ontos, function(onto){
            result = lapply(datasets, function(dataset){
                result = enrich_list[[dataset]][[onto]]
                if(!is.null(result)){
                    result = result@result
                    result[['dataset']] = dataset
                    result[['ontology']] = onto
                }
                return(result)
            })
            result[sapply(result, is.null)] = NULL
            result = do.call(rbind,result)
            ## filter by p.adjusted
            result = result %>% filter(p.adjust<thresh)
        }, simplify=FALSE)
        results = do.call(rbind,results)
    }else{
        results = data.frame(
            ID=NA, Description=NA, GeneRatio=NA, BgRatio=NA, pvalue=NA,
            p.adjust=NA, qvalue=NA, geneID=NA, Count=NA, dataset=NA, ontology=NA
        ) %>% drop_na()
    }
    
    return(results)
}


parseargs = function(){
    
    option_list = list( 
        make_option("--S_file", type="character"),
        make_option("--msigdb_dir", type="character"),
        make_option("--output_file", type="character"),
        make_option("--features_col", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    S_file = args[["S_file"]]
    msigdb_dir = args[["msigdb_dir"]]
    output_file = args[["output_file"]]
    features_col = args[["features_col"]]
    
    # load
    S = read_tsv(S_file)
    ontologies = list(
        "hallmarks" = read.gmt(file.path(msigdb_dir,'h.all.v7.4.symbols.gmt')),
        "oncogenic_signatures" = read.gmt(file.path(msigdb_dir,'c6.all.v7.4.symbols.gmt')),
        "GO_BP" = read.gmt(file.path(msigdb_dir,'c5.go.bp.v7.4.symbols.gmt'))
    )
    
    # run enrichments
    genes_lists = get_genes_lists(S, features_col)
    universe = get_universe(S, features_col)
    enrichments = run_enrichments(genes_lists, universe, ontologies)
    results_enrich = get_enrichment_result(enrichments, THRESH_FDR)
    
    # save
    write_tsv(results_enrich, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
