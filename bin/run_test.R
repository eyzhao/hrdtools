#!/usr/bin/env Rscript

print('loading packages...')
library(optparse)
print('done')

if (args[['debug']]) {
    library(devtools)
    load_all('scripts/hrdtools')
    n <- readLines('scripts/hrdtools/NAMESPACE')
    for (p in unique(gsub('.*?\\((.*?)[,)].*', '\\1', n[grepl('^import', n)]))) library(p, character.only = TRUE)
} else {
    suppressMessages(library(hrdtools))
}

prepare_optparser <- function(){
#Prepare optparser object. New options will be added in this function first.
    option_list <- list(
        make_option(c("-l","--loh"), type="character", help="loh segs file from APOLLOH", dest="loh_file"),
        make_option(c("-g","--genome"), type="character", help="Genome version (default: hg19)", dest="genome", default='hg19')
        )

    optparser <- OptionParser(option_list=option_list,
                              usage="usage: %prog -i <input_file> -s <pi0_file> -c <normal_contamination>",
                              description="Imports necessary data for kaltimer")

    return(optparser)
}

verify_option <- function(options){
# check if mandatory options are present in the options object
    argpass <- TRUE

    if(!file.exists(options$loh_file)){
        argpass <- FALSE
        print("Please provide a valid loh file")
    }

    return(argpass)
}

get_args <- function() {
# Run the argument parser and return arguments
    optparser <- prepare_optparser()
    opt <- parse_args(optparser)
    argpass <- verify_option(opt)
    
    if(!argpass){
       print("Non valid arguments.")
       quit("no", status=1)
    }

    return(opt)
}

if(getOption('run_test.from.rscript', default=TRUE)) {
    args = get_args()
    run_test(args$loh_file, genome=args$genome)
}
