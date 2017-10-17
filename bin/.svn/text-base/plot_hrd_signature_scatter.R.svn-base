#!/usr/bin/env Rscript

library('argparse')
library('RSvgDevice')
library('stringr')

plot_hrd_signature_scatter <- function(hrd_file, sig_file, scores_file, output_file, figure_path, cancerType) {
    scoresTable <- read.table(scores_file, header=TRUE)
    hrdScoreValue = as.numeric(scoresTable$hrdscore[1])
    hrdExposureValue = as.numeric(scoresTable$exposure[1])

    hrdTable <- read.table(hrd_file, header=FALSE, stringsAsFactors=FALSE, sep='\t')
    sigTable <- read.table(sig_file, header=TRUE, stringsAsFactors=FALSE, sep='\t')

    sigTable[, 1] <- str_replace(sigTable[, 1], '.sig.txt', '')
    hrdTable[, 1] <- str_replace(hrdTable[, 1], 'data/Breast/loh/', '')
    hrdTable[, 1] <- str_replace(hrdTable[, 1], 'data/Ovarian/loh/', '')
    hrdTable[, 1] <- str_replace(hrdTable[, 1], 'data/Sarcoma/loh/', '')
    hrdTable[, 1] <- str_replace(hrdTable[, 1], 'data/Melanoma/loh/', '')
    hrdTable[, 1] <- str_replace(hrdTable[, 1], '.loh.txt', '')

    colnames(hrdTable) <- c('sample', 'loh', 'tai', 'lst', 'total')

    merged <- merge(sigTable[, c(1,4)], hrdTable)
     
    x <- merged$total
    y <- merged$Signature.3

    fileConn <- file(output_file)
    writeLines(c(sprintf("HRD Score: %s - %s / %s = %s",
                   hrdScoreValue,
                   length(x[x <= hrdScoreValue]), 
                   length(x),
                   length(x[x <= hrdScoreValue]) / length(x)
                ), sprintf("HRD Exposure: %s - %s / %s = %s",
                   hrdExposureValue,
                   length(y[y <= hrdExposureValue]), 
                   length(y),
                   length(y[y <= hrdExposureValue]) / length(y)
                )), fileConn)

    close(fileConn)

    devSVG(figure_path, height=4, width=4)
    do.call(plot, list(x=x, y=y, ylim=c(0, 1), pch=20, xlab='HRD Score', ylab='Signature 3 Exposure Fraction', main=cancerType))
    points(c(hrdScoreValue), c(hrdExposureValue), pch=19, cex=2, col='red')
    dev.off()
}

prepare_optparser <- function() {
    parser <- ArgumentParser()
    parser$add_argument("-i", dest="hrd_file", help="HRD output table")
    parser$add_argument("-s", dest="sig_file", help="Mutation signature output")
    parser$add_argument("-d", dest="scores_file", help="The file with scores to highlight")
    parser$add_argument("-o", dest="output_file", help="Output Information File")
    parser$add_argument("-f", dest="figure_path", help="Output SVG Figure")
    parser$add_argument("-c", dest="cancer_type", help="Cancer Type")
    args <- parser$parse_args()
    return(args)
}

verify_option <- function(args){
# check if mandatory options are present in the options object
    argpass <- TRUE

    if(!file.exists(args$hrd_file)){
        argpass <- FALSE
        print(sprintf('File does not exist: $s', args$hrd_file))
    }

    if(!file.exists(args$sig_file)){
        argpass <- FALSE
        print(sprintf('File does not exist: $s', args$sig_file))
    }

    return(argpass)
}

get_args <- function() {
# Run the argument parser and return arguments
    args <- prepare_optparser()
    argpass <- verify_option(args)
    
    if(!argpass){
       print("Non valid arguments.")
       quit("no", status=1)
    }

    return(args)
}

if(getOption('plot.from.rscript', default=TRUE)) {
    args = get_args()
    plot_hrd_signature_scatter(args$hrd_file, args$sig_file, args$scores_file, args$output_file, args$figure_path, args$cancer_type)
}
