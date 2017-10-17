#' Default sanitization function for GRange objects
#' 
#' This function returns a sanitized GRanges object. It can be passed as an argument to import_ranges.
#' @import GenomeInfoDb
#' @param gr GRanges object which is to be sanitized

sanitize_granges <- function(gr) {
    seqlevels(gr) <- gsub('hs', '', seqlevels(gr))
    seqlevels(gr) <- gsub('chr', '', seqlevels(gr))
    seqlevels(gr) <- gsub('23', 'X', seqlevels(gr))
    seqlevels(gr) <- gsub('24', 'Y', seqlevels(gr))

    return(gr)
}

#' Imports ranges output from APOLLOH
#' 
#' This function returns a GRanges object based on the data contained in a TSV.
#' This TSV can be output from APOLLOH or similar program.
#' @param path Path to LOH ranges file
#' @param col.names A vector providing the column names of the LOH ranges file
#' @param sanitization_function A callback function which takes a GRanges object as input and returns a sanitized version
#' @import GenomicRanges
#' @import IRanges
#' @export

import_ranges <- function(path, col.names, sanitization_function=sanitize_granges) {
    ranges_table <- read.table(path, header=FALSE, stringsAsFactors=FALSE, sep='\t')
    colnames(ranges_table) <- col.names
    gr <- sanitization_function(GRanges(ranges_table))

    return(gr)
}
