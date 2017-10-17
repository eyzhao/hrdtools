#' SNV Importer
#' This function imports SNVs from a tab-delimited file. The file must be formatted
#' in at least 4 tab-separated columns. The first four columns of the file must be
#' chromosome number (without 'chr'), position, reference base, mutant base.
#' @param path Path to the tab-delimited SNV data file
#' @importFrom IRanges IRanges
#' @importFrom VariantAnnotation VRanges 
#' @export

import_snvs <- function(path, genomeVersion='hg19', hasHeader=FALSE) {
    table = read.table(path, header=hasHeader, stringsAsFactors=FALSE, sep='\t')
    table = table[table[, 1] != 'MT', ]
    
    chrNum = as.character(table[,1])

    chr = paste0('chr', chrNum)
    ranges = IRanges(as.numeric(table[, 2]), as.numeric(table[, 2]))
    ref = as.character(table[, 3])
    alt = as.character(table[, 4])
    sampleNames = rep('sample', length(chr))
    seqinfo = 'Genomes for HRD SNV signature analysis'

    vr <- VRanges(
        seqnames = as.character(chr),
        ranges = ranges,
        ref = ref,
        alt = alt,
        sampleNames = sampleNames,
        study = rep('none', length(chr))
    )

    return(vr)
}


