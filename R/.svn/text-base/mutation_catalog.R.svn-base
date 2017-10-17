#' Obtain the Mutational Catalog / Spectrum
#'
#' This function returns the mutational catalogue / mutational spectrum of a cancer genome.
#' @param vr VRanges file, obtained from import_snvs()
#' @importFrom SomaticSignatures mutationContext motifMatrix
#' @importFrom BSgenome getBSgenome
#' @export

mutation_catalog <- function(vr) {
    genome = getBSgenome('BSgenome.Hsapiens.UCSC.hg19')
    context <- mutationContext(vr, genome)
    catalog <- motifMatrix(context, group='sampleNames', normalize=FALSE)
    return(catalog)
}
