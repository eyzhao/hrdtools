#' Function that imports VCFs into VariantAnnotation format
#' 
#' @param vcfPath Path to VCF file
#' @param genomeVersion String representing the genome version (default: 'hg19')
#' @return A VRanges object with variants and allelic depths
#' @importFrom VariantAnnotation readVcf makeVRangesFromGRanges geno
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom IRanges CharacterList
#' @importFrom S4Vectors elementNROWS
#' @export

import_vcf <- function(vcfPath, genomeVersion='hg19') {
    snvRanges <- readVcf(vcfPath, 'hg19')

    vrObj = rowRanges(snvRanges)

    # Remove rows where ALT has more than one character
    rowsToKeep = elementNROWS(vrObj$ALT) == 1
    vrObj <- vrObj[rowsToKeep, ]

    vrObj$REF = as.character(vrObj$REF)
    vrObj$ALT = unlist(CharacterList(vrObj$ALT))
    vrObj$totalDepth = geno(snvRanges)$DP[rowsToKeep, 2]

    depths = data.frame(A=geno(snvRanges)$AU[rowsToKeep, 2, 1],
                        C=geno(snvRanges)$CU[rowsToKeep, 2, 1],
                        G=geno(snvRanges)$GU[rowsToKeep, 2, 1],
                        T=geno(snvRanges)$TU[rowsToKeep, 2, 1]
                        )

    vrObj$refDepth = sapply(1:length(vrObj$REF), function(i) {depths[i, vrObj$REF[i]]} )
    vrObj$altDepth = sapply(1:length(vrObj$ALT), function(i) {depths[i, vrObj$ALT[i]]} )
    vrObj$sampleNames = rep('sample', length(vrObj$REF))

    return(makeVRangesFromGRanges(vrObj))
}

#' Modification of seqlevels to match those of the genome
#' 
#' @param vr VRanges object
#' @return VRanges object with seqlevels changed to match those of the genome
#' @importFrom VariantAnnotation readVcf makeVRangesFromGRanges geno
#' @importFrom BSgenome getBSgenome
#' @importFrom GenomeInfoDb seqlevels seqlevels<-

sanitize_vranges <- function(vr) {
    genome = getBSgenome('BSgenome.Hsapiens.UCSC.hg19')
    seqlevels(vr, force=TRUE) <- seqlevels(vr)[1:24]
    if (identical(seqlevels(vr), seqlevels(genome)[1:24])) {
    } else if (identical(paste0('chr', seqlevels(vr)), seqlevels(genome)[1:24])) {
        seqlevels(vr) <- paste0('chr', seqlevels(vr))
    } else {
        stop('Cannot match VCF and reference genome seqlevels')
    }

    return(vr)
}
