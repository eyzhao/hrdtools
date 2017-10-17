#' Retrieving Reference Mutation Signatures
#' 
#' This function imports reference mutation signatures to serve as a comparison point when
#' performing non-negative least squares decomposition.
#' @param path String indicating path to reference signatures table file. Defaults to bundled 30-signatures.
#' @export

get_reference_signatures <- function(path = NULL) {
    if (is.null(path)) {
        data(reference.signatures)
    } else {
        reference.signatures <- read.table(path, header=TRUE, sep='\t', stringsAsFactors=FALSE)
    }

    referenceAlteration = gsub('>', '', reference.signatures$Substitution.Type)
    referenceContext = gsub('(.).(.)', '\\1\\.\\2', reference.signatures$Trinucleotide)
    reference.signatures = reference.signatures[order(paste(referenceAlteration, referenceContext)), ]

    alteration = referenceAlteration[order(paste(referenceAlteration, referenceContext))]
    context = referenceContext[order(paste(referenceAlteration, referenceContext))]

    reference.signatures = cbind(alteration,
                           context,
                           reference.signatures[, 4:dim(reference.signatures)[2]]
                           )

    return(reference.signatures)
}
