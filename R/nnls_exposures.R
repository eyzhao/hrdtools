#' Calculate signature exposures
#'
#' This function determines the signature exposures using non-negative least squares.
#' @param subjectMotifs The mutation catalog of a subject derived from mutation_catalog()
#' @param refSignature The reference signature table, which can be imported by get_reference_signatures()
#' @param fractions If TRUE, will return results as a fraction of all mutations rather than number of mutations.
#' @param montecarlo If TRUE, will perform Monte Carlo simulation to obtain 95 percent confidence limits
#' @param iterations The number of iterations of Monte Carlo to run
#' @export

nnls_exposures <- function(subjectMotifs, refSignature, fractions=FALSE, montecarlo=FALSE, iterations=1000) {
    if (montecarlo) {
        return(nnls_montecarlo(subjectMotifs, refSignature, fractions, iterations))
    } else {
        return(nnls_exposures_worker(subjectMotifs, refSignature, fractions))
    }
}


#' Monte Carlo Simulation of NNLS
#'
#' Performs Monte Carlo simulation to obtain NNLS exposure confidence intervals
#' @param subjectMotifs Mutation Catalog from get_mutation_catalog()
#' @param refSignature Reference signatures from get_reference_signatures()
#' @param fractions If TRUE, will return exposures as a fraction of total mutation burden
#' @param iterations Number of iterations to run
#' @param alpha Confidence level to retrieve confidence intervals at (default 0.05)
#' @export

nnls_montecarlo <- function(subjectMotifs, refSignature, fractions, iterations, alpha=0.05) {
    naiveSolution <- nnls_exposures_worker(subjectMotifs, refSignature, fractions)

    types = rownames(subjectMotifs)
    n = sum(subjectMotifs)
    p = subjectMotifs / n
    
    nnlsMC <- sapply(1:iterations, function(iterationNumber) {
        data <- sample(types, n, replace=TRUE, prob=p)
        vec <- table(data)
        if (! identical(names(vec), rownames(subjectMotifs))) {
            if (all(names(vec) %in% rownames(subjectMotifs))) {
                # Monte carlo motifs are a subset of subjectmotifs.
                # This happens when some subject motifs have frequency of 0.
                motifNames = rownames(subjectMotifs)
                newVec <- rep(0, length(motifNames))
                temp <- sapply(motifNames, function(z) {newVec[z] = vec[z]})
                vec <- temp
                vec[is.na(vec)] = 0
                names(vec) = motifNames
            } else {
                print('Error: Monte Carlo has different motifs')
                print('Monte Carlo motifs names:')
                print(names(vec))
                print('Subject Motif Names:')
                print(rownames(subjectMotifs))
                stop('Signature names are different from Monte Carlo signature names')
            }
        }

        nnlsSolution <- nnls_exposures_worker(vec, refSignature, fractions)
        return(nnlsSolution$x)
    })

    nnlsMC <- t(nnlsMC)

    exposureMeans <- apply(nnlsMC, 2, mean)
    exposureCI <- t(apply(nnlsMC, 2, function(col) {quantile(col, c(alpha/2, 1-(alpha/2)))}))

    signatureNames = names(refSignature[3:dim(refSignature)[2]])

    return(list(x=naiveSolution$x, means=exposureMeans,
        lCI=exposureCI[, 1], uCI=exposureCI[, 2], names=signatureNames))
}

#' NNLS Worker Function
#'
#' Performs the NNLS calculation itself.
#' @param subjectMotifs Mutation Catalog from get_mutation_catalog()
#' @param refSignature Reference signatures from get_reference_signatures()
#' @param fractions If TRUE, will return exposures as a fraction of total mutation burden
#' @import nnls
#' @export

nnls_exposures_worker <- function(subjectMotifs, refSignature, fractions) {
    refSignature = refSignature[order(paste(refSignature$alteration, refSignature$context)), ]

    refSignatureMotifNames <- paste(refSignature$alteration, refSignature$context)
    if ((!identical(rownames(subjectMotifs), refSignatureMotifNames)) && !identical(names(subjectMotifs), refSignatureMotifNames)) {
        print('Subject Motif Names:')
        print(rownames(subjectMotifs))
        print(names(subjectMotifs))
        print('Reference Motif Names:')
        print(refSignatureMotifNames)
        stop('Signatures are different.')
    }

    A = as.matrix(refSignature[3:dim(refSignature)[2]])
    subjectMotifs = as.vector(subjectMotifs)
    nnlsEstimateObj = nnls(A, subjectMotifs)
    nnlsEstimate = nnlsEstimateObj$x
    estimationError = A %*% nnlsEstimate - subjectMotifs

    if (fractions) {
        nnlsEstimate = nnlsEstimate / sum(nnlsEstimate)
    } 

    signatureNames = names(refSignature[3:dim(refSignature)[2]])

    return(list(x=nnlsEstimate, names=signatureNames))
}
