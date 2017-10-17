#' The HRD-TAI Test
#'
#' This function runs the Telomeric Allelic Imbalance test (HRD-TAI).
#' @param gr GRanges object obtained from import_ranges()
#' @export

tai_test <- function(gr) {
    return(sum(tai_test_bool(gr)))
}

#' The HRD-TAI Test Boolean Calculator
#'
#' Returns logical vector indicating whether each region has TAI
#' @param gr GRanges object obtained from import_ranges()
#' @importFrom GenomicRanges GRanges ranges seqnames findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom plyr ddply
#' @export

tai_test_bool <- function(gr) {
    subtelomeres <- GRanges(ddply(get_subtelomere_regions(), 'chromosome', function(z) {
        return(data.frame(start = c(0, z$end), end = c(z$start, .Machine$integer.max)))
    }))

    centromeres <- GRanges(get_centromere_regions())

    overlaps_with_subtelomere <- rep(FALSE, length(gr))
    overlaps_with_subtelomere[queryHits(findOverlaps(gr, subtelomeres))] <- TRUE

    overlaps_with_centromere <- rep(FALSE, length(gr))
    overlaps_with_centromere[queryHits(findOverlaps(gr, centromeres))] <- TRUE

    is_balanced <- gr$lohtype %in% c('HET', 'BCNA')

    return(overlaps_with_subtelomere & ! overlaps_with_centromere & ! is_balanced)
}

#' Get Subtelomere data
#' 
#' This function returns subtelomere regions as a table.
#' @export

get_subtelomere_regions <- function(chr_label = F) {
    data('subtelomeres')
    if (! chr_label) {
        subtelomeres$chromosome <- gsub('chr', '', subtelomeres$chromosome)
    }
    return(subtelomeres)
}

#' Get Centromere regions
#'
#' This function returns centromere regions as a table
#' @export

get_centromere_regions <- function() {
    data(centromeres)
    return(centromeres)
}

#' Does it Cross the Centromere?
#'
#' This function checks whether a region crosses the centromere
#' @param chr The chromosome that the region belongs to
#' @param start The start position of the region
#' @param end The end position of the region
#' @param centromeres Centromere positions from get_centromere_regions()
#' @export

crosses_centromere <- function(chr, start, end, centromeres) {
    rowIndex <- which(centromeres[, 2] == chr)
    centroStart = centromeres[rowIndex, 3]
    centroEnd = centromeres[rowIndex, 4]

    if (start > centroEnd && end > centroEnd) {
        return(FALSE)
    } else if (start < centroStart && end < centroStart) {
        return(FALSE)
    } else {
        return(TRUE)
    }
}

#' Get Subtelomere Coordinate
#'
#' Convenience function that retrieves the start or end coordinate of a specific subtelomere
#' @param chr The chromosome of interest
#' @param subtelomere A subtelomeres object from get_subtelomere_regions()
#' @param startOrEnd Takes on value either "start" or "end" depending on which position to return
#' @export

get_subtelomere_coordinate <- function(chr, subtelomere, startOrEnd) {
    chr = paste0('chr', chr)
    rowIndex <- which(subtelomere[, 1] == chr)

    if (startOrEnd == 'start') {
        colIndex = 2
    } else if (startOrEnd == 'end') {
        colIndex = 3
    } else {
        stop("get_subtelomere_coordinate: startOrEnd must be either 'start' or 'end'")
    }
        
    coordinate = subtelomere[rowIndex, colIndex]

    return(coordinate)
}


