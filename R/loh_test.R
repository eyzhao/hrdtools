#' The HRD-LOH Test
#'
#' This function runs the Telomeric Allelic Imbalance test (HRD-LOH).
#' @param gr GRanges object obtained from import_ranges()
#' @importFrom GenomicRanges ranges
#' @export

loh_test <- function(gr) {
    bool <- loh_test_bool(gr)

    return(sum(bool))
}


#' LOH test boolean calculator
#'
#' Determines boolean vector indicating whether each range is a large LOH event
#' @param gr GRanges object obtained from import_ranges()
#' @export

loh_test_bool <- function(gr) {
    SIZE_THRESHOLD = 15 * 1000 * 1000
    df <- as.data.frame(gr)
    is_long <- df$width > SIZE_THRESHOLD
    is_loh <- df$lohtype %in% c('NLOH', 'DLOH', 'ALOH')

    return(is_long & is_loh)
}
