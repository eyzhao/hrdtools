#' The HRD-LST Test
#' 
#' This function carries out the large scale transition (HRD-LST) test on a GRanges object.
#' @param gr GRanges object, obtained from import_ranges().
#' @export

lst_test <- function(gr, loh.col = 'lohtype') {
    return(sum(lst_test_bool(gr, loh.col)))
}

#' The HRD-LST Test Boolean Calculator
#' 
#' This function returns a boolean vector indicating whether each region has an LST with the next adjacent region.
#' @param gr GRanges object, obtained from import_ranges().
#' @param loh.col Name of the column containing LOH information
#' @export

lst_test_bool <- function(gr, loh.col) {
    df <- as.data.frame(gr)
    SIZE_THRESHOLD = 10 * 1000 * 1000
    is_long <- df$width > SIZE_THRESHOLD 
    next_is_long <- c(is_long[2:length(is_long)], FALSE)
    same_lohtype_as_next <- c(df[1:(dim(df)[1] - 1), loh.col] == df[2:dim(df)[1], loh.col], FALSE)
    same_chr_as_next <- as.logical(c(df[1:(dim(df)[1] - 1), 'seqnames'] == df[2:dim(df)[1], 'seqnames'], TRUE))
    
    return(is_long & next_is_long & !same_lohtype_as_next & same_chr_as_next)
}

