#' Processes a GRanges object, producing ranges which fill in all gaps between adjacent entries
#'
#' @param gr GRanges object containing the ranges to process
#' @param loh_colname String representing the column name containing LOH classification
#' @param copynumber_colname String representing the column name containing copy number classification
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @export

contig <- function(gr, loh_colname, copynumber_colname) {
    col_names <- c(loh_colname, copynumber_colname)
    values(gr)[, 'lohtype'] <- as.character(values(gr)[, 'lohtype'])

    ## Contig the ranges
    gr <- gr[, col_names]; gr <- sortSeqlevels(gr); gr <- sort(gr)

    g <- gaps(gr)
    g <- attach_values(g, col_names,
        c(rep('TEMP_CONTIG', length(loh_colname)),
          rep(-1, length(copynumber_colname))
        )
    )

    combined <- sort(c(gr, g))
    temp_contig_indices <- which(values(combined)[, loh_colname] == 'TEMP_CONTIG')
    upstream_indices <- temp_contig_indices - 1; downstream_indices <- temp_contig_indices + 1;
    upstream_length <- sapply(upstream_indices, function(z) {
        if (z == 0) {return(0)} else {(return(width(combined)[z]))}
    })
    downstream_length <- width(combined)[downstream_indices]

    # The following code block gives the new gap-filling rows the values of the
    # adjacent segment which is longest. It effectively fills the gaps in the range
    # by extending the longer segment to fill them.

    g <- do.call('c', lapply(1:length(g), function(g_idx) { 
        indices <- c(upstream_indices[g_idx], downstream_indices[g_idx])
        lengths <- c(upstream_length[g_idx], downstream_length[g_idx])
        new_value_idx <- indices[which(lengths == max(lengths))[1]]
        new_values <- as.data.frame(values(combined)[new_value_idx, ])
        return(attach_values(g[g_idx, 0], col_names, new_values))
    }))

    gr <- stitch(sort(c(gr, g)), col_names)

    return(gr)
}

#' Filters a GRanges object to remove short events between identical longer events
#'
#' @param gr GRanges object containing the ranges to process
#' @param loh_col String representing the column name containing LOH classification
#' @param cnv_col String representing the column name containing copy number classification
#' @param threshold Numeric value indicating the ratio of longer to shorter adjacent regions which result in removal of short region
#' @export

filter_ranges <- function(gr, loh_col, cnv_col, threshold = 100) {
    col_names <- c(loh_col, cnv_col)
    print(length(gr))
    gr_new <- run_filter(gr, col_names, threshold)
    print(length(gr_new))

    while (!identical(length(gr_new), length(gr))) {
        gr <- gr_new
        gr_new <- run_filter(gr, col_names, threshold)
        print(length(gr_new))
    }
    return(gr_new)
}

#' Worker function which runs a single iteration of the filtering algorithm
#'
#' @param gr GRanges object containing the ranges to process
#' @param col_names Column names on which basis to check whether two adjacent ranges are different (typically c(loh_col, cnv_col))
#' @param threshold Numeric value indicating the ratio of longer to shorter adjacent regions which result in removal of short region
#' @import GenomicRanges

run_filter <- function(gr, col_names, threshold = 100) {
    ## Identify the short ranges with equivalent neighbours
    gr <- add_neighbours_column(gr, col_names)
    gr <- add_length_ratio_column(gr)
    
    indices_to_remove <- which(gr$neighbours_equivalent & (gr$length_ratio > threshold))
    new_gr_list <- lapply(indices_to_remove, function(index) {
        if (index == 1) {
            ref_index <- 2
            index_range <- index : (index+1)
        } else {
            ref_index <- index - 1
            if (index == length(gr)) {
                index_range <- (index-1) : index
            } else {
                index_range <- (index-1) : (index+1)
            }
        }

        values <- as.data.frame(values(gr[ref_index, col_names]))
        new <- attach_values(reduce(gr[index_range]), col_names, values)

        return(new)
    })
    new_gr <- do.call('c', new_gr_list)
    if (! is.null(new_gr)) {
        gr <- sort(c(gr[! overlapsAny(gr, new_gr)][, col_names], new_gr))
        gr <- stitch(gr, col_names)
    }
    return(gr)
}


#' Adds a column to the data which indicates how much longer surrounding ranges are than the given range
#' Helps facilitate the filtering functions.
#'
#' @param gr GRanges object containing the ranges to process
#' @import GenomicRanges
#' @importFrom plyr ddply

add_length_ratio_column <- function(gr) {
    length_ratio <- ddply(as.data.frame(gr), 'seqnames', function(chr_data) {
        data.frame(length_ratio = sapply(1:dim(chr_data)[1], function(i) {
            if (i == 1) {
                ratio = chr_data[i + 1, 'width'] / chr_data[i, 'width']
            } else if (i == dim(chr_data)[1]) {
                ratio = chr_data[i - 1, 'width'] / chr_data[i, 'width'] 
            } else {
                ratio = ( chr_data[i - 1, 'width'] + chr_data[i + 1, 'width'] ) / chr_data[i, 'width']
            }
            return(ratio)
        }))
    })
    values(gr) <- cbind(values(gr), length_ratio[2])
    return(gr)
}

#' Adds a column to the data which indicates whether neighbouring ranges are identical to each other.
#' Assists in the filtering functions.
#'
#' @param gr GRanges object containing the ranges to process
#' @param col_names Column names on which basis to determine whether neighbouring regions are identical
#' @import GenomicRanges 
#' @importFrom plyr ddply

add_neighbours_column <- function(gr, col_names) {
    neighbours_equivalent <- ddply(as.data.frame(gr), 'seqnames', function(chr_data) {
        data.frame(neighbours_equivalent = sapply(1:dim(chr_data)[1], function(i) {
            return( i != 1 && i != dim(chr_data)[1] && chr_data[i-1, col_names] == chr_data[i+1, col_names] )
        }))
    })
    values(gr) <- cbind(values(gr), neighbours_equivalent[2])
    return(gr)
}

#' Combines identical neighbours into a single range - performed at the end of contig and filter steps
#'
#' @param gr GRanges object containing the ranges to process
#' @param col_names Column names on which basis to determine whether neighbouring regions are identical
#' @import GenomicRanges

stitch <- function(gr, col_names) {
    unique_types <- unique(as.data.frame(values(gr))[, col_names])
    reduced_gr_list <- apply(unique_types, 1, function(u) {
          matching_rows <- apply(as.data.frame(values(gr))[, col_names], 1, function(z) { all(z == u) })
          gr_subset <- reduce(gr[matching_rows])
          gr_subset <- attach_values(gr_subset, col_names, u)
          return(gr_subset)
    })
    names(reduced_gr_list) <- NULL
    gr <- do.call('c', reduced_gr_list)
    gr <- sort(gr)
    return(gr)
}

#' Attaches new columns to a GRanges object
#'
#' @param gr GRanges object containing the ranges to process
#' @param col_names The names of new columns, in order
#' @param values The values to fill new columns with
#' @import GenomicRanges

attach_values <- function(gr, col_names, values) {
    if (is.null(dim(values))) {
        values <- as.data.frame(t(values))
    }
    
    if (dim(values)[2] == 1) { 
        values = rep(values[1, 1], length(gr))
        new_values <- as.data.frame(values)
    } else {
        new_values <- values[rep(1, length(gr)), ]
    }

    colnames(new_values) <- col_names
    rownames(new_values) <- NULL
    new_values <- as.data.frame(new_values)
    values(gr) <- cbind(values(gr), new_values)
    return(gr)
}


