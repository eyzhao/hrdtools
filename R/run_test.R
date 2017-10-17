#' Runs the LOH, TAI, and LST Tests
#' 
#' This function loads an LOH segments file and performs the three HRD tests on it.
#' @param loh_file Path to an LOH file suitable for import_ranges()
#' @param genome ID of the genome being used (default: 'hg19')
#' @param silent If TRUE, does not print HRD scores, just returns them (default: 'FALSE')
#' @export

run_test <- function(loh_file,
                     file.col=NULL,
                     loh.col='lohtype', cnv.col='copy_number',
                     genome='hg19', silent=FALSE) {

    if (is.null(file.col)) {
        file.col <- c('chr', 'start', 'end', 'width', 'vars', 'copy_number', 'lohtype', 'allele1', 'allele2')
    }

    print('Importing data...')

    # Import data and sort

    if (typeof(loh_file) == "character") {
        loh_ranges = sort(import_ranges(loh_file, col.names=file.col))
    } else if (typeof(loh_file) == "S4") {
        loh_ranges = sort(loh_file)
    } else {
        stop('loh_file must be of type character or GRanges S4')
    }

    loh_ranges <- contig(loh_ranges, loh.col, cnv.col)    # Fill in seg gaps
    loh_ranges <- filter_ranges(loh_ranges, loh.col, cnv.col)          # Filter out very short segs

    lohValue = loh_test(loh_ranges)
    taiValue = tai_test(loh_ranges)
    lstValue = lst_test(loh_ranges, loh.col)

    testObj <- list(
        loh = lohValue,
        tai = taiValue,
        lst = lstValue,
        total = lohValue + taiValue + lstValue
    )

    if (!silent) {
        print(sprintf('LOH Score: %s', lohValue))
        print(sprintf('TAI Score: %s', taiValue))
        print(sprintf('LST Score: %s', lstValue))
        print(sprintf('Total HRD Score: %s', lohValue + taiValue + lstValue))
    }

    return(testObj)
}
