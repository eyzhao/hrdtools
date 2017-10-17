#' Plots the output of NNLS signature exposures
#' 
#' This function outputs a ggplot plot object based on output from the run_snv function.
#' @param exposures Output from run_snv function.
#' @param confidence_intervals If TRUE, will include confidence intervals in the plot
#' @import ggplot2
#' @export

plot_nnls_exposures <- function(exposures, confidence_intervals=TRUE, sampleName=NULL) {
    exposures_df <- as.data.frame(exposures); exposures_df <- exposures_df[rev(rownames(exposures_df)), ]
    exposures_df$names <- factor(gsub('Signature.', '', exposures_df$names),
                                 levels=gsub('Signature.', '', exposures_df$names))
    print(exposures_df$names)

    plot <- ggplot(exposures_df, aes_string(x = 'means', y = 'names')) +
        ylab('Signature') +
        geom_point() +
        geom_errorbarh(aes_string(xmin = 'lCI', xmax = 'uCI', height=0))

    if (!is.null(sampleName)) {
        plot <- plot + labs(title = sampleName)
    }

    if (sum(exposures_df$x) == 1) {
        plot <- plot + xlab('Exposure Fraction') + xlim(0, 1)
    } else {
        plot <- plot + xlab('Exposure (Mutation Count)')
    }

    return(plot)
}
