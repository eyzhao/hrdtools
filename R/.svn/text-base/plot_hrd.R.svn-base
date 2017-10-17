#' Plots pre-filtered and post-filtered states of LOH 
#' 
#' @param pre_filtered_gr GRanges object, prior to filtering
#' @param post_filtered_gr GRanges object representing pre_filtered_gr after contig and filtering steps
#' @param copy.col Name of the column with copy number information, defaults to 'copy_number'
#' @param loh.col Name of the column with LOH information, defaults to 'lohtype'
#'
#' @import ggplot2

get_filtering_plot <- function(pre_filtered_gr, post_filtered_gr, copy.col = 'copy_number', loh.col = 'lohtype') {
    pre_filtered <- as.data.frame(pre_filtered_gr)[, c('seqnames', 'start', 'end', copy.col, loh.col)]
    post_filtered <- as.data.frame(post_filtered_gr)[, c('seqnames', 'start', 'end', copy.col, loh.col)]
    pre_filtered$copy_number <- pre_filtered[, copy.col]; post_filtered$copy_number <- post_filtered[, copy.col]
    pre_filtered$lohtype <- pre_filtered[, loh.col]; post_filtered$lohtype <- post_filtered[, loh.col]

    pre_filtered$type <- 'LOH calls'
    post_filtered$type <- 'filtered'
    combined <- rbind(pre_filtered, post_filtered)
    combined$type <- factor(as.character(combined$type), levels = c('LOH calls', 'filtered'))

    combined$lohtype <- factor(as.character(combined$lohtype), levels = c('HET', 'BCNA', 'ASCNA', 'DLOH', 'NLOH', 'ALOH', 'HOMD'))

    p <- ggplot(combined) +
        geom_rect(aes(xmin = start, xmax = end, 
                      ymin = as.numeric(copy_number) - 0.5, ymax = as.numeric(copy_number) + 0.5,
                      fill = lohtype)) +
        scale_fill_brewer(palette = 'Set1') +
        facet_grid(type ~ seqnames, space = 'free_x', scales = 'free_x') +
        theme(panel.spacing.x = unit(0, 'lines'),
              axis.text.x=element_text(angle = 45, hjust = 1, size=6))
    return(p)
}

#' Returns a ggplot2 object illustrating the LOH, TAI, and LST calls
#' 
#' @param post_filtered_gr GRanges object representing LOH calls after contig and filter steps
#' @param loh.col Name of the column with LOH information, defaults to 'lohtype'
#'
#' @import ggplot2

get_hrd_score_plot <- function(post_filtered_gr, loh.col = 'lohtype') {
    df <- as.data.frame(post_filtered_gr)
    df$is_loh <- loh_test_bool(post_filtered_gr)
    df$is_tai <- tai_test_bool(post_filtered_gr)
    df$is_lst <- lst_test_bool(post_filtered_gr, loh.col)
    df$lohtype <- factor(as.character(df[, loh.col]), levels = c('HET', 'BCNA', 'ASCNA', 'DLOH', 'NLOH', 'ALOH', 'HOMD'))

    df$region_call <- rep('none', dim(df)[1])
    df$region_call[df$is_loh] <- 'loh'
    df$region_call[df$is_tai] <- 'tai'
    df$region_call[df$is_loh & df$is_tai] <- 'both'
    df$region_call <- factor(df$region_call, levels = c('none', 'loh', 'tai', 'both'))

    p <- ggplot(data = df) +
         geom_segment(aes(x = start, xend = end, colour = lohtype, y = region_call, yend = region_call),
                      size=3, lineend = 'butt') + 
         scale_colour_brewer(palette = 'Set1') +
         geom_segment(data = df[df$is_lst, ], aes(x = end, xend = end, y = 1, yend=4)) +
         facet_grid(. ~ seqnames, space = 'free_x', scales = 'free_x') +
         theme(panel.spacing.x = unit(0, 'lines'),
               axis.text.x=element_text(angle = 45, hjust = 1, size=5),
               axis.title.y=element_blank()) +
         xlab(NULL)

    return(p)
}

#' Outputs a plot combining the filtering plots and HRD score plots
#' 
#' @param pre_filtered_gr GRanges object, prior to filtering
#' @param post_filtered_gr GRanges object representing pre_filtered_gr after contig and filtering steps
#' @param output Path to output plot file (pdf)
#' @param copy.col Name of the column with copy number information, defaults to 'copy_number'
#' @param loh.col Name of the column with LOH information, defaults to 'lohtype'
#'
#' @import ggplot2
#' @export

export_process_plot <- function(pre_filtered_gr, post_filtered_gr, output, copy.col = 'copy_number', loh.col = 'lohtype') {
    p1 <- get_filtering_plot(pre_filtered_gr, post_filtered_gr)
    p2 <- get_hrd_score_plot(post_filtered_gr)

    pdf(output, width = 15, height = 5)
    print(plot_grid(p1, p2, labels = c('A', 'B'), align = 'v', ncol = 1))
    dev.off()
}
