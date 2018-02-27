#' Beta-diversity plot
#'
#' @param data
#' @param distmeasure
#' @param label_by
#' @param show_values
#' @param values_size
#'
#' @return
#' @export
#'
#' @import ggplot2
#' @importFrom tidyr gather
#'
#' @examples
#' amp_betadiv(AalborgWWTPs, show_values = FALSE)
amp_betadiv <- function(data,
                        distmeasure = "bray",
                        label_by = NULL,
                        show_values = TRUE,
                        values_size = 3) {
  distcoeffs <- as.matrix(vegan::vegdist(t(data$abund), method = distmeasure))
  distcoeffs[upper.tri(distcoeffs)] <- NA
  if(!is.null(label_by)) {
    labels <- apply(data$metadata[,c(1, which(colnames(data$metadata) %in% label_by))], 1, paste0, collapse = "; ")
  } else if (is.null(label_by)) {
    labels <- rownames(distcoeffs)
  }
  distcoeffs <- data.frame(label = labels, distcoeffs, check.names = FALSE)
  levels <- as.factor(distcoeffs[["label"]])
  distcoeffs <- tidyr::gather(distcoeffs,
                              key = "y",
                              value = "coeff",
                              -label,
                              na.rm = TRUE,
                              factor_key = TRUE)
  distcoeffs$label <- factor(distcoeffs$label, levels = levels)

  p <- ggplot(distcoeffs,
              aes(x = label,
                  y = y,
                  fill = coeff,
                  label = round(1-coeff, 2))) +
    geom_tile(stat = "identity") +
    xlab("") +
    ylab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90,
                                     vjust = 0.5,
                                     hjust = 1)) +
    scale_x_discrete(labels = levels) +
    scale_y_discrete(labels = levels) +
    scale_fill_continuous(low = "White", high = "red")
  if(isTRUE(show_values)) {
    p <- p +
      geom_text(size = values_size) +
      theme(legend.position = "none")
  }
  return(p)
}
