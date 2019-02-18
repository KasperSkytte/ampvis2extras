#' Merge replicate samples
#'
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' @param merge_var Variable in the sample metadata that defines the sample groups
#' @param merge_fun Name of the function used to merge the samples, fx \code{mean} or \code{sum}. (\emph{default:} \code{"mean"})
#' @param round If the read counts have not been normalised, read counts are integers, so if any decimals after merging they must be rounded either \code{"up"} or \code{"down"}. Make sure this makes sense if the read counts have been normalised, as it may result in 0's, 1's, and 2's everywhere. (\emph{default:} \code{NULL})
#'
#' @return An object of class \code{ampvis2}
#'
#' @importFrom magrittr %>%
#' @importFrom rlang enquo syms !!!
#' @importFrom dplyr group_by summarise_at bind_rows
#' @importFrom tibble column_to_rownames
#' @export
#'
#' @examples
#' library(ampvis2extras)
#' d <- amp_load(ampvis2::example_otutable, ampvis2::example_metadata)
#' d$metadata
#' d$metadata$group <- c("group1", "group1", "group2", "group2", "group2", "group3", "group4", "group4")
#' d$metadata
#' dmerged <- amp_mergereplicates(d,
#'   merge_var = "group",
#'   merge_fun = "mean",
#'   round = "up"
#' )
#' dmerged
amp_mergereplicates <- function(data,
                                merge_var,
                                merge_fun = "mean",
                                round = NULL) {
  ### Data must be in ampvis2 format
  if (class(data) != "ampvis2") {
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)", call. = FALSE)
  }

  ### merge_var must be a length 1 string
  if (!is.character(merge_var) || length(merge_var) != 1) {
    stop("\"merge_var\" must be a character string of length 1", call. = FALSE)
  }

  # warn when user rounds up or down if data is normalised
  if (isTRUE(attributes(data)$normalised) & any(tolower(round) %in% c("down", "up"))) {
    warning("The data has been normalised, rounding up or down likely does not make sense, at least if merge_fun is \"mean\"",
      call. = FALSE
    )
  }

  # capture merge function
  merge_fun <- rlang::enquo(merge_fun)

  # find unique group names excluding than NA and empty strings ("")
  groups <- data$metadata[which(data$metadata[, merge_var] != ""), merge_var] %>%
    na.omit() %>%
    unique()

  # retain name of first column
  nameoffirstcol <- colnames(data$metadata)[1]

  # add a group column to the abundance table to define the groups
  tempabund <- data$metadata[, c(1, which(colnames(data$metadata) == merge_var))] %>%
    merge(as.data.frame(t(data$abund)), by.x = 1, by.y = 0) %>%
    {
      .[, 1] <- ifelse(.[, 2] %in% groups, .[, 2], .[, 1])
      colnames(.)[1] <- nameoffirstcol
      return(.[, -2])
    }

  # merge samples based on the groups
  abund_merged <- tempabund[which(tempabund[, 1] %in% groups), ] %>%
    dplyr::group_by(!!!syms(nameoffirstcol)) %>%
    dplyr::summarise_all(.vars = list(2:ncol(.)), .funs = merge_fun)

  # combine the merged samples and unmerged samples if any
  if (any(!data$metadata[, merge_var] %in% groups)) {
    newabund <- dplyr::bind_rows(abund_merged, tempabund[which(!tempabund[, 1] %in% groups), ])
  } else {
    newabund <- abund_merged
  }

  # transpose the new abundance table, remove first column and round up or down
  out <- data
  out$abund <- newabund %>%
    as.data.frame() %>%
    tibble::column_to_rownames(nameoffirstcol) %>%
    t() %>%
    {
      if (is.null(round)) {
        return(.)
      } else if (tolower(round) == "up") {
        return(ceiling(.))
      } else if (tolower(round) == "down") {
        return(floor(.))
      } else {
        return(.)
      }
    } %>%
    as.data.frame()

  # make new metadata keeping only the first rowof each sample group
  merge_var_id <- which(colnames(out$metadata) %in% merge_var)
  out$metadata <- out$metadata[!duplicated(out$metadata[, merge_var_id]), -merge_var_id]
  out$metadata[, 1] <- colnames(out$abund)
  rownames(out$metadata) <- out$metadata[, 1]
  return(out)
}
