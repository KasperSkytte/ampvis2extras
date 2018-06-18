#' Merge replicate samples
#'
#' @param data (\emph{required}) Data list as loaded with \code{amp_load()}.
#' @param merge_var Variable in the sample metadata that defines the sample groups
#' @param merge_fun Name of the function used to merge the samples, fx \code{mean} or \code{sum}. (\emph{default:} \code{"mean"})
#' @param round Read counts are whole numbers, so if any decimals after merging they must be rounded either \code{"up"} or \code{"down"}. (\emph{default:} \code{"up"})
#'
#' @return An object of class \code{ampvis2}
#'
#' @importFrom magrittr %>%
#' @importFrom rlang enquo
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
#'                                merge_var = "group",
#'                                merge_fun = "mean",
#'                                round = "up")
#' dmerged
amp_mergereplicates <- function(data,
                                merge_var,
                                merge_fun = "mean",
                                round = "up") {
  ### Data must be in ampvis2 format
  if(class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)", call. = FALSE)
  ### merge_var must be a length 1 string
  if(!is.character(merge_var) || length(merge_var) != 1)
    stop("\"merge_var\" must be a character string of length 1", call. = FALSE)
  #capture merge function
  merge_fun <- rlang::enquo(merge_fun)
  #find unique group names excluding than NA and empty strings ("")
  groups <- data$metadata[which(data$metadata[,merge_var] != ""), merge_var] %>%
    na.omit() %>%
    unique()
  #retain name of first column
  nameoffirstcol <- colnames(data$metadata)[1]
  #add a group column to the abundance table to define the groups
  tempabund <- data$metadata[,c(1,which(colnames(data$metadata) == merge_var))] %>%
    merge(as.data.frame(t(data$abund)), by.x = 1, by.y = 0) %>% {
      .[,1] <- ifelse(.[,2] %in% groups, .[,2], .[,1])
      colnames(.)[1] <- nameoffirstcol
      return(.[,-2])
    }
  #merge samples based on the groups
  abund_merged <- tempabund[which(tempabund[,1] %in% groups),] %>%
    dplyr::group_by(!!! syms(nameoffirstcol)) %>%
    dplyr::summarise_at(.vars = 2:ncol(.), .funs = merge_fun)
  #combine the merged samples and unmerged samples if any
  if(any(!data$metadata[,merge_var] %in% groups)) {
    newabund <- dplyr::bind_rows(abund_merged, tempabund[which(!tempabund[,1] %in% groups),])
  } else
    newabund <- abund_merged
  #transpose the new abundance table, remove first column and round up or down
  out <- data
  out$abund <- newabund %>%
    as.data.frame() %>%
    tibble::column_to_rownames(nameoffirstcol) %>%
    t() %>% {
      if(round == "up") {
        return(ceiling(.))
      } else if(round == "down") {
        return(floor(.))
      } else
        return(.)
    }
  #make new metadata keeping only the first row of each sample group
  out$metadata <- out$metadata[!duplicated(out$metadata$group),]
  out$metadata[,1] <- colnames(out$abund)
  rownames(out$metadata) <- out$metadata[,1]
  return(out)
}
