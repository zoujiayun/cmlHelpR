#' Calculate LRT statistic for models from each group (nsm/bcm)
#'
#' Calcualtes LRT statistic between two models, where one belongs to Branch/Clade/Branch-site and the other
#' to Null/Site. Name derived from "half-n-half".
#' @param df Long-format dataframe to parse
#' @param mdl Vector of models for LRT statistic
#' @keywords internal
#' @importFrom rlang .data
#' @importFrom rlang :=
#' @export
.lrt_hnh <- function(df, mdl){
  o <- dplyr::filter(.data = df, stringr::str_detect(string = .data$model, pattern = paste0(mdl, "$", collapse = "|")))

  ## Cleaning null/site
  r <- dplyr::filter(.data = o, !stringr::str_detect(string = .data$id, pattern = "_")) ## Selecting row without tree
  m <- dplyr::pull(.data = r, .data$model) ## Model ID
  m <- unique(x = m)

  newNP <- paste0("np_", m)   ## New column names
  newlnL <- paste0("lnL_", m) ## New column names

  r <- dplyr::rename(.data = r,
                     gene = .data$id, ## Renaming columns
                     !! newNP := .data$np,
                     !! newlnL := .data$lnL)
  r <- dplyr::select(.data = r, -.data$model)    ## Removing unneeded column

  ## Cleaning branch/-site/clade
  o <- dplyr::filter(.data = o, stringr::str_detect(string = .data$id, pattern = "_")) ## Row with tree information
  o <- tibble::as_tibble(data.table::dcast(data.table::setDT(o), id ~ model, value.var = c("np", "lnL")))
  o <- tidyr::separate(data = o, col = .data$id, into = c("gene", "tree"), sep = "_")

  ## Joining to form comparison dataframe
  o <- dplyr::left_join(x = o, y = r, by = "gene")
  o <- dplyr::select(.data = o, .data$gene, .data$tree, dplyr::starts_with("np"), dplyr::starts_with("lnL"))

  ## Model comparisons
  o <- dplyr::mutate(.data = o,
                     delta = (2*(abs(o[[6]] - o[[5]]))),
                     degFree = abs(o[[4]] - o[[3]]),
                     pval = stats::pchisq(.data$delta, .data$degFree, lower.tail=FALSE))

  return(o)
}
