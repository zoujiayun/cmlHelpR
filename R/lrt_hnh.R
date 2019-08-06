#' Calculate LRT statistic for models from each group (nsm/bcm)
#'
#' Calcualtes LRT statistic between two models, where one belongs to Branch/Clade/Branch-site and the other
#' to Null/Site. Name derived from "half-n-half".
#' @param df Long-format dataframe to parse
#' @param mdl Vector of models for LRT statistic
#' @keywords internal
#' @export
.lrt_hnh <- function(df, mdl){
  o <- dplyr::filter(.data = df, stringr::str_detect(string = model, pattern = paste0(mdl, "$", collapse = "|")))

  ## Cleaning null/site
  r <- dplyr::filter(.data = o, !stringr::str_detect(string = id, pattern = "_")) ## Selecting row without tree
  m <- dplyr::pull(.data = r, model) ## Model ID
  m <- unique(m)

  newNP <- paste0("np_", m)   ## New column names
  newlnL <- paste0("lnL_", m) ## New column names

  r <- dplyr::rename(.data = r,
                     gene = id, ## Renaming columns
                     !! newNP := np,
                     !! newlnL := lnL)
  r <- dplyr::select(.data = r, -model)    ## Removing unneeded column

  ## Cleaning branch/-site/clade
  o <- dplyr::filter(.data = o, stringr::str_detect(string = id, pattern = "_")) ## Row with tree information
  o <- tibble::as_tibble(data.table::dcast(data.table::setDT(o), id ~ model, value.var = c("np", "lnL")))
  o <- tidyr::separate(data = o, col = id, into = c("gene", "tree"), sep = "_")

  ## Joining to form comparison dataframe
  o <- dplyr::left_join(x = o, y = r, by = "gene")
  o <- dplyr::select(.data = o, gene, tree, dplyr::starts_with("np"), dplyr::starts_with("lnL"))

  ## Model comparisons
  o <- dplyr::mutate(.data = o,
                     delta = (2*(abs(o[[6]] - o[[5]]))),
                     degFree = abs(o[[4]] - o[[3]]),
                     pval = pchisq(delta, degFree, lower.tail=FALSE))

  return(o)
}
