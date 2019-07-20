#' Calculate LRT statistic for Null/Site models
#'
#' Builds a table from the np/lnL data and runs an LRT analysis between two models.
#' @param df Long form dataframe
#' @param mdl Vector of two null/site models for comparison
#' @keywords internal
#' @export
#' @examples
#' .lrt_nsm(df = df.long, mdl = c("M1a", "M2a"))
.lrt_nsm <- function(df, mdl){

  ## Parsing input dataframe
  o <- dplyr::filter(.data = df, stringr::str_detect(string = model, pattern = paste0(mdl, "$", collapse = "|")))
  o <- tibble::as_tibble(data.table::dcast(data.table::setDT(o), id ~ model, value.var = c("np", "lnL")))

  ## LRT calculation
  o <- dplyr::mutate(.data = o,
                      delta = (2*(abs(o[[5]] - o[[4]]))),
                      degFree = abs(o[[3]] - o[[2]]),
                      pval = pchisq(delta, degFree, lower.tail=FALSE))

  return(o)
}
