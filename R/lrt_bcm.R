#' Calculate LRT statistic for Branch/clade/branch-site models
#'
#' Builds a table from the np/lnL data and runs an LRT analysis between two models.
#' @param df Long-format dataframe to parse for model information
#' @param mdl Branch/clade/branch-site models to compare
#' @keywords internal
#' @export
#' @examples
#' .lrt_bcm(df = df.long, mdl = c("CmC", "M2a_Rel"))
.lrt_bcm <- function(df, mdl){

  ## Subset long-format dataframe for clade/brach/-site models
  o <- dplyr::filter(.data = df, stringr::str_detect(string = model, pattern = paste0(mdl, "$", collapse = "|")))

  ## Cast the data to wide format and separate on pre-determined column names
  o <- tibble::as.tibble(data.table::dcast(data.table::setDT(o), id ~ model, value.var = c("np", "lnL")))
  o <- tidyr::separate(data = o,
                        col = id,
                        into = c("gene","tree"),
                        sep = "_")

  ## LRT calculation
  o <- dplyr::mutate(.data = o,
                      delta = (2*(abs(o[[6]] - o[[5]]))),
                      degFree = abs(o[[4]] - o[[3]]),
                      pval = pchisq(delta, degFree, lower.tail=FALSE))
}
