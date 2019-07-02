#' Create a frequency table from the parsed Bayes Empiracle Bayes data frame
#'
#' Creates a frequency table based on the parsed output from the parse_BEB() function
#' @param df_long data frame from parse_BEB() function
#' @param only_signif Set to TRUE to select only significant sites
#' @param max_freq Keep genes who have n-sites under selection up to or equal to this
#' @keywords frequency table
#' @export
#' @examples
#' beb2freq(df_long = df.beb, only_signif = TRUE)
beb2freq <- function(df_long, significant = NULL) {

  ## Filtering on significance
  if (isTRUE(significant)) {
    df_long <- dplyr::filter(.data = df_long, !is.na(signif))
  }

  ## Getting frequency
  df <- dplyr::mutate(.data = df_long, id = dplyr::if_else(rowSums(is.na(df_long[4:10])) == 7, 0, 1))
  df <- dplyr::group_by(.data = df, model, gene, tree, id)
  df <- tidyr::nest(data = df, .key = beb)
  df <- dplyr::mutate(.data = df, freq = unlist(purrr::map(beb, dplyr::tally)))
  df <- dplyr::mutate(.data = df, freq = id * freq)
  df <- dplyr::select(.data = df, model, gene, tree, freq)
  df <- tidyr::replace_na(data = df, replace = 0)

  ## Building plotting format
  plt <- tidyr::unite(data = df, model_tree, c("model", "tree"))
  plt <- tidyr::spread(data = plt, model_tree,  freq)
  plt <- tibble::column_to_rownames(.data = plt, "gene")
  plt <- tibble::as_tibble(x = plt, rownames = NA)

  ## Object to return
  lst <- list(long = df, heatmap = plt)
  return(lst)

}
