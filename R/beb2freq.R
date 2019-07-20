#' Create a frequency table from the parsed Bayes Empiracle Bayes data
#'
#' Creates a frequency table based on the parsed output from the getBEB() function in the cmlHelpR package.
#' @param df_long Long form data frame from the getBEB() function
#' @param only_signif Set to TRUE to select only significant sites
#' @param site_as_prop Return frequencies values as a proportion of the genes length instead of simple counts
#' @keywords frequency table
#' @export
#' @examples
#' getFreq(df_long = df.beb, significant = TRUE, site_as_prop = TRUE)
getFreq <- function(df_long, significant = NULL, site_as_prop = NULL) {

  ## Filtering on significance
  if (isTRUE(significant)) {
    df_long <- dplyr::filter(.data = df_long, !is.na(signif))
  }

  ## Getting frequency
  df <- dplyr::mutate(.data = df_long, id = dplyr::if_else(rowSums(is.na(df_long[4:10])) == 7, 0, 1))
  df <- dplyr::group_by(.data = df, model, gene, tree, id, seqLen)
  df <- tidyr::nest(data = df, .key = beb)
  df <- dplyr::mutate(.data = df, freq = unlist(purrr::map(beb, dplyr::tally)))
  df <- dplyr::mutate(.data = df, freq = id * freq)
  df <- dplyr::select(.data = df, model, gene, tree, freq, seqLen)

  if(!is.null(site_as_prop)){
    tmp <- dplyr::mutate(.data = df, prop = freq/seqLen)
    tmp <- dplyr::select(.data = tmp, model, gene, tree, prop)

    ## Building plotting format
    plt <- tidyr::unite(data = tmp, model_tree, c("model", "tree"))
    plt <- tidyr::spread(data = plt, model_tree,  prop)
    plt[is.na(plt)] <- 0
    plt <- tibble::column_to_rownames(.data = plt, "gene")
    plt <- tibble::as_tibble(x = plt, rownames = NA)

  } else {

    ## Selecting necessary columns
    tmp <- dplyr::select(.data = df, model, gene, tree, freq)

    ## Building plotting format
    plt <- tidyr::unite(data = tmp, model_tree, c("model", "tree"))
    plt <- tidyr::spread(data = plt, model_tree,  freq)
    plt[is.na(plt)] <- 0
    plt <- tibble::column_to_rownames(.data = plt, "gene")
    plt <- tibble::as_tibble(x = plt, rownames = NA)

  }

  ## Object to return
  lst <- list(long = df, heatmap = plt)
  return(lst)

}

