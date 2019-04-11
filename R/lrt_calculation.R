#' Description
#'
#' General run down of what it does
#' @param dir_path Path to codeml_parallel output directory.
#' @param lst_models A list of models to use in the LRT calculations
#' @param lst_comparisons Combinations/comparisons of the above models
#' @keywords Likelihood ratio test, LRT
#' @export
#' @examples
#' lrt_statistic(dir_path = "path/to/codeml_out", lst_models = c("Model1Neutral", "Model2Selection"), lst_comparisons = c("Model1Neutral", "Model2Selection"))
lrt_statistic <- function(dir_path, lst_models, lst_comparisons) {

  ## Importing all log files
  print("Getting list of files")
  vec.logs <- list.files(path = dir_path, pattern = ".log", full.names = TRUE, recursive = TRUE)
  names(vec.logs) <- gsub(".+/(.*)/(.*)/codeml.log", "\\1::\\2", vec.logs)

  print("Reading file lines")
  lst.logs <- purrr::map(.x = vec.logs, .f = readr::read_lines)

  ## Splitting list by models that have been run
  print("Splitting list of files by model")
  lst.logs.models <- purrr::map(lst_models, ~{
    lst.logs[grep(paste0(.x, "_[a-zA-Z]+$"), names(lst.logs))]
  })
  lst.logs.models <- magrittr::set_names(x = lst.logs.models, value = lst_models)

  ## Getting np and lnL values from output files
  df.out <- purrr::map(lst.logs.models, .parse_np_lnL)

  print("cleaning tree names")
  df.out <- magrittr::set_names(x = df.out, value = lst_models)
  df.out <- dplyr::bind_rows(df.out, .id = "model")
  df.out <- dplyr::mutate(.data = df.out, gene_tree = gsub("(.*)-.+-.+_(.*)", "\\1_\\2", gene_tree))

  ## Long to wide format on multiple variables
  print("Long to wide format")
  df.out <- data.table::dcast(data.table::setDT(df.out), gene_tree ~ model, value.var = c("np", "lnL"))
  df.out <- tibble::as_tibble(x = df.out)
  df.out <- tidyr::separate(data = df.out,
                     col = gene_tree,
                     into = c("gene","tree"),
                     sep = "_")

  ## Model comparisons
  print("Generating list comparisons")
  lst.comparisons <- purrr::map(lst_comparisons, ~{

    tmp.out <- dplyr::select(.data = df.out, gene, tree, dplyr::matches(paste(.x, collapse = "|")))
    tmp.out <- dplyr::mutate(.data = tmp.out,
                  delta = (2*(abs(tmp.out[[6]] - tmp.out[[5]]))),
                  df = abs(tmp.out[[4]] - tmp.out[[3]]),
                  pval = pchisq(delta, df, lower.tail=FALSE))

  })

  ## Assigning names to list object
  names(lst.comparisons) <- purrr::map(lst_comparisons, ~{paste(.x, collapse = "_")})
  return(lst.comparisons)

}
