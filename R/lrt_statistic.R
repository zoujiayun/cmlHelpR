#' Description
#'
#' General run down of what it does
#' @param dir_path Path to codeml_parallel output directory.
#' @param models A list of models to use in the LRT calculations
#' @param lst_comparisons Combinations/comparisons of the above models
#' @keywords Likelihood ratio test, LRT
#' @export
#' @examples
#' lrt_statistic(dir_path = "path/to/codeml_out", models = c("M1a", "M2a"), lst_comparisons = c("M1a", "M2a"))
lrt_statistic <- function(dir_path, models, lst_comparisons) {

  ## Importing all log files
  dirs <- list.dirs(path = dir_path, full.names = TRUE)[stringr::str_detect(string = list.dirs(path = dir_path, full.names = TRUE), pattern = paste(models, "$", sep = "", collapse = "|"))]
  dirs <- magrittr::set_names(x = dirs, value = models)

  print("Reading output files")
  files <- purrr::map(dirs, ~{
    fl <- list.files(path = .x, pattern = ".out", full.names = TRUE, recursive = TRUE)
    fl <- magrittr::set_names(x = fl, value = sub(".out", "", basename(fl)))
    lapply(fl, readr::read_lines)
  })

  ## Extracting np + lnL values
  print("Extracting np + lnL")
  np_lnL <- purrr::map(files, .parse_np_lnL)
  np_lnL <- dplyr::bind_rows(np_lnL, .id = "model")

  ## Long to wide format on multiple variables
  print("Long to wide format")
  np_lnL <- data.table::dcast(data.table::setDT(np_lnL), gene_tree ~ model, value.var = c("np", "lnL"))
  np_lnL <- tibble::as_tibble(x = np_lnL)
  np_lnL <- tidyr::separate(data = np_lnL,
                            col = gene_tree,
                            into = c("gene","tree"),
                            sep = "_")

  ## Model comparisons
  print("Generating list comparisons")
  lst.comparisons <- purrr::map(lst_comparisons, ~{
    tmp.out <- dplyr::select(.data = np_lnL, gene, tree, dplyr::matches(paste(.x, "$", collapse = "|", sep = "")))
    tmp.out <- dplyr::mutate(.data = tmp.out,
                             delta = (2*(abs(tmp.out[[6]] - tmp.out[[5]]))),
                             df = abs(tmp.out[[4]] - tmp.out[[3]]),
                             pval = pchisq(delta, df, lower.tail=FALSE))

  })

  ## Assigning names to list object
  names(lst.comparisons) <- purrr::map(lst_comparisons, ~{paste(.x, collapse = "_")})
  return(lst.comparisons)

}
