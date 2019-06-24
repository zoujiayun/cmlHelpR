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

comp <- readr::read_csv(file = "~/Documents/pipelines/pararllel_codeml/ctl_parameter_files/190410_model_comparisons.csv") %>%
  select(3:4) %>%
  filter(comparison != "NA") %>%
  group_split(id) %>%
  map(~flatten_chr(.), . == "")

comp <- comp[4]

t <- lrtStatistic(dir_path = "~/Documents/pipelines/pararllel_codeml/test_data/test_output", lst_comparisons = comp)

lrtStatistic <- function(dir_path, lst_comparisons) {

  ## Models
  m <- unlist(lst_comparisons)

  ## Selecting output dirs for analysis
  d <- list.dirs(path = dir_path, full.names = TRUE)
  d <- d[stringr::str_detect(string = d, pattern = paste(m, "$", sep = "", collapse = "|"))]

  ## Naming directory vectors
  n <- sub(pattern = ".*/", "", d)
  d <- magrittr::set_names(x = d, value = n)

  ## Reading codeml output files
  print("Reading output")
  f <- purrr::map(d, ~{
    fl <- list.files(path = .x, pattern = ".out", full.names = TRUE)
    fl <- magrittr::set_names(x = fl, value = sub(".out", "", basename(fl)))
    fl <- purrr::map(fl, readr::read_lines)
  })

  ## Extracting np + lnL values
  print("Extracting np + lnL")
  np_lnL <- purrr::map(f, .parse_np_lnL)
  np_lnL <- dplyr::bind_rows(np_lnL, .id = "model")

  ## Long to wide format on multiple variables
  np_lnL <- data.table::dcast(data.table::setDT(np_lnL), gene_tree ~ model, value.var = c("np", "lnL"))
  np_lnL <- tibble::as_tibble(x = np_lnL)
  np_lnL <- tidyr::separate(data = np_lnL,
                            col = gene_tree,
                            into = c("gene","tree"),
                            sep = "_")

  ## Model comparisons
  print("Generating list comparisons")
  comp <- purrr::map(lst_comparisons, ~{
    tmp <- dplyr::select(.data = np_lnL, gene, tree, dplyr::matches(paste(.x, "$", collapse = "|", sep = "")))
    tmp <- dplyr::mutate(.data = tmp,
                         delta = (2*(abs(tmp[[6]] - tmp[[5]]))),
                         df = abs(tmp[[4]] - tmp[[3]]),
                         pval = pchisq(delta, df, lower.tail=FALSE))

  })

  ## Assigning names to list object
  names(comp) <- purrr::map(lst_comparisons, ~{paste(.x, collapse = "_")})
  return(comp)

}
