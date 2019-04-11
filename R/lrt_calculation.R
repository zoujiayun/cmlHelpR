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

  # ## Subsetting all lines by what is needed
  # df.out <- purrr::map(names(lst.logs.models), ~{
  #
  #   if(.x == "Model1Neutral") {
  #
  #     print("Subsetting Model1Neutral")
  #     out.m1a <- purrr::map(lst.logs.models[[.x]], ~{
  #
  #       ## Extracting df
  #       int.np <- .x[grepl("np =", .x)]
  #       int.np <- gsub("np =\\s+", "", int.np)
  #       int.np <- as.numeric(x = int.np)
  #
  #       ## Extracting log-likelihood
  #       int.lnL <- .x[grepl("^lnL\\s+=", .x)]
  #       int.lnL <- gsub("lnL\\s+=\\s+", "", int.lnL)
  #       int.lnL <- as.numeric(int.lnL)
  #
  #       df.out <- tibble::tibble(np = int.np, lnL = int.lnL)
  #       return(df.out)
  #
  #     })
  #     dplyr::bind_rows(out.m1a, .id = "gene_tree")
  #
  #   } else if (.x == "Model2Selection") {
  #
  #     print("Subsetting Model2Selection")
  #     out.m2a <- purrr::map(lst.logs.models[[.x]], ~{
  #
  #       int.np <- .x[grepl("np =", .x)]
  #       int.np <- gsub("np = \\s+", "", int.np)
  #       int.np <- as.numeric(int.np)
  #
  #       int.lnL <- .x[grepl("lnL\\s+=", .x)]
  #       int.lnL <- gsub("lnL\\s+=\\s+", "", int.lnL)
  #       int.lnL <- as.numeric(int.lnL)
  #
  #       df.out <- tibble::tibble(np = int.np, lnL = int.lnL)
  #       return(df.out)
  #
  #     })
  #     dplyr::bind_rows(out.m2a, .id = "gene_tree")
  #
  #   }
  #
  # })
  #
  # print("cleaning tree names")
  # df.out <- magrittr::set_names(x = df.out, value = lst_models)
  # df.out <- dplyr::bind_rows(df.out, .id = "model")
  # # df.out <- dplyr::mutate(.data = df.out, gene_tree = gsub("/.*_", "_", gene_tree))
  # df.out <- dplyr::mutate(.data = df.out, gene_tree = gsub(".+/(.*)-.+-.+/Model.+_(.*)/codeml.log", "\\1_\\2", gene_tree))
  #
  # ## Long to wide format on multiple variables
  # print("Long to wide format")
  # df.out <- data.table::dcast(data.table::setDT(df.out), gene_tree ~ model, value.var = c("np", "lnL"))
  # df.out <- tibble::as_tibble(x = df.out)
  # df.out <- tidyr::separate(data = df.out,
  #                    col = gene_tree,
  #                    into = c("gene","tree"),
  #                    sep = "_")
  #
  # ## Model comparisons
  # print("Generating list comparisons")
  # lst.comparisons <- purrr::map(lst_comparisons, ~{
  #
  #   tmp.out <- dplyr::select(.data = df.out, gene, tree, dplyr::matches(paste(.x, collapse = "|")))
  #   tmp.out <- dplyr::mutate(.data = tmp.out,
  #                 delta = (2*(abs(tmp.out[[5]] - tmp.out[[6]]))),
  #                 df = 2, #abs(tmp.out[[3]] - tmp.out[[4]]), ## Hard-coded as was told this is always the value...?
  #                 pval = pchisq(delta, df, lower.tail=FALSE))
  #
  # })
  #
  # ## Assigning names to list object
  # names(lst.comparisons) <- purrr::map(lst_comparisons, ~{paste(.x, collapse = "_")})
  # return(lst.comparisons)

}

t <- lrt_statistic(dir_path = "~/Desktop/test_zone/codeml_out", lst_models = c("M1a", "M2a"), lst_comparisons = list(c("M2a", "M1a")))
