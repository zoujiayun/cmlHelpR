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
  print("Gettting list of files")
  vec.logs <- list.files(path = dir_path, pattern = ".log", full.names = TRUE, recursive = TRUE)
  names(vec.logs) <- gsub(".+02_codeML/(.*)/codeml.log", "\\1", vec.logs)

  print("Reading file lines")
  lst.logs <- map(.x = vec.logs, .f = read_lines)

  ## Splitting list by models that have been run
  print("Splitting list of files by model")
  lst.logs.models <- map(lst_models, ~{
    lst.logs[grep(.x, names(lst.logs))]
  })
  lst.logs.models <- set_names(x = lst.logs.models, value = lst_models)

  ## Subsetting all lines by what is needed
  df.out <- map(names(lst.logs.models), ~{

    if(.x == "Model1Neutral") {

      print("Subsetting Model1Neutral")
      out.m1a <- map(lst.logs.models[[.x]], ~{

        ## Extracting df
        int.np <- .x[grepl("np =", .x)]
        int.np <- gsub("np =\\s+", "", int.np)
        int.np <- as.numeric(x = int.np)

        ## Extracting log-likelihood
        int.lnL <- .x[grepl("^lnL\\s+=", .x)]
        int.lnL <- gsub("lnL\\s+=\\s+", "", int.lnL)
        int.lnL <- as.numeric(int.lnL)

        df.out <- tibble(np = int.np, lnL = int.lnL)
        return(df.out)

      })
      bind_rows(out.m1a, .id = "gene_tree")

    } else if (.x == "Model2Selection") {

      print("Subsetting Model2Selection")
      out.m2a <- map(lst.logs.models[[.x]], ~{

        int.np <- .x[grepl("np =", .x)]
        int.np <- gsub("np = \\s+", "", int.np)
        int.np <- as.numeric(int.np)

        int.lnL <- .x[grepl("lnL\\s+=", .x)]
        int.lnL <- gsub("lnL\\s+=\\s+", "", int.lnL)
        int.lnL <- as.numeric(int.lnL)

        df.out <- tibble(np = int.np, lnL = int.lnL)
        return(df.out)

      })
      bind_rows(out.m2a, .id = "gene_tree")

    }

  })

  print("cleaning tree names")
  df.out <- set_names(x = df.out, value = lst_models)
  df.out <- bind_rows(df.out, .id = "model")
  df.out <- mutate(.data = df.out, gene_tree = gsub("/.*_", "_", gene_tree))

  ## Long to wide format on multiple variables
  print("Long to wide format")
  df.out <- data.table::dcast(data.table::setDT(df.out), gene_tree ~ model, value.var = c("np", "lnL"))
  df.out <- as_tibble(x = df.out)
  df.out <- separate(data = df.out,
                     col = gene_tree,
                     into = c("gene","tree"),
                     sep = "_")

  ## Model comparisons
  print("Generating list comparisons")
  lst.comparisons <- map(lst_comparisons, ~{

    tmp.out <- select(.data = df.out, gene, tree, matches(paste(.x, collapse = "|")))
    tmp.out <- mutate(.data = tmp.out,
                  delta = (2*(abs(tmp.out[[5]] - tmp.out[[6]]))),
                  df = 2, #abs(tmp.out[[3]] - tmp.out[[4]]),
                  pval = pchisq(delta, df, lower.tail=FALSE))

  })

  ## Assigning names to list object
  names(lst.comparisons) <- map(lst_comparisons, ~{paste(.x, collapse = "_")})
  return(lst.comparisons)

}
