#' Parse BEB output from ModelA codeml output
#'
#' Extracts the BEB information from ModelA codeml output files.
#' Uses anchors specific to Model A to parse data.
#' @param lines List of codeml output files
#' @param c number of cores
#' @keywords internal
#' @importFrom rlang .data
.parse_ModelA <- function(lines, c) {

  obj <- parallel::mclapply(lines, mc.cores = c, function(y) {

    s <- grep(pattern = "Bayes Empirical Bayes", x = y)
    e <- grep(pattern = "^The grid", x = y)

    if (y[s + 2] == "") {

      o <- c("NA NA NA NA NA NA")

    } else {

      ## Need to figure what the offset is below
      o <- y[(s+2):(e-3)]
      o <- trimws(x = o)
    }

  })

  df <- purrr::map(.x = obj, .f = tibble::enframe, name = NULL)
  df <- dplyr::bind_rows(df, .id = "id")
  df <- dplyr::mutate(.data = df, value = stringr::str_replace_all(string = .data$value, pattern = "\\s+", replacement = "_"))
  df <- dplyr::filter(.data = df, .data$value != "")
  df <- suppressWarnings(tidyr::separate(data = df, col = .data$value, into = c("pos", "aa", "val", "postMean", "plusMinus", "SE"), sep = "_"))
  df <- dplyr::mutate(.data = df,
                      signif = stringr::str_extract(string = .data$val, pattern = "\\*+"),
                      val = suppressWarnings(as.numeric(stringr::str_remove_all(string = .data$val, pattern = "\\*"))),
                      pval = 1 - .data$val,
                      aa = dplyr::na_if(x = .data$aa, "NA"))
  df <- suppressWarnings(dplyr::mutate_at(.tbl = df, .vars = c("pos", "postMean", "SE"), .funs = as.numeric))
  # df <- tidyr::separate(data = df, col = id, into = c("gene", "tree"), sep = "_")
  # df <- dplyr::select(.data = df, gene, tree, pos, aa, val, pval, signif, postMean, SE)
  df <- dplyr::select(.data = df, .data$id, .data$pos, .data$aa, .data$val, .data$pval, .data$signif, .data$postMean, .data$SE)

  return(df)

}
