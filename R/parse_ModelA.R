#' Parse BEB output from ModelA codeml output
#'
#' General run down of what it does
#' @param lines List of codeml output files
#' @param c number of cores
#' @keywords internal
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
  df <- dplyr::mutate(.data = df, value = stringr::str_replace_all(string = value, pattern = "\\s+", replacement = "_"))
  df <- dplyr::filter(.data = df, value != "")
  df <- suppressWarnings(tidyr::separate(data = df, col = value, into = c("pos", "aa", "val", "postMean", "plusMinus", "SE"), sep = "_"))
  df <- dplyr::mutate(.data = df,
                      signif = stringr::str_extract(string = val, pattern = "\\*+"),
                      val = suppressWarnings(as.numeric(stringr::str_remove_all(string = val, pattern = "\\*"))),
                      pval = 1 - val,
                      aa = dplyr::na_if(x = aa, "NA"))
  df <- suppressWarnings(dplyr::mutate_at(.tbl = df, .vars = c("pos", "postMean", "SE"), .funs = as.numeric))
  df <- tidyr::separate(data = df, col = id, into = c("gene", "tree"), sep = "_")
  df <- dplyr::select(.data = df, gene, tree, pos, aa, val, pval, signif, postMean, SE)

  return(df)

}
