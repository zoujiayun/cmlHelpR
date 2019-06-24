#' Create a frequency table from the parsed Bayes Empiracle Bayes data frame
#'
#' Creates a frequency table based on the parsed output from the parse_BEB() function
#' @param df_beb data frame from parse_BEB() function
#' @param only_signif Set to TRUE to select only significant sites
#' @param max_freq Keep genes who have n-sites under selection up to or equal to this
#' @keywords frequency table
#' @export
#' @examples
#' BEBDF2Freq(df_beb = df.beb, only_signif = TRUE, max_freq = 100)

beb2freq <- function(df_beb, only_signif = NULL, max_freq = NULL) {

  ## Split on tree
  lst.byTree <- split(df_beb, df_beb$tree)

  ## Filter to significant sites only
  if (isTRUE(only_signif)) {

    lst.byTree <- map(lst.byTree, filter, signif != "NA")

  }

  ## Get frequencies
  lst.freq <- map(names(lst.byTree), ~{

    tbl.consistent <- table(lst.byTree[[.x]]$gene)             ## Table of frequencies
    df.consistent <- as_tibble(as.data.frame(tbl.consistent))  ## Convert to df
    df.consistent <- rename(.data = df.consistent,             ## Renaming columns
                            gene = Var1,
                            !! .x := Freq)
  })

  df.freq <- reduce(lst.freq, full_join)
  df.freq <- replace(df.freq, is.na(df.freq), 0)

  ## Filter by maximum value
  if (!is.null(max_freq)) {

    df.freq <- filter_if(.tbl = df.freq, is.numeric, all_vars(. <= max_freq))

  }

  ## Convert to data.frame with rownames
  df.freq <- column_to_rownames(.data = df.freq, var = "gene")

  return(df.freq)

}
