#' Parse dN/dS branch values from output files
#'
#' Given a list of files, this function will parse the branch dN/dS value from the
#' CODEML output. It returns a list and nested dataframe object of dN/dS values for each file.
#' @param file_list Vector of file names with path
#' @keywords helper
#' @export
#' @examples
#' getBranchdNdS(file_list = list.files(path = "path/to/files"))
getBranchdNdS <- function(file_list) {

  ## Accepted models
  g <- c("FreeRatio" ,"M0" ,"M1a" ,"M2a_Rel" ,"M3" ,"M7" ,"M8" , "TwoRatio")
  g <- paste0("/", g, "/", collapse = "|")

  ## Filter out models not in g
  f <- file_list[str_detect(string = file_list, pattern = g)]
  f <- set_names(x = f, value = str_remove(string = sub(".*/(.*)/(.*)$", "\\1_\\2", f), pattern = ".out"))
  f <- map(f, read_lines)

  ## Iterate through each file and get values
  d <- map(f, ~{
    p1 <- grep(pattern = "^ branch", x = .x)
    p2 <- grep(pattern = "^Naive\ Empirical\ Bayes\ \\(NEB\\)\ analysis", .x) - 1
    v <- trimws(x = .x[p1:p2], which = "left")
    v <- v[v != ""]
    v <- str_replace_all(string = v, pattern = "\\s+", replacement = "\t")
    read_tsv(v)
  })

  ## Build dataframe - nested
  df <- dplyr::bind_rows(d, .id = "id")
  df <- tidyr::separate(data = df, col = id, into = c("model", "gene", "tree"), sep = "_")
  df <- dplyr::group_by(.data = df, model, gene, tree)
  df <- tidyr::nest(data = df, .key = "branch")

  ## Return list + dataframe
  lst <- list(list = d, dataframe = df)
  return(lst)

}
