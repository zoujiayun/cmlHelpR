#' Parse dN/dS branch values from output files
#'
#' Given a list of files, this function will parse the branch dN/dS value from the
#' CODEML output. It returns a list and nested dataframe object of dN/dS values for each file.
#' @param directory_path Vector of file names with path
#' @param models Models to get branch dN/dS values for
#' @param ext Extension of CODEML output files
#' @keywords helper
#' @export
#' @examples
#' getBranchdNdS(directory_path = list.files(path = "path/to/codeml/outDir", models = c("M1a", "M8")))
getBranchDNDS <- function(directory_path, models, ext = ".out") {

  ## Check: do the models provided match all contain branch dN/dS values
  g <- c("FreeRatio" ,"M0" ,"M1a" ,"M2a_Rel" ,"M3" ,"M7" ,"M8" , "TwoRatio")
  if(!all(is.element(models, g))){
    stop(call. = TRUE, "One of the specified models will not contain branch dN/dS values")
  }
  m <- paste0("/", models, collapse = "|")

  ## Listing all directories in directory path
  d <- list.dirs(path = directory_path, full.names = TRUE, recursive = TRUE)
  d <- d[stringr::str_detect(string = d, pattern = m)]

  ## Listing files + reading them in
  f <- purrr::map(d, list.files, pattern = ext, full.names = TRUE)
  f <- unlist(x = f)
  f <- magrittr::set_names(x = f, value = stringr::str_remove(string = sub(".*/(.*)/(.*)$", "\\1_\\2", f), pattern = ext))
  f <- purrr::map(f, readr::read_lines)

  ## Iterate through each file and get values
  d <- purrr::map(f, ~{
    p1 <- grep(pattern = "^ branch", x = .x)
    p2 <- grep(pattern = "^Naive\ Empirical\ Bayes\ \\(NEB\\)\ analysis", .x) - 1
    v <- trimws(x = .x[p1:p2], which = "left")
    v <- v[v != ""]
    v <- stringr::str_replace_all(string = v, pattern = "\\s+", replacement = "\t")
    readr::read_tsv(v)
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
