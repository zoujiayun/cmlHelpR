#' Parse dN/dS branch values from output files
#'
#' Given the path to the parallel_codeml output directory, this function will parse the branch dN/dS values from the
#' CODEML output. It returns a list and nested dataframe object of dN/dS values for each file.
#' @param dir_path Path to parallel_codeml output directory
#' @param models Models to get branch dN/dS values for
#' @param ext Extension of CODEML output files
#' @keywords helper
#' @importFrom rlang .data
#' @export
#' @examples
#' getBranchdNdS(dir_path = list.files(path = "path/to/codeml/outDir", models = c("M1a", "M8")))
getBranchDNDS <- function(dir_path, models, ext = ".out") {

  ## Check: do the models provided match all contain branch dN/dS values
  g <- c("FreeRatio" ,"M0" ,"M1a" ,"M2a_Rel" ,"M3" ,"M7" ,"M8" , "TwoRatio")
  if(!all(is.element(models, g))){
    stop(call. = TRUE, "One of the specified models will not contain branch dN/dS values")
  }
  m <- paste0("/", models, collapse = "|")

  ## Listing all directories in directory path
  d <- normalizePath(path = list.dirs(path = dir_path, full.names = TRUE, recursive = TRUE))
  d <- d[stringr::str_detect(string = d, pattern = m)]

  ## Listing files + reading them in
  f <- purrr::map(d, list.files, pattern = ext, full.names = TRUE)
  f <- unlist(x = f)
  f <- magrittr::set_names(x = f, value = stringr::str_remove(string = sub(".*/(.*)/(.*)$", "\\1::\\2", f), pattern = ext))
  f <- purrr::map(f, readr::read_lines)

  ## Iterate through each file and get values
  d <- purrr::map(f, ~{

    subString <- .x[grep(pattern = "^ branch", x = .x):length(.x)]

    end <- stringr::str_locate(string = subString, pattern = "")[,1]
    end <- which(is.na(end))[2]

    out <- trimws(x = subString[1:end], which = "left")
    out <- out[out != ""]
    out <- stringr::str_replace_all(string = out, pattern = "\\s+", replacement = "\t")
    readr::read_tsv(out)
  })

  ## Build dataframe - nested
  df <- dplyr::bind_rows(d, .id = "id")
  df <- tidyr::separate(data = df, col = .data$id, into = c("model", "gene_tree"), sep = "::")
  df <- tidyr::separate(data = df, col = .data$gene_tree, into = c("gene", "tree"), sep = "_")
  df <- dplyr::mutate(.data = df, tree = dplyr::if_else(is.na(.data$tree), "base", .data$tree))
  df <- dplyr::group_by(.data = df, .data$model, .data$gene, .data$tree)
  df <- tidyr::nest(data = df, .key = "branch")

  ## Return list + dataframe
  lst <- list(list = d, dataframe = df)
  return(lst)

}
