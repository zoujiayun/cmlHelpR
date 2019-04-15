#' Parse lnL and np data from codeml output files
#'
#' Internal function for the lrt_calculation() function.
#' @param vec.list list of vectors from readr::read_lines() command (corresponds to single codeml.log file)
#' @keywords internal
.parse_np_lnL <- function(vec.list){

  lst <- purrr::map(vec.list, ~{

    r <- .x[grep(pattern = "^lnL\\(", x = .x)]
    np <- as.numeric(trimws(sub("lnL\\(ntime:.+np:(.*)\\).+", "\\1", r)))
    lnL <- as.numeric(trimws(sub(".*\\):(.*)\\s+.+$", "\\1", r)))
    out <- tibble::tibble(np = np, lnL = lnL)
    return(out)
  })

  dplyr::bind_rows(lst, .id = "gene_tree")

}
