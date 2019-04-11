#' Parse lnL and np data from codeml output files
#'
#' Internal function for the lrt_calculation() function.
#' @param vec.list list of vectors from readr::read_lines() command (corresponds to single codeml.log file)
#' @keywords internal
.parse_np_lnL <- function(vec.list){

  lst <- purrr::map(vec.list, ~{
    ## Number of processes
    int.np <- .x[grepl("^np\\s+=", .x)]
    int.np <- gsub("np =\\s+", "", int.np)
    int.np <- as.numeric(x = int.np)

    ## log-likelihood
    int.lnL <- .x[grepl("^lnL  =", .x)]
    int.lnL <- gsub("lnL +=\\s+", "", int.lnL)
    int.lnL <- as.numeric(int.lnL)

    df.out <- tibble::tibble(np = int.np, lnL = int.lnL)
    return(df.out)
  })

  dplyr::bind_rows(lst, .id = "gene_tree")

}
