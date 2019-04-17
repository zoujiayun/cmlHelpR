#' Identify gene alignments that have more than one stop codon in the sequence
#'
#' Parses the codeml output for the number of columns converted to '???' in the sequence alignment for each gene
#' @param dir_path Directory path to codeml_parallel output directory.
#' @param model_string Model to import information for. This is arbitrary as all comparisons will use the same alignment.
#' @param tree_string If there are more than one tree, specify which tree to read in. Again, this is arbitrary as all gene-models will have the same alignment
#' @param nSites Integer threshold. Extract genes that have this many or greater stop codons in their alignment.
#' @keywords stop, column, converted
#' @export
#' @examples
#' internal_stops(dir_path = "path/to/codeml_out", tree_string = "brownForeground", nSites = 2)
internal_stops <- function(dir_path, tree_string, nSites = NULL) {

  fl <- list.files(path = dir_path, pattern = ".output", full.names = TRUE, recursive = TRUE)
  fl <- fl[str_detect(string = fl, pattern = tree_string)]
  fl <- magrittr::set_names(x = fl, value = basename(fl))
  lst.fl <- map(fl, read_lines)

  ## Parsing number of columns converted to '???'
  lst.nCols <- map(names(lst.fl), ~{

    s <- grep(pattern = "Before deleting alignment gaps", x = lst.fl[[.x]]) + 3

    v <- lst.fl[[.x]][s]
    c <- stringr::str_count(string = v, pattern = "\\?\\?\\?")

    tibble(gene_tree = .x,
           nCols = c)
  })

  df.nCols <-  bind_rows(lst.nCols)
  df.nCols <- separate(df.nCols, gene_tree, c("gene", "tree"), "::")

  if (!is.null(nSites)) {
    df.filt <- filter(df.nCols, nCols <= nSites)
    return(df.filt)
  } else {
    return(df.nCols)
  }

}
