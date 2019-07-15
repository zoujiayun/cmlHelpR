#' Return a wide and nested form othologue table
#'
#' This function is designed to operate on THREE samples only! It will fail with more!
#' The function reads in Conditional Reciprocal Best Blast (CRBB) tables and relables
#' the data using the file name. For each comparison table, a best hit is selected (if
#' multiple hits exist for each query) using the percentage-identity and bitscore.
#'
#' This function requires the filenames contain the following; Samples in the comparison
#' are in the filename separated by a delimiter (specified in the function). That the gene
#' symbol is at the end of the fasta header separated by an underscore e.g. '>someIdentifier_SYMBOL'
#' @param crbb_path Path to the directory containing CRBB files
#' @param crbb_ext Extension of the CRBB files
#' @param sample_separator Delimiter between sample IDs in the filename
#' @param id_pct Percentage identity between BLAST hits
#' @param aln_pct Proportion of query/target that needs to be accounted for in BLAST alignment
#' @keywords helper
#' @export
#' @examples
#' getOrthologHeaders(crbb_path = "path/to/directory", crbb_ext = ".tsv", sample_separator = "__", id_pct = 95, aln_pct = 90)
getOrthologHeaders <- function(crbb_path, crbb_ext, sample_separator, id_pct, aln_pct) {

  ## Column names of CRBB output
  colNames <- c("query", "target", "id", "alnlen", "evalue",
                "bitscore", "qstart", "tstart", "qlen", "tlen")

  ## Importing files in CRBB output directory
  crbb_files <- list.files(path = crbb_path, pattern = crbb_ext, full.names = TRUE)
  crbb_files <- magrittr::set_names(x = crbb_files, value = stringr::str_remove_all(string = basename(crbb_files), pattern = crbb_ext))
  crbb_list <- purrr::map(crbb_files, readr::read_tsv, col_names = colNames)

  ## Selecting BEST Reciprocal best hit: arrange by %id --> evalue --> take top value of group (by gene) == Overall most likely ortholog candidate
  bestHits_list <- purrr::map(names(crbb_list), ~{

    q <- sub(paste0(sample_separator, ".*") , "", .x)   ## Getting query_sequence name
    t <- sub(paste0(".*", sample_separator) , "", .x)   ## Gettign target sequence name

    g_q <- paste("gene", q, sep = "_")  ## Column names
    g_t <- paste("gene", t, sep = "_")  ## Column names

    df <- dplyr::rename(.data = crbb_list[[.x]], !! q := query, !! t := target)
    df <- dplyr::mutate(.data = df,
                        !! g_q := sub(".+_(.*)", "\\1", df[[q]]),
                        !! g_t := sub(".+_(.*)", "\\1", df[[t]]),
                        prop_qlen = alnlen/qlen * 100,
                        prop_tlen = alnlen/tlen * 100)
    df <- dplyr::select(.data = df, 1, 2, 11, 12, 3, 5, 6, 13, 14)
    df <- dplyr::arrange(.data = df, df[[g_q]],
                         dplyr::desc(bitscore),
                         dplyr::desc(id))
    df <- dplyr::group_by_at(.tbl = df, g_q)
    df <- dplyr::slice(.data = df, 1)
    df <- dplyr::ungroup(x = df)
    df <- dplyr::group_by_at(.tbl = df, t) ## Arranging by target and slicing first hit incase duplicate.
    df <- dplyr::arrange(.data = df, df[[t]],
                         dplyr::desc(bitscore),
                         dplyr::desc(id))
    df <- dplyr::slice(.data = df, 1)
    df <- dplyr::ungroup(x = df)
    df <- dplyr::filter(.data = df, id >= id_pct)
    df <- dplyr::filter(.data = df, prop_qlen >= aln_pct & prop_tlen >= aln_pct)
    df <- dplyr::select(df, q, t, g_q, g_t)
  })

  ## Matching between the different comparisons
  orth_df <- purrr::reduce(.x = bestHits_list, dplyr::full_join)   ## Reduce command applies left join to all dataframes in the list
  orth_df <- tidyr::drop_na(data = orth_df)                        ## Removing rows where information isn't present for ALL samples
  orth_df <- tidyr::unite(data = orth_df,                          ## Uniting gene symbol columns to generate file name column
                        col = "file_name",
                        dplyr::starts_with(match = "gene_"), sep = "-")
  orth_df <- dplyr::select(.data = orth_df, file_name, dplyr::everything()) ## Ordering the output

  ## Building list output
  nested <- tidyr::gather(data = orth_df, key = "sample", value = "header", 2:ncol(orth_df))
  nested <- dplyr::group_by(.data = nested, file_name)
  nested <- tidyr::nest(data = nested)

  out <- list(wide_format = orth_df, nested_format = nested)
  return(out)
}
