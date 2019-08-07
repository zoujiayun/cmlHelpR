#' Returns a key-value table of the gene symbols (filename) and fasta headers corresponding to those genes
#'
#' Given a path to a directory containing conditional reciprocal best blast (CRBB) files, this function will
#' parse the data and find sequences that are considered a CRBB across multiple comparisons. This function
#' relies on one sample being the reference for the remaining comparisons to be joined to.
#'
#' CRBB files need to have the filenaming structure as follows; sampleA_sampleB.ext. Both sample identifiers need
#' to be present in the file name, separated by some kind of delimeter. The delimiter does not matter, as you can
#' specify this at the function call.
#'
#' @param crbb_path Path to the directory containing CRBB files
#' @param crbb_ext Extension of the CRBB files
#' @param sample_separator Delimiter between sample IDs in the filename
#' @param id_pct Percentage identity between hits
#' @param aln_pct Proportion of query/target that needs to be accounted for in BLAST alignment portion
#' @keywords helper
#' @importFrom rlang .data
#' @export
getOrthologHeaders <- function(crbb_path, crbb_ext, sample_separator, id_pct, aln_pct) {

  ## Column names of CRBB output
  colNames <- c("query", "target", "id", "alnlen", "evalue",
                "bitscore", "qstart", "tstart", "qlen", "tlen")

  ## Importing files in CRBB output directory
  crbb_files <- list.files(path = crbb_path, pattern = crbb_ext, full.names = TRUE)
  crbb_files <- magrittr::set_names(x = crbb_files, value = stringr::str_remove_all(string = basename(crbb_files), pattern = crbb_ext))
  crbb_list <- purrr::map(crbb_files, readr::read_tsv, col_names = colNames)

  ## Selecting BEST Reciprocal best hit: arrange by descending bitscore then %id - take top value of group (by gene) == Overall most likely ortholog candidate
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
                        prop_tlen = alnlen/tlen * 100)               ## Getting alignment proportion values
    df <- dplyr::select(.data = df, 1, 2, 11, 12, 3, 5, 6, 13, 14)   ## Organising to be cleaner
    df <- dplyr::group_by_at(.tbl = df, g_q)
    df <- dplyr::arrange(.data = df, df[[g_q]],
                         dplyr::desc(bitscore),
                         dplyr::desc(id))                            ## Arranging grouped data by bitscore + %id = best hit
    df <- dplyr::slice(.data = df, 1)                                ## Selecting best hit for group
    df <- dplyr::ungroup(x = df)
    df <- dplyr::group_by_at(.tbl = df, t) ## Arranging by target and slicing first hit incase duplicate.
    df <- dplyr::arrange(.data = df, df[[t]],
                         dplyr::desc(bitscore),
                         dplyr::desc(id))
    df <- dplyr::slice(.data = df, 1)                                ## Same as above but for the query (just in case)
    df <- dplyr::ungroup(x = df)
    df <- dplyr::filter(.data = df, id >= id_pct)
    df <- dplyr::filter(.data = df, prop_qlen >= aln_pct & prop_tlen >= aln_pct)   ## Filtering on thresholds
    df <- dplyr::select(df, q, t, g_q, g_t)
  })

  ## Matching between the different comparisons
  orth_df <- purrr::reduce(.x = bestHits_list, dplyr::full_join)   ## Reduce command applies left join to all dataframes in the list
  orth_df <- tidyr::drop_na(data = orth_df)                        ## Removing rows where information isn't present for ALL samples
  orth_df <- tidyr::unite(data = orth_df,                          ## Uniting gene symbol columns to generate file name column
                        col = "file_name",
                        dplyr::starts_with(match = "gene_"), sep = "-")
  orth_df <- dplyr::select(.data = orth_df, .data$file_name, dplyr::everything()) ## Ordering the output

  ## Building list output
  nested <- tidyr::gather(data = orth_df, key = "sample", value = "header", 2:ncol(orth_df))
  nested <- dplyr::group_by(.data = nested, .data$file_name)
  nested <- tidyr::nest(data = nested)

  out <- list(wide_format = orth_df, nested_format = nested)
  return(out)
}
