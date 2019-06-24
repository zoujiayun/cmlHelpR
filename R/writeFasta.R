#' Writes a multi-fasta file by gene
#'
#' This function creates fasta files for each gene sequence from the orthologue table (get_ortholog_headers()).
#' @param orthologs Orthologue list object from get_ortholog_headers
#' @param fasta_dir Directory where sample fasta files are located (peptide sequences)
#' @param fasta_ext Extension of fasta files
#' @param out_dir Directory where multi-fasta files should be written
#' @keywords write fasta, by-gene
#' @export
#' @examples
#' write_fasta(orthologs = listObj, fasta_dir = "path/to/fasta/directory", fasta_ext = ".fasta", out_dir = "path/to/output/dir")
writeFasta <- function(orthologs, fasta_dir, pep_ext, nuc_ext, pep_out, nuc_out){

  ## Create output directories
  dir.create(path = pep_out, recursive = TRUE)
  dir.create(path = nuc_out, recursive = TRUE)

  ## Import peptide files
  pep_seq <- list.files(path = fasta_dir, pattern = pep_ext, full.names = TRUE)
  pep_seq <- magrittr::set_names(x = pep_seq, stringr::str_remove(string = basename(pep_seq), pattern = pep_ext))
  pep_seq <- purrr::map(.x = pep_seq, .f = Biostrings::readAAStringSet, format = "fasta")

  ## Importing nucleotide files
  nuc_seq <- list.files(path = fasta_dir, pattern = nuc_ext, full.names = TRUE)
  nuc_seq <- magrittr::set_names(x = nuc_seq, stringr::str_remove(string = basename(nuc_seq), pattern = nuc_ext))
  nuc_seq <- purrr::map(.x = nuc_seq, .f = Biostrings::readAAStringSet, format = "fasta")

  ## Iterate through table and write mutli-fasta files
  genes_long <- orthologs[[2]]

  ##
  out <- purrr::pmap(.l = genes_long, .f = function(file_name, data){

    pep_temp <- purrr::pmap(.l = data, .f = function(sample, header) {
      ## Extract sequence - pep_seq
      p <- pep_seq[[sample]][header] ## Single brackets to retain sequence name field
      names(p) <- sample
      return(p)
    })

    ## Assigning names and returning biostrings set
    pep_temp <- set_names(x = pep_temp, nm = data[["file_name"]])
    pep_temp <- Biostrings::AAStringSetList(pep_temp)@unlistData

    ## Any stop codons?
    stop_df <- tibble::as_tibble(Biostrings::vmatchPattern(pattern = ".", subject = pep_temp))
    stop_df <- dim(stop_df)[[1]]

    ## If no internal stop codons
    if(stop_df == 0){

      nuc_temp <- purrr::pmap(.l = data, .f = function(sample, header) {
        ## Extracting sequence - nuc_seq
        n <- nuc_seq[[sample]][header]
        names(n) <- sample
        return(n)
      })

      ## Assigning names and returning biostrings set
      nuc_temp <- set_names(x = nuc_temp, nm = data[["file_name"]])
      nuc_temp <- Biostrings::AAStringSetList(nuc_temp)@unlistData

      ## Check if file already exists - delete if it does
      if(file.exists(paste0(pep_out, "/", file_name, ".fasta"))) file.remove(paste0(pep_out, "/", file_name, ".fasta"))
      if(file.exists(paste0(nuc_out, "/", file_name, ".fasta"))) file.remove(paste0(nuc_out, "/", file_name, ".fasta"))

      ## Create file
      file.create(paste0(pep_out, "/", file_name, ".fasta"))
      file.create(paste0(nuc_out, "/", file_name, ".fasta"))

      ## Write file
      Biostrings::writeXStringSet(x = pep_temp, filepath = paste0(pep_out, "/", file_name, ".fasta"), append = TRUE)
      Biostrings::writeXStringSet(x = nuc_temp, filepath = paste0(nuc_out, "/", file_name, ".fasta"), append = TRUE)

      ## Preparing return object
      list(peptide = pep_temp, nucleotide = nuc_temp, internalStop = FALSE)

    } else {

      list(peptide = pep_temp, internalStop = TRUE)

    }


  })

  ## Naming output list
  names(out) <- genes_long[["file_name"]]
  return(out)

}