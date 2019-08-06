#' Writes a multi-fasta file by gene for n-species
#'
#' This function creates fasta files for each gene sequence from the orthologue table (`getOrthologueHeaders``).
#' @param orthologs Output object from `getOrthologueHeaders()` function
#' @param fasta_dir Directory where fasta files (peptide and nucleotide) are located for each sample in table
#' @param pep_ext Extension of peptide files
#' @param nuc_ext Extension of nucleotide files
#' @param pep_out Path to output directory for peptide files
#' @param nuc_out Path to output directory for nucleotide files
#' @param stop_codon Character that is used to specify what a stop codon is in your sequence. Defaul = "."
#' @keywords write fasta, by-gene
#' @export
#' @examples
#' write_fasta(orthologs = listObj, fasta_dir = "path/to/fasta/directory", nuc_ext = ".fna", pep_ext = ".pep", nuc_out = "path/to/nuc_out/dir", pep_out = "/path/to/pep_out/dir", stop_codon = "*")
writeFasta <- function(orthologs, fasta_dir, pep_ext, nuc_ext, pep_out, nuc_out, stop_codon = ".", write_file = TRUE){

  ## Create output directories
  print("Creating output dirs")
  dir.create(path = pep_out, recursive = TRUE)
  dir.create(path = nuc_out, recursive = TRUE)

  ## Import peptide files
  print("Importing peptide sequence")
  pep_seq <- list.files(path = fasta_dir, pattern = pep_ext, full.names = TRUE)
  pep_seq <- magrittr::set_names(x = pep_seq, stringr::str_remove(string = basename(pep_seq), pattern = pep_ext))
  pep_seq <- purrr::map(.x = pep_seq, .f = Biostrings::readAAStringSet, format = "fasta")

  ## Importing nucleotide files
  print("Importing nucleotide sequence")
  nuc_seq <- list.files(path = fasta_dir, pattern = nuc_ext, full.names = TRUE)
  nuc_seq <- magrittr::set_names(x = nuc_seq, stringr::str_remove(string = basename(nuc_seq), pattern = nuc_ext))
  nuc_seq <- purrr::map(.x = nuc_seq, .f = Biostrings::readAAStringSet, format = "fasta")

  ## Iterate through table and write mutli-fasta files
  genes_long <- orthologs[[2]]

  ##
  print("Writing fasta files")
  out <- purrr::pmap(.l = genes_long, .f = function(file_name, data){

    pep_temp <- purrr::pmap(.l = data, .f = function(sample, header) {
      ## Extract sequence - pep_seq
      p <- pep_seq[[sample]][header] ## Single brackets to retain sequence name field
      names(p) <- sample
      return(p)
    })

    ## Assigning names and returning biostrings set
    pep_temp <- magrittr::set_names(x = pep_temp, nm = data[["file_name"]])
    pep_temp <- Biostrings::AAStringSetList(pep_temp)@unlistData

    ## Any stop codons?
    stop_df <- tibble::as_tibble(Biostrings::vmatchPattern(pattern = stop_codon, subject = pep_temp))
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
      nuc_temp <- magrittr::set_names(x = nuc_temp, nm = data[["file_name"]])
      nuc_temp <- Biostrings::AAStringSetList(nuc_temp)@unlistData

      ## Preparing return object
      lst <- list(peptide = pep_temp, nucleotide = nuc_temp, internalStop = FALSE)

      ## If write_file == TRUE - go ahead and write the file
      if(isTRUE(write_file)){

        ## Check if file already exists - delete if it does
        if(file.exists(paste0(pep_out, "/", file_name, ".fasta"))) file.remove(paste0(pep_out, "/", file_name, ".fasta"))
        if(file.exists(paste0(nuc_out, "/", file_name, ".fasta"))) file.remove(paste0(nuc_out, "/", file_name, ".fasta"))

        ## Create file
        file.create(paste0(pep_out, "/", file_name, ".fasta"))
        file.create(paste0(nuc_out, "/", file_name, ".fasta"))

        ## Write file
        Biostrings::writeXStringSet(x = pep_temp, filepath = paste0(pep_out, "/", file_name, ".fasta"), append = TRUE)
        Biostrings::writeXStringSet(x = nuc_temp, filepath = paste0(nuc_out, "/", file_name, ".fasta"), append = TRUE)

        ## Returning summary of sequence
        return(lst)

      } else{
        return(lst)
      }

    } else {

      lst <- list(peptide = pep_temp, internalStop = TRUE)
      return(lst)

    }


  })

  ## Naming output list
  names(out) <- genes_long[["file_name"]]
  return(out)

}







