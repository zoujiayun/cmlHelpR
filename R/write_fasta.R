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
write_fasta <- function(orthologs, fasta_dir, fasta_ext, out_dir){

  dir.create(path = out_dir, recursive = TRUE)

  ## Import fasta files
  fastas <- list.files(path = fasta_dir, pattern = fasta_ext, full.names = TRUE)
  fastas <- magrittr::set_names(x = fastas, stringr::str_remove(string = basename(fastas), pattern = fasta_ext))
  fastas <- purrr::map(.x = fastas, .f = Biostrings::readAAStringSet, format = "fasta")

  ## Iterate through table and write mutli-fasta files
  df <- orthologs[[2]]

  out <- purrr::pmap(.l = df, .f = function(file_name, data){

    ## Output file preparation
    print(paste0(out_dir, "/", file_name, ".fasta"))
    if(file.exists(paste0(out_dir, "/", file_name, ".fasta"))) file.remove(paste0(out_dir, "/", file_name, ".fasta"))
    file.create(paste0(out_dir, "/", file_name, ".fasta"))

    ## Iterate over data-frame lists for genes
    tmp <- purrr::pmap(.l = data, .f = function(sample, header) {

      ## Extract sequence
      s <- fastas[[sample]][header] ## Single brackets to retain sequence name field
      names(s) <- sample

      ## Writing to file
      Biostrings::writeXStringSet(x = s, filepath = paste0(out_dir, "/", file_name, ".fasta"), append = TRUE)
      s ## Returning object

    })

    ## Assigning names and returning biostrings set
    tmp <- set_names(x = tmp, nm = data[["file_name"]])
    tmp <- Biostrings::AAStringSetList(tmp)@unlistData

  })

  ## Naming output list
  names(out) <- df[["file_name"]]
  return(out)

}


