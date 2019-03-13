#' Accessory function to read in fasta files to R
#'
#' Imports fasta files to R using the seqinR package.
#' @param cds_file_path Path to the fasta files directory
#' @param ext_cds Extension of fasta files
#' @param custom_names Names for the fasta list elements. NEED to be in the same order as the files are read in.
#' @keywords helper
#' @export
#' @examples
#' import_fasta(cds_file_path = "/path/to/directory", ext_cds = ".phy", custom_names = c("sample_1", "sample_2", "sample_3"))
import_fasta <- function(cds_file_path, ext_cds, custom_names){
  ## Fasta files
  vec.fileList <- list.files(path = cds_file_path, pattern = ext_cds, full.names = TRUE)
  lst.fasta <- purrr::map(vec.fileList,
                          seqinr::read.fasta, seqtype = "DNA",
                          as.string = FALSE,
                          strip.desc = TRUE,
                          set.attributes = FALSE,
                          forceDNAtolower = FALSE)
  names(lst.fasta) <- custom_names

  return(lst.fasta)

}
