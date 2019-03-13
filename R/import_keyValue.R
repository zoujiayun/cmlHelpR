#' Accessory function to import key-value pair file
#'
#' Imports the key-value pair file for use with write_fasta()
#' @param kv_file_path Path to the key-value pair file
#' @keywords helper
#' @export
#' @examples
#' import_keyValue(kv_file_path = "path/to/kv.tsv")
import_keyValue <- function(kv_file_path) {

  ## Reading in key-value pairs
  df.key_val <- readr::read_tsv(file = kv_file_path, col_names = TRUE)

  ##  trimming column names for phylip spec - NAMES NEED TO BE THE SAME AS TREE!
  colnames(df.key_val) <- sub("^.*?([A-Z])", "\\1", names(df.key_val))
  return(df.key_val)
}
