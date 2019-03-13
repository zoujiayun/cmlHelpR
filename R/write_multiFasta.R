#' Writes a multi-fasta file by gene
#'
#' Writes a multi-fasta by gene, where each sequence in a file corresponds to the gene sequence from a different sample.
#' @param df_kevValue Key-value data frame
#' @param lst_fastas List of fasta files
#' @param out_dir Output directory path
#' @keywords write fasta, by-gene
#' @export
#' @examples
#' write_multiFasta(df_kevValue = df.kv, lst_fastas = lst.fa, out_dir = "path/to/output/dir")
write_multiFasta <- function(df_kevValue, lst_fastas, out_dir) {

  ## Iterating through each row of a dataframe
  out <- apply(X = df_kevValue, MARGIN = 1, FUN = function(current_row){

    ## Vector of headers to extract from list objects
    lst <- current_row[1:3]

    ## Each row is a character vector - index to get header inforamtion
    lst.out <- purrr::map2(.x = lst, .y = lst_fastas, .f = ~{

      .y[[.x]]

    })

    ## Didn't solve the issue - needed to implement at the MSA level - but like the code
    ## Checking if longest sequence is a multiple of 3 - if not, append empty sequence to make it so
    # num.val <- length(lst.out[which.max(lengths(lst.out))][[1]])/3
    #
    # ## if the number is NOT divisable by 3
    # if ((3*num.val)%%3 != 0) {
    #
    #   ## Number of bases to add to fasta
    #   num.val <- ceiling(num.val/3) * 3 - num.val
    #   chr.name <- names(lst.out[which.max(lengths(lst.out))])
    #
    #   ## Appending n bases to vector string
    #   lst.out[[chr.name]] <- c(lst.out[[chr.name]], rep("N", num.val))
    #
    # }

    ## Building file names
    dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)

    chr.filename <- paste0(current_row[4], "-", current_row[5], "-", current_row[6],"_dnds",".fa")
    chr.fileout <- paste(out_dir, chr.filename, sep = "/")
    chr.fileout

    # print(chr.fileout)
    seqinr::write.fasta(sequences = lst.out, file.out = chr.fileout, names = names(lst.out))

  })

}
