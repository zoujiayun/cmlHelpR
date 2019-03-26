#' Trim phylip files to make them a multiple of three
#'
#' Checks if the sequence in a phylip file is a multiple of three. If it is not, it trims the minimum number of bases needed to make it so.
#' @param infile File path to phylip file
#' @param headers A vector of the sequence headers
#' @param file_dir An output file directory (doesn't need to exist)
#' @keywords trim, phylip, clean, seqret
#' @export
#' @examples
#' clean_seqret(infile = "/path/to/phylip.phy", headers = c("sample1", "sample2", "sample3"), file_dir = "/path/to/output_directory")
clean_seqret <- function (infile, headers, file_dir = NULL) {

  ## Getting file basename
  chr.basename <- gsub(".+/(.*).phy", "\\1", infile)

  ## Reading header
  dat <- readr::read_lines(file = infile, skip = 1)
  dat <- trimws(x = dat, which = "left")
  dat <- paste(dat, collapse = " ")
  dat <- gsub(" ","", dat)

  ## Splitting on header information
  dat <- stringr::str_split(string = dat, pattern = paste(headers, collapse = "|"))
  dat <- dat[[1]][-1] ## Makes an empty list object at the beginning always

  ## Assigning sequences names
  names(dat) <- headers

  ## Getting length + making multiple of three
  mod <- nchar(dat[[1]]) %% 3

  ## Removing n bases from end of string (where n = mod)
  dat <- purrr::map(dat, ~{
    substr(x = .x, start = 1, stop =  nchar(.x) - mod)
  })

  ## Generating NEW meta information for phylip file
  meta <- paste(length(headers),nchar(dat[[1]]), sep = " ")

  ## Adding sample headers to phylip sequences
  dat <- purrr::map2(.x = dat, .y = headers, ~{
    paste(.y, .x, sep = "  ")
  })

  ## Preappending meta inforamtion to list object
  dat <- append(x = meta, values = dat, after = 1)

  ## Writing to file
  if (!is.null(file_dir)) {

    dir.create(path = file_dir, showWarnings = FALSE, recursive = TRUE)

    lapply(dat, write, paste0(file_dir,"/",chr.basename,".phy"), append = TRUE)
  }

  ## Returning list
  return(dat)

}

