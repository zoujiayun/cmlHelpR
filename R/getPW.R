#' Parse the `p` and `w` fields from codeml output files
#'
#' Extracts the p and w lines from codeml output and stores it in a list/data frame object
#' @param dir_path Path to codeml_parallel output directory.
#' @param df_out Set at TRUE to output a long form data frame. Default is a list object; one element per gene-tree combination.
#' @keywords parse
#' @export
#' @examples
#' parse_PW(parse_PW = "path/to/codeml_out", df_out = TRUE)
getPW <- function(file_list) {

  ## Accepted model
  m <- "ModelA"
  m <- paste0("/", m, "/", collapse = "|")

  ## Filter out models not in g
  f <- file_list[stringr::str_detect(string = file_list, pattern = m)]
  f <- magrittr::set_names(x = f, value = stringr::str_remove(string = sub(".*/(.*)/(.*)$", "\\1_\\2", f), pattern = ".out"))
  f <- purrr::map(f, readr::read_lines)

  ## Pattern matching rows
  lst <- purrr::map(f, ~{

    v <- .x[grep(pattern = "^site class", x = .x)]
    v <- stringr::str_replace(string = v, pattern = "\\s+", replacement = "-")
    v <- unlist(stringr::str_split(string = v, pattern = "\\s+"))

    if (!is.null(v)) {

      prop <- .x[grep(pattern = "^proportion", x = .x)]
      prop <- unlist(stringr::str_split(string = prop, pattern = "\\s+"))

      bkg <- .x[grep(pattern = "^background", x = .x)]
      bkg <- stringr::str_replace(string = bkg, pattern = "\\s+", replacement = "-")
      bkg <- unlist(stringr::str_split(string = bkg, pattern = "\\s+"))

      fg <- .x[grep(pattern = "^foreground", x = .x)]
      fg <- stringr::str_replace(string = fg, pattern = "\\s+", replacement = "-")
      fg <- unlist(stringr::str_split(string = fg, pattern = "\\s+"))

      ## Converting to data frame
      df <- tibble::tibble(
        !! prop[1] := prop[2:length(prop)],
        !! bkg[1] := bkg[2:length(bkg)],
        !! fg[1] := fg[2:length(fg)]
      )

      ## Transposing to have correct orientation
      df <- t(df)
      df <- tibble::as_tibble(x = df, rownames = "temp")
      colnames(df) <- v
      df

    } else {

      df <- tibble::tibble(`site-class` = NA,
                       `0` = NA,
                       `1` = NA,
                       `2a` = NA,
                       `2b` = NA)
    }

  })

  ## Binding rows
  d <- dplyr::bind_rows(lst, .id = "temp")
  d <- tidyr::separate(data = d, col = `temp`, into = c("model", "gene", "tree"), sep = "_")
  d <- dplyr::group_by(.data = d, model, gene, tree)
  d <- tidyr::nest(data = d, .key = "values")

  lst <- list(list = lst, dataframe = d)
  return(lst)

}
