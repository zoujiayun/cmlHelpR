#' Parse the `p` and `w` fields from codeml output files
#'
#' Extracts the p and w lines from codeml output and stores it in a list/data frame object
#' @param dir_path Path to codeml_parallel output directory.
#' @param df_out Set at TRUE to output a long form data frame. Default is a list object; one element per gene-tree combination.
#' @keywords parse
#' @export
#' @examples
#' parse_PW(parse_PW = "path/to/codeml_out", df_out = TRUE)

parse_PW <- function(dir_path, df_out = NULL) {

  fl <- list.files(path = dir_path, pattern = "Model2Selection_", full.names = TRUE, recursive = TRUE)
  fl <- set_names(x = fl, value = basename(str_remove(string = fl, pattern = ".output")))
  lst.fl <- map(fl, read_lines)

  ## Pattern matching rows
  lst.df <- map(lst.fl, ~{

    vec.sc <- .x[grep(pattern = "^site class", x = .x, perl = TRUE)]
    vec.sc <- str_replace(string = vec.sc, pattern = "\\s+", replacement = "-")
    vec.sc <- unlist(str_split(string = vec.sc, pattern = "\\s+"))

    if (!is.null(vec.sc)) {

      vec.prop <- .x[grep(pattern = "^proportion", x = .x, perl = TRUE)]
      vec.prop <- unlist(str_split(string = vec.prop, pattern = "\\s+"))

      vec.bkg <- .x[grep(pattern = "^background", x = .x, perl = TRUE)]
      vec.bkg <- str_replace(string = vec.bkg, pattern = "\\s+", replacement = "-")
      vec.bkg <- unlist(str_split(string = vec.bkg, pattern = "\\s+"))

      vec.fg <- .x[grep(pattern = "^foreground", x = .x, perl = TRUE)]
      vec.fg <- str_replace(string = vec.fg, pattern = "\\s+", replacement = "-")
      vec.fg <- unlist( str_split(string = vec.fg, pattern = "\\s+"))

      ## Converting to data frame
      df.out <- tibble(
        !! vec.prop[1] := vec.prop[2:length(vec.prop)],
        !! vec.bkg[1] := vec.bkg[2:length(vec.bkg)],
        !! vec.fg[1] := vec.fg[2:length(vec.fg)]
      )

      ## Transposing to have correct orientation
      df.out <- t(df.out)
      df.out <- as_tibble(x = df.out, rownames = "temp")
      colnames(df.out) <- vec.sc
      df.out

    } else {

      df.out <- tibble(`site-class` = NA,
                       `0` = NA,
                       `1` = NA,
                       `2a` = NA,
                       `2b` = NA)
      df.out
    }

  })

  if (!is_null(df_out)) {
    return(bind_rows(lst.df, .id = "condition"))
  } else {

    return(lst.df)
  }


}
