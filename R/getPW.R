#' Parse the `p` and `w` fields from codeml output files
#'
#' Extracts the p and w lines from codeml output and stores it in a list/data frame object
#' @param dir_path Path to codeml_parallel output directory.
#' @param models Vector of models to extract p/w information from
#' @param ext Extension of output files
#' @keywords parse
#' @export
#' @examples
#' parse_PW(parse_PW = "path/to/codeml_out")
getPW <- function(dir_path, models, ext=".out") {

  ## Models that have p/w rows in some form
  m <- c("CmC", "CmD", "M1a", "M2a", "M2a_Rel", "M3", "M8", "ModelA", "ModelA1")

  ## Groupings based on structure of output
  g1 <- c("CmC", "CmD")
  g2 <- c("M1a", "M2a", "M2a_Rel", "M3")
  g3 <- c("M8")
  # g4 <- c("ModelA", "ModelA1")

  ## Check if all provided models match list above
  if(!all(is.element(models, m))){
    stop(call. = TRUE, "One of the specified models will not contain p/w values that can be parsed")
  }

  ## Listing directories
  d <- list.dirs(path = dir_path, full.names = TRUE, recursive = TRUE)
  d <- d[stringr::str_detect(string = d, pattern = paste0("/", models, "$", collapse = "|"))]
  d <- magrittr::set_names(x = d, value = stringr::str_remove(string = d, pattern = ".*/"))

  ## Get list of output files
  f <- purrr::map(d, list.files, pattern = ext, full.names = TRUE)
  f <- purrr::map(f, ~{
    n <- stringr::str_remove(string = sub(".*/(.*)/(.*)$", "\\1_\\2", .x), pattern = ext)
    names(.x) <- n
    purrr::map(.x, readr::read_lines)
  })

  ## Iterate through models names
  lst <- lapply(names(f), function(y){

    ## Model determines how file is parsed
    if(any(stringr::str_detect(string = y, pattern = g1))){

      ## Iterate through samples
      ret <- purrr::map(f[[y]], ~{

        ## Column values
        v <- .x[grep(pattern = "^site class", x = .x)]
        v <- stringr::str_replace(string = v, pattern = "\\s+", replacement = "-")
        v <- unlist(stringr::str_split(string = v, pattern = "\\s+"))

        prop <- .x[grep(pattern = "^proportion", x = .x)]
        prop <- unlist(stringr::str_split(string = prop, pattern = "\\s+"))

        bt0 <- .x[grep(pattern = "^branch type 0", x = .x)]
        bt0 <- stringr::str_replace(string = bt0, pattern = "\\s+", replacement = "-")
        bt0 <- stringr::str_replace(string = bt0, pattern = "\\s+", replacement = "-")
        bt0 <- unlist(stringr::str_split(string = bt0, pattern = "\\s+"))

        bt1 <- .x[grep(pattern = "^branch type 1", x = .x)]
        bt1 <- stringr::str_replace(string = bt1, pattern = "\\s+", replacement = "-")
        bt1 <- stringr::str_replace(string = bt1, pattern = "\\s+", replacement = "-")
        bt1 <- unlist(stringr::str_split(string = bt1, pattern = "\\s+"))

        ## Converting to data frame
        df <- tibble::tibble(
          !! prop[1] := prop[2:length(prop)],
          !! bt0[1] := bt0[2:length(bt0)],
          !! bt1[1] := bt1[2:length(bt1)]
        )

        ## Transposing to have correct orientation
        df <- t(df)
        df <- tibble::as_tibble(x = df, rownames = "temp")
        colnames(df) <- v
        return(df)

      })

      return(ret)

    } else if(any(stringr::str_detect(string = y, pattern = g2))){

      ret <- purrr::map(f[[y]], ~{

        p <- .x[grep(pattern = "^p:", x = .x)]
        p <- unlist(stringr::str_split(string = p, pattern = "\\s+"))
        p <- stringr::str_remove(string = p, pattern = ":.*")

        w <- .x[grep(pattern = "^w:", x = .x)]
        w <- unlist(stringr::str_split(string = w, pattern = "\\s+"))
        w <- stringr::str_remove(string = w, pattern = ":.*")

        df <- tibble::tibble(
          var = c(p[1], w[1]),
          K1 = c(p[2], w[2]),
          K2 = c(p[3], w[3])
        )
        return(df)

      })
      return(ret)

    } else if(any(stringr::str_detect(string = y, pattern = g3))){

      ret <- purrr::map(f[[y]],~{

        p <- .x[grep(pattern = "^p:", x = .x)]
        p <- unlist(stringr::str_split(string = p, pattern = "\\s+"))
        p <- stringr::str_remove(string = p, pattern = ":.*")

        w <- .x[grep(pattern = "^w:", x = .x)]
        w <- unlist(stringr::str_split(string = w, pattern = "\\s+"))
        w <- stringr::str_remove(string = w, pattern = ":.*")

        cn <- paste0("K", seq(1:11))

        df <- tibble::tibble(
          !! p[1] := p[2:length(p)],
          !! w[1] := w[2:length(w)]
        )

        df <- t(df)
        colnames(df) <- cn
        df <- tibble::as_tibble(x = df, rownames = "var")
        return(df)
      })

      return(ret)

    } else {

      ret <- purrr::map(f[[y]], ~{

        ## Column values
        v <- .x[grep(pattern = "^site class", x = .x)]
        v <- stringr::str_replace(string = v, pattern = "\\s+", replacement = "-")
        v <- unlist(stringr::str_split(string = v, pattern = "\\s+"))

        ## dataframe
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
        return(df)
      })

      return(ret)

    }

  })

  ## Assigning names back to list object
  names(lst) <- names(f)

  ## Formatting output
  lst <- purrr::map(lst, dplyr::bind_rows, .id = "model_gene_tree")

  ## Need to do this if M2a_Rel is included - underscore stuffs up seperators
  nest <- purrr::map(names(lst), ~{
    if(.x == "M2a_Rel"){

      t <- tidyr::separate(data = lst[[.x]], col = "model_gene_tree", into = c("model", "gene", "tree"), sep = "_", extra = "merge")
      t <- tidyr::unite(t, col = model, model, gene)
      t <- tidyr::separate(data = t, col = "tree", into = c("gene", "tree"), sep = "_")
      return(t)

    } else {
      tidyr::separate(data = lst[[.x]], col = "model_gene_tree", into = c("model", "gene", "tree"), sep = "_")
    }
  })
  names(nest) <- names(f)

  ## Group and nest for easy access
  nest <- purrr::map(nest, dplyr::group_by, model, gene, tree)
  nest <- purrr::map(nest, tidyr::nest)

  out <- list(long_list = lst, nested_list = nest)
  return(out)

}
