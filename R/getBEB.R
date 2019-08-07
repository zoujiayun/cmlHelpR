#' Extract Bayes Empiracle Bayes information from codeml output
#'
#' Extracts all sites reported in the Bayes Empiracle Bayes section of the codeml output as a dataframe object
#' @param dir_path Directory path to codeml_parallel output directory
#' @param models Vector of models to run. Valid models include `M2a`, `M8` and `ModelA`
#' @param sig P-value threshold to filter BEB output by e.g. `sig = 0.05`
#' @param cores Number of cores to parallelise over
#' @param ext Extension of CODEML output files
#' @keywords Bayes, empiracle, Bayes, parse
#' @importFrom rlang .data
#' @export
getBEB <- function(dir_path, models, sig = NULL, cores = 1, ext = ".out") {

  ## Models with BEB output
  m <- c("M2a", "M8", "ModelA")

  ## Check model inputs
  if (!all(models %in% m)) {
    print(paste("Accepted models:", m, sep = " "))
    print(paste("Your models:", models, sep = " "))
    stop("One of the models above does not contain BEB values")
  }

  ## Importing all log files
  d <- list.dirs(path = dir_path, full.names = TRUE)
  d <- d[stringr::str_detect(string = d, pattern = paste(models, "$", sep = "", collapse = "|"))]

  ## Setting names of dirs vector
  n <- sub(".*/", "", d)
  d <- magrittr::set_names(x = d, value = n)

  ## Importing data
  print("Reading output files")
  f <- purrr::map(d, ~{
    fl <- list.files(path = .x, pattern = ext, full.names = TRUE)
    fl <- magrittr::set_names(x = fl, value = sub(ext, "", basename(fl)))
    purrr::map(fl, readr::read_lines)
  })

  ## Parsing BEB from models
  beb <- purrr::map(.x = names(f), ~{

    # print(paste("Parsing BEB from:", .x, sep = " "))
    if (.x == "M2a" | .x == "M8") {

      .parse_m2a_M8(lines = f[[.x]], c = cores)

    } else if (.x == "ModelA") {

      .parse_ModelA(lines = f[[.x]], c = cores)

    }

  })

  ## Single tibble
  beb <- magrittr::set_names(x = beb, value = names(f))
  beb <- dplyr::bind_rows(beb, .id = "model")
  beb <- tidyr::separate(data = beb, col = .data$id, into = c("gene", "tree"), sep = "_")
  beb <- dplyr::mutate(.data = beb, tree = stringr::str_replace_na(string = .data$tree, replacement = "base"))

  print("Getting no-gap-length value")
  alnLength <- purrr::map(f[[1]], ~{

    r <- head(x = .x, n = 3)
    r <- trimws(x = r)
    r <- r[stringr::str_detect(string = r, pattern = "^3")]
    r <- stringr::str_split(string = r, pattern = "\\s+")
    r <- unlist(x = r)
    tibble::tibble(seqLen = as.numeric(r[[2]]))

  })

  alnLength <- dplyr::bind_rows(alnLength, .id = "gene")
  beb <- dplyr::left_join(beb, alnLength)

  ## Significance filtering
  if(!is.null(sig)){
    beb <- dplyr::filter(.data = beb, .data$pval <= sig)
  }

  # ## Return a list object of long form and nested
  lst <- split(x = beb, beb[["model"]])
  nst <- dplyr::group_by(.data = beb, .data$model)
  nst <- tidyr::nest(data = nst, .key = "BEB")

  o <- list(long = beb, list = lst, nested = nst)

  return(o)
}
