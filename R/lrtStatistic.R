#' Description
#'
#' General run down of what it does
#' @param dir_path Path to codeml_parallel output directory.
#' @param lst_comparisons Combinations/comparisons of the above models
#' @param ext Extension of CODEML output files. Defaults to ".out"
#' @export
#' @examples
#' lrt_statistic(dir_path = "path/to/codeml_out", models = c("M1a", "M2a"), lst_comparisons = c("M1a", "M2a"))
lrtStatistic <- function(dir_path, lst_comparisons, ext = ".out") {

  ## Models
  m <- unlist(lst_comparisons)
  m <- unique(m)

  ## Predefined groupings
  nsm <- c("M0", "M1a", "M2a", "M3", "M7", "M8")
  bsm <- c("FreeRatio", "TwoRatio", "ModelA", "ModelA1", "M2a_Rel", "CmC", "CmD")

  ## Selecting output dirs for analysis
  d <- list.dirs(path = dir_path, full.names = TRUE)
  d <- d[stringr::str_detect(string = d, pattern = paste(m, "$", sep = "", collapse = "|"))]

  ## Naming directory vectors
  n <- sub(pattern = ".*/", "", d)
  d <- magrittr::set_names(x = d, value = n)

  ## Reading codeml output files
  print("Reading output")
  f <- purrr::map(d, ~{
    fl <- list.files(path = .x, pattern = ext, full.names = TRUE)
    fl <- magrittr::set_names(x = fl, value = sub(ext, "", basename(fl)))
    fl <- purrr::map(fl, readr::read_lines)
  })

  ## Extracting np + lnL values
  print("Extracting np + lnL")
  np_lnL <- purrr::map(f, .parse_np_lnL)
  np_lnL <- dplyr::bind_rows(np_lnL, .id = "model")

  ## Calculating LRT statistic
  out_lrt <- purrr::map(lst_comparisons, ~{

    if(all(is.element(.x, nsm))) {       ## Both models are null/site models

      df <- .lrt_nsm(df = np_lnL, mdl = .x)
      return(df)

    } else if(all(is.element(.x, bsm))){ ## Both models are branch/branch-site/clade models

      df <- .lrt_bcm(df = np_lnL, mdl = .x)
      return(df)

    } else {                             ## One model from each camp

      df <- .lrt_hnh(df = np_lnL, mdl = .x)
      return(df)

    }

  })

  ## Assigning names to list object
  names(out_lrt) <- purrr::map(lst_comparisons, ~{paste(.x, collapse = "-")})
  return(out_lrt)

}
