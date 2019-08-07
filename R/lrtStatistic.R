#' Conduct Likelihood-Ratio-Test (LRT) between two evolutionary models
#'
#' This function takes the path to the parallel_codeml output directory and a list of comparisons as input.
#' Based on these variables, it will go and run LRTs for each comparison, generating a list of data frames
#' containing all parameters used to run the statistic, as well as the result.
#' @param dir_path Path to `parallel_codeml` output directory.
#' @param comp_list Model combinations as a list of vectors
#' @param ext Extension of CODEML output files. Defaults to ".out"
#' @export
lrtStatistic <- function(dir_path, comp_list, ext = ".out") {

  ## Models
  m <- unlist(comp_list)
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
  out_lrt <- purrr::map(comp_list, ~{

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
  names(out_lrt) <- purrr::map(comp_list, ~{paste(.x, collapse = "-")})
  return(out_lrt)

}
