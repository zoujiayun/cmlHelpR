#' Identify sites under selection that are in close proximity to other sites under selection
#'
#' This function checks the distance between sites under selection and reports those that
#' fall within a specified distance of each other. For example, if you have two sites reported as
#' under selection, but they are the two amino acids next to each other, both sites would be reported
#' if the gap_threshold was set to '1'.
#' @param beb_long Long form BEB table from getBEB() function
#' @param gap_threshold Report sites that are within this proximity to other reported sites under selection
#' @keywords helper
#' @export
#' @examples
#' isSequential(beb_long = df.beb, gap_threshold = 1)
isSequential <- function(beb_long, gap_threshold){

  ## Nesting long-form data
  nst <- dplyr::group_by(.data = beb_long, model, gene, tree)
  nst <- tidyr::nest(data = nst)

  ## New column with sites of interest
  nst <- dplyr::mutate(.data = nst,
                       within_threshold = purrr::map(data, ~{

                         pos <- .x[["pos"]]
                         pos_gap <- which(diff(pos) <= gap_threshold)
                         pos_len <- pos_len <- length(pos_gap)

                         if(pos_len > 0){ ## When there is at least one site that is too close to another
                           idx <- pos_gap + 1
                           idx <- sort(unique(c(pos_gap, idx)))
                           .x[idx, ]
                         } else {
                           return(NULL)
                         }
                       })
  )

  ## Naming output
  nm <- paste(nst[["model"]], nst[["gene"]], nst[["tree"]], sep = "_")
  names(nst[["within_threshold"]]) <- nm

  return(nst)

}
