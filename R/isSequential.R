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

  ## Sites that fall within proximity threshold (close to other sites under slection)
  nst <- dplyr::mutate(.data = nst,
                       within_threshold = purrr::map(data, ~{

                         pos <- .x[["pos"]] ## Positions
                         rn <- 1:nrow(.x) ## Row number used as position name
                         names(pos) <- rn

                         pos_gap <- which(diff(pos) <= gap_threshold) ## Sites within the proximity
                         pos_len <- length(pos_gap) ## N-sites that fail proximity threshold

                         if(pos_len > 0){ ## When there is at least one site that is too close to another
                           gap_idx <- as.numeric(names(pos_gap))
                           idx <- gap_idx - 1 ## staggered id now (matches gap between vector values) - need to subtract 1 to align ids back to actual rows
                           idx <- sort(unique(c(gap_idx, idx))) ## Joining vectors to get all sites to extract

                           o <- .x[1:nrow(.x) %in% idx, ] ## Subsetting dataframe for sites that are within proximity threshold
                           return(o)
                         } else {
                           return(NULL)
                         }
                       }))

  ## Sites outside of proximity distance (dispersed)
  nst <- dplyr::mutate(.data = nst,
                       confident_sites = purrr::map(data, ~{

                         pos <- .x[["pos"]]
                         rn <- 1:nrow(.x)
                         names(pos) <- rn

                         pos_gap <- which(diff(pos) <= gap_threshold)
                         pos_len <- length(pos_gap)

                         if(pos_len > 0){
                           gap_idx <- as.numeric(names(pos_gap))
                           idx <- gap_idx - 1
                           idx <- sort(unique(c(gap_idx, idx)))

                           o <- .x[!1:nrow(.x) %in% idx, ]
                           return(o)
                         } else {
                           .x
                         }
                       }))

  ## Naming output
  nm <- paste(nst[["model"]], nst[["gene"]], nst[["tree"]], sep = "_")
  names(nst[["within_threshold"]]) <- nm
  names(nst[["confident_sites"]]) <- nm

  return(nst)

}
