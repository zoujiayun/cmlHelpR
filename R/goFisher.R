#' Significance testing of GO terms for grouped data
#'
#' Conducts a Fisher's exact test on a 2x2 contigency table of presence/absence of GO terms for grouped data.
#' @param gene_2_go Dataframe that contains two columns. Column 1 is gene ID and column 2 is GO terms associated with that gene ID.
#' @param gene_2_group Long form dataframe mapping which group a gene belongs to. Genes can appear multiple times belonging to different groups.
#' Column 1 should be the gene ID and column two is the grouping variable.
#' @param pAdj_method P-value adjustment method. Takes any value from the `p.adjust()`` function
#' @param p_cutoff P-value threshold that values must be less than or equal to
#' @keywords helper
#' @importFrom rlang .data
#' @export
goFisher <- function(gene_2_go, gene_2_group, pAdj_method, p_cutoff){

  ## Improvement: Join both inputs into one dataframe that is then parsed. Should
  ## save on a few of the left joins later on in the code.

  ## Total number of genes in analysis
  nGenes_total <- length(unique(gene_2_go[[1]]))

  ## Frequency of GO term in total set
  nGO_total <- dplyr::group_by_at(.tbl = gene_2_go, 2) ## Grouping by GO ids
  nGO_total <- dplyr::summarise(.data = nGO_total, nGO_total = dplyr::n())

  ## Number of unique genes in each GROUP
  df_nGenes_group <- dplyr::tally(dplyr::group_by_at(.tbl = gene_2_group, 2), name = "nGene_group")
  df_nGenes_group <- dplyr::rename(.data = df_nGenes_group, group = 1)

  ## Frequency of GO term WITHIN CONDITION
  df_cluster_go <- tidyr::nest(data = dplyr::group_by_at(.tbl = gene_2_group, 2), .key = "genes_in_condition") ## Genes in group
  df_cluster_go <- dplyr::mutate(.data = df_cluster_go,
                                 genes_cluster_GOID = purrr::map(.data$genes_in_condition, ~{     ## Iterate over genes in each condition
                                   t <- dplyr::left_join(x = .x, gene_2_go)                 ## Getting GO terms for each gene in group
                                   t <- dplyr::group_by_at(.tbl = t, 2)                     ## Grouping data by GO term
                                   t <- dplyr::tally(x = t, name = "inClust_nGene_withGO")    ## Frequency of GO term within grouping - More terms than genes as genes can have MULTIPLE GO terms associated with it
                                 }))

  df_cluster_go <- magrittr::set_names(x = df_cluster_go[["genes_cluster_GOID"]],
                                       value = df_cluster_go[[1]])                  ## Extracting frequency table for each group as list - naming list
  df_cluster_go <- dplyr::bind_rows(df_cluster_go, .id = "group")                   ## Binding each groups frequency table - using list names in group

  ## Joining all information together
  df_build_from <- dplyr::left_join(x = df_nGenes_group, df_cluster_go)
  df_build_from <- dplyr::left_join(x = df_build_from, nGO_total)

  ## Dealing with doubling up: Subtracting gene in group from total
  df_build_from <- dplyr::mutate(.data = df_build_from,                                        ## Subtracting gene groupings from total
                                 nGene_remain = nGenes_total - .data$nGene_group,                    ## This should remove the genes in the group from the total
                                 outClust_nGene_withGO = nGO_total - .data$inClust_nGene_withGO,     ## meaning we're not doubling up on genes in the group AND
                                 inClust_nGene_outGO = .data$nGene_group - .data$inClust_nGene_withGO,     ## in the total.
                                 outClust_nGene_outGO = .data$nGene_remain - .data$outClust_nGene_withGO)
  df_build_from <- tidyr::drop_na(data = df_build_from)                                        ## Don't care about NA as this represents no GO term
  df_build_from <- dplyr::select(df_build_from,
                                 c(1,3), -.data$nGene_group,
                                 .data$inClust_nGene_withGO, .data$outClust_nGene_withGO,
                                 .data$inClust_nGene_outGO, .data$outClust_nGene_outGO)

  ## Fishers Exact Test
  out <- dplyr::group_by_at(.tbl = df_build_from, c(1,2))
  out <- tidyr::nest(data = out)
  out <- dplyr::mutate(.data = out,
                       fisher = purrr::map(.x = .data$data, ~{fisher.test(x = matrix(data = unlist(.x), nrow = 2), alternative = "two.sided")}),  ## Fisher test
                       pVal = unlist(purrr::map(.data$fisher, ~{.x["p.value"]})),                                                                 ## Extracting p-value
                       conf_low = unlist(purrr::map(.data$fisher, ~{.x[["conf.int"]][1]})),                                                       ## Confidence intervals
                       conf_high = unlist(purrr::map(.data$fisher, ~{.x[["conf.int"]][2]})),
                       OR = unlist(purrr::map(.data$fisher, ~{.x["estimate"]})))                                                                  ## Odds ratio

  ## Statistical correction
  out <- dplyr::group_by(.data = out, .data$group)            ## Grouping to limit correction to group - not to all samples
  out <- dplyr::mutate(.data = out, adjP = stats::p.adjust(p = .data$pVal, method = pAdj_method))  ## Correcting using p.adjust
  out <- dplyr::select(out, 1, 2, 5, 9, 6, 7, 8, dplyr::everything())     ## Arrange columns
  out <- dplyr::arrange(.data = out, .data$adjP)                         ## Sort by group and adjusted-p
  out <- dplyr::ungroup(out)

  ## Adjusted p-value cut-off
  if(!is.null(p_cutoff)){
    out <- dplyr::filter(.data = out, .data$adjP <= p_cutoff)
  }

  return(out)

}
