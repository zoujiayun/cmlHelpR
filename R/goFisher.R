#' Conduct Gene Ontology within gene clusters using Fisher Test
#'
#' Conducts a Fisher's Exact Test on all GO terms within a gene cluster.
#' @param df_clusters Data frame that contains
#' @param pAdj_method
#' @keywords helper
#' @export
#' @examples
#' goFisher(df_clusters = df.clusters.goTerms, pAdj_method = "bonferroni")
goFisher <- function(gene_GOID_mapping, gene_cluster_mapping, pAdj_method){

  ## Total number of genes assessed (excluding those with internal stop codons)
  nGenes_total <- length(unique(gene_GOID_mapping[["gene"]]))

  ## Total frequency of GO terms
  nGO_total <- dplyr::group_by(.data = gene_GOID_mapping, GO_ID)
  nGO_total <- dplyr::summarise(.data = nGO_total, nGO_total = dplyr::n())

  ## Number of unique genes in each CONDITION
  df.nGenes.cluster <- dplyr::tally(dplyr::group_by(.data = gene_cluster_mapping, cluster), name = "nGene_cluster")

  ## Frequency of GO term WITHIN CONDITION
  df.cluster.go <- tidyr::nest(data = dplyr::group_by(.data = gene_cluster_mapping, cluster), .key = "genes_in_condition")
  df.cluster.go <- dplyr::mutate(.data = df.cluster.go,
                                 genes_cluster_GOID = purrr::map(genes_in_condition, ~{
                                   t <- dplyr::left_join(x = .x, gene_GOID_mapping)
                                   t <- dplyr::group_by(.data = t, GO_ID)
                                   t <- dplyr::tally(x = t, name = "inClust_nGene_inGO")
                                 }))

  df.cluster.go <- magrittr::set_names(x = df.cluster.go[["genes_cluster_GOID"]], value = df.cluster.go[["cluster"]])
  df.cluster.go <- dplyr::bind_rows(df.cluster.go, .id = "cluster")

  ## Joining all information together
  df.build.from <- dplyr::left_join(x = df.nGenes.cluster, df.cluster.go)
  df.build.from <- dplyr::left_join(x = df.build.from, nGO_total)
  df.build.from <- dplyr::mutate(.data = df.build.from,
                                 nGene_remain = nGenes_total - nGene_cluster,
                                 outClust_nGene_inGO = nGO_total - inClust_nGene_inGO,
                                 inClust_nGene_outGO = nGene_cluster - inClust_nGene_inGO,
                                 outClust_nGene_outGO = nGene_remain - outClust_nGene_inGO)
  df.build.from <- tidyr::drop_na(data = df.build.from)
  df.build.from <- dplyr::select(df.build.from,
                                 cluster, GO_ID, -nGene_cluster,
                                 inClust_nGene_inGO, outClust_nGene_inGO,
                                 inClust_nGene_outGO, outClust_nGene_outGO)

  ## Fishers Exact Test
  out <- dplyr::group_by(.data = df.build.from, cluster, GO_ID)
  out <- tidyr::nest(data = out)
  out <- dplyr::mutate(.data = out,
                       fisher = purrr::map(.x = data, ~{fisher.test(x = matrix(data = unlist(.x), nrow = 2), alternative = "two.sided")}),
                       pVal = unlist(purrr::map(fisher, ~{.x["p.value"]})),
                       conf_low = unlist(purrr::map(fisher, ~{.x[["conf.int"]][1]})),
                       conf_high = unlist(purrr::map(fisher, ~{.x[["conf.int"]][2]})),
                       OR = unlist(purrr::map(fisher, ~{.x["estimate"]}))
  )
  out <- dplyr::group_by(.data = out, cluster)
  out <- dplyr::mutate(.data = out, adjP = p.adjust(p = pVal, method = pAdj_method))
  out <- dplyr::select(out, cluster, GO_ID, pVal, adjP, conf_low, conf_high, OR, dplyr::everything())
  out <- dplyr::arrange(out, cluster, adjP)
  out <- dplyr::ungroup(out)

  return(out)

}

