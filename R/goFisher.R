#' Conduct Gene Ontology within gene clusters using Fisher Test
#'
#' Conducts a Fisher's Exact Test on all GO terms within a gene cluster.
#' @param df_clusters Data frame that contains
#' @param pAdj_method
#' @keywords helper
#' @export
#' @examples
#' goFisher(df_clusters = df.clusters.goTerms, pAdj_method = "bonferroni")
goFisher <- function(df_clusters, pAdj_method){

  ## Total number of genes in analysis
  nGenes_total <- df_clusters %>%
    dplyr::pull(gene) %>%
    unique() %>%
    length()

  ## Frequency counts of GO terms
  nGO_total <- df_clusters %>%
    dplyr::group_by(GO_ID) %>%
    dplyr::summarise(nGO_total = n())

  ## Number of unique genes per cluster
  df.nGenes.cluster <- df_clusters %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(nGene_cluster = length(unique(gene)))

  ## Number of genes with each GOterms per cluster
  df.cluster.go <- df_clusters %>%
    dplyr::group_by(cluster, GO_ID) %>%
    dplyr::summarise(inClust_nGene_inGO = n())

  ## Joining all information together
  df.build.from <- df.nGenes.cluster %>%
    dplyr::left_join(df.cluster.go) %>%
    dplyr::left_join(nGO_total) %>%
    dplyr::mutate(nGene_remain = nGenes_total - nGene_cluster,
                  outClust_nGene_inGO = nGO_total - inClust_nGene_inGO,
                  inClust_nGene_outGO = nGene_cluster - inClust_nGene_inGO,
                  outClust_nGene_outGO = nGene_remain - outClust_nGene_inGO,
                  GO_ID = dplyr::na_if(x = GO_ID, "")) %>%
    tidyr::drop_na() %>%
    dplyr::select(cluster, GO_ID, -nGene_cluster,
                  inClust_nGene_inGO, outClust_nGene_inGO,
                  inClust_nGene_outGO, outClust_nGene_outGO)

  ## Dropping NA above:
  ## Genes with no GO term have been used to build the data frame for the analyses.
  ## However, NA row is just a grouping and does not need to have a fisher's statistic
  ## calculated.

  ## Fishers Exact Test
  df.build.from %>%
    dplyr::group_by(cluster, GO_ID) %>%
    tidyr::nest() %>%
    dplyr::mutate(fisher = purrr::map(data, ~{fisher.test(x = matrix(data = unlist(.x),
                                                                     nrow = 2),alternative = "two.sided")}),
                  pVal = unlist(purrr::map(fisher, ~{.x["p.value"]})),
                  conf_low = unlist(purrr::map(fisher, ~{.x[["conf.int"]][1]})),
                  conf_high = unlist(purrr::map(fisher, ~{.x[["conf.int"]][2]})),
                  OR = unlist(purrr::map(fisher, ~{.x["estimate"]}))) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(adjP = p.adjust(p = pVal, method = pAdj_method)) %>%
    dplyr::select(cluster, GO_ID, pVal, adjP, conf_low, conf_high, OR, dplyr::everything()) %>%
    dplyr::arrange(cluster, adjP) %>%
    dplyr::ungroup()

}
