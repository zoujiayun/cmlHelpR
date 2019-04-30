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
  nGenes_total <- gene_GOID_mapping %>%
    dplyr::pull(gene) %>%
    unique() %>%
    length()

  ## Total frequency of GO terms
  nGO_total <- gene_GOID_mapping %>%
    dplyr::group_by(GO_ID) %>%
    dplyr::summarise(nGO_total = n())

  ## Number of unique genes in each CONDITION
  df.nGenes.cluster <- gene_cluster_mapping %>%
    dplyr::group_by(cluster) %>%
    dplyr::tally(name = "nGene_cluster")

  ## Frequency of GO term WITHIN CONDITION
  df.cluster.go <- gene_cluster_mapping %>%
    dplyr::group_by(cluster) %>%
    tidyr::nest(.key = "genes_in_condition") %>%
    mutate(genes_cluster_GOID = purrr::map(genes_in_condition, ~{
      dplyr::left_join(.x, go_mapped) %>%
        dplyr::group_by(GO_ID) %>%
        dplyr::tally(name = "inClust_nGene_inGO")
    }))

  df.cluster.go <- df.cluster.go$genes_cluster_GOID %>%
    magrittr::set_names(df.cluster.go$cluster) %>%
    dplyr::bind_rows(.id = "cluster")

  ## Joining all information together
  df.build.from <- df.nGenes.cluster %>%
    dplyr::left_join(df.cluster.go) %>%
    dplyr::left_join(nGO_total) %>%
    dplyr::mutate(nGene_remain = nGenes_total - nGene_cluster,
                  outClust_nGene_inGO = nGO_total - inClust_nGene_inGO,
                  inClust_nGene_outGO = nGene_cluster - inClust_nGene_inGO,
                  outClust_nGene_outGO = nGene_remain - outClust_nGene_inGO) %>%
    tidyr::drop_na() %>%
    dplyr::select(cluster, GO_ID, -nGene_cluster,
                  inClust_nGene_inGO, outClust_nGene_inGO,
                  inClust_nGene_outGO, outClust_nGene_outGO)

  ## Dropping NA above:
  ## Genes with no GO term have been used to build the data frame for the analyses.
  ## They need to be included as they are still part of the dataset, they just don't have a GO annotation.
  ## However, NA row is just a grouping and does not need to have a fisher's statistic
  ## calculated. It is not of interest to find out if NA is significant.

  ## Fishers Exact Test
  df.build.from %>%
    dplyr::group_by(cluster, GO_ID) %>%
    tidyr::nest() %>%
    dplyr::mutate(fisher = purrr::map(data, ~{fisher.test(x = matrix(data = unlist(.x),
                                                                     nrow = 2),
                                                          alternative = "two.sided")}),
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
