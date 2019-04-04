#' Build directory structure and control files for codeml in parallel
#'
#' Creates sequence directory with nested directories relating to model and tree combinations.
#' Control files for each of these comparisons are then made according to the specification the user provides in the
#' codeml model parameters file.
#' @param phylipDir Path to directory that contains phylip files
#' @param phylipExt Extension of phylip files
#' @param treeDir Path to directory that contains tree files
#' @param treeExt Extension of tree files
#' @param modelsVec Vector of models that you wish to run
#' @param pathCTLParam File path to CTL parameter file
#' @param outDir Path to output directory (will create if it doesn't exist)
#' @keywords control files, build, directory structure
#' @export
#' @examples
#' build_codemlCTL(phylipDir = "path/to/directory", phylipExt = ".phy",
#' treeDir = "path/to/tree/directory", treeExt = ".tre", modelsVec = c("Model1Neutral", "Model2Selection"),
#' pathCTLParam = "path/to/ctl_params.csv", outDir = "path/to/output/directory")
build_codemlCTL <- function(phylipDir, phylipExt, treeDir, treeExt, modelsVec, pathCTLParam, outDir) {

  ## listing phylip files
  lst.phy <- list.files(path = phylipDir, pattern = phylipExt, full.names = TRUE)
  lst.phy <- magrittr::set_names(x = lst.phy,
                                 value = stringr::str_remove(string = basename(lst.phy),pattern = phylipExt))

  ## Listing tree files
  lst.tree <- list.files(path = treeDir, pattern = treeExt, full.names = TRUE)
  lst.tree <- magrittr::set_names(x = lst.tree, value = stringr::str_remove(string = basename(lst.tree), pattern = treeExt))
  lst.tree <- magrittr::set_names(x = lst.tree, value = stringr::str_remove(string = names(lst.tree), pattern = ".*_"))

  ## Importing ctl parameters + subsetting
  df.ctl.param <- readr::read_csv(file = pathCTLParam, col_names = TRUE)
  # df.ctl.param <- df.ctl.param[df.ctl.param[["models"]] %in% modelsVec, ]

  ## All directories that need to be created
  vec.model_tree <- outer(X = modelsVec, Y = names(lst.tree), FUN = "paste", sep = "_")
  vec.gene.model_tree <- outer(X = names(lst.phy), Y = vec.model_tree, FUN = "paste", sep = "/")
  vec.out.gene.model_tree <- paste(outDir, vec.gene.model_tree, sep = "/")

  ## Creating all output directories
  map(vec.out.gene.model_tree, dir.create, recursive = TRUE)

  ## Copying phylip file to corresponding directory
  purrr::map(names(lst.phy), ~{

    ## Output dirs split by gene
    lst.dirs.gene <- vec.out.gene.model_tree[grepl(pattern = .x, x = vec.out.gene.model_tree, perl = TRUE)]

    ## Copying tree file to subset of dirs
    lapply(lst.dirs.gene, function(y){
      file.copy(from = lst.phy[.x], to = y, overwrite = TRUE)
      write(x = paste0("seqfile = ", basename(path = lst.phy[.x])), file = paste(y,"codeml.ctl", sep = "/"))
    })
  })

  ## Copying tree files to corresponding directory + making first lines of codeml
  purrr::map(names(lst.tree), ~{
    ## output dirs split by tree
    lst.dirs.tree <- vec.out.gene.model_tree[stringr::str_detect(string = vec.out.gene.model_tree, pattern = .x)]

    ## Copying tree file to subset of dirs
    lapply(lst.dirs.tree, function(y){
      file.copy(from = lst.tree[.x], to = y, overwrite = TRUE)
      write(x = paste0("treefile = ", basename(path = lst.tree[.x])), file = paste(y, "codeml.ctl", sep = "/"), append = TRUE)
    })

  })

  ## Appending rest of codeml to ctl files
  lst.cml <- list.files(path = outDir, pattern = "codeml.ctl", full.names = TRUE, recursive = TRUE)
  lst.cml <- map(modelsVec, ~{
    lst.cml[grepl(pattern = .x, x = lst.cml, perl = TRUE)]
  })
  lst.cml <- magrittr::set_names(x = lst.cml, value = modelsVec)

  ## Iterating through models
  map(names(lst.cml), ~{

    df <- df.ctl.param[df.ctl.param[["models"]] %in% .x, ]

    ## Iterating through vectors
    lapply(lst.cml[[.x]], function(y) {

      df %>%
        glue::glue_data("outfile = codeml.output

                        noisy = {noisy}  * 0,1,2,3,9: how much rubbish on the screen
                        verbose = {verbose}  * 0: concise; 1: detailed, 2: too much
                        runmode = {runmode}  * 0: user tree;  1: semi-automatic;  2: automatic
                        * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

                        seqtype = {seqtype}  * 1:codons; 2:AAs; 3:codons-->AAs
                        CodonFreq = {CodonFreq}  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
                        aaDist = {aaDist}  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
                        aaRatefile = {aaRatefile} * only used for aa seqs with model=empirical(_F)
                        * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own
                        clock = {clock}  * 0:no clock, 1:global clock; 2:local clock; 3:TipDate

                        model = {model}
                        * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                        * models for AAs or codon-translated AAs:
                        * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                        * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

                        NSsites = {NSsites}  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                        * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                        * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                        * 13:3normal>0

                        icode = {icode}  * 0:universal code; 1:mammalian mt; 2-10:see below
                        fix_kappa = {fix_kappa}  * 1: kappa fixed, 0: kappa to be estimated
                        kappa = {kappa}  * initial or fixed kappa
                        fix_omega = {fix_omega}  * 1: omega or omega_1 fixed, 0: estimate
                        omega = {omega}  * initial or fixed omega, for codons or codon-based AAs
                        ncatG = {ncatG}  * # of categories in dG of NSsites models
                        getSE = {getSE}  * 0: don't want them, 1: want S.E.s of estimates
                        RateAncestor = {RateAncestor}  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
                        fix_alpha = {fix_alpha}  * 0: estimate gamma shape parameter; 1: fix it at alpha
                        alpha = {alpha}  * initial or fixed alpha, 0:infinity (constant rate)
                        Small_Diff = {Small_Diff}
                        cleandata = {cleandata}") %>%
        write(file = y, append = TRUE)

    })

  })

  lst.ctl <- list.files(path = outDir, pattern = ".ctl", full.names = TRUE, recursive = TRUE)
  write(x = lst.ctl, file = "ctl_list.txt")

}
