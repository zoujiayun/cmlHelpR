# cmlHelpR

This package and repository is a simple way for me to collect functions I've written to analyses `CODEML` output. 

The functions are continually being tweaked, and in many instances serve a specific purpose for what I am trying to do. Nevertheless, I've endeavoured to make the functions as extensible as possible.

In the sections below I'll demonstrate some use cases of the major functions of interest.

# Installation
In the `R-studio` console, execute the following code.

```
> devtools::github_install("a-lud/cmlHelpR")
> library(cmlHelpR)
```

# Tutorial

I have included below a few examples of how to use the main package functions. I'll try and include what the starting data looks like and what the output is.

## Function: getOrthologHeaders()

This function builds the sequencee header - gene filename table. Using the Conditional Reciprocal Best Blast (CRBB) output, it identifies CRBB-hits across all sample comparisons and joins consistent CRBBs. 

For example, let's say we have three samples and are dealing with one sequence; *Gene-A*. A CRBB consistent across three samples would look like so

`Sample 1 == Sample 1 == sample 3`

That is to say, the same sequences in each of the species hits the equivalent sequence in the other two samples.

**Usage**:
```
library(cmlHelpR)
out <- getOrthologueHeaders(crbb_path = "path/to/crbb/out", 
                                      crbb_ext = ".tsv", 
                                      sample_separator = "_", 
                                      id_pct = 95, 
                                      aln_pct = 90)
```

**Output**:

The output is a list object that stores a wide format dataframe and a nested dataframe.

The wide format dataframe is a good example what the explanation above. The first column is the gene symbols parsed from the end of each samples header line (see `parallel_codeml` README for building file this way). Columns 2-4 are the actual headers from each samples respective fasta file which will be used to parse the sequences in later steps.

```
$wide_format
# A tibble: 3,119 x 4
   file_name                       Laevis                          Scutatus                             Textilis                            
   <chr>                           <chr>                           <chr>                                <chr>                               
 1 IGHM-id342200-id336538          evm.model.scaffold6153.1_IGHM   id342200                             id336538                            
 2 HV333-id379141-id368235         evm.model.scaffold91631.1_HV333 id379141                             id368235                            
 3 EDRF1-EDRF1-EDRF1               evm.model.scaffold1962.2_EDRF1  rna10_XM_026668936.1_EDRF1           rna13700_XM_026704949.1_EDRF1       
 4 TM14C-LOC113413631-LOC113445701 evm.model.scaffold133.1_TM14C   rna10007_XM_026670054.1_LOC113413631 rna23543_XM_026715005.1_LOC113445701
 5 MYNN-MYNN-MYNN                  evm.model.scaffold16556.1_MYNN  rna1001_XM_026666067.1_MYNN          rna14532_XM_026705790.1_MYNN        
 6 KCRB-CKB-CKB                    evm.model.scaffold1822.3_KCRB   rna10012_XM_026670080.1_CKB          rna7197_XM_026698556.1_CKB          
 7 BAG5-BAG5-BAG5                  evm.model.scaffold1369.1_BAG5   rna10016_XM_026670087.1_BAG5         rna7199_XM_026698559.1_BAG5         
 8 KLC1-KLC1-KLC1                  evm.model.scaffold1369.3_KLC1   rna10023_XM_026670084.1_KLC1         rna7205_XM_026698564.1_KLC1         
 9 ZFY21-ZFYVE21-ZFYVE21           evm.model.scaffold7335.2_ZFY21  rna10030_XM_026670104.1_ZFYVE21      rna7211_XM_026698576.1_ZFYVE21      
10 SIVA-SIVA1-SIVA1                evm.model.scaffold18146.1_SIVA  rna10054_XM_026670116.1_SIVA1        rna7235_XM_026698589.1_SIVA1        
# … with 3,109 more rows

$nested_format
# A tibble: 3,119 x 2
   file_name                       data            
   <chr>                           <list>          
 1 IGHM-id342200-id336538          <tibble [3 × 2]>
 2 HV333-id379141-id368235         <tibble [3 × 2]>
 3 EDRF1-EDRF1-EDRF1               <tibble [3 × 2]>
 4 TM14C-LOC113413631-LOC113445701 <tibble [3 × 2]>
 5 MYNN-MYNN-MYNN                  <tibble [3 × 2]>
 6 KCRB-CKB-CKB                    <tibble [3 × 2]>
 7 BAG5-BAG5-BAG5                  <tibble [3 × 2]>
 8 KLC1-KLC1-KLC1                  <tibble [3 × 2]>
 9 ZFY21-ZFYVE21-ZFYVE21           <tibble [3 × 2]>
10 SIVA-SIVA1-SIVA1                <tibble [3 × 2]>
# … with 3,109 more rows

```

## Function: writeFasta()

As the name suggests, this function writes multi-fasta files using the table from above. For each gene, it parses the sequence from each respective fasta file and writes them all to one multi-fasta file.

**Usage**:

```
library(cmlHelpR)
writeFasta(orthologies = tbl_lst, 
                         fasta_dir = "/path/to/sequence_fastas", 
                         pep_ext = ".pep", 
                         nuc_ext = ".fa", 
                         pep_out = "/path/to/peptide/outDir", 
                         nuc_out = "/path/to/nucleotide/outDir")
```

**Output**

Here, the output is both the fasta files being written in peptide/nucelotide format at the specified directories, but also a list object that contains the following for each gene:

- Peptide sequences
- Nucleotide sequences
- Logical indicating if an internal stop codon was found

```
$peptide
  A AAStringSet instance of length 3
    width seq                                                                                                                                                                                     names               
[1]   457 SPKAPSLFPLIPSGDNSETIDITIGCLAKNFLPDSIDFSWDNQQNQSIGNQNYIKFPSILSSGTYTAVSQAKVPRSTWDGFQLFYCKATH...WLQNEQPVSESVYFTSKAILESKIQSKGYFAYSMLNISEQEWSAGDSFTCVVGHEAFPYNSIQKTVNKNTGKPSIVNVSLVLSDTSTPCY Laevis
[2]   480 MTIVHCDNWFDYWGKGTSVTVTAESPKAPSLFPLIPSGDNSETIDITIGCLAKNFLPDSIDFSWDNQRNQSIGNQNYIKFPSILSSGTYT...WLQNEQPVSESVYFTSKAILESKIQSKEYFAYSMLNINEQEWSAGDSFTCVVGHEALHYNSIQKTVNKNTGKPSIVNVSLVLSDTSTPCY Scutatus
[3]   465 MVTVSSETPKAPSLFPLIPSGDNSETIDITIGCLAKNFLPDSIDFSWDNQRNESIGNQNYIKFPSILSSGTYTAVSQAKVPRSTWDEFQL...WLQNEQPVSESVYFTSKAISESKIQSKEYFAYSMLNISEQEWSAGDSFTCVVGHEALYYGSIQKTVNKNTGKPSIVNVSLVLSDTSTPCY Textilis

$nucleotide
  A AAStringSet instance of length 3
    width seq                                                                                                                                                                                     names               
[1]  1374 TCCCCAAAGGCCCCTTCTCTTTTCCCACTCATCCCATCTGGGGACAACTCAGAGACCATAGATATCACCATCGGATGTCTTGCCAAGAAT...ATTCAGAAGACTGTAAACAAGAACACGGGTAAACCATCCATAGTCAACGTCTCCTTAGTCCTCTCCGACACTTCCACCCCTTGCTATTAA Laevis
[2]  1443 ATGACTATTGTTCATTGTGACAACTGGTTTGACTATTGGGGAAAAGGCACTTCAGTCACCGTTACTGCAGAATCCCCAAAGGCCccttct...ATTCAGAAGACTGTAAACAAGAACACGGGTAAACCATCCATAGTCAACGTCTCCTTAGTCCTCTCCGACACTTCCACCCCTTGCTATTAA Scutatus
[3]  1398 ATGGTCACAGTCAGCTCAGAAACCCCAAAGGCCCCTTCTCTTTTCCCACTCATCCCATCTGGGGACAACTCAGAGACCATAGATATCACC...ATTCAGAAGACTGTAAACAAGAACACGGGTAAACCATCCATAGTCAACGTCTCCTTAGTCCTCTCCGACACTTCCACCCCTTGCTATTAA Textilis

$internalStop
[1] FALSE
```

## Function: lrtStatistic()

The `lrtStatistic()` function calculates likelihood ratio tests betweeen two different evolutionary models. In the `Parallel Codeml` repository is an excel file describing which models should be compared (i.e. which models make sense to compare).

**Usage**:

The input for this command is a list of model comparisons and the file path to the output directory that houses the `CODEML` output.

```
## Libraries needed
> library(cmlHelpR)
> library(tidyverse)
> library(magrittr)

## Build all model comparisons from excel file from paralle codeml
> comp <- readr::read_csv(file = "/path/to/pararllel_codeml/ctl_parameter_files/190410_model_comparisons.csv") %>%
                          select(3:4) %>%
                          filter(comparison != "NA") %>%
                          group_split(id) %>%
                          map(~flatten_chr(.), . == "")
                          
## An alternative way to make the comparisons
> comp <- list(c("M1a", "M2a"), c("ModelA", "M1a"))

## Run comand
> out <- lrtStatistic(dir_path = "/path/to/codeml/outputDir", lst_comparisons = comp)
```

**Output**:

The output is a list of dataframes, where each data frame corresponds to each model comparison.

```
$M2a_M1a
# A tibble: 9 x 9
  gene              tree             np_M1a np_M2a lnL_M1a lnL_M2a      delta    df         pval
  <chr>             <chr>             <dbl>  <dbl>   <dbl>   <dbl>      <dbl> <dbl>        <dbl>
1 1433E-YWHAE-YWHAE brownForeground       7      9  -1026.  -1026.  0.0000660     2        1.000       
2 1433E-YWHAE-YWHAE laevisForeground      7      9  -1026.  -1026.  0.0000660     2        1.000       
3 1433E-YWHAE-YWHAE tigerForeground       7      9  -1026.  -1026.  0.0000660     2        1.000       
4 1433T-YWHAQ-YWHAQ brownForeground       7      9  -1024.  -1006.  35.7          2        0.0000000179
5 1433T-YWHAQ-YWHAQ laevisForeground      7      9  -1024.  -1006.  35.7          2        0.0000000179
6 1433T-YWHAQ-YWHAQ tigerForeground       7      9  -1024.  -1006.  35.7          2        0.0000000179
7 1433Z-YWHAZ-YWHAZ brownForeground       7      9   -967.   -967.  0.0000340     2        1.000       
8 1433Z-YWHAZ-YWHAZ laevisForeground      7      9   -967.   -967.  0.0000340     2        1.000       
9 1433Z-YWHAZ-YWHAZ tigerForeground       7      9   -967.   -967.  0.0000340     2        1.000
```

## Function: parseBEB()

Bayes Empirical Bayes (BEB) data is reported for some models (e.g. ModelA, M2a and M8) which represents putative sites under selection. This function parses that information into an easy to use dataframe object that can be manipulated for down-stream analysis.

**Usage**

Input is simply the directory path to where the `CODEML` output is, the models you want to extract the BEB data for and the number of cores to parse the outputs with. It's important to note here that this function operates off the output directories having the name of their model. If you change the output directory names, e.g. to `M2a_first_output`, this function will not work.

I have written these scripts to work directly with the output of the `parallel codeml` pipeline. If you want to update them, by all means do.

```
out <- parseBEB(dir_path = "/path/to/codeml/output", 
              models = c("M2a", "ModelA"), 
              cores = 4)
```

**Output**:

The output is a list object housing a list split by model and a nested dataframe object, where the data are grouped by their models. Rows that are all `NA` are those that did not have any BEB information to parse.

```
$long
   # A tibble: 62,594 x 11
   model gene              tree               pos aa       val   pval signif postMean    SE no_gap_length
   <chr> <chr>             <chr>            <dbl> <chr>  <dbl>  <dbl> <chr>     <dbl> <dbl>         <dbl>
    1 M2a   1433E-YWHAE-YWHAE brownForeground     NA NA    NA     NA     NA        NA    NA              765
    2 M2a   1433E-YWHAE-YWHAE laevisForeground    NA NA    NA     NA     NA        NA    NA              765
    3 M2a   1433E-YWHAE-YWHAE tigerForeground     NA NA    NA     NA     NA        NA    NA              765
    4 M2a   1433T-YWHAQ-YWHAQ brownForeground    227 E      0.991  0.009 **         9.63  1.37           735
    5 M2a   1433T-YWHAQ-YWHAQ brownForeground    228 Y      0.978  0.022 *          9.51  1.69           735
    6 M2a   1433T-YWHAQ-YWHAQ brownForeground    229 G      0.973  0.027 *          9.47  1.79           735
    7 M2a   1433T-YWHAQ-YWHAQ brownForeground    230 F      0.941  0.059 NA         9.19  2.37           735
    8 M2a   1433T-YWHAQ-YWHAQ brownForeground    231 I      0.967  0.033 *          9.42  1.91           735
    9 M2a   1433T-YWHAQ-YWHAQ brownForeground    232 I      0.854  0.146 NA         8.40  3.39           735
   10 M2a   1433T-YWHAQ-YWHAQ laevisForeground   227 E      0.991  0.009 **         9.63  1.37           735
$list
   $list$M2a
   # A tibble: 23,967 x 11
      model gene              tree               pos aa       val   pval signif postMean    SE no_gap_length
      <chr> <chr>             <chr>            <dbl> <chr>  <dbl>  <dbl> <chr>     <dbl> <dbl>         <dbl>
    1 M2a   1433E-YWHAE-YWHAE brownForeground     NA NA    NA     NA     NA        NA    NA              765
    2 M2a   1433E-YWHAE-YWHAE laevisForeground    NA NA    NA     NA     NA        NA    NA              765
    3 M2a   1433E-YWHAE-YWHAE tigerForeground     NA NA    NA     NA     NA        NA    NA              765
    4 M2a   1433T-YWHAQ-YWHAQ brownForeground    227 E      0.991  0.009 **         9.63  1.37           735
    5 M2a   1433T-YWHAQ-YWHAQ brownForeground    228 Y      0.978  0.022 *          9.51  1.69           735
    6 M2a   1433T-YWHAQ-YWHAQ brownForeground    229 G      0.973  0.027 *          9.47  1.79           735
    7 M2a   1433T-YWHAQ-YWHAQ brownForeground    230 F      0.941  0.059 NA         9.19  2.37           735
    8 M2a   1433T-YWHAQ-YWHAQ brownForeground    231 I      0.967  0.033 *          9.42  1.91           735
    9 M2a   1433T-YWHAQ-YWHAQ brownForeground    232 I      0.854  0.146 NA         8.40  3.39           735
   10 M2a   1433T-YWHAQ-YWHAQ laevisForeground   227 E      0.991  0.009 **         9.63  1.37           735
   # … with 23,957 more rows

$list$ModelA
   # A tibble: 38,627 x 11
      model  gene              tree               pos aa       val    pval signif postMean    SE no_gap_length
      <chr>  <chr>             <chr>            <dbl> <chr>  <dbl>   <dbl> <chr>     <dbl> <dbl>         <dbl>
    1 ModelA 1433E-YWHAE-YWHAE brownForeground     NA NA    NA     NA      NA           NA    NA           765
    2 ModelA 1433E-YWHAE-YWHAE laevisForeground    NA NA    NA     NA      NA           NA    NA           765
    3 ModelA 1433E-YWHAE-YWHAE tigerForeground     NA NA    NA     NA      NA           NA    NA           765
    4 ModelA 1433T-YWHAQ-YWHAQ brownForeground     NA NA    NA     NA      NA           NA    NA           735
    5 ModelA 1433T-YWHAQ-YWHAQ laevisForeground   227 E      0.995  0.005  **           NA    NA           735
    6 ModelA 1433T-YWHAQ-YWHAQ laevisForeground   228 Y      0.987  0.013  *            NA    NA           735
    7 ModelA 1433T-YWHAQ-YWHAQ laevisForeground   229 G      0.986  0.014  *            NA    NA           735
    8 ModelA 1433T-YWHAQ-YWHAQ laevisForeground   230 F      0.97   0.03   *            NA    NA           735
    9 ModelA 1433T-YWHAQ-YWHAQ laevisForeground   231 I      0.985  0.015  *            NA    NA           735
   10 ModelA 1433T-YWHAQ-YWHAQ laevisForeground   232 I      0.909  0.0910 NA           NA    NA           735
   # … with 38,617 more rows

$nested
# A tibble: 2 x 2
  model  BEB                   
  <chr>  <list>                
1 M2a    <tibble [23,967 × 10]>
2 ModelA <tibble [38,627 × 10]>
```

**BEBDF2Freq()**

**getBranchdNdS()**

**getPW()**

**goFisher()**

## Author 

- **Alastair Ludington**: alastair.ludington@adelaide.edu.au
