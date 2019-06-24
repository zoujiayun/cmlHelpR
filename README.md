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
out <- cmlHelpR::getOrthologueHeaders(crbb_path = "path/to/crbb/out", 
                                      crbb_ext = ".tsv", 
                                      sample_separator = "_", 
                                      id_pct = 95, 
                                      aln_pct = 90)
```

**Output**:

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


**Input** 

**write_fasta()**

**lrt_statistics()**

**parse_BEB()**

**BEBDF2Freq()**

**getBranchdNdS()**

**getPW()**

**goFisher()**

## Author 

- **Alastair Ludington**: alastair.ludington@adelaide.edu.au
