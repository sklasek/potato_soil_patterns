better upset plots
================
Scott Klasek
2023-01-03

#### Purpose

Make better upset plots for figure 2.

#### load libraries

``` r
packages <- c("tidyverse", "phyloseq", "speedyseq", "patchwork", "NatParksPalettes", "UpSetR", "scales")
t <- lapply(packages, require, character.only = TRUE) # load all packages at once
```

    ## Loading required package: tidyverse

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.4.0     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.1
    ## ✔ readr   2.1.2     ✔ forcats 0.5.2
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## Loading required package: phyloseq
    ## 
    ## Loading required package: speedyseq
    ## 
    ## 
    ## Attaching package: 'speedyseq'
    ## 
    ## 
    ## The following objects are masked from 'package:phyloseq':
    ## 
    ##     filter_taxa, plot_bar, plot_heatmap, plot_tree, psmelt, tax_glom,
    ##     tip_glom, transform_sample_counts
    ## 
    ## 
    ## Loading required package: patchwork
    ## 
    ## Loading required package: NatParksPalettes
    ## 
    ## Loading required package: UpSetR
    ## 
    ## Loading required package: scales
    ## 
    ## 
    ## Attaching package: 'scales'
    ## 
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     discard
    ## 
    ## 
    ## The following object is masked from 'package:readr':
    ## 
    ##     col_factor

#### load and clean up data

``` r
# ITS phyloseq object
all.fung.1920.c.ps <- readRDS(file = "/Users/scottklasek/Desktop/UMN/phyloseqs/fung1920.cleaned.ps") # cleaned fungal phyloseq object
filter <- phyloseq::genefilter_sample(all.fung.1920.c.ps, filterfun_sample(function(x) x > 0)) # for each ASV get a true/false whether it has over 0 reads
table(filter) # 68 taxa were present only in blanks or omitted samples < 10k reads
```

    ## filter
    ## FALSE  TRUE 
    ##    68 22035

``` r
all.fung.1920.c.ps <- prune_taxa(filter, all.fung.1920.c.ps) # remove these zero-count taxa from the ps object
all.fung.1920.c.ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:          [ 22035 taxa and 1355 samples ]:
    ## sample_data() Sample Data:        [ 1355 samples by 28 sample variables ]:
    ## tax_table()   Taxonomy Table:     [ 22035 taxa by 7 taxonomic ranks ]:
    ## refseq()      DNAStringSet:       [ 22035 reference sequences ]
    ## taxa are columns

``` r
# bact venn object binary
bact.venn_obj.binary <- read.csv(file = "/Users/scottklasek/Desktop/UMN/figures/from_server/16s.venn_obj.binary.csv", row.names = 1)
```

#### make a venn binary dataframe for ITS

``` r
# merge samples by state
ps.state.upset <- merge_samples(all.fung.1920.c.ps, sample_data(all.fung.1920.c.ps)$state, fun = sum)
```

    ## Warning in asMethod(object): NAs introduced by coercion

    ## Warning in asMethod(object): NAs introduced by coercion

    ## Warning in asMethod(object): NAs introduced by coercion

    ## Warning in asMethod(object): NAs introduced by coercion

    ## Warning in asMethod(object): NAs introduced by coercion

    ## Warning in asMethod(object): NAs introduced by coercion

    ## Warning in asMethod(object): NAs introduced by coercion

    ## Warning in asMethod(object): NAs introduced by coercion

    ## Warning in asMethod(object): NAs introduced by coercion

    ## Warning in asMethod(object): NAs introduced by coercion

    ## Warning in asMethod(object): NAs introduced by coercion

    ## Warning in asMethod(object): NAs introduced by coercion

    ## Warning in asMethod(object): NAs introduced by coercion

    ## Warning in asMethod(object): NAs introduced by coercion

    ## Warning in asMethod(object): NAs introduced by coercion

``` r
# transposes the otu table and writes it as a dataframe
venn_obj <- as.data.frame(t(otu_table(ps.state.upset)))
venn_obj.binary <- sapply(venn_obj, function(x) ifelse(x > 0, 1, 0),
                          USE.NAMES = T) # makes the dataframe binary 
rownames(venn_obj.binary) <- rownames(venn_obj) # assigns the rownames as the same
venn_obj.binary <- as.data.frame(venn_obj.binary) # creates a dataframe where rows are ASVs, columns are sites, and entries are either 1 or 0
```

#### plot ITS

``` r
cols.west.east <- c("#DC4405", "#B3A369", "#236192", "#7A0019", "#FFC72A", "#C5050C", "#18453B", "#B0D7FF", "gray") # west-east color palette
colnames(venn_obj.binary) <- c("CO", "ID", "ME1", "MI", "MN1", "MN2", "OR", "ME2", "WI") # fix the names of the sites

its.state.upset <- upset(venn_obj.binary, 
      sets = rev(c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2")),
      sets.bar.color = rev(cols.west.east),
      mainbar.y.label = 'ITS ASV count',
      sets.x.label = 'Observed ASVs by site',
      show.numbers = FALSE,
      order.by = 'freq', nintersects = 25, # can change order.by to "degree", or adjust number of intersections to show 
      keep.order = T)
its.state.upset # 
```

![](19_better_upset_plots_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

#### plot bacteria

``` r
colnames(bact.venn_obj.binary) <- c("CO", "ID", "ME1", "MI", "MN1", "MN2", "OR", "ME2", "WI") # fix the names of the sites
bact.state.upset <- upset(bact.venn_obj.binary, 
      sets = rev(c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2")),
      sets.bar.color = rev(cols.west.east),
      mainbar.y.label = '16S ASV count',
      sets.x.label = 'Observed ASVs by site',
      show.numbers = FALSE,
      order.by = 'freq', nintersects = 39, # can change order.by to "degree", or adjust number of intersections to show 
      keep.order = T)
bact.state.upset
```

![](19_better_upset_plots_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

#### session info

``` r
sessionInfo()
```

    ## R version 4.2.2 (2022-10-31)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] scales_1.2.1           UpSetR_1.4.0           NatParksPalettes_0.2.0
    ##  [4] patchwork_1.1.2        speedyseq_0.5.3.9018   phyloseq_1.40.0       
    ##  [7] forcats_0.5.2          stringr_1.4.1          dplyr_1.0.9           
    ## [10] purrr_0.3.4            readr_2.1.2            tidyr_1.2.0           
    ## [13] tibble_3.1.8           ggplot2_3.4.0          tidyverse_1.3.2       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-160           bitops_1.0-7           fs_1.5.2              
    ##  [4] lubridate_1.8.0        httr_1.4.4             GenomeInfoDb_1.32.3   
    ##  [7] tools_4.2.2            backports_1.4.1        vegan_2.6-2           
    ## [10] utf8_1.2.2             R6_2.5.1               mgcv_1.8-41           
    ## [13] DBI_1.1.3              BiocGenerics_0.42.0    colorspace_2.0-3      
    ## [16] permute_0.9-7          rhdf5filters_1.8.0     ade4_1.7-19           
    ## [19] withr_2.5.0            gridExtra_2.3          tidyselect_1.1.2      
    ## [22] compiler_4.2.2         cli_3.4.1              rvest_1.0.3           
    ## [25] Biobase_2.56.0         xml2_1.3.3             labeling_0.4.2        
    ## [28] digest_0.6.29          rmarkdown_2.15         XVector_0.36.0        
    ## [31] pkgconfig_2.0.3        htmltools_0.5.3        highr_0.9             
    ## [34] dbplyr_2.2.1           fastmap_1.1.0          rlang_1.0.6           
    ## [37] readxl_1.4.1           rstudioapi_0.14        farver_2.1.1          
    ## [40] generics_0.1.3         jsonlite_1.8.0         googlesheets4_1.0.1   
    ## [43] RCurl_1.98-1.8         magrittr_2.0.3         GenomeInfoDbData_1.2.8
    ## [46] biomformat_1.24.0      Matrix_1.5-3           Rcpp_1.0.9            
    ## [49] munsell_0.5.0          S4Vectors_0.34.0       Rhdf5lib_1.18.2       
    ## [52] fansi_1.0.3            ape_5.6-2              lifecycle_1.0.3       
    ## [55] stringi_1.7.8          yaml_2.3.5             MASS_7.3-58.1         
    ## [58] zlibbioc_1.42.0        rhdf5_2.40.0           plyr_1.8.7            
    ## [61] grid_4.2.2             parallel_4.2.2         crayon_1.5.1          
    ## [64] lattice_0.20-45        splines_4.2.2          Biostrings_2.64.1     
    ## [67] haven_2.5.1            multtest_2.52.0        hms_1.1.2             
    ## [70] knitr_1.40             pillar_1.8.1           igraph_1.3.4          
    ## [73] reshape2_1.4.4         codetools_0.2-18       stats4_4.2.2          
    ## [76] reprex_2.0.2           glue_1.6.2             evaluate_0.16         
    ## [79] data.table_1.14.2      modelr_0.1.9           vctrs_0.5.1           
    ## [82] tzdb_0.3.0             foreach_1.5.2          cellranger_1.1.0      
    ## [85] gtable_0.3.0           assertthat_0.2.1       xfun_0.32             
    ## [88] broom_1.0.0            survival_3.4-0         googledrive_2.0.0     
    ## [91] gargle_1.2.0           iterators_1.0.14       IRanges_2.30.1        
    ## [94] cluster_2.1.4          ellipsis_0.3.2
