table info
================
Scott Klasek
2023-01-03

#### Purpose

To make pretty tables. JK, the first table will show numbers of
microbiome samples from each site x year x season x amplicon, and the
second will show cultivar info by site.

#### libraries

``` r
packages <- c("tidyverse", "phyloseq", "plyr")
t <- lapply(packages, require, character.only = TRUE) # load all packages at once
```

    ## Loading required package: tidyverse

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.4.1     ✔ purrr   1.0.1
    ## ✔ tibble  3.1.8     ✔ dplyr   1.1.0
    ## ✔ tidyr   1.3.0     ✔ stringr 1.5.0
    ## ✔ readr   2.1.4     ✔ forcats 1.0.0
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## Loading required package: phyloseq
    ## 
    ## Loading required package: plyr
    ## 
    ## ------------------------------------------------------------------------------
    ## 
    ## You have loaded plyr after dplyr - this is likely to cause problems.
    ## If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
    ## library(plyr); library(dplyr)
    ## 
    ## ------------------------------------------------------------------------------
    ## 
    ## 
    ## Attaching package: 'plyr'
    ## 
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     arrange, count, desc, failwith, id, mutate, rename, summarise,
    ##     summarize
    ## 
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     compact

#### load data

``` r
all.fung.1920.c.ps <- readRDS(file = "/Users/klas0061/Desktop/UMN/phyloseqs/fung1920.cleaned.ps") # cleaned fungal phyloseq object
core.bact.ps <- readRDS(file = "/Users/klas0061/Desktop/UMN/phyloseqs/core/core.bact.asvs.ps") # core bacterial phyloseq object to save memory
```

#### not a table anymore, but counts of ITS and 16S communities by site.

``` r
its.site <- table(data.frame(all.fung.1920.c.ps@sam_data)$state) # ITS microbiome counts by state
bact.site <- table(data.frame(core.bact.ps@sam_data)$state) # 16S microbiome counts by state
site.counts <- rbind("ITS count" = its.site, "16S count" = bact.site) # bind its and 16s counts
colnames(site.counts) <- c("CO", "ID", "ME1", "MI", "MN1", "MN2", "OR", "ME2", "WI") # correct site labels
site.counts <- site.counts[,c(7,2,1,5,6,9,4,3,8)] # order west to east 
site.counts
```

    ##            OR  ID  CO MN1 MN2  WI  MI ME1 ME2
    ## ITS count 138 109 172 199 168 127 140 160 142
    ## 16S count 144 120 180 220 179 150 150 180 150

#### table S1

Counts of ITS and 16S communities by site, year, and season.

``` r
its.counts <- ddply(data.frame(all.fung.1920.c.ps@sam_data), .(year, season, state),nrow) # sum up ITS counts by year, season, and state
colnames(its.counts)[4] <- "ITS count"

bact.counts <- ddply(data.frame(core.bact.ps@sam_data), .(year, season, state),nrow) # sum up 16S counts by year, season, and state
colnames(bact.counts)[4] <- "16S count"

# merge them together and make the table neat
both.counts <- left_join(bact.counts, its.counts, by = c("year", "season", "state")) %>% drop_na(state) # remove na states (blanks and tech reps)
both.counts$year <- paste(20, both.counts$year, sep = "") # add 20 to year
both.counts <- both.counts %>% arrange(state, year, factor(season, levels = c("Spring", "Summer", "Fall"))) # order by season
both.counts[,2] <- ifelse(both.counts[,2] == "Spring", "Pre-plant",
                          ifelse(both.counts[,2] == "Summer", "60 days after planting", 
                                 ifelse(both.counts[,2] == "Fall", "After harvest", NA))) # change season labels to harvest descriptions
both.counts[which(both.counts[,3] == "ME"),3] <- "ME1" # update site labels
both.counts[which(both.counts[,3] == "US"),3] <- "ME2"
both.counts[which(both.counts[,3] == "MN"),3] <- "MN1"
both.counts[which(both.counts[,3] == "ND"),3] <- "MN2"
colnames(both.counts)[3] <- "site"
both.counts
```

    ##    year                 season site 16S count ITS count
    ## 1  2019              Pre-plant   CO        30        29
    ## 2  2019 60 days after planting   CO        30        30
    ## 3  2019          After harvest   CO        30        27
    ## 4  2020              Pre-plant   CO        30        30
    ## 5  2020 60 days after planting   CO        30        30
    ## 6  2020          After harvest   CO        30        26
    ## 7  2019              Pre-plant   ID        24        21
    ## 8  2019 60 days after planting   ID        24        24
    ## 9  2020              Pre-plant   ID        24        22
    ## 10 2020 60 days after planting   ID        24        20
    ## 11 2020          After harvest   ID        24        22
    ## 12 2019              Pre-plant  ME1        30        29
    ## 13 2019 60 days after planting  ME1        30        30
    ## 14 2019          After harvest  ME1        30        25
    ## 15 2020              Pre-plant  ME1        30        18
    ## 16 2020 60 days after planting  ME1        30        30
    ## 17 2020          After harvest  ME1        30        28
    ## 18 2019              Pre-plant   MI        30        22
    ## 19 2019 60 days after planting   MI        30        29
    ## 20 2020              Pre-plant   MI        30        29
    ## 21 2020 60 days after planting   MI        30        30
    ## 22 2020          After harvest   MI        30        30
    ## 23 2018          After harvest  MN1        40        39
    ## 24 2019              Pre-plant  MN1        30        30
    ## 25 2019 60 days after planting  MN1        30        30
    ## 26 2019          After harvest  MN1        30        22
    ## 27 2020              Pre-plant  MN1        29        27
    ## 28 2020 60 days after planting  MN1        31        25
    ## 29 2020          After harvest  MN1        30        26
    ## 30 2019              Pre-plant  MN2        30        30
    ## 31 2019 60 days after planting  MN2        30        25
    ## 32 2019          After harvest  MN2        30        29
    ## 33 2020              Pre-plant  MN2        30        30
    ## 34 2020 60 days after planting  MN2        29        25
    ## 35 2020          After harvest  MN2        30        29
    ## 36 2019              Pre-plant   OR        24        24
    ## 37 2019 60 days after planting   OR        24        24
    ## 38 2019          After harvest   OR        24        23
    ## 39 2020              Pre-plant   OR        24        23
    ## 40 2020 60 days after planting   OR        24        24
    ## 41 2020          After harvest   OR        24        20
    ## 42 2019              Pre-plant  ME2        25        24
    ## 43 2019 60 days after planting  ME2        25        25
    ## 44 2019          After harvest  ME2        25        24
    ## 45 2020              Pre-plant  ME2        25        25
    ## 46 2020 60 days after planting  ME2        25        20
    ## 47 2020          After harvest  ME2        25        24
    ## 48 2019              Pre-plant   WI        30        28
    ## 49 2019 60 days after planting   WI        30        22
    ## 50 2020              Pre-plant   WI        29        28
    ## 51 2020 60 days after planting   WI        29        20
    ## 52 2020          After harvest   WI        32        29

We are missing: ID, MI, and WI fall 2019. And we have an extra MN fall
(2018).

#### table S2

Counts of 16S communities from each potato cultivar at each site.

``` r
cultivar.count <- ddply(data.frame(core.bact.ps@sam_data), .(cultivar, state),nrow) # get counts of cultivars
colnames(cultivar.count)[c(2,3)] <- c("site","count") # fix colnames

cultivar.count[which(cultivar.count[,2] == "ME"),2] <- "ME1" # update site labels
cultivar.count[which(cultivar.count[,2] == "US"),2] <- "ME2"
cultivar.count[which(cultivar.count[,2] == "MN"),2] <- "MN1"
cultivar.count[which(cultivar.count[,2] == "ND"),2] <- "MN2"

cultivar.count <- pivot_wider(cultivar.count, names_from = "cultivar", values_from = "count") # pivot wider, cultivars as columns
cultivar.count <- cultivar.count[1:9,1:10] # cut out NAs
cultivar.count <- data.frame(t(column_to_rownames(cultivar.count, "site"))) # put site as row name, transpose, and convert to df
cultivar.count <- cultivar.count[,c(7,4,3,1,2,8,6,5,9)] # place sites in west-east order
cultivar.count[is.na(cultivar.count)] <- 0 # turn NAs into zeros
cultivar.count
```

    ##                      OR  ID CO MN1 MN2 WI  MI ME1 ME2
    ## Bannock               0   0  0  37  29  0   0   0   0
    ## Burbank              22  20 15 113 150 76  25  90   0
    ## Burbank and Caribou   0   0  0   0   0  0   0   0 150
    ## Canela                0   0 75   0   0  0   0   0   0
    ## Caribou               0   0  0   0   0  0   0  90   0
    ## Norkotah            122 100 90  69   0  0   0   0   0
    ## Snowden               0   0  0   0   0  0  25   0   0
    ## Superior              0   0  0   0   0  0 100   0   0
    ## Yukon Gold            0   0  0   0   0 74   0   0   0

All sites used Burbank. ME2 samples were between rows of different
varieties. Norkotah is popular in the west, Bannock in MN, and Caribou
in ME. All other varieties were grown in one site only.

#### session info

``` r
sessionInfo()
```

    ## R version 4.2.2 (2022-10-31)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Ventura 13.2.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] plyr_1.8.8      phyloseq_1.42.0 forcats_1.0.0   stringr_1.5.0  
    ##  [5] dplyr_1.1.0     purrr_1.0.1     readr_2.1.4     tidyr_1.3.0    
    ##  [9] tibble_3.1.8    ggplot2_3.4.1   tidyverse_1.3.2
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-160           bitops_1.0-7           fs_1.6.1              
    ##  [4] lubridate_1.9.2        httr_1.4.4             GenomeInfoDb_1.34.9   
    ##  [7] tools_4.2.2            backports_1.4.1        vegan_2.6-4           
    ## [10] utf8_1.2.3             R6_2.5.1               mgcv_1.8-41           
    ## [13] DBI_1.1.3              BiocGenerics_0.44.0    colorspace_2.1-0      
    ## [16] permute_0.9-7          rhdf5filters_1.10.0    ade4_1.7-22           
    ## [19] withr_2.5.0            tidyselect_1.2.0       compiler_4.2.2        
    ## [22] cli_3.6.0              rvest_1.0.3            Biobase_2.58.0        
    ## [25] xml2_1.3.3             scales_1.2.1           digest_0.6.31         
    ## [28] rmarkdown_2.20         XVector_0.38.0         pkgconfig_2.0.3       
    ## [31] htmltools_0.5.4        dbplyr_2.3.0           fastmap_1.1.0         
    ## [34] rlang_1.0.6            readxl_1.4.2           rstudioapi_0.14       
    ## [37] generics_0.1.3         jsonlite_1.8.4         googlesheets4_1.0.1   
    ## [40] RCurl_1.98-1.10        magrittr_2.0.3         GenomeInfoDbData_1.2.9
    ## [43] biomformat_1.26.0      Matrix_1.5-1           Rcpp_1.0.10           
    ## [46] munsell_0.5.0          S4Vectors_0.36.1       Rhdf5lib_1.20.0       
    ## [49] fansi_1.0.4            ape_5.7                lifecycle_1.0.3       
    ## [52] stringi_1.7.12         yaml_2.3.7             MASS_7.3-58.1         
    ## [55] zlibbioc_1.44.0        rhdf5_2.42.0           grid_4.2.2            
    ## [58] parallel_4.2.2         crayon_1.5.2           lattice_0.20-45       
    ## [61] splines_4.2.2          Biostrings_2.66.0      haven_2.5.1           
    ## [64] multtest_2.54.0        hms_1.1.2              knitr_1.42            
    ## [67] pillar_1.8.1           igraph_1.4.0           reshape2_1.4.4        
    ## [70] codetools_0.2-18       stats4_4.2.2           reprex_2.0.2          
    ## [73] glue_1.6.2             evaluate_0.20          data.table_1.14.8     
    ## [76] modelr_0.1.10          vctrs_0.5.2            tzdb_0.3.0            
    ## [79] foreach_1.5.2          cellranger_1.1.0       gtable_0.3.1          
    ## [82] assertthat_0.2.1       xfun_0.37              broom_1.0.3           
    ## [85] survival_3.4-0         googledrive_2.0.0      gargle_1.3.0          
    ## [88] iterators_1.0.14       IRanges_2.32.0         cluster_2.1.4         
    ## [91] timechange_0.2.0       ellipsis_0.3.2
