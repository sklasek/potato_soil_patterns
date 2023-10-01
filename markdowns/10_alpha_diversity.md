alpha diversity II
================
2023-02-27

#### purpose

to make some better alpha diversity plots now that we’re going to split
the paper up into three.

#### load libraries

``` r
packages <- c("tidyverse", "phyloseq", "speedyseq", "patchwork")
t <- lapply(packages, require, character.only = TRUE) # load all packages at once
```

    ## Loading required package: tidyverse

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.2     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.2     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.1     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
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

#### load phyloseq data and clean it

Bacterial data didn’t have ASVs named, or samples pruned to 10k reads.
This below gives me a cleaned bacterial phyloseq object that I didn’t
have before.

``` r
fung.ps <- readRDS(file = "/Users/klas0061/Desktop/UMN/phyloseqs/first_two_years_all_sites/fung1920.cleaned.ps") # cleaned fungal phyloseq object
bact.ps <- readRDS(file = "/Users/klas0061/Desktop/UMN/phyloseqs/first_two_years_all_sites/all.bact1920.ps") # "unclean" bacterial phyloseq object

#### clean our bacterial object
# assign ASV numbers cuz i'm sick of the sequences as the rownames
number_asvs <- function(ps){
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
  return(ps)
}
bact.ps <- number_asvs(bact.ps)

# prune samples below 10k reads
bact.ps <- subset_samples(bact.ps, sample_sums(bact.ps) > 10000) 

# write out phyloseq object
#saveRDS(bact.ps, "/Users/klas0061/Desktop/UMN/phyloseqs/bact1920.cleaned.ps") # cleaned fungal phyloseq object
```

#### rarefy to 10k reads

We’re doing this for alpha diversity only, because shannon and
invsimpson indices had small (but significant) increases with library
size.

``` r
bact.rps <- rarefy_even_depth(bact.ps, sample.size = 10000, rngseed = 111) # setting a seed for reproducibility
```

    ## `set.seed(111)` was used to initialize repeatable random subsampling.

    ## Please record this for your records so others can reproduce.

    ## Try `set.seed(111); .Random.seed` for the full vector

    ## ...

    ## 128883OTUs were removed because they are no longer 
    ## present in any sample after random subsampling

    ## ...

``` r
euks.rps <- rarefy_even_depth(fung.ps, sample.size = 10000, rngseed = 111) 
```

    ## `set.seed(111)` was used to initialize repeatable random subsampling.

    ## Please record this for your records so others can reproduce.

    ## Try `set.seed(111); .Random.seed` for the full vector

    ## ...

    ## 1207OTUs were removed because they are no longer 
    ## present in any sample after random subsampling

    ## ...

#### calculate alpha diversity

``` r
# eukaryotes first
adiv.euks <- estimate_richness(euks.rps, measures = c("Observed", "Shannon", "InvSimpson")) # calculate these three metrics
adiv.euks$sample_name <- str_sub(rownames(adiv.euks), 1, 15) # get the sample name from the first 15 chars of the rowname
colnames(adiv.euks)[1:3] <- c("Observed.euk", "Shannon.euk", "InvSimpson.euk") # specify domain in column name

# bacteria next
adiv.bact <- estimate_richness(bact.rps, measures = c("Observed", "Shannon", "InvSimpson")) # calculate these three metrics
adiv.bact$sample_name <- str_sub(rownames(adiv.bact), 1, 15) # get the sample name from the first 15 chars of the rowname
colnames(adiv.bact)[1:3] <- c("Observed.bact", "Shannon.bact", "InvSimpson.bact") # specify domain in column name

# join them together
adiv.r <- inner_join(adiv.bact, adiv.euks) %>% drop_na() # bind bacterial and eukaryotic data by sample names
```

    ## Joining with `by = join_by(sample_name)`

    ## Warning in inner_join(adiv.bact, adiv.euks): Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 1025 of `x` matches multiple rows in `y`.
    ## ℹ Row 974 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

``` r
adiv.r$state <- str_sub(adiv.r$sample_name, 1, 2) # add state info
adiv.r[which(adiv.r$state == "Co"),"state"] <- "CO"
adiv.r[which(adiv.r$state == "ME"),"state"] <- "ME1"
adiv.r[which(adiv.r$state == "US"),"state"] <- "ME2"
adiv.r[which(adiv.r$state == "MN"),"state"] <- "MN1"
adiv.r[which(adiv.r$state == "ND"),"state"] <- "MN2"
```

#### plot bacterial and eukaryotic alpha-diversity together

``` r
# plot 
cols.alphabetical <- c("#236192", "#B3A369", "#B0D7FF", "gray", "#18453B", "#7A0019", "#FFC72A", "#DC4405", "#C5050C")

# plotting bacterial and eukaryotic data together
both.observed <- ggplot(adiv.r, aes(Observed.bact, Observed.euk, color = state))+
  geom_point()+
  scale_color_manual("Site", values = cols.alphabetical)+
  scale_x_continuous("Observed Bacterial ASVs")+scale_y_continuous("Observed Eukaryotic ASVs")+
  facet_wrap(~state)+
  theme_bw()

both.shannon <- ggplot(adiv.r, aes(Shannon.bact, Shannon.euk, color = state))+
  geom_point()+
  scale_x_continuous("Shannon, Bacteria")+scale_y_continuous("Shannon, Eukaryotes")+
  scale_color_manual("Site", values = cols.alphabetical)+
  facet_wrap(~state)+
  theme_bw()

both.invsimpson <- ggplot(adiv.r, aes(InvSimpson.bact, InvSimpson.euk, color = state))+
  geom_point()+
  scale_x_continuous("Inverse Simpson, Bacteria")+scale_y_continuous("Inverse Simpson, Eukaryotes")+
  scale_color_manual("Site", values = cols.alphabetical)+
  facet_wrap(~state)+
  theme_bw()

both.observed + both.shannon + both.invsimpson + plot_layout(guides = "collect")
```

![](25_alpha_diversity_II_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
# linear models
summary(lm(Observed.euk ~ state + Observed.bact, adiv.r))
```

    ## 
    ## Call:
    ## lm(formula = Observed.euk ~ state + Observed.bact, data = adiv.r)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -294.789  -54.381    2.118   56.542  311.469 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   -2.517e+02  2.984e+01  -8.435  < 2e-16 ***
    ## stateID        9.020e+01  1.144e+01   7.882 6.52e-15 ***
    ## stateME1       4.151e+02  9.635e+00  43.083  < 2e-16 ***
    ## stateME2       3.017e+02  1.088e+01  27.722  < 2e-16 ***
    ## stateMI        2.793e+02  1.043e+01  26.786  < 2e-16 ***
    ## stateMN1       9.531e+01  1.031e+01   9.243  < 2e-16 ***
    ## stateMN2       1.274e+02  9.600e+00  13.273  < 2e-16 ***
    ## stateOR        2.901e+02  9.103e+00  31.866  < 2e-16 ***
    ## stateWI        8.096e+01  1.036e+01   7.816 1.07e-14 ***
    ## Observed.bact  1.699e-01  8.227e-03  20.647  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 85.67 on 1378 degrees of freedom
    ## Multiple R-squared:  0.7111, Adjusted R-squared:  0.7093 
    ## F-statistic: 376.9 on 9 and 1378 DF,  p-value: < 2.2e-16

``` r
summary(lm(Shannon.euk ~ state + Shannon.bact, adiv.r))
```

    ## 
    ## Call:
    ## lm(formula = Shannon.euk ~ state + Shannon.bact, data = adiv.r)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.65053 -0.15366  0.03308  0.18669  1.22700 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  -6.88101    0.56251 -12.233  < 2e-16 ***
    ## stateID       0.42602    0.04793   8.889  < 2e-16 ***
    ## stateME1      1.21646    0.03987  30.512  < 2e-16 ***
    ## stateME2      0.79696    0.04526  17.609  < 2e-16 ***
    ## stateMI       0.56285    0.04319  13.031  < 2e-16 ***
    ## stateMN1      0.28454    0.04243   6.707 2.90e-11 ***
    ## stateMN2      0.26858    0.04107   6.540 8.65e-11 ***
    ## stateOR       0.86900    0.03821  22.743  < 2e-16 ***
    ## stateWI      -0.10167    0.04280  -2.376   0.0177 *  
    ## Shannon.bact  1.44067    0.07398  19.473  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3596 on 1378 degrees of freedom
    ## Multiple R-squared:  0.601,  Adjusted R-squared:  0.5984 
    ## F-statistic: 230.6 on 9 and 1378 DF,  p-value: < 2.2e-16

``` r
summary(lm(InvSimpson.euk ~ state + InvSimpson.bact, adiv.r))
```

    ## 
    ## Call:
    ## lm(formula = InvSimpson.euk ~ state + InvSimpson.bact, data = adiv.r)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -53.308 -11.798  -0.956   8.981  86.247 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      4.973392   2.150714   2.312 0.020900 *  
    ## stateID         21.336331   2.418458   8.822  < 2e-16 ***
    ## stateME1        42.540421   2.073054  20.521  < 2e-16 ***
    ## stateME2        33.462780   2.326203  14.385  < 2e-16 ***
    ## stateMI         10.641234   2.231140   4.769 2.04e-06 ***
    ## stateMN1         7.693821   2.135946   3.602 0.000327 ***
    ## stateMN2         3.557290   2.149181   1.655 0.098115 .  
    ## stateOR         34.908557   1.980914  17.622  < 2e-16 ***
    ## stateWI          1.613369   2.195699   0.735 0.462595    
    ## InvSimpson.bact  0.025410   0.002058  12.347  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 18.63 on 1378 degrees of freedom
    ## Multiple R-squared:  0.4651, Adjusted R-squared:  0.4616 
    ## F-statistic: 133.1 on 9 and 1378 DF,  p-value: < 2.2e-16

Interestingly, there are positive correlations

#### Plot them separately

``` r
# bacteria- pivot longer
adiv.bact2 <- adiv.r %>% select(Observed.bact, Shannon.bact, InvSimpson.bact, state) 
colnames(adiv.bact2)[1:3] <- c("Observed", "Shannon", "InvSimpson")
adiv.bact2 <- adiv.bact2 %>% pivot_longer(!state, names_to = "index", values_to = "value")
adiv.bact2$state <- factor(adiv.bact2$state, levels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))

# euks- same
adiv.euks2 <- adiv.r %>% select(Observed.euk, Shannon.euk, InvSimpson.euk, state) 
colnames(adiv.euks2)[1:3] <- c("Observed", "Shannon", "InvSimpson")
adiv.euks2 <- adiv.euks2 %>% pivot_longer(!state, names_to = "index", values_to = "value")
adiv.euks2$state <- factor(adiv.euks2$state, levels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))

# plot
cols.factor.order <- c("#DC4405", "#B3A369", "#236192", "#7A0019", "#FFC72A", "#C5050C","#18453B", "#B0D7FF", "gray")
allbact <- ggplot(adiv.bact2, aes(state, value, color = state))+
  geom_jitter(width = 0.3, alpha = 0.5)+
  facet_wrap(~factor(index, levels = c("Observed", "Shannon", "InvSimpson")), scales = "free")+
  scale_x_discrete("")+scale_y_continuous("Bacterial alpha diversity")+
  scale_color_manual("Site", values = cols.factor.order)+
  theme_bw()+ggtitle("A")+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

alleuks <- ggplot(adiv.euks2, aes(state, value, color = state))+
  geom_jitter(width = 0.3, alpha = 0.5)+
  facet_wrap(~factor(index, levels = c("Observed", "Shannon", "InvSimpson")), scales = "free")+
  scale_x_discrete("")+scale_y_continuous("Eukaryotic alpha diversity")+
  scale_color_manual("Site", values = cols.factor.order)+
  theme_bw()+ggtitle("B")+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

allbact / alleuks + plot_layout(guides = "collect")
```

![](25_alpha_diversity_II_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Takeaways from alpha diversity: Both bacterial and eukaryotic alpha
diversity metrics vary by state, but not in the same way for both
communities (For example, ID has high bacterial diversity but medium
eukaryotic diversity). There is a lot of spread within communities of
the same state. There are positive relationships between bacterial and
eukaryotic diversity across Observed and Shannon indices, not as much
the case in Inv Simpson.

#### anovas

``` r
# bacterial diversity by state
summary(aov(value ~ state, data = adiv.bact2 %>% filter(index == "Observed")))
```

    ##               Df    Sum Sq  Mean Sq F value Pr(>F)    
    ## state          8 114486481 14310810     182 <2e-16 ***
    ## Residuals   1379 108445081    78640                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(value ~ state, data = adiv.bact2 %>% filter(index == "Shannon")))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)    
    ## state          8  19.18  2.3971     140 <2e-16 ***
    ## Residuals   1379  23.62  0.0171                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(value ~ state, data = adiv.bact2 %>% filter(index == "InvSimpson")))
```

    ##               Df   Sum Sq Mean Sq F value Pr(>F)    
    ## state          8 44005675 5500709   92.52 <2e-16 ***
    ## Residuals   1379 81990727   59457                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# eukaryotic
summary(aov(value ~ state, data = adiv.euks2 %>% filter(index == "Observed")))
```

    ##               Df   Sum Sq Mean Sq F value Pr(>F)    
    ## state          8 21772822 2721603   283.4 <2e-16 ***
    ## Residuals   1379 13243680    9604                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(value ~ state, data = adiv.euks2 %>% filter(index == "Shannon")))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)    
    ## state          8  219.3  27.413   166.4 <2e-16 ***
    ## Residuals   1379  227.2   0.165                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(value ~ state, data = adiv.euks2 %>% filter(index == "InvSimpson")))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)    
    ## state          8 363039   45380   117.8 <2e-16 ***
    ## Residuals   1379 531441     385                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

While we see considerable within-site variation in alpha diversity,
site-specific variance is much higher.

#### session info

``` r
sessionInfo()
```

    ## R version 4.3.0 (2023-04-21)
    ## Platform: x86_64-apple-darwin20 (64-bit)
    ## Running under: macOS Ventura 13.3.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] patchwork_1.1.2      speedyseq_0.5.3.9018 phyloseq_1.44.0     
    ##  [4] lubridate_1.9.2      forcats_1.0.0        stringr_1.5.0       
    ##  [7] dplyr_1.1.2          purrr_1.0.1          readr_2.1.4         
    ## [10] tidyr_1.3.0          tibble_3.2.1         ggplot2_3.4.2       
    ## [13] tidyverse_2.0.0     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] ade4_1.7-22             tidyselect_1.2.0        farver_2.1.1           
    ##  [4] Biostrings_2.68.1       bitops_1.0-7            fastmap_1.1.1          
    ##  [7] RCurl_1.98-1.12         digest_0.6.31           timechange_0.2.0       
    ## [10] lifecycle_1.0.3         cluster_2.1.4           survival_3.5-5         
    ## [13] magrittr_2.0.3          compiler_4.3.0          rlang_1.1.1            
    ## [16] tools_4.3.0             igraph_1.4.3            utf8_1.2.3             
    ## [19] yaml_2.3.7              data.table_1.14.8       knitr_1.42             
    ## [22] labeling_0.4.2          plyr_1.8.8              withr_2.5.0            
    ## [25] BiocGenerics_0.46.0     grid_4.3.0              stats4_4.3.0           
    ## [28] fansi_1.0.4             multtest_2.56.0         biomformat_1.28.0      
    ## [31] colorspace_2.1-0        Rhdf5lib_1.22.0         scales_1.2.1           
    ## [34] iterators_1.0.14        MASS_7.3-58.4           cli_3.6.1              
    ## [37] rmarkdown_2.21          vegan_2.6-4             crayon_1.5.2           
    ## [40] generics_0.1.3          rstudioapi_0.14         reshape2_1.4.4         
    ## [43] tzdb_0.4.0              ape_5.7-1               rhdf5_2.44.0           
    ## [46] zlibbioc_1.46.0         splines_4.3.0           parallel_4.3.0         
    ## [49] XVector_0.40.0          vctrs_0.6.2             Matrix_1.5-4           
    ## [52] jsonlite_1.8.4          IRanges_2.34.0          hms_1.1.3              
    ## [55] S4Vectors_0.38.1        foreach_1.5.2           glue_1.6.2             
    ## [58] codetools_0.2-19        stringi_1.7.12          gtable_0.3.3           
    ## [61] GenomeInfoDb_1.36.0     munsell_0.5.0           pillar_1.9.0           
    ## [64] htmltools_0.5.5         rhdf5filters_1.12.1     GenomeInfoDbData_1.2.10
    ## [67] R6_2.5.1                evaluate_0.21           lattice_0.21-8         
    ## [70] Biobase_2.60.0          highr_0.10              Rcpp_1.0.10            
    ## [73] nlme_3.1-162            permute_0.9-7           mgcv_1.8-42            
    ## [76] xfun_0.39               pkgconfig_2.0.3
