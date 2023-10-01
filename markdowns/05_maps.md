maps and figure 1
================
Scott Klasek
2022-10-04

Purpose: To make Figure 1. Map of locations (A), and ordinations of
total communities (16S, B; ITS, C).

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.8     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.1
    ## ✔ readr   2.1.2     ✔ forcats 0.5.2
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(sf)
```

    ## Linking to GEOS 3.10.2, GDAL 3.4.2, PROJ 8.2.1; sf_use_s2() is TRUE

``` r
library(mapproj)
```

    ## Loading required package: maps
    ## 
    ## Attaching package: 'maps'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     map

``` r
library(maps)
library(ggrepel)
```

#### Map the field sites for Objective 1

Make a dataframe of the site names and their approximate lat-long
coordinates (not terribly precise but I’m hiding it anyway)

``` r
us.map <- map_data("state") # get the map of states
usmap <- ggplot(us.map, aes(long, lat, group = group)) +
  geom_polygon(fill = "white", colour = "grey50") # make a dataframe to map

cols.alphabetical <- c("#236192", "#B3A369", "#B0D7FF", "#18453B", "#7A0019", "#FFC72A", "#DC4405", "gray", "#C5050C")

# plotting 
# could use this instead for a mercator-type projection: coord_map("azequalarea", orientation = c(40, -100, 0), xlim = c(-123, -68), ylim = c(35, 49)) 
potato.map <- usmap + coord_map("azequalarea", orientation = c(40, -100, 0), xlim = c(-123, -68), ylim = c(37, 48))+
  scale_x_continuous("")+
  scale_y_continuous("")+
  geom_point(size = 2, fields, mapping = aes(x = long, y = lat, group = site), color = cols.alphabetical)+
  geom_label_repel(fields, mapping = aes(x = long, y = lat, group = site, label = site))+
  theme_void() # theme_void removes lat-long 
potato.map
```

![](14_maps_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

et voila, un map  
yes, ND site is actually in MN

``` r
sessionInfo()
```

    ## R version 4.2.0 (2022-04-22)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur/Monterey 10.16
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
    ##  [1] ggrepel_0.9.1   mapproj_1.2.8   maps_3.4.0      sf_1.0-8       
    ##  [5] forcats_0.5.2   stringr_1.4.1   dplyr_1.0.9     purrr_0.3.4    
    ##  [9] readr_2.1.2     tidyr_1.2.0     tibble_3.1.8    ggplot2_3.3.6  
    ## [13] tidyverse_1.3.2
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.9          lubridate_1.8.0     class_7.3-20       
    ##  [4] assertthat_0.2.1    digest_0.6.29       utf8_1.2.2         
    ##  [7] R6_2.5.1            cellranger_1.1.0    backports_1.4.1    
    ## [10] reprex_2.0.2        evaluate_0.16       e1071_1.7-11       
    ## [13] highr_0.9           httr_1.4.4          pillar_1.8.1       
    ## [16] rlang_1.0.4         googlesheets4_1.0.1 readxl_1.4.1       
    ## [19] rstudioapi_0.14     rmarkdown_2.15      labeling_0.4.2     
    ## [22] googledrive_2.0.0   munsell_0.5.0       proxy_0.4-27       
    ## [25] broom_1.0.0         compiler_4.2.0      modelr_0.1.9       
    ## [28] xfun_0.32           pkgconfig_2.0.3     htmltools_0.5.3    
    ## [31] tidyselect_1.1.2    fansi_1.0.3         crayon_1.5.1       
    ## [34] tzdb_0.3.0          dbplyr_2.2.1        withr_2.5.0        
    ## [37] grid_4.2.0          jsonlite_1.8.0      gtable_0.3.0       
    ## [40] lifecycle_1.0.1     DBI_1.1.3           magrittr_2.0.3     
    ## [43] units_0.8-0         scales_1.2.1        KernSmooth_2.23-20 
    ## [46] cli_3.3.0           stringi_1.7.8       farver_2.1.1       
    ## [49] fs_1.5.2            xml2_1.3.3          ellipsis_0.3.2     
    ## [52] generics_0.1.3      vctrs_0.4.1         tools_4.2.0        
    ## [55] glue_1.6.2          hms_1.1.2           fastmap_1.1.0      
    ## [58] yaml_2.3.5          colorspace_2.0-3    gargle_1.2.0       
    ## [61] classInt_0.4-8      rvest_1.0.3         knitr_1.40         
    ## [64] haven_2.5.1
