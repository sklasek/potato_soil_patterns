Soil data and how it varies across sites
================
Scott Klasek
2022-12-19

#### purpose

To see how different soil data varies across states and develop
hypotheses about how it could be driving trends in microbiomes across
sites.

#### load libraries

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.1     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.2     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.1     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
# library(soiltexture)
library(ggtern)
```

    ## Registered S3 methods overwritten by 'ggtern':
    ##   method           from   
    ##   grid.draw.ggplot ggplot2
    ##   plot.ggplot      ggplot2
    ##   print.ggplot     ggplot2
    ## --
    ## Remember to cite, run citation(package = 'ggtern') for further info.
    ## --
    ## 
    ## Attaching package: 'ggtern'
    ## 
    ## The following objects are masked from 'package:ggplot2':
    ## 
    ##     aes, annotate, ggplot, ggplot_build, ggplot_gtable, ggplotGrob,
    ##     ggsave, layer_data, theme_bw, theme_classic, theme_dark,
    ##     theme_gray, theme_light, theme_linedraw, theme_minimal, theme_void

``` r
library(patchwork)
```

#### load and combine data

``` r
soil <- read.csv(file = "/Users/klas0061/Desktop/UMN/jim_info/soil_characteristics_jim_2023_06_13.csv", na.strings = c(".")) # convert . to NA
table(soil$State) # no USDA ME2 site
```

    ## 
    ##  CO  ID  ME  MI  MN  ND  OR  WI 
    ## 360 288 360 360 360 360 288 360

``` r
larkin.data <- read.csv(file = "/Users/klas0061/Desktop/UMN/larkin_info/chem_data_first_two_years.csv") # ME2 data that I could find (this csv was manually curated)

soil <- bind_rows(soil, larkin.data) # bind in ME2 data with the rest of the site data
```

#### plotting edaphic factors that vary noticeably across sites

Note: I tried to write these ggplots as functions, using embracing to
specify each variable to plot. But it didn’t work, because REASONS.

``` r
cols.alphabetical.east.west <- c("#DC4405", "#B3A369", "#236192", "#7A0019", "#FFC72A", "#C5050C", "#18453B", "#B0D7FF", "gray") # color scale

# plotting pH in samples from the first two years
ggph <- ggplot(soil %>% filter(Year < 21), aes(State, pH, color = factor(State, levels = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))))+
  geom_jitter()+
  scale_color_manual("Site", values = cols.alphabetical.east.west)+
  scale_x_discrete("", limits = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))+
  theme_bw()+ggtitle("A")+theme(legend.position = "none")

# plotting organic matter and C in samples from the first two years
ggom <- ggplot(soil %>% filter(Year < 21), aes(State, OM...., color = factor(State, levels = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))))+
  geom_jitter()+
  scale_color_manual("Site", values = cols.alphabetical.east.west)+
  scale_x_discrete("", limits = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))+
  theme_bw()+scale_y_continuous("% Organic matter")+ggtitle("B")+theme(legend.position = "none")

ggtc <- ggplot(soil %>% filter(Year < 21), aes(State, Total.C...., color = factor(State, levels = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))))+
  geom_jitter()+
  scale_color_manual("Site", values = cols.alphabetical.east.west)+
  scale_x_discrete("", limits = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))+
  theme_bw()+scale_y_continuous("% total C")+ggtitle("C")+theme(legend.position = "none")

# ggom + ggtc + plot_layout(guides = "collect") 

# NPK
ggn <- ggplot(soil %>% filter(Year < 21), aes(State, Nitrate.N..ppm., color = factor(State, levels = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))))+
  geom_jitter()+
  scale_color_manual("Site", values = cols.alphabetical.east.west)+
  scale_x_discrete("", limits = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))+
  theme_bw()+scale_y_continuous("Nitrate (ppm)")+ggtitle("D")+theme(legend.position = "none")

ggp <- ggplot(soil %>% filter(Year < 21), aes(State, P.Bray..ppm., color = factor(State, levels = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))))+
  geom_jitter()+
  scale_color_manual("Site", values = cols.alphabetical.east.west)+
  scale_x_discrete("", limits = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))+
  theme_bw()+scale_y_continuous("P (ppm, Bray)")+ggtitle("E")+theme(legend.position = "none")

ggk <- ggplot(soil %>% filter(Year < 21), aes(State, K..ppm., color = factor(State, levels = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))))+
  geom_jitter()+
  scale_color_manual("Site", values = cols.alphabetical.east.west)+
  scale_x_discrete("", limits = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))+
  theme_bw()+scale_y_continuous("K (ppm)")+ggtitle("F")+theme(legend.position = "none")

# CEC 
ggcec <- ggplot(soil %>% filter(Year < 21), aes(State, CEC..meq.100g., color = factor(State, levels = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))))+
  geom_jitter()+
  scale_color_manual("Site", values = cols.alphabetical.east.west)+
  scale_x_discrete("", limits = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))+
  theme_bw()+scale_y_continuous("Cation Exchange Capacity \n(mEq/100g)")+ggtitle("G")+theme(legend.position = "none")


ggph + ggom + ggtc + ggn + ggp + ggk + ggcec
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

    ## Warning: Removed 481 rows containing missing values (`geom_point()`).
    ## Removed 481 rows containing missing values (`geom_point()`).

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).
    ## Removed 1 rows containing missing values (`geom_point()`).
    ## Removed 1 rows containing missing values (`geom_point()`).
    ## Removed 1 rows containing missing values (`geom_point()`).

![](16_edaphics_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

CO and ID clearly have higher pH than the rest, ME has highest OM.

Other noteworthy trends:

ME is a bit of an outlier, with highest OM and total C, also high
nitrate and P.

CO has high K and nitrate.

OR has high K.

MN and ND (both MN really) all have way higher total C than WI does.
Microbial communities from these sites all cluster together.

#### Make a ternary plot of soil texture classifications

Help from [this resource
here](https://saryace.github.io/flipbook_soiltexture_en/#37).

``` r
soil.texture <- soil %>% select(State, Sand.., Silt.., Clay..) %>% drop_na() # trim dataframe and remove NAs 
colnames(soil.texture)[2:4] <- c("Sand", "Silt", "Clay") # clean up column names

theme_set(theme_bw()) # call a theme somehow
data(USDA) # load USDA soil classification polygon boundaries
USDA_text <- USDA  %>% group_by(Label) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) # 

# plot with USDA classifications
tplot <- ggplot(data = USDA, aes(y = Clay, x = Sand, z = Silt)) +
  coord_tern(L = "x", T = "y", R = "z") +
  geom_polygon(
    aes(fill = Label), alpha = 0.0, size = 0.5,
    color = "black") +
  geom_text(data = USDA_text, aes(label = Label), color = 'black', size = 2) +
  geom_point(data = soil.texture, aes(x = Sand, y = Clay, z = Silt, color = State, alpha = 0.2))+
  scale_color_manual("Site", limits = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"), values = cols.alphabetical.east.west)+guides(fill="none", alpha = "none")
```

    ## Coordinate system already present. Adding new coordinate system, which will
    ## replace the existing one.

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning in geom_point(data = soil.texture, aes(x = Sand, y = Clay, z = Silt, :
    ## Ignoring unknown aesthetics: z

``` r
tplot
```

![](16_edaphics_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

MI has a lot of textural variability. It doesn’t seem to impact fungal
beta-diversity, but I wonder if this explains the high bacterial
beta-diversity.

#### session info

``` r
sessionInfo()
```

    ## R version 4.2.3 (2023-03-15)
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
    ##  [1] patchwork_1.1.2 ggtern_3.4.2    lubridate_1.9.2 forcats_1.0.0  
    ##  [5] stringr_1.5.0   dplyr_1.1.1     purrr_1.0.1     readr_2.1.4    
    ##  [9] tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.2   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.0   xfun_0.38          lattice_0.20-45    latex2exp_0.9.6   
    ##  [5] colorspace_2.1-0   vctrs_0.6.1        generics_0.1.3     htmltools_0.5.5   
    ##  [9] yaml_2.3.7         compositions_2.0-6 utf8_1.2.3         rlang_1.1.0       
    ## [13] hexbin_1.28.3      pillar_1.9.0       glue_1.6.2         withr_2.5.0       
    ## [17] lifecycle_1.0.3    plyr_1.8.8         robustbase_0.95-1  munsell_0.5.0     
    ## [21] gtable_0.3.3       evaluate_0.20      labeling_0.4.2     knitr_1.42        
    ## [25] tzdb_0.3.0         fastmap_1.1.1      fansi_1.0.4        highr_0.10        
    ## [29] DEoptimR_1.0-13    proto_1.0.0        Rcpp_1.0.10        scales_1.2.1      
    ## [33] farver_2.1.1       gridExtra_2.3      tensorA_0.36.2     hms_1.1.3         
    ## [37] digest_0.6.31      stringi_1.7.12     grid_4.2.3         cli_3.6.1         
    ## [41] tools_4.2.3        magrittr_2.0.3     crayon_1.5.2       pkgconfig_2.0.3   
    ## [45] MASS_7.3-58.2      bayesm_3.1-5       timechange_0.2.0   rmarkdown_2.21    
    ## [49] rstudioapi_0.14    R6_2.5.1           compiler_4.2.3
