Reviewables
================
Scott Klasek
2023-09-13

### Purpose

Revisit some data relating to reviewer comments on Potato Soil Core
Microbiomes manuscript (first two years).

#### load libraries

``` r
packages <- c("tidyverse", "phyloseq", "speedyseq", "patchwork", "NatParksPalettes", "UpSetR", "scales", "grid", "gridExtra", "cowplot")
t <- lapply(packages, require, character.only = TRUE) # load all packages at once
```

    ## Loading required package: tidyverse

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
    ## 
    ## 
    ## Loading required package: grid
    ## 
    ## Loading required package: gridExtra
    ## 
    ## 
    ## Attaching package: 'gridExtra'
    ## 
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine
    ## 
    ## 
    ## Loading required package: cowplot
    ## 
    ## 
    ## Attaching package: 'cowplot'
    ## 
    ## 
    ## The following object is masked from 'package:patchwork':
    ## 
    ##     align_plots
    ## 
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     stamp

#### load data

``` r
jimdata <- read_csv(file = "/Users/klas0061/Desktop/UMN/jim_info/PSHP_ALL_obj_1_all_data_2023_06_13.csv")
```

    ## New names:
    ## Rows: 2736 Columns: 142
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (133): State, Season, Block, Cultivar, End market, Seed piece size (oz),... dbl
    ## (9): Objective, Rotation, Plot, Year...5, Year within rotation, Month,...
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `Year` -> `Year...5`
    ## • `Year` -> `Year...102`
    ## • `Green plant cover before vine kill (%)` -> `Green plant cover before vine
    ##   kill (%)...103`
    ## • `Culls` -> `Culls...104`
    ## • `0 - 4 oz` -> `0 - 4 oz...105`
    ## • `4 - 6 oz.` -> `4 - 6 oz....106`
    ## • `6 - 10 oz` -> `6 - 10 oz...107`
    ## • `10 - 14 oz` -> `10 - 14 oz...108`
    ## • `>14 oz.` -> `>14 oz....109`
    ## • `Total yield` -> `Total yield...110`
    ## • `U.S. No. 1` -> `U.S. No. 1...111`
    ## • `U.S. No.  2` -> `U.S. No.  2...112`
    ## • `Marketable yield` -> `Marketable yield...113`
    ## • `Greater than 6 oz.` -> `Greater than 6 oz....114`
    ## • `Greater than 10 oz.` -> `Greater than 10 oz....115`
    ## • `Specific gravity` -> `Specific gravity...116`
    ## • `3- to 5-type scab (% of tubers)` -> `3- to 5-type scab (% of tubers)...117`
    ## • `Verticillium vascular ring (% of tubers)` -> `Verticillium vascular ring (%
    ##   of tubers)...118`
    ## • `Disqualifying hollow heart (% of tubers)` -> `Disqualifying hollow heart (%
    ##   of tubers)...119`
    ## • `Mean scab cover overall (% of surface)` -> `Mean scab cover overall (% of
    ##   surface)...120`
    ## • `Mean scab cover among infected (%)` -> `Mean scab cover among infected
    ##   (%)...121`
    ## • `Year` -> `Year...122`
    ## • `Green plant cover before vine kill (%)` -> `Green plant cover before vine
    ##   kill (%)...123`
    ## • `Culls` -> `Culls...124`
    ## • `0 - 4 oz` -> `0 - 4 oz...125`
    ## • `4 - 6 oz.` -> `4 - 6 oz....126`
    ## • `6 - 10 oz` -> `6 - 10 oz...127`
    ## • `10 - 14 oz` -> `10 - 14 oz...128`
    ## • `>14 oz.` -> `>14 oz....129`
    ## • `Total yield` -> `Total yield...130`
    ## • `U.S. No. 1` -> `U.S. No. 1...131`
    ## • `U.S. No.  2` -> `U.S. No.  2...132`
    ## • `Marketable yield` -> `Marketable yield...133`
    ## • `Greater than 6 oz.` -> `Greater than 6 oz....134`
    ## • `Greater than 10 oz.` -> `Greater than 10 oz....135`
    ## • `Specific gravity` -> `Specific gravity...136`
    ## • `3- to 5-type scab (% of tubers)` -> `3- to 5-type scab (% of tubers)...137`
    ## • `Verticillium vascular ring (% of tubers)` -> `Verticillium vascular ring (%
    ##   of tubers)...138`
    ## • `Disqualifying hollow heart (% of tubers)` -> `Disqualifying hollow heart (%
    ##   of tubers)...139`
    ## • `Mean scab cover overall (% of surface)` -> `Mean scab cover overall (% of
    ##   surface)...140`
    ## • `Mean scab cover among infected (%)` -> `Mean scab cover among infected
    ##   (%)...141`

``` r
view(jimdata)

# entire communities
bact.ps <- readRDS(file = "/Users/klas0061/Desktop/UMN/phyloseqs/first_two_years_all_sites/bact1920.cleaned.ps")
euks.ps <- readRDS(file = "/Users/klas0061/Desktop/UMN/phyloseqs/first_two_years_all_sites/fung1920.cleaned.ps")

# core communities - last col of otu table is "summing", all its tax categories are also "summing"
bact.core.ps <- readRDS(file = "/Users/klas0061/Desktop/UMN/phyloseqs/core/core.bact.asvs.ps")
euks.core.ps <- readRDS(file = "/Users/klas0061/Desktop/UMN/phyloseqs/core/core.its.asvs.ps")

# asvcumulate results with number of states the ASVs were observed in 
bact.asvc <- read_csv(file = "/Users/klas0061/Desktop/UMN/figures/from_server/16s.asvcumulate.csv")
```

    ## New names:
    ## Rows: 402273 Columns: 8
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (1): Row.names dbl (7): ...1, reads, pctreads, pct_prevalence, csum, rank,
    ## num_states_per_a...
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
euks.asvc <- read_csv(file = "/Users/klas0061/Desktop/UMN/figures/from_server/ITS.ASVs.asvc")
```

    ## New names:
    ## Rows: 22035 Columns: 8
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (1): Row.names dbl (7): ...1, reads, pctreads, pct_prevalence, csum, rank,
    ## num_states_per_a...
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

#### define functions

``` r
# merge.and.tidy.samples takes a ps object, groups samples by replicates according to phyloseq::merge_samples, 
# then calculate percent abundances, and fix state labels
merge.and.tidy.samples <- function(ps){

  # assign repgroup for each set of communities that can be considered a "replicate"
  sample.df <- data.frame(ps@sam_data)
  sample.df.us <- sample.df %>% filter(state=="US")
  sample.df.rest <- sample.df %>% filter(state!="US")
  
  # previously added sample.df.us$treatment_new to separate US samples by treatment, but found it didn't matter much and could be an unnecessary distraction/tangent. 
  sample.df.us$repgroup <- paste(sample.df.us$state, sample.df.us$year, sample.df.us$season, sep = "") 
  sample.df.rest$repgroup <- paste(sample.df.rest$state, sample.df.rest$year, sample.df.rest$season, sep = "")
  sample.df <- rbind(sample.df.rest, sample.df.us)
  sample_data(ps) <- sample.df
  
  # merge samples by repgroup
  ps.grouped <- merge_samples(ps, "repgroup")
  
  # re-add state names because merge_samples decided to omit them all
  sample_data(ps.grouped)$state <- str_sub(rownames(sample_data(ps.grouped)), 1, 2)
  
  # transform to percent abundances
  ps.grouped <- transform_sample_counts(ps.grouped, function(OTU) 100*OTU/sum(OTU))
  
  # redo state factoring
  ps.grouped@sam_data$state <- factor(ps.grouped@sam_data$state, levels = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"),
                                      labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))
  return(ps.grouped)
}

# remove summing columns from core ps objects
prune_sums <- function(ps){
  keep.taxa.names <- taxa_names(ps)[which(taxa_names(ps) != "summing")]
  prune.ps <- prune_taxa(keep.taxa.names, ps)
  return(prune.ps)
}

# prune all but most abundant n classes
take.top.classes <- function(ps, nclasses){
  asvs.of.classes.descending <- names(sort(colSums(ps@otu_table) / nrow(ps@otu_table), decreasing = TRUE)) # the representative ASVs of classes
  cla <- data.frame(ps@tax_table)[asvs.of.classes.descending,]$Class # their corresponding classes, in order of descending abundance
  cat("these classes are selected: ", cla[1:nclasses], "\n") # lists the classes 
  ps.pruned <- prune_taxa(asvs.of.classes.descending[1:nclasses], ps) # prunes all but the top n classes from the ps object
  return(ps.pruned)
}
```

### Do other soil chemical, physical measurements vary by state, and if so, are they worth mentioning?

``` r
jimdata %>% 
  select(State, `Ca (ppm)`, `Mg (ppm)`, `Na (ppm)`, `Zn (ppm)`) %>% 
  mutate_at(c("Ca (ppm)", "Mg (ppm)", "Na (ppm)", "Zn (ppm)"), as.numeric) %>% 
  pivot_longer(!State, values_to = "value", names_to = "chem") %>% 
  ggplot(aes(State, value))+geom_violin()+facet_grid(chem~., scales = "free")
```

    ## Warning: There were 4 warnings in `mutate()`.
    ## The first warning was:
    ## ℹ In argument: `Ca (ppm) = .Primitive("as.double")(`Ca (ppm)`)`.
    ## Caused by warning:
    ## ! NAs introduced by coercion
    ## ℹ Run `dplyr::last_dplyr_warnings()` to see the 3 remaining warnings.

    ## Warning: Removed 28 rows containing non-finite values (`stat_ydensity()`).

![](32_reviewables_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### You want Class? I’ll give you Class

As in, plotting the top Class-level taxonomy for the full and the core
16S and ITS microbiomes.

``` r
# make a list of phyloseq objects
ps.list <- list(bact.ps, euks.ps, bact.core.ps, euks.core.ps)
names(ps.list) <- c("all.bact", "all.euks", "core.bact", "core.euks")

# tax-glom to class, transform to percent abundance, omit low-abundance phyla
classglom.ps.list <- lapply(ps.list, tax_glom, taxrank = "Class")

# run this to merge samples by replicates, calculate percent abundances, and fix state labels
mtcglom.ps.list <- suppressWarnings(lapply(classglom.ps.list, merge.and.tidy.samples)) # Warning: NAs introduced by coercion but its ok

# remove "summing" columns from core phyloseq objects
core.bact.pruned <- prune_sums(mtcglom.ps.list[["core.bact"]])
core.euks.pruned <- prune_sums(mtcglom.ps.list[["core.euks"]])

# prune all but most n abundant Classes
top.c.all.bact <- take.top.classes(mtcglom.ps.list[["all.bact"]], 10)
```

    ## these classes are selected:  Alphaproteobacteria Actinobacteria Gammaproteobacteria Thermoleophilia Acidobacteriae Bacteroidia Bacilli Vicinamibacteria Gemmatimonadetes Verrucomicrobiae

``` r
top.c.core.bact <- take.top.classes(core.bact.pruned, 10)
```

    ## these classes are selected:  Actinobacteria Alphaproteobacteria Bacilli Gammaproteobacteria Acidimicrobiia Thermoleophilia Gemmatimonadetes Bacteroidia KD4-96 Acidobacteriae

``` r
top.c.all.euks <- take.top.classes(mtcglom.ps.list[["all.euks"]], 10)
```

    ## these classes are selected:  Sordariomycetes Leotiomycetes Mortierellomycetes Dothideomycetes Eurotiomycetes Tremellomycetes Chlorophyceae Agaricomycetes Pezizomycetes Chromadorea

``` r
top.c.core.euks <- take.top.classes(core.euks.pruned, 10)
```

    ## these classes are selected:  Sordariomycetes Leotiomycetes Mortierellomycetes Dothideomycetes Tremellomycetes Eurotiomycetes Agaricomycetes Mucoromycetes Chlorophyceae Pezizomycetes

#### Bacteria at class-level

``` r
# all bacteria
keep.these <- sample_names(top.c.all.bact)[which(!is.na(data.frame(top.c.all.bact@sam_data)$state))] # prune a few state NAs
ca <- prune_samples(keep.these, top.c.all.bact) %>% 
        plot_bar(fill = "Class")+ 
         geom_bar(stat="identity", position="stack")+ 
         scale_y_continuous("")+
         scale_x_discrete("")+
         scale_fill_manual("Class", values = natparks.pals("Torres", 12)[c(2:8,10:12)])+
         facet_grid(~factor(state, levels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2")), scales = "free", space = "free")+
         ggtitle("A")+theme_bw()+
         theme(axis.text.x = element_blank(), legend.position = "none", axis.ticks.x = element_blank())

# core bacteria
keep.these <- sample_names(top.c.core.bact)[which(!is.na(data.frame(top.c.core.bact@sam_data)$state))] # prune a few state NAs
cb <- prune_samples(keep.these, top.c.core.bact) %>% 
        plot_bar(fill = "Class")+ 
         geom_bar(stat="identity", position="stack")+
          scale_y_continuous("", limits = c(0,20))+ 
          scale_x_discrete("")+
          scale_fill_manual("Class", values = natparks.pals("Torres", 12)[c(1:10)])+
          facet_grid(~factor(state, levels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2")), scales = "free", space = "free")+
          ggtitle("B")+theme_bw()+
          theme(axis.text.x = element_blank(), legend.position = "none", axis.ticks.x = element_blank())

# combine bacterial legend
bact.for.legend <- take.top.classes(mtcglom.ps.list[["all.bact"]], 15) # want all these classes, except Chloroflexia, Ktedo, and Polyangia
```

    ## these classes are selected:  Alphaproteobacteria Actinobacteria Gammaproteobacteria Thermoleophilia Acidobacteriae Bacteroidia Bacilli Vicinamibacteria Gemmatimonadetes Verrucomicrobiae Chloroflexia KD4-96 Acidimicrobiia Ktedonobacteria Polyangia

``` r
keep.taxa.names <- taxa_names(bact.for.legend)[!(taxa_names(bact.for.legend) %in% c("ASV105058", "ASV146159", "ASV198"))] # which correspond to these ASV names

bact.for.legend <- prune_taxa(keep.taxa.names, bact.for.legend) # remove the offending taxa
view(bact.for.legend@tax_table)

bact.legend <- plot_bar(bact.for.legend, fill="Class")+ # make a dummy plot
   geom_bar(stat="identity", position="stack")+ 
   scale_fill_manual("Class", values = natparks.pals("Torres", 12))+
   theme(axis.text.x = element_blank())
bact.legend <- get_legend(bact.legend) # extract its legend

# plot bacteria 
bact.class <- ((ca / cb + plot_layout(heights = c(2,1))) | bact.legend) + plot_layout(widths = c(3,1)) # specify proportions
bact.class
grid::grid.draw(grid::textGrob("Percent abundance", x = 0.02, rot = 90)) # add in a common y-axis label
```

![](32_reviewables_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

#### Eukaryotes at class-level

``` r
# all euks
keep.these <- sample_names(top.c.all.euks)[which(!is.na(data.frame(top.c.all.euks@sam_data)$state))] # prune a few state NAs
cc <- prune_samples(keep.these, top.c.all.euks) %>% 
        plot_bar(fill = "Class")+ 
         geom_bar(stat="identity", position="stack")+ 
         scale_y_continuous("")+
         scale_x_discrete("")+
         scale_fill_manual("Class", values = natparks.pals("DeathValley", 11)[c(1:7,9:11)])+
         facet_grid(~factor(state, levels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2")), scales = "free", space = "free")+
         ggtitle("C")+theme_bw()+
         theme(axis.text.x = element_blank(), legend.position = "none", axis.ticks.x = element_blank())

# core euks
keep.these <- sample_names(top.c.core.euks)[which(!is.na(data.frame(top.c.core.euks@sam_data)$state))] # prune a few state NAs
cd <- prune_samples(keep.these, top.c.core.euks) %>% 
        plot_bar(fill = "Class")+ 
         geom_bar(stat="identity", position="stack")+
          scale_y_continuous("", limits = c(0,60))+ 
          scale_x_discrete("")+
          scale_fill_manual("Class", values = natparks.pals("DeathValley", 11)[c(1:2,4:11)])+
          facet_grid(~factor(state, levels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2")), scales = "free", space = "free")+
          ggtitle("D")+theme_bw()+
          theme(axis.text.x = element_blank(), legend.position = "none", axis.ticks.x = element_blank())

# combine euk legend
euks.for.legend <- take.top.classes(mtcglom.ps.list[["all.euks"]], 15) # want all these classes, except Treboux, Spizello, Lobulo, and Spiro.
```

    ## these classes are selected:  Sordariomycetes Leotiomycetes Mortierellomycetes Dothideomycetes Eurotiomycetes Tremellomycetes Chlorophyceae Agaricomycetes Pezizomycetes Chromadorea Trebouxiophyceae Spizellomycetes Lobulomycetes Spirotrichea Mucoromycetes

``` r
keep.taxa.names <- taxa_names(euks.for.legend)[!(taxa_names(euks.for.legend) %in% c("ASV5313", "ASV1449", "ASV3508", "ASV317"))] # which correspond to these ASV names
euks.for.legend <- prune_taxa(keep.taxa.names, euks.for.legend) # remove the offending taxa
euks.legend <- plot_bar(euks.for.legend, fill="Class")+ # make a dummy plot
   geom_bar(stat="identity", position="stack")+ 
   scale_fill_manual("Class", values = natparks.pals("DeathValley", 11))+
   theme(axis.text.x = element_blank())
euks.legend <- get_legend(euks.legend) # extract its legend

# plot euks
euks.class <- ((cc / cd + plot_layout(heights = c(1,0.7))) | euks.legend) + plot_layout(widths = c(3, 1)) # specify proportions
euks.class
grid::grid.draw(grid::textGrob("Percent abundance", x = 0.02, rot = 90)) # add in a common y-axis label
```

![](32_reviewables_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### Occupancy (num states) and abundance

Just to scratch the surface enough to address reviewer concerns about
sequence depth.

``` r
# plot mean occupancy and number of states each ASV was detected in 
bact.state <- ggplot(bact.asvc, aes(num_states_per_asv_test, log10(pctreads), 
                       fill = as.factor(as.numeric(num_states_per_asv_test))))+geom_violin()+
  scale_y_continuous("log10 mean % abundance", limits = c(-5,1))+
  scale_x_continuous("", breaks = c(1:9))+
  scale_fill_manual(values = natparks.pals("Acadia", 9))+
  theme_bw()+theme(legend.position = "none")+ggtitle("A")

euks.state <- ggplot(euks.asvc, aes(num_states_per_asv_test, log10(pctreads), 
                       fill = as.factor(as.numeric(num_states_per_asv_test))))+geom_violin()+
  scale_y_continuous("", limits = c(-5,1))+
  scale_x_continuous("", breaks = c(1:9))+
  scale_fill_manual(values = natparks.pals("Acadia", 9))+
  theme_bw()+theme(legend.position = "none")+ggtitle("B")

# numbers of ASVs across each state since I didn't even show those earlier
bact.count <- ggplot(bact.asvc %>% count(num_states_per_asv_test), aes(num_states_per_asv_test, log10(n), 
                       color = as.factor(as.numeric(num_states_per_asv_test))))+
  geom_point()+
  scale_y_continuous("log10 ASV count", limits = c(2,6))+
  scale_x_continuous("", breaks = c(1:9))+
  scale_color_manual(values = natparks.pals("Acadia", 9))+
  theme_bw()+theme(legend.position = "none")+ggtitle("C")

euks.count <- ggplot(euks.asvc %>% count(num_states_per_asv_test), aes(num_states_per_asv_test, log10(n), 
                       color = as.factor(as.numeric(num_states_per_asv_test))))+
  geom_point()+
  scale_y_continuous("", limits = c(1,5))+
  scale_x_continuous("", breaks = c(1:9))+
  scale_color_manual(values = natparks.pals("Acadia", 9))+
  theme_bw()+theme(legend.position = "none")+ggtitle("D")

(bact.state + euks.state) / (bact.count + euks.count) + plot_layout(heights = c(2, 1))
```

    ## Warning: Removed 214515 rows containing non-finite values (`stat_ydensity()`).

    ## Warning: Removed 1827 rows containing non-finite values (`stat_ydensity()`).

``` r
grid::grid.draw(grid::textGrob("Number of sites an ASV was found in", x = 0.5, y = 0.05, rot = 0))
```

![](32_reviewables_files/figure-gfm/unnamed-chunk-8-1.png)<!-- --> A-
Bacteria, B- Eukaryotes. ASVs with three reads or less in the entire
dataset are omitted from this plot, as they are cut off at y \< -5. 34M
and 31M reads in the bacterial and eukaryotic phyloseqs respectively.

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
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] cowplot_1.1.1          gridExtra_2.3          scales_1.2.1          
    ##  [4] UpSetR_1.4.0           NatParksPalettes_0.2.0 patchwork_1.1.2       
    ##  [7] speedyseq_0.5.3.9018   phyloseq_1.42.0        lubridate_1.9.2       
    ## [10] forcats_1.0.0          stringr_1.5.0          dplyr_1.1.1           
    ## [13] purrr_1.0.1            readr_2.1.4            tidyr_1.3.0           
    ## [16] tibble_3.2.1           ggplot2_3.4.2          tidyverse_2.0.0       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Biobase_2.58.0         bit64_4.0.5            vroom_1.6.1           
    ##  [4] jsonlite_1.8.4         splines_4.2.3          foreach_1.5.2         
    ##  [7] highr_0.10             stats4_4.2.3           GenomeInfoDbData_1.2.9
    ## [10] yaml_2.3.7             pillar_1.9.0           lattice_0.20-45       
    ## [13] glue_1.6.2             digest_0.6.31          XVector_0.38.0        
    ## [16] colorspace_2.1-0       htmltools_0.5.5        Matrix_1.5-3          
    ## [19] plyr_1.8.8             pkgconfig_2.0.3        zlibbioc_1.44.0       
    ## [22] tzdb_0.4.0             timechange_0.2.0       mgcv_1.8-42           
    ## [25] farver_2.1.1           generics_0.1.3         IRanges_2.32.0        
    ## [28] withr_2.5.0            BiocGenerics_0.44.0    cli_3.6.1             
    ## [31] survival_3.5-3         magrittr_2.0.3         crayon_1.5.2          
    ## [34] evaluate_0.20          fansi_1.0.4            nlme_3.1-162          
    ## [37] MASS_7.3-58.2          vegan_2.6-4            tools_4.2.3           
    ## [40] data.table_1.14.8      hms_1.1.3              lifecycle_1.0.3       
    ## [43] Rhdf5lib_1.20.0        S4Vectors_0.36.2       munsell_0.5.0         
    ## [46] cluster_2.1.4          Biostrings_2.66.0      ade4_1.7-22           
    ## [49] compiler_4.2.3         GenomeInfoDb_1.34.9    rlang_1.1.0           
    ## [52] rhdf5_2.42.0           RCurl_1.98-1.12        iterators_1.0.14      
    ## [55] rhdf5filters_1.10.1    biomformat_1.26.0      rstudioapi_0.14       
    ## [58] igraph_1.4.1           labeling_0.4.2         bitops_1.0-7          
    ## [61] rmarkdown_2.21         gtable_0.3.3           codetools_0.2-19      
    ## [64] multtest_2.54.0        reshape2_1.4.4         R6_2.5.1              
    ## [67] knitr_1.42             bit_4.0.5              fastmap_1.1.1         
    ## [70] utf8_1.2.3             permute_0.9-7          ape_5.7-1             
    ## [73] stringi_1.7.12         parallel_4.2.3         Rcpp_1.0.10           
    ## [76] vctrs_0.6.1            tidyselect_1.2.0       xfun_0.38
