24_taxonomy_II
================
2023-02-24

#### Purpose

Make better plots of taxa and core taxa for 16S and ITS at the Phylum
level, and plot variation of top core taxa across states.

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

#### load previously made phyloseq objects, number the ASVs, and dump them into a list

``` r
# communities-all
bact.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/first_two_years_all_sites/bact1920.cleaned.ps")
euks.ps <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/first_two_years_all_sites/fung1920.cleaned.ps")

# number_asvs assigns ASV numbers to replace the sequences as the rownames. takes only a ps object
number_asvs <- function(ps){
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
  return(ps)
}
# number bact and euk ASVs IF NECESSARY
# bact.ps <- number_asvs(bact.ps) 
# euks.ps <- number_asvs(euks.ps)

# core communities - last col of otu table is "summing", all its tax categories are also "summing"
bact.core.ps <- readRDS(file = "/Users/klas0061/Desktop/UMN/phyloseqs/core/core.bact.asvs.ps")
euks.core.ps <- readRDS(file = "/Users/klas0061/Desktop/UMN/phyloseqs/core/core.its.asvs.ps")

# make a list
ps.list <- list(bact.ps, euks.ps, bact.core.ps, euks.core.ps)
names(ps.list) <- c("all.bact", "all.euks", "core.bact", "core.euks")
```

#### tax-glom to Phylum, transform to percent abundance, omit low-abundance phyla

``` r
# glom at phylum level
phyglom.ps.list <- lapply(ps.list, tax_glom, taxrank = "Phylum")

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

# run this to merge samples by replicates, calculate percent abundances, and fix state labels
mtphyglom.ps.list <- lapply(phyglom.ps.list, merge.and.tidy.samples) 
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

    ## Warning in asMethod(object): NAs introduced by coercion

    ## Warning in asMethod(object): NAs introduced by coercion

``` r
# remove summing columns from core ps objects
prune_sums <- function(ps){
  keep.taxa.names <- taxa_names(ps)[which(taxa_names(ps) != "summing")]
  prune.ps <- prune_taxa(keep.taxa.names, ps)
  return(prune.ps)
}
core.bact.pruned <- prune_sums(mtphyglom.ps.list[["core.bact"]])
core.euks.pruned <- prune_sums(mtphyglom.ps.list[["core.euks"]])

# prune all but most abundant n phyla
take.top.phyla <- function(ps, nphyla){
  asvs.of.phyla.descending <- names(sort(colSums(ps@otu_table) / nrow(ps@otu_table), decreasing = TRUE)) # the representative ASVs of phyla
  phy <- data.frame(ps@tax_table)[asvs.of.phyla.descending,]$Phylum # their corresponding phyla, in order of descending abundance
  cat("these phyla are selected: ", phy[1:nphyla], "\n") # lists the phyla 
  ps.pruned <- prune_taxa(asvs.of.phyla.descending[1:nphyla], ps) # prunes all but the top n phyla from the ps object
  return(ps.pruned)
}
top.phy.all.bact <- take.top.phyla(mtphyglom.ps.list[["all.bact"]], 10)
```

    ## these phyla are selected:  Proteobacteria Actinobacteriota Chloroflexi Acidobacteriota Firmicutes Bacteroidota Gemmatimonadota Verrucomicrobiota Myxococcota Planctomycetota

``` r
top.phy.core.bact <- take.top.phyla(core.bact.pruned, 10)
```

    ## these phyla are selected:  Actinobacteriota Proteobacteria Firmicutes Acidobacteriota Chloroflexi Gemmatimonadota Bacteroidota Myxococcota Nitrospirota Methylomirabilota

``` r
top.phy.all.euks <- take.top.phyla(mtphyglom.ps.list[["all.euks"]], 10)
```

    ## these phyla are selected:  Ascomycota Basidiomycota Mortierellomycota Chlorophyta Chytridiomycota Nematoda Ciliophora Mucoromycota Bryophyta Arthropoda

``` r
top.phy.core.euks <- take.top.phyla(core.euks.pruned, 9) # there were only 9 core eukaryotic phyla
```

    ## these phyla are selected:  Ascomycota Mortierellomycota Basidiomycota Mucoromycota Chlorophyta Monoblepharomycota Nematoda Chytridiomycota Ciliophora

``` r
length(data.frame(core.euks.pruned@tax_table)$Phylum)
```

    ## [1] 9

#### barplots of total communities, by phylum

``` r
# all bacteria
keep.these <- sample_names(top.phy.all.bact)[which(!is.na(data.frame(top.phy.all.bact@sam_data)$state))] # prune a few state NAs
top.phy.all.bact <- prune_samples(keep.these, top.phy.all.bact)
all.bact.bar <- plot_bar(top.phy.all.bact, fill="Phylum")+
   geom_bar(stat="identity", position="stack")+ 
   scale_y_continuous("")+
   scale_x_discrete("")+
   scale_fill_manual("Phylum", values = natparks.pals("Torres", 12)[c(1:6,8,10:12)])+
   facet_grid(~factor(state, levels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2")), scales = "free", space = "free")+
   ggtitle("A")+theme_bw()+
   theme(axis.text.x = element_blank(), legend.position = "none", axis.ticks.x = element_blank())

# all euks
keep.these <- sample_names(top.phy.all.euks)[which(!is.na(data.frame(top.phy.all.euks@sam_data)$state))] # prune a few state NAs
top.phy.all.euks <- prune_samples(keep.these, top.phy.all.euks)
all.euks.bar <- plot_bar(top.phy.all.euks, fill="Phylum")+
   geom_bar(stat="identity", position="stack")+ 
   scale_y_continuous("")+
   scale_x_discrete("")+
   scale_fill_manual("Phylum", values = natparks.pals("DeathValley", 11)[c(1:7, 9:11)])+
   facet_grid(~factor(state, levels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2")), scales = "free", space = "free")+
   ggtitle("C")+theme_bw()+
   theme(axis.text.x = element_blank(), legend.position = "none", axis.ticks.x = element_blank())
```

#### barplots of core communities, by phylum

``` r
# core bacteria
keep.these <- sample_names(top.phy.core.bact)[which(!is.na(data.frame(top.phy.core.bact@sam_data)$state))] # prune a few state NAs
top.phy.core.bact <- prune_samples(keep.these, top.phy.core.bact)
core.bact.bar <- plot_bar(top.phy.core.bact, fill="Phylum")+
   geom_bar(stat="identity", position="stack")+ 
   scale_y_continuous("", limits = c(0,20))+ 
   scale_x_discrete("")+
   scale_fill_manual("Phylum", values = natparks.pals("Torres", 12)[c(1:9, 11)])+
   facet_grid(~factor(state, levels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2")), scales = "free", space = "free")+
   ggtitle("B")+theme_bw()+
   theme(axis.text.x = element_blank(), legend.position = "none", axis.ticks.x = element_blank())

# core euks
keep.these <- sample_names(top.phy.core.euks)[which(!is.na(data.frame(top.phy.core.euks@sam_data)$state))] # prune a few state NAs
top.phy.core.euks <- prune_samples(keep.these, top.phy.core.euks)
core.euks.bar <- plot_bar(top.phy.core.euks, fill="Phylum")+
   geom_bar(stat="identity", position="stack")+ 
   scale_y_continuous("")+
   scale_x_discrete("")+
   scale_fill_manual("Phylum", values = natparks.pals("DeathValley", 11)[c(2,3,5:11)])+
   facet_grid(~factor(state, levels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2")), scales = "free", space = "free")+
   ggtitle("D")+theme_bw()+
   theme(axis.text.x = element_blank(), legend.position = "none", axis.ticks.x = element_blank())
```

#### Attempt to smush all plots together

combining the color legend in patchwork as shown below gives two
legends, when I want one combined one. Solution is probably to make a
dummy legend and add it as a separate plot element…

``` r
# get bacterial combined legend
take.top.phyla(mtphyglom.ps.list[["all.bact"]], 15)@tax_table
```

    ## these phyla are selected:  Proteobacteria Actinobacteriota Chloroflexi Acidobacteriota Firmicutes Bacteroidota Gemmatimonadota Verrucomicrobiota Myxococcota Planctomycetota Patescibacteria Nitrospirota Bdellovibrionota Methylomirabilota Cyanobacteria

    ## Taxonomy Table:     [ 15 taxa by 7 taxonomic ranks ]:
    ##              Kingdom  Phylum            Class Order Family Genus Species
    ##              <chr>    <chr>             <chr> <chr> <chr>  <chr> <chr>  
    ##  1 ASV2      Bacteria Actinobacteriota  <NA>  <NA>  <NA>   <NA>  <NA>   
    ##  2 ASV109340 Bacteria Chloroflexi       <NA>  <NA>  <NA>   <NA>  <NA>   
    ##  3 ASV4      Bacteria Firmicutes        <NA>  <NA>  <NA>   <NA>  <NA>   
    ##  4 ASV7      Bacteria Proteobacteria    <NA>  <NA>  <NA>   <NA>  <NA>   
    ##  5 ASV98903  Bacteria Acidobacteriota   <NA>  <NA>  <NA>   <NA>  <NA>   
    ##  6 ASV329    Bacteria Gemmatimonadota   <NA>  <NA>  <NA>   <NA>  <NA>   
    ##  7 ASV34     Bacteria Nitrospirota      <NA>  <NA>  <NA>   <NA>  <NA>   
    ##  8 ASV198    Bacteria Myxococcota       <NA>  <NA>  <NA>   <NA>  <NA>   
    ##  9 ASV907    Bacteria Bacteroidota      <NA>  <NA>  <NA>   <NA>  <NA>   
    ## 10 ASV2178   Bacteria Cyanobacteria     <NA>  <NA>  <NA>   <NA>  <NA>   
    ## 11 ASV8411   Bacteria Verrucomicrobiota <NA>  <NA>  <NA>   <NA>  <NA>   
    ## 12 ASV6228   Bacteria Methylomirabilota <NA>  <NA>  <NA>   <NA>  <NA>   
    ## 13 ASV101177 Bacteria Planctomycetota   <NA>  <NA>  <NA>   <NA>  <NA>   
    ## 14 ASV1791   Bacteria Bdellovibrionota  <NA>  <NA>  <NA>   <NA>  <NA>   
    ## 15 ASV40003  Bacteria Patescibacteria   <NA>  <NA>  <NA>   <NA>  <NA>

``` r
bact.for.legend <- take.top.phyla(mtphyglom.ps.list[["all.bact"]], 15) # want all these phyla, except Patescibacteria, Bdellovibrionota, and Latescibacterota
```

    ## these phyla are selected:  Proteobacteria Actinobacteriota Chloroflexi Acidobacteriota Firmicutes Bacteroidota Gemmatimonadota Verrucomicrobiota Myxococcota Planctomycetota Patescibacteria Nitrospirota Bdellovibrionota Methylomirabilota Cyanobacteria

``` r
keep.taxa.names <- taxa_names(bact.for.legend)[!(taxa_names(bact.for.legend) %in% c("ASV40003", "ASV1791", "ASV2178"))] # which correspond to these ASV names

bact.for.legend <- prune_taxa(keep.taxa.names, bact.for.legend) # remove the offending taxa
bact.legend <- plot_bar(bact.for.legend, fill="Phylum")+ # make a dummy plot
   geom_bar(stat="identity", position="stack")+ 
   scale_fill_manual("Phylum", values = natparks.pals("Torres", 12))+
   theme(axis.text.x = element_blank())
bact.legend <- get_legend(bact.legend) # extract its legend

# plot bacteria 
bact.phyla <- ((all.bact.bar / core.bact.bar + plot_layout(heights = c(2,1))) | bact.legend) + plot_layout(widths = c(3,1)) # specify proportions
bact.phyla
grid::grid.draw(grid::textGrob("Percent abundance", x = 0.02, rot = 90)) # add in a common y-axis label
```

![](24_taxonomy_II_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# get euk combined legend
take.top.phyla(mtphyglom.ps.list[["all.euks"]], 14)@tax_table
```

    ## these phyla are selected:  Ascomycota Basidiomycota Mortierellomycota Chlorophyta Chytridiomycota Nematoda Ciliophora Mucoromycota Bryophyta Arthropoda Blastocladiomycota Heterolobosa_phy_Incertae_sedis Monoblepharomycota Cercozoa

    ## Taxonomy Table:     [ 14 taxa by 7 taxonomic ranks ]:
    ##            Kingdom       Phylum                 Class Order Family Genus Species
    ##            <chr>         <chr>                  <chr> <chr> <chr>  <chr> <chr>  
    ##  1 ASV578  Fungi         Ascomycota             <NA>  <NA>  <NA>   <NA>  <NA>   
    ##  2 ASV4    Fungi         Basidiomycota          <NA>  <NA>  <NA>   <NA>  <NA>   
    ##  3 ASV17   Fungi         Mortierellomycota      <NA>  <NA>  <NA>   <NA>  <NA>   
    ##  4 ASV76   Fungi         Blastocladiomycota     <NA>  <NA>  <NA>   <NA>  <NA>   
    ##  5 ASV3346 Heterolobosa  Heterolobosa_phy_Ince… <NA>  <NA>  <NA>   <NA>  <NA>   
    ##  6 ASV3508 Fungi         Chytridiomycota        <NA>  <NA>  <NA>   <NA>  <NA>   
    ##  7 ASV668  Metazoa       Nematoda               <NA>  <NA>  <NA>   <NA>  <NA>   
    ##  8 ASV317  Alveolata     Ciliophora             <NA>  <NA>  <NA>   <NA>  <NA>   
    ##  9 ASV136  Fungi         Mucoromycota           <NA>  <NA>  <NA>   <NA>  <NA>   
    ## 10 ASV139  Rhizaria      Cercozoa               <NA>  <NA>  <NA>   <NA>  <NA>   
    ## 11 ASV9829 Metazoa       Arthropoda             <NA>  <NA>  <NA>   <NA>  <NA>   
    ## 12 ASV5313 Viridiplantae Chlorophyta            <NA>  <NA>  <NA>   <NA>  <NA>   
    ## 13 ASV936  Fungi         Monoblepharomycota     <NA>  <NA>  <NA>   <NA>  <NA>   
    ## 14 ASV5320 Viridiplantae Bryophyta              <NA>  <NA>  <NA>   <NA>  <NA>

``` r
euks.for.legend <- take.top.phyla(mtphyglom.ps.list[["all.euks"]], 14) # want all these phyla, except three: Heterolobosa_phy_Incertae_sedis, Blastocladiomycota, Cercozoa
```

    ## these phyla are selected:  Ascomycota Basidiomycota Mortierellomycota Chlorophyta Chytridiomycota Nematoda Ciliophora Mucoromycota Bryophyta Arthropoda Blastocladiomycota Heterolobosa_phy_Incertae_sedis Monoblepharomycota Cercozoa

``` r
keep.taxa.names <- taxa_names(euks.for.legend)[!(taxa_names(euks.for.legend) %in% c("ASV139", "ASV3346", "ASV76"))]
euks.for.legend <- prune_taxa(keep.taxa.names, euks.for.legend) # remove the offending taxa
euks.legend <- plot_bar(euks.for.legend, fill="Phylum")+ # make a dummy plot
   geom_bar(stat="identity", position="stack")+ 
   scale_fill_manual("Phylum", values = natparks.pals("DeathValley", 11))+
   theme(axis.text.x = element_blank())
euks.legend <- get_legend(euks.legend) # extract its legend

# plot euks
euks.phyla <- ((all.euks.bar / core.euks.bar + plot_layout(heights = c(1,0.7))) | euks.legend) + plot_layout(widths = c(3,1)) # specify proportions
euks.phyla 
grid::grid.draw(grid::textGrob("Percent abundance", x = 0.02, rot = 90)) # add in a common y-axis label
```

![](24_taxonomy_II_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
# ALL together: too unstable, can't adjust proportions with plot_layout for some reason
# (all.bact.bar / core.bact.bar | bact.legend) / (all.euks.bar / core.euks.bar | euks.legend)
# grid::grid.draw(grid::textGrob("Percent abundance", x = 0.015, rot = 90)) # add in a common y-axis label
```

I also had issues with the pdf barplots made as outputs: the colors of
the final figures were way darker than they looked within RStudio.

#### Plotting a few examples of differentially abundant ASVs from differentially abundant genera across states

I ran aldex and the results basically said: everything is differentially
abundant. Well ok, all core ITS ASVs and genera they belong to are, but
only 93% of core bacterial ASVs and 96% of genera they belong to. The
most illustrative examples of this are Mortierella in the fungi, and
Bacillus in the bacteria, as they are both highly abundant genera with
many core ASVs that all happen to be differentially abundant across
sites. They’re also implicated in plant health and growth promotion.

I also realized what was wrong with my barplots. Since you can’t reduce
the thin black outlines of the bars any more than I already have, better
to merge_samples that are considered “replicates”. This is like 24-30
communities per barplot in all sites except for ME2, where there are 5
communities per barplot.

``` r
ps.euks <- merge.and.tidy.samples(euks.core.ps) # merge and tidy the samples by site/year/season & transform into percent abundances
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

    ## Warning in asMethod(object): NAs introduced by coercion

``` r
biomarkerlist <- rownames(data.frame(ps.euks@tax_table) %>% filter(Genus == "Mortierella")) # get the core Mortierella ASVs 
colSums(ps.euks@otu_table[,biomarkerlist])/nsamples(ps.euks) # look at their mean percent abundances
```

    ##    ASV1350      ASV17     ASV516     ASV878 
    ## 0.01781564 2.19391533 1.02369428 1.27923308

``` r
biomarkerlist <- c("ASV17", "ASV516", "ASV878") # that M. polygonia is too low percent abundance to be visible

# stuff from the function
num <- vector("numeric") # define a numeric vector
  for (i in biomarkerlist) {num[i] <- mean(otu_table(ps.euks)[,i])} # calculate mean percent abundances for all biomarker ASVs
  asvs.to.subset <- vector("integer") # define an integer vector
  for (i in names(num)) {asvs.to.subset[i] <- which(rownames(tax_table(ps.euks))==i)} # obtains the row numbers corresponding to biomarker ASVs
  bmtt <- tax_table(ps.euks)[asvs.to.subset,] # subsets the tax table
  bmtt[,"Species"] <- paste("M.", bmtt[,"Species"])
  bmasvt <- otu_table(ps.euks)[,names(num)] # subsets the ASV table
  ps.bm <- phyloseq(tax_table(bmtt), 
                 sample_data(ps.euks),
                 otu_table(bmasvt, taxa_are_rows = FALSE)) # writes subsetted files into a phyloseq object
mort.barplot <- plot_bar(ps.bm, fill="Species")+
  geom_bar(stat="identity", position="stack", color = NA)+
  scale_y_continuous("% Abundance")+
  facet_grid(~state, scales = "free", space = "free")+
  scale_fill_manual("", values = c('#fc8d62','#8da0cb','#66c2a5'))+ # M. polygonia was '#e78ac3'  
  scale_x_discrete("")+ggtitle("B")+theme_bw()+
  theme(axis.text.x = element_blank(), legend.position = "bottom", axis.ticks.x = element_blank()) # bar graphs % abundances across all samples

### next Bacillus
ps.bact <- merge.and.tidy.samples(bact.core.ps)
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
biomarkerlist <- rownames(data.frame(ps.bact@tax_table) %>% filter(Genus == "Bacillus"))
omit.samples <- c("tech_rep_1_16S_S1460", "tech_rep_2_16S_S1461", "tech_rep_3_16S_S1462", "OR_1_Compost_10_19_16S_S1299", "OR_1_2_blank_20_4_16S_S1039", "OR_1_2_blank_20_6_16S_S1040", "OR_1_3_blank_20_8_16S_S1065")
ps <- prune_samples(setdiff(sample_names(ps.bact), omit.samples), ps.bact)

bacillus.mean.pa <- colSums(ps@otu_table[,biomarkerlist])/nsamples(ps) # look at their mean percent abundances
sort(bacillus.mean.pa, decreasing = TRUE)
```

    ##         ASV4         ASV5       ASV249        ASV36       ASV579       ASV255 
    ## 0.3274190847 0.3007835711 0.1984539712 0.0817172036 0.0286459274 0.0237375330 
    ##      ASV4472       ASV557     ASV14131       ASV971     ASV37738 
    ## 0.0152871168 0.0143340453 0.0138929383 0.0129416299 0.0004817457

``` r
biomarkerlist <- names(sort(bacillus.mean.pa, decreasing = TRUE)[1:4]) # just take the top n

# stuff from the function
num <- vector("numeric") # define a numeric vector
  for (i in biomarkerlist) {num[i] <- mean(otu_table(ps)[,i])} # calculate mean percent abundances for all biomarker ASVs
  asvs.to.subset <- vector("integer") # define an integer vector
  for (i in names(num)) {asvs.to.subset[i] <- which(rownames(tax_table(ps))==i)} # obtains the row numbers corresponding to biomarker ASVs
  bmtt <- tax_table(ps)[asvs.to.subset,] # subsets the tax table
  bmtt.df <- data.frame(bmtt)
  bmtt.df$ASV <- rownames(bmtt.df)
  bmtt.df$Label <- paste(bmtt.df$Genus, bmtt.df$ASV)
  # bmtt.df$Label[1] <- "Bacillus niacini ASV14131" # now I don't have B. niacini anymore
  bmtt <- as.matrix(bmtt.df)
  bmasvt <- otu_table(ps)[,names(num)] # subsets the ASV table
  ps.bm <- phyloseq(tax_table(bmtt), 
                 sample_data(ps),
                 otu_table(bmasvt, taxa_are_rows = FALSE)) # writes subsetted files into a phyloseq object
bac.barplot <- plot_bar(ps.bm, fill="Label")+
  geom_bar(stat="identity", position="stack", color = NA)+
  scale_y_continuous("% Abundance", breaks = c(0, 1, 2))+
  facet_grid(~state, scales = "free", space = "free")+
  scale_fill_manual("", values = natparks.pals("Volcanoes", 6)[c(2,5,3,4)])+  
  scale_x_discrete("")+ggtitle("A")+theme_bw()+
  theme(axis.text.x = element_blank(), legend.position = "none", axis.ticks.x = element_blank())

bac.barplot / mort.barplot
```

![](24_taxonomy_II_files/figure-gfm/unnamed-chunk-7-1.png)<!-- --> Notes
on color legend for A, Bacillus: I used only the top four ASVs, none of
which had species names. Lowest percent abundance of this one was 0.08%,
and the top one omitted was 0.03%. You can now actually see the colors
because there aren’t over a thousand barplots, and the
replicate-averaged bars remove any outlier communities. For the most
part these spatially-variable populations are pretty stable over
year/season with maybe a couple exceptions (like OR, WI M. elongata vs
hyalina).

#### notes about aldex diff abundance

I ran the tax-glommed differential abundance calculations by state
interactively, using the code below. Not running here cuz it took a long
time and once was enough.

``` r
# core its asvs genus-glommed- started 5:18, went on run, finished within an hour
core.its.g <- tax_glom(core.its.asvs.ps, taxrank = "Genus")
table(core.its.g@sam_data$state) # there should just be the nine sites
ald.its <- aldex.clr(data.frame(t(core.its.g@otu_table)), core.its.g@sam_data$state, mc.samples = 128, verbose = FALSE, useMC=TRUE) # do NOT print output
ald.its.kw <- aldex.kw(ald.its, useMC = TRUE, verbose = FALSE) # this takes a while 
write.csv(ald.its.kw, file = "/Users/scottklasek/Desktop/aldex.kw.core.its.genus.glom.csv")
# core.its.g@tax_table

# core bact asvs genus-glommed- started 6:19 
core.bact.g <- tax_glom(core.bact.asvs.ps, taxrank = "Genus")
table(core.bact.g@sam_data$state) # there should just be the nine sites
ald.bact <- aldex.clr(data.frame(t(core.bact.g@otu_table)), core.bact.g@sam_data$state, mc.samples = 128, verbose = FALSE, useMC=TRUE) # do NOT print output
ald.bact.kw <- aldex.kw(ald.bact, useMC = TRUE, verbose = FALSE) # this takes a while 
write.c
```

#### bonus

How much of the reads are from core members?

``` r
sum(colSums(bact.core.ps@otu_table[,1:606])) / sum(colSums(bact.core.ps@otu_table)) * 100
```

    ## [1] 13.70002

``` r
sum(colSums(euks.core.ps@otu_table[,1:74])) / sum(colSums(euks.core.ps@otu_table)) * 100
```

    ## [1] 27.40905

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
    ##  [1] Biobase_2.58.0         jsonlite_1.8.4         splines_4.2.3         
    ##  [4] foreach_1.5.2          highr_0.10             stats4_4.2.3          
    ##  [7] GenomeInfoDbData_1.2.9 yaml_2.3.7             pillar_1.9.0          
    ## [10] lattice_0.20-45        glue_1.6.2             digest_0.6.31         
    ## [13] XVector_0.38.0         colorspace_2.1-0       htmltools_0.5.5       
    ## [16] Matrix_1.5-3           plyr_1.8.8             pkgconfig_2.0.3       
    ## [19] zlibbioc_1.44.0        tzdb_0.3.0             timechange_0.2.0      
    ## [22] mgcv_1.8-42            farver_2.1.1           generics_0.1.3        
    ## [25] IRanges_2.32.0         withr_2.5.0            BiocGenerics_0.44.0   
    ## [28] cli_3.6.1              survival_3.5-3         magrittr_2.0.3        
    ## [31] crayon_1.5.2           evaluate_0.20          fansi_1.0.4           
    ## [34] nlme_3.1-162           MASS_7.3-58.2          vegan_2.6-4           
    ## [37] tools_4.2.3            data.table_1.14.8      hms_1.1.3             
    ## [40] lifecycle_1.0.3        Rhdf5lib_1.20.0        S4Vectors_0.36.2      
    ## [43] munsell_0.5.0          cluster_2.1.4          Biostrings_2.66.0     
    ## [46] ade4_1.7-22            compiler_4.2.3         GenomeInfoDb_1.34.9   
    ## [49] rlang_1.1.0            rhdf5_2.42.0           RCurl_1.98-1.12       
    ## [52] iterators_1.0.14       rhdf5filters_1.10.1    biomformat_1.26.0     
    ## [55] rstudioapi_0.14        igraph_1.4.1           labeling_0.4.2        
    ## [58] bitops_1.0-7           rmarkdown_2.21         gtable_0.3.3          
    ## [61] codetools_0.2-19       multtest_2.54.0        reshape2_1.4.4        
    ## [64] R6_2.5.1               knitr_1.42             fastmap_1.1.1         
    ## [67] utf8_1.2.3             permute_0.9-7          ape_5.7-1             
    ## [70] stringi_1.7.12         parallel_4.2.3         Rcpp_1.0.10           
    ## [73] vctrs_0.6.1            tidyselect_1.2.0       xfun_0.38
