Fungal Beginnings
================
Scott Klasek
2022-08-12

Purpose: fungal analysis for all 2018-2020 samples across the entire
platform. Include Larkin’s samples here, because they’re all in potato.
Lots of this is built off the previous document, 04_add_potato_data, but
this one aims to be more a final, presentable product.

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
library(phyloseq)
library(patchwork)
library(vegan)
```

    ## Loading required package: permute
    ## Loading required package: lattice
    ## This is vegan 2.6-4

``` r
library(DESeq2)
```

    ## Loading required package: S4Vectors
    ## Loading required package: stats4
    ## Loading required package: BiocGenerics
    ## 
    ## Attaching package: 'BiocGenerics'
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     intersect, setdiff, union
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min
    ## 
    ## 
    ## Attaching package: 'S4Vectors'
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     second, second<-
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename
    ## 
    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname
    ## 
    ## Loading required package: IRanges
    ## 
    ## Attaching package: 'IRanges'
    ## 
    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     %within%
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce
    ## 
    ## Loading required package: GenomicRanges
    ## Loading required package: GenomeInfoDb
    ## Loading required package: SummarizedExperiment
    ## Loading required package: MatrixGenerics
    ## Loading required package: matrixStats
    ## 
    ## Attaching package: 'matrixStats'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count
    ## 
    ## 
    ## Attaching package: 'MatrixGenerics'
    ## 
    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars
    ## 
    ## Loading required package: Biobase
    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.
    ## 
    ## 
    ## Attaching package: 'Biobase'
    ## 
    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians
    ## 
    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians
    ## 
    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     sampleNames

``` r
library(metagenomeSeq)
```

    ## Loading required package: limma
    ## 
    ## Attaching package: 'limma'
    ## 
    ## The following object is masked from 'package:DESeq2':
    ## 
    ##     plotMA
    ## 
    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA
    ## 
    ## Loading required package: glmnet
    ## Loading required package: Matrix
    ## 
    ## Attaching package: 'Matrix'
    ## 
    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     expand
    ## 
    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack
    ## 
    ## Loaded glmnet 4.1-7
    ## Loading required package: RColorBrewer

``` r
library(NatParksPalettes)
library(UpSetR)
```

    ## 
    ## Attaching package: 'UpSetR'
    ## 
    ## The following object is masked from 'package:lattice':
    ## 
    ##     histogram

#### define functions

put functions before data import because some functions are used *for*
data import

``` r
# merge_with_jim_data takes as input a phyloseq object and a dataframe of objective 1 data that Jim has curated, which can be subsetted
# it writes data from Jim's spreadsheet into the corresponding phyloseq object, based on combinations
# make sure jim's data is cleaned up according to the steps above!
merge_with_jim_data <- function(ps){
  df <- left_join(data.frame(ps@sam_data), jim.info.s, 
                by = c("state" = "State", "objective" = "Objective", "rotation" = "Rotation",
                       "plot" = "Plot", "year" = "Year", "month" = "Month")) # merges 
  #the unique combo of columns from phyloseq sample data and jim's data
  rownames(df) <- rownames(ps@sam_data) # ok as long as the rows are in the same order, which they are
  sample_data(ps) <- df # put df into the sample_data
  return(ps)
}

# assign ASV numbers cuz i'm sick of the sequences as the rownames
number_asvs <- function(ps){
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
  return(ps)
}

# investigate library sizes
plot_libsize <- function(ps){
  libsizes <- as.data.frame(sample_data(ps)) # make data frame
  libsizes$LibrarySize <- sample_sums(ps) # obtain library size
  libsizes <- libsizes[order(libsizes$LibrarySize),] # order data frame
  libsizes$Index <- seq(nrow(libsizes)) # index samples for plotting
  gg <- ggplot(libsizes, aes(x=Index, y=LibrarySize)) + # make plot
    geom_point() + 
    geom_hline(yintercept=10000,linetype=2)
  return(gg)
}

# compare_taxa_abunds makes a dataframe of the percent abundances of a certain taxon at whatever taxonomic level
# and uses sample info to return data to test differences between groups of samples
compare_taxa_abunds <- function(ps, taxlevel, taxon){
  tax.df <- data.frame(ps@tax_table) # make taxonomy table a dataframe
  asvs <- rownames(tax.df[which(tax.df[taxlevel]==taxon),]) # get the ASVs that correspond to the taxon of interest
  percent.per.sample <- rowSums(ps@otu_table[,asvs]) / sample_sums(ps) * 100 # percentages per sample
  df <- cbind(data.frame(ps@sam_data) %>% select(state, year, season, block, in.potato, post.fumigation), percent.per.sample)
  return(df)
}

# hellinger transformation
helltransform <- function(ps){
  otu.ps3hell <- otu_table(decostand(otu_table(ps), method = "hellinger"), taxa_are_rows=FALSE)
  ps.hel <- phyloseq(tax_table(ps),
                    sample_data(ps),
                    otu_table(otu.ps3hell),
                    refseq(ps))
  return(ps.hel)
}

# asvcumulate calculates prevalence and cumulative abundance of ASVs, ranking by most common to least common, returning a dataframe
# input: a phyloseq object
asvcumulate <- function(ps){
  ps.cumul.abund <- data.frame(taxa_sums(ps)) # make a dataframe of the number of reads for each ASV across the dataset
  colnames(ps.cumul.abund)[1] <- "reads" # relabel the column
  ps.cumul.abund$pctreads <- 100*ps.cumul.abund$reads/sum(taxa_sums(ps)) # this column is the percent of reads in the entire dataset from each ASV
  asv.binary.ps <- data.frame(otu_table(ps)) # write ASV table to data frame
  asv.binary.ps[asv.binary.ps > 0] <- 1 # convert all reads to presence/absence
  asv.prevalence <- colSums(asv.binary.ps) # count the number of communities each ASV appears in 
  numsamples <- nrow(sample_data(ps)) # nsamples decided not to work so i had to add this in here (phyloseq v 1.34)
  asv.pct.prevalence <- (asv.prevalence/numsamples)*100 # convert this into percent
  ps.cumul.abund$pct_prevalence <- asv.pct.prevalence # add in prevalence information (what % of communities have this ASV in them)
  ps.cumul.abund <- ps.cumul.abund[order(-ps.cumul.abund$reads),] # this orders ASVs by most abundant to least (make sure you have more than one row in your dataframe when you do this)
  ps.cumul.abund$csum <- cumsum(ps.cumul.abund$pctreads) # calculate cumulative sum
  ps.cumul.abund$rank <- c(1:nrow(ps.cumul.abund)) # define a rank to plot with
  return(ps.cumul.abund)
}

# count the number of states in which an ASV was detected
num_states_per_asv <- function(ps){
  filter <- phyloseq::genefilter_sample(ps, filterfun_sample(function(x) x > 0)) # for each ASV get a true/false whether it has over 0 reads
  ps <- prune_taxa(filter, ps) # remove the zero-count taxa from the ps object 
  num_states <- vector("numeric") # define our numeric vector, and below: count the number of states in which an ASV is found in 
  for (i in rownames(ps@tax_table)){num_states[i] <- length(unique(ps@sam_data[rownames(ps@otu_table[which(ps@otu_table[,i] > 0),]),]$state))} 
  return(num_states)
}
```

#### import data

1)  Import Jim’s yield and soil info data, do some clean up  
2)  The fungal phyloseq object
3)  Add Jim’s data to the phyloseq object using the custom function,
    merge_with_jim_data At the end of this, we have only one input file,
    the phyloseq object, to work from.

``` r
### 
# Jim's yield and soil data to merge with phyloseq object metadata 
jim.info <- read.csv(file="/Users/klas0061/Desktop/UMN/jim_info/PSHP_ALL_obj_1_all_data_2022_07_08.csv") # import

# basic clean-up stuff
colnames(jim.info) <- jim.info[1,] # fix column names (there are two header cols)
jim.info <- jim.info[2:nrow(jim.info),] # omit redundant column name 
jim.info <- jim.info[,c(1:101,103:ncol(jim.info))] # omit the redundant year in column 102
for (i in 2:7){jim.info[,i] <- as.numeric(jim.info[,i])} # convert columns that should be numeric into numeric

# change OR Spring 2020 month to 4 so Cultivar (whether field is in potato or not) is not NA
jim.info[which(jim.info$State=="OR" & jim.info$Year==20 & jim.info$Month=="3"),"Month"] <- "4"

# select data columns to merge by, and to add (an example, can always select more later)
jim.info.s <- jim.info %>% select(State, Objective, Rotation, Plot, Year, Month, Cultivar, VPPG,
                                  `OM (%)`, `Total yield`, `Marketable yield`, `Greater than 6 oz.`) 

for (i in 8:12){jim.info.s[,i] <- as.numeric(jim.info.s[,i])} # again, convert columns that should be numeric into numeric
```

    ## Warning: NAs introduced by coercion

    ## Warning: NAs introduced by coercion

    ## Warning: NAs introduced by coercion

    ## Warning: NAs introduced by coercion

    ## Warning: NAs introduced by coercion

``` r
### import our phyloseq object of all fungi
all.fung.1920.ps  <- readRDS(file="/Users/klas0061/Desktop/UMN/phyloseqs/first_two_years_all_sites/all.fung1920.ps")

### add jim's data to the sample_data of phyloseq object
all.fung.1920.ps <- merge_with_jim_data(all.fung.1920.ps)
all.fung.1920.ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 23377 taxa and 1476 samples ]
    ## sample_data() Sample Data:       [ 1476 samples by 27 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 23377 taxa by 7 taxonomic ranks ]

#### more data cleaning

define column “in.potato”, whether samples come from potato or rotation
years

``` r
# rewrite Cultivar column of sample data as in.potato
colnames(all.fung.1920.ps@sam_data)[which(colnames(all.fung.1920.ps@sam_data) == "Cultivar")] <- "in.potato"
# change the cultivar to "potato" because we have it in another column
all.fung.1920.ps@sam_data[which(all.fung.1920.ps@sam_data$in.potato != "Rotation"), "in.potato"] <- "Potato"
# Minnesota 2018 samples were NOT in potato
all.fung.1920.ps@sam_data[which(all.fung.1920.ps@sam_data$year==18), "in.potato"] <- "Rotation"
# larkin samples are ALL potato, except after harvest
all.fung.1920.ps@sam_data[which(all.fung.1920.ps@sam_data$state=="US"), "in.potato"] <- "Potato"
all.fung.1920.ps@sam_data[which(all.fung.1920.ps@sam_data$state=="US" & all.fung.1920.ps@sam_data$season=="Fall"), "in.potato"] <- "Rotation"
# make sure larkin treatment general_categories 3A and 3B are "amended", not "control"
all.fung.1920.ps@sam_data[which(all.fung.1920.ps@sam_data$state=="US" & all.fung.1920.ps@sam_data$treatment_new==3), "general_category"] <- "Amended"
```

For treatments that were fumigated, fumigation would have happened fall
2018 for the 3-year rotations, and fall 2019 for the two-year rotations.
Add in another category “post.fumigation” that applies to in.potato
samples from spring and summer that were fumigated the previous fall.

``` r
table(all.fung.1920.ps@sam_data %>% pull(general_category)) # breakdown of our general categories
```

    ## 
    ##           Amended Amended/Fumigated           Control         Fumigated 
    ##               410               124               656               282

``` r
# samples considered post.fumigation TRUE are 2019 3-yr rotation spring/summer, or 2020 2-yr spring/summer, and general_category = Fumigated or Amended/Fumigated
all.fung.1920.ps@sam_data$post.fumigation <- FALSE # write a new sample data column, post.fumigation is FALSE

all.fung.1920.ps@sam_data[which(all.fung.1920.ps@sam_data$rotation == 3 & all.fung.1920.ps@sam_data$year == 19 &
                                grepl("umigated", all.fung.1920.ps@sam_data$general_category) &
                                startsWith(all.fung.1920.ps@sam_data$season, "S")), "post.fumigation"] <- TRUE # post fumigation true for 3-yr rotations

all.fung.1920.ps@sam_data[which(all.fung.1920.ps@sam_data$rotation == 2 & all.fung.1920.ps@sam_data$year == 20 &
                                grepl("umigated", all.fung.1920.ps@sam_data$general_category) &
                                startsWith(all.fung.1920.ps@sam_data$season, "S")), "post.fumigation"] <- TRUE # post fumigation true for 2-yr rotations
```

*We now have several variables by which to examine community
variation:*  
state, year, season, block, in.potato, and post.fumigation  
Keep in mind ALL post.fumigation samples are Potato, and that only
Potato samples should be compared for fumigation effects.

``` r
summary(data.frame(all.fung.1920.ps@sam_data) %>% select(state, year, season, block, in.potato, post.fumigation)) # our variables summary
```

    ##     state                year          season              block      
    ##  Length:1476        Min.   :18.00   Length:1476        Min.   :1.000  
    ##  Class :character   1st Qu.:19.00   Class :character   1st Qu.:2.000  
    ##  Mode  :character   Median :20.00   Mode  :character   Median :3.000  
    ##                     Mean   :19.49                      Mean   :2.923  
    ##                     3rd Qu.:20.00                      3rd Qu.:4.000  
    ##                     Max.   :20.00                      Max.   :5.000  
    ##                     NA's   :3                          NA's   :183    
    ##   in.potato         post.fumigation
    ##  Length:1476        Mode :logical  
    ##  Class :character   FALSE:1198     
    ##  Mode  :character   TRUE :278      
    ##                                    
    ##                                    
    ##                                    
    ## 

``` r
table(all.fung.1920.ps@sam_data$state) 
```

    ## 
    ##  CO  ID  ME  MI  MN  ND  OR  US  WI 
    ## 180 120 180 150 220 179 144 150 150

``` r
table(all.fung.1920.ps@sam_data$year) 
```

    ## 
    ##  18  19  20 
    ##  40 675 758

``` r
table(all.fung.1920.ps@sam_data$season)
```

    ## 
    ##   Fall Spring Summer 
    ##    464    504    505

``` r
table(all.fung.1920.ps@sam_data$block)
```

    ## 
    ##   1   2   3   4   5 
    ## 268 257 286 271 211

``` r
table(all.fung.1920.ps@sam_data$in.potato)
```

    ## 
    ##   Potato Rotation 
    ##     1009      463

``` r
table(all.fung.1920.ps@sam_data$post.fumigation)
```

    ## 
    ## FALSE  TRUE 
    ##  1198   278

Different numbers of samples per state and block are roughly comparable
and are probably trustworthy with permanovas. Season is ok too if
comparing Spring vs Summer, which are in potato. (in.potato itself is
confounded by season, Spring/Summer potato vs Fall rotation samples).
For year, make sure to omit 2018 Minnesota samples. Be more careful with
comparisons made between fumigated/unfumigated treatments, because
sample sizes are definitely not balanced.

#### Viridiplantae: to remove or not to remove? Some plots

``` r
all.fung.1920.ps <- number_asvs(all.fung.1920.ps) # store the sequences as a biostrings object, and assign ASV numbers for ease

table(data.frame(all.fung.1920.ps@tax_table) %>% filter(Kingdom=="Viridiplantae") %>% pull(Phylum)) # the major ASVs are Anthophyta (real plants) but there are many ASVs belonging to Chlorophyta (green algae)
```

    ## 
    ##                       Anthophyta                        Bryophyta 
    ##                              464                              151 
    ##                      Chlorophyta                   Lycopodiophyta 
    ##                             2479                                1 
    ##                  Marchantiophyta Viridiplantae_phy_Incertae_sedis 
    ##                                1                                7

``` r
### plotting percent abundances of Viridiplantae by Phylum
# tax glom at phylum
all.fung.1920.phy <- tax_glom(all.fung.1920.ps, taxrank = "Phylum") # 

# transform to percent abundance
all.fung.1920.phy.pa <- transform_sample_counts(all.fung.1920.phy, function(OTU) OTU/sum(OTU) * 100)

# subset only the kingdom Viridiplantae
all.fung.1920.phy.pa <- subset_taxa(all.fung.1920.phy.pa, Kingdom == "Viridiplantae")

# plot percent abundances of phyla from kingdom Viridiplantae
plot_bar(all.fung.1920.phy.pa, fill="Phylum")+
   geom_bar(stat="identity", position="stack")+ 
   scale_fill_manual("Phylum", values = natparks.pals("Acadia", length(unique(all.fung.1920.phy.pa@tax_table[,"Phylum"]))))+
   facet_grid(~state, scales = "free", space = "free")+
   scale_y_continuous(limits = c(0,70))+
   ggtitle("")+
   theme(axis.text.x = element_blank())
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
### plotting percent abundances of different taxonomic levels across samples

# first the major phyla of Viridiplantae
gg.anth <- ggplot(compare_taxa_abunds(all.fung.1920.ps, "Phylum", "Anthophyta"), aes(state, percent.per.sample))+
  geom_jitter()+ggtitle("Anthophyta percent abundances")

gg.chlo <- ggplot(compare_taxa_abunds(all.fung.1920.ps, "Phylum", "Chlorophyta"), aes(state, percent.per.sample))+
  geom_jitter()+ggtitle("Chlorophyta percent abundances")

gg.anth + gg.chlo
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
gg.solan <- ggplot(compare_taxa_abunds(all.fung.1920.ps, "Genus", "Solanum"), aes(state, percent.per.sample, color=in.potato))+
  geom_jitter()+ggtitle("Solanum percent abundances")

gg.brass <- ggplot(compare_taxa_abunds(all.fung.1920.ps, "Genus", "Brassica"), aes(state, percent.per.sample, color=in.potato))+
  geom_jitter()+ggtitle("Brassica percent abundances")

gg.trit <- ggplot(compare_taxa_abunds(all.fung.1920.ps, "Genus", "Triticum"), aes(state, percent.per.sample, color=in.potato))+
  geom_jitter()+ggtitle("Triticum percent abundances")

gg.solan + gg.brass + gg.trit + plot_layout(guides = "collect")
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->
Most of the Viridiplantae actually are Chlorophyta, green algae. Is this
biologically relevant?  
The most abundant orders are: Solanales (no surprise), Caryophyllales,
Poales, and Brassicales. Solanales: pretty much all potatoes as far as I
can tell (mainly tuberosum or NA, some bukasovii)  
Caryophyllales: Amaranths, like Bassia scoparia (a weed), Chenopodium
quinoa/acuminatum (quinoa and a wild relative?), Amaranthus deflexus (an
amaranth), Salsola collina (saltwort), Fagopyrum (buckwheat), Portulaca
oleraca (purslane)… so many of them may be weeds Poales: Triticum
(wheat), Hordeum (barley) Brassicales: Raphanus (radish), Brassica
carinata/napus/juncea (rape/mustard). So basically they represent target
crops (potatoes), weeds, and rotation crops.

From the plots above you can guess that Colorado is rotating with
Brassica, and Michigan with wheat.

#### Chlorophyta: big orders are Chlymadomonadales, Sphaeropleales

``` r
gg.chlam <- ggplot(compare_taxa_abunds(all.fung.1920.ps, "Order", "Chlamydomonadales"), aes(state, percent.per.sample))+
  geom_jitter()+ggtitle("Chlamydomonadales percent abundances")
gg.sphae <- ggplot(compare_taxa_abunds(all.fung.1920.ps, "Order", "Sphaeropleales"), aes(state, percent.per.sample))+
  geom_jitter()+ggtitle("Sphaeropleales percent abundances")
gg.chlam + gg.sphae
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->
\#### Library size (before pruning any plant ITS sequences)

``` r
plot_libsize(all.fung.1920.ps)+ggtitle("ITS library sizes, 2018-2020")
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->
\#### Cleanup Remove NA samples, prune at some threshold, omit Phylum
Anthophyta

``` r
# remove OR blanks
all.fung.1920.c.ps <- subset_samples(all.fung.1920.ps, sample_names(all.fung.1920.ps) != 
                                       c("OR_1_2_blank_20_4_ITS_S1039", "OR_1_2_blank_20_6_ITS_S1040", "OR_1_3_blank_20_8_ITS_S1065")) 
# omit sequences belonging to true plants 
all.fung.1920.c.ps <- subset_taxa(all.fung.1920.c.ps, (Phylum!="Anthophyta")) 

plot_libsize(all.fung.1920.c.ps)+ggtitle("ITS library sizes from cleaned samples, 2018-2020") # looks nearly identical as previous
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
# prune samples below 10k reads
all.fung.1920.c.ps <- subset_samples(all.fung.1920.c.ps, sample_sums(all.fung.1920.c.ps) > 10000) # throws away 118 samples yikes

# write out phyloseq object
#saveRDS(all.fung.1920.c.ps, "/Users/scottklasek/Desktop/UMN/phyloseqs/fung1920.cleaned.ps") # cleaned fungal phyloseq object

# check the balance of sample numbers
table(all.fung.1920.c.ps@sam_data$state) / table(all.fung.1920.ps@sam_data$state) # proportions of libraries > 10k reads, by state
```

    ## 
    ##        CO        ID        ME        MI        MN        ND        OR        US 
    ## 0.9555556 0.9083333 0.8888889 0.9333333 0.9045455 0.9385475 0.9583333 0.9466667 
    ##        WI 
    ## 0.8466667

``` r
table(all.fung.1920.c.ps@sam_data$year) / table(all.fung.1920.ps@sam_data$year)
```

    ## 
    ##        18        19        20 
    ## 0.9750000 0.9274074 0.9102902

``` r
table(all.fung.1920.c.ps@sam_data$season) / table(all.fung.1920.ps@sam_data$season)
```

    ## 
    ##      Fall    Spring    Summer 
    ## 0.9116379 0.9305556 0.9168317

``` r
table(all.fung.1920.c.ps@sam_data$block) / table(all.fung.1920.ps@sam_data$block)
```

    ## 
    ##         1         2         3         4         5 
    ## 0.9291045 0.8949416 0.8811189 0.9261993 0.9526066

``` r
table(all.fung.1920.c.ps@sam_data$in.potato) / table(all.fung.1920.ps@sam_data$in.potato)
```

    ## 
    ##    Potato  Rotation 
    ## 0.9236868 0.9136069

``` r
table(all.fung.1920.c.ps@sam_data$post.fumigation) / table(all.fung.1920.ps@sam_data$post.fumigation)
```

    ## 
    ##     FALSE      TRUE 
    ## 0.9173623 0.9208633

Doesn’t seem to be huge bias in the libraries removed due to low reads,
except mayyyybe by state?

### ANALYSIS TIME

#### ITS alpha diversity

``` r
# calculate the alpha diversity (three metrics)
fung.div <- estimate_richness(all.fung.1920.c.ps, measures = c("Observed", "Shannon", "InvSimpson")) 

# merge select sample data columns with the alpha diversity data, and clean up the dataframe
d <- merge(fung.div, data.frame(all.fung.1920.c.ps@sam_data)[,c("state", "year", "month", "season", "block", "cultivar", "in.potato", "OM....", "post.fumigation")], by=0)
rownames(d) <- d$Row.names
d <- d %>% select(-Row.names)

# plot the three metrics I've calculated
cols.threediscrete <- c("#1b9e77", "#d95f02", "#7570b3")
f1920obs <- ggplot(d, aes(state, Observed, color=as.character(year)))+geom_jitter(width=0.3, size = 0.9)+
  scale_x_discrete("", limits = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))+
  scale_y_continuous("")+
  scale_color_manual("Year", values = cols.threediscrete)+
  theme_bw()#+ggtitle("", subtitle = "B")

f1920sha <- ggplot(d, aes(state, Shannon, color=as.character(year)))+geom_jitter(width=0.3, size = 0.9)+
  scale_x_discrete("", limits = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))+
  scale_y_continuous("")+
  scale_color_manual("Year", values = cols.threediscrete)+
  theme_bw()+theme(legend.position = "none")#+ggtitle("", subtitle = "D")

f1920sim <- ggplot(d, aes(state, InvSimpson, color=as.character(year)))+geom_jitter(width=0.3, size = 0.9)+
  scale_x_discrete("", limits = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2"))+
  scale_y_continuous("")+
  scale_color_manual("Year", values = cols.threediscrete)+
  theme_bw()+theme(legend.position = "none")#+ggtitle("", subtitle = "F")

f1920obs / f1920sha / f1920sim + plot_layout(guides = "collect")
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
Now for some follow-up stats on the alpha diversity

``` r
summary(aov(Observed~state, data=d)) # state is significant across all
```

    ##               Df   Sum Sq Mean Sq F value Pr(>F)    
    ## state          8 39459660 4932457   220.7 <2e-16 ***
    ## Residuals   1346 30085947   22352                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(Shannon~state, data=d))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)    
    ## state          8  221.8  27.725   165.1 <2e-16 ***
    ## Residuals   1346  226.0   0.168                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(InvSimpson~state, data=d))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)    
    ## state          8 347850   43481   120.8 <2e-16 ***
    ## Residuals   1346 484590     360                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(Observed~year, data=d)) # year is significant across all
```

    ##               Df   Sum Sq Mean Sq F value Pr(>F)    
    ## year           1  5792137 5792137   122.9 <2e-16 ***
    ## Residuals   1353 63753470   47120                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(Shannon~year, data=d))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)    
    ## year           1   33.7   33.73   110.2 <2e-16 ***
    ## Residuals   1353  414.1    0.31                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(InvSimpson~year, data=d))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)    
    ## year           1  59853   59853   104.8 <2e-16 ***
    ## Residuals   1353 772586     571                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# for all POTATO soils in states that fumigated, are there differences across fumigations?
unique(d %>% filter(post.fumigation==TRUE) %>% pull(state)) # all states except CO and US fumigated
```

    ## [1] "ID" "ME" "MI" "MN" "ND" "OR" "WI"

``` r
alpha.fum.d <- d %>% filter(state != "US" & state != "CO" & in.potato == "Potato") # write them in a separate dataframe

# apparently WI fumigated the entire field before planting in 2020, change this:
alpha.fum.d[which(alpha.fum.d$state=="WI" & alpha.fum.d$year==20),"post.fumigation"] <- TRUE
alpha.fum.d[which(alpha.fum.d$state=="WI" & alpha.fum.d$year==19),"post.fumigation"] <- FALSE

# and ND did the reverse
alpha.fum.d[which(alpha.fum.d$state=="ND" & alpha.fum.d$year==19),"post.fumigation"] <- TRUE
alpha.fum.d[which(alpha.fum.d$state=="ND" & alpha.fum.d$year==20),"post.fumigation"] <- FALSE

table(alpha.fum.d$post.fumigation) # remember nearly twice as many communities are unfumigated
```

    ## 
    ## FALSE  TRUE 
    ##   411   308

``` r
t.test(alpha.fum.d %>% filter(post.fumigation == TRUE) %>% pull(Observed),
       alpha.fum.d %>% filter(post.fumigation == FALSE) %>% pull(Observed))
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  alpha.fum.d %>% filter(post.fumigation == TRUE) %>% pull(Observed) and alpha.fum.d %>% filter(post.fumigation == FALSE) %>% pull(Observed)
    ## t = -6.027, df = 635.2, p-value = 2.832e-09
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -129.60470  -65.90402
    ## sample estimates:
    ## mean of x mean of y 
    ##  720.8734  818.6277

``` r
t.test(alpha.fum.d %>% filter(post.fumigation == TRUE) %>% pull(Shannon),
       alpha.fum.d %>% filter(post.fumigation == FALSE) %>% pull(Shannon))
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  alpha.fum.d %>% filter(post.fumigation == TRUE) %>% pull(Shannon) and alpha.fum.d %>% filter(post.fumigation == FALSE) %>% pull(Shannon)
    ## t = -10.279, df = 431.38, p-value < 2.2e-16
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.5703920 -0.3872763
    ## sample estimates:
    ## mean of x mean of y 
    ##  4.500652  4.979486

``` r
t.test(alpha.fum.d %>% filter(post.fumigation == TRUE) %>% pull(InvSimpson),
       alpha.fum.d %>% filter(post.fumigation == FALSE) %>% pull(InvSimpson))
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  alpha.fum.d %>% filter(post.fumigation == TRUE) %>% pull(InvSimpson) and alpha.fum.d %>% filter(post.fumigation == FALSE) %>% pull(InvSimpson)
    ## t = -13.065, df = 675.53, p-value < 2.2e-16
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -25.35359 -18.72861
    ## sample estimates:
    ## mean of x mean of y 
    ##  35.73100  57.77211

``` r
gg.fum.alpha <- ggplot(alpha.fum.d %>% pivot_longer(cols = c("Observed", "Shannon", "InvSimpson"), names_to = "diversity"), 
       aes(state, value, color=post.fumigation))+
  geom_jitter(width=0.3, size = 0.9)+
  scale_x_discrete("State", limits = c("OR", "ID", "ND", "MN", "WI", "MI", "ME"))+
  scale_y_continuous("Diversity metric")+
  scale_color_discrete("Fumigation")+
  facet_grid(diversity~factor(year, levels = c(19,20), labels = c("2019", "2020")), scales = "free")+
  theme_bw()+ggtitle("Eukaryotic alpha diversity")
gg.fum.alpha
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
# Michigan
# t.test(alpha.fum.d %>% filter(state=="MI" & year==19) %>% pull(InvSimpson), # generally little difference in alpha diversity by year
#        alpha.fum.d %>% filter(state=="MI" & year==20) %>% pull(InvSimpson))

t.test(alpha.fum.d %>% filter(state=="MI" & post.fumigation == FALSE) %>% pull(Observed),
       alpha.fum.d %>% filter(state=="MI" & post.fumigation == TRUE) %>% pull(Observed))
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  alpha.fum.d %>% filter(state == "MI" & post.fumigation == FALSE) %>% pull(Observed) and alpha.fum.d %>% filter(state == "MI" & post.fumigation == TRUE) %>% pull(Observed)
    ## t = 4.3278, df = 35.515, p-value = 0.000117
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##   75.53206 208.87280
    ## sample estimates:
    ## mean of x mean of y 
    ## 1042.8947  900.6923

``` r
t.test(alpha.fum.d %>% filter(state=="MI" & post.fumigation == FALSE) %>% pull(Shannon),
       alpha.fum.d %>% filter(state=="MI" & post.fumigation == TRUE) %>% pull(Shannon))
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  alpha.fum.d %>% filter(state == "MI" & post.fumigation == FALSE) %>% pull(Shannon) and alpha.fum.d %>% filter(state == "MI" & post.fumigation == TRUE) %>% pull(Shannon)
    ## t = 5.5781, df = 57.389, p-value = 6.864e-07
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  0.2183240 0.4628047
    ## sample estimates:
    ## mean of x mean of y 
    ##  5.206113  4.865549

``` r
t.test(alpha.fum.d %>% filter(state=="MI" & post.fumigation == FALSE) %>% pull(InvSimpson),
       alpha.fum.d %>% filter(state=="MI" & post.fumigation == TRUE) %>% pull(InvSimpson))
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  alpha.fum.d %>% filter(state == "MI" & post.fumigation == FALSE) %>% pull(InvSimpson) and alpha.fum.d %>% filter(state == "MI" & post.fumigation == TRUE) %>% pull(InvSimpson)
    ## t = 3.2999, df = 32.672, p-value = 0.002344
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##   5.673404 23.935308
    ## sample estimates:
    ## mean of x mean of y 
    ##  56.62755  41.82320

``` r
# Minnesota
# t.test(alpha.fum.d %>% filter(state=="MN" & year==19) %>% pull(Shannon), # alpha diversity is different by year
#        alpha.fum.d %>% filter(state=="MN" & year==20) %>% pull(Shannon))

t.test(alpha.fum.d %>% filter(state=="MN" & post.fumigation == FALSE) %>% pull(Observed),
       alpha.fum.d %>% filter(state=="MN" & post.fumigation == TRUE) %>% pull(Observed))
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  alpha.fum.d %>% filter(state == "MN" & post.fumigation == FALSE) %>% pull(Observed) and alpha.fum.d %>% filter(state == "MN" & post.fumigation == TRUE) %>% pull(Observed)
    ## t = 2.7033, df = 74.406, p-value = 0.008503
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##   23.9343 158.0891
    ## sample estimates:
    ## mean of x mean of y 
    ##  770.9722  679.9605

``` r
t.test(alpha.fum.d %>% filter(state=="MN" & post.fumigation == FALSE) %>% pull(Shannon),
       alpha.fum.d %>% filter(state=="MN" & post.fumigation == TRUE) %>% pull(Shannon))
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  alpha.fum.d %>% filter(state == "MN" & post.fumigation == FALSE) %>% pull(Shannon) and alpha.fum.d %>% filter(state == "MN" & post.fumigation == TRUE) %>% pull(Shannon)
    ## t = 5.1218, df = 87.064, p-value = 1.796e-06
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  0.2083147 0.4725248
    ## sample estimates:
    ## mean of x mean of y 
    ##  4.981675  4.641255

``` r
t.test(alpha.fum.d %>% filter(state=="MN" & post.fumigation == FALSE) %>% pull(InvSimpson),
       alpha.fum.d %>% filter(state=="MN" & post.fumigation == TRUE) %>% pull(InvSimpson))
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  alpha.fum.d %>% filter(state == "MN" & post.fumigation == FALSE) %>% pull(InvSimpson) and alpha.fum.d %>% filter(state == "MN" & post.fumigation == TRUE) %>% pull(InvSimpson)
    ## t = 5.8295, df = 66.169, p-value = 1.81e-07
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  13.74683 28.06699
    ## sample estimates:
    ## mean of x mean of y 
    ##  56.66019  35.75328

``` r
# ND
t.test(alpha.fum.d %>% filter(state=="ND" & post.fumigation == FALSE) %>% pull(Observed),
       alpha.fum.d %>% filter(state=="ND" & post.fumigation == TRUE) %>% pull(Observed))
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  alpha.fum.d %>% filter(state == "ND" & post.fumigation == FALSE) %>% pull(Observed) and alpha.fum.d %>% filter(state == "ND" & post.fumigation == TRUE) %>% pull(Observed)
    ## t = -8.4423, df = 105.24, p-value = 1.831e-13
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -183.7921 -113.8806
    ## sample estimates:
    ## mean of x mean of y 
    ##  608.1273  756.9636

``` r
t.test(alpha.fum.d %>% filter(state=="ND" & post.fumigation == FALSE) %>% pull(Shannon),
       alpha.fum.d %>% filter(state=="ND" & post.fumigation == TRUE) %>% pull(Shannon))
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  alpha.fum.d %>% filter(state == "ND" & post.fumigation == FALSE) %>% pull(Shannon) and alpha.fum.d %>% filter(state == "ND" & post.fumigation == TRUE) %>% pull(Shannon)
    ## t = -6.2553, df = 101.04, p-value = 9.628e-09
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.3456759 -0.1792183
    ## sample estimates:
    ## mean of x mean of y 
    ##  4.488124  4.750571

``` r
t.test(alpha.fum.d %>% filter(state=="ND" & post.fumigation == FALSE) %>% pull(InvSimpson),
       alpha.fum.d %>% filter(state=="ND" & post.fumigation == TRUE) %>% pull(InvSimpson))
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  alpha.fum.d %>% filter(state == "ND" & post.fumigation == FALSE) %>% pull(InvSimpson) and alpha.fum.d %>% filter(state == "ND" & post.fumigation == TRUE) %>% pull(InvSimpson)
    ## t = -3.2529, df = 105.5, p-value = 0.001535
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -8.909603 -2.161565
    ## sample estimates:
    ## mean of x mean of y 
    ##  36.10948  41.64507

Alpha diversity varies across states and years, and does not across
block (fortunately). Whether or not a plot was in potato was less clear,
with different metrics giving different results. In states that had
fumigation treatments and were in potato, t-tests say that fumigated
communities unambiguously (ie all three metrics agreed) had lower
measurements of alpha diversity, but wait a second– when plotting by
state and year, it seems like fumigation only reduces alpha diversity in
Minnesota and Michigan.

In Michigan, alpha diversity was largely not different across years, so
I t-tested all samples together and for each metric, fumigated samples
had lower diversity.  
In Minnesota, diversity was higher in 2019, so I tested each year
separately. For each metric, fumigated samples had lower diversity.

#### ordinations

``` r
# transform
all.fung.1920.c.ps.hel <- helltransform(all.fung.1920.c.ps) # hellinger transformation

# make distance matrices
dm.fung1920.bjac <- phyloseq::distance(all.fung.1920.c.ps, method = "jaccard", binary = TRUE) # binary jaccard
dm.fung1920.hel.bc <- phyloseq::distance(all.fung.1920.c.ps.hel, method = "bray") # bray-curtis with hellinger

# ordinate
ord.fung1920.hel.bc.pcoa <- ordinate(all.fung.1920.c.ps.hel, method="PCoA", distance=dm.fung1920.hel.bc) # hellinger
ord.fung1920.bjac.pcoa <- ordinate(all.fung.1920.c.ps, method="PCoA", distance=dm.fung1920.bjac) # binary jaccard

# plots
cols.alphabetical <- c("#236192", "#B3A369", "#B0D7FF", "#18453B", "#7A0019", "#FFC72A", "#DC4405", "gray", "#C5050C")
gg.bc.ord <- plot_ordination(all.fung.1920.c.ps.hel, ord.fung1920.hel.bc.pcoa, color="state", title="ITS PCoA, Bray-Curtis, Hellinger transformation")+
  scale_color_manual(values = cols.alphabetical)
gg.bjac.ord <- plot_ordination(all.fung.1920.c.ps, ord.fung1920.bjac.pcoa, color="state", title="ITS PCoA, Binary Jaccard")+
  scale_color_manual(values = cols.alphabetical)
gg.bc.ord + gg.bjac.ord + plot_layout(guides = "collect")
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->
Similar patterns seen as before, and the Larkin USDA samples do plot
closely with U of Maine samples as I suspected. Other variables don’t
show much, if any, separation.  
Now let’s see which variables explain variance among communities (this
chunk takes a while):

``` r
# first by state
adonis.state <- adonis2(formula = dm.fung1920.hel.bc ~ state, data = data.frame(all.fung.1920.c.ps@sam_data)) # 59.7% of variance is by state!
adonis.state
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = dm.fung1920.hel.bc ~ state, data = data.frame(all.fung.1920.c.ps@sam_data))
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## state       8   264.21 0.59769 249.96  0.001 ***
    ## Residual 1346   177.84 0.40231                  
    ## Total    1354   442.05 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# drop certain samples for certain comparisons

# blocks, because Colorado
subset.ps <- subset_samples(all.fung.1920.c.ps, state !="CO") # removing CO because it has no blocks
subset.ps <- helltransform(subset.ps) # hellinger transformation
dm.subset <- phyloseq::distance(subset.ps, method = "bray") # bray-curtis distance matrix
adonis.block <- adonis2(formula = dm.subset ~ block, data = data.frame(subset.ps@sam_data)) # compare blocks
adonis.block # block is 0.4% of variance, that's good
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = dm.subset ~ block, data = data.frame(subset.ps@sam_data))
    ##            Df SumOfSqs    R2      F Pr(>F)    
    ## block       1     1.47 0.004 4.7375  0.001 ***
    ## Residual 1181   365.36 0.996                  
    ## Total    1182   366.83 1.000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# year, because Minnesota 2018
subset.ps <- subset_samples(all.fung.1920.c.ps, year !=18) # removing MN 2018 samples
subset.ps <- helltransform(subset.ps) # hellinger transformation
dm.subset <- phyloseq::distance(subset.ps, method = "bray") # bray-curtis distance matrix
adonis.year <- adonis2(formula = dm.subset ~ year, data = data.frame(subset.ps@sam_data)) # compare years
adonis.year # 0.9% of variance by year
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = dm.subset ~ year, data = data.frame(subset.ps@sam_data))
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## year        1     3.88 0.00899 11.921  0.001 ***
    ## Residual 1314   427.53 0.99101                  
    ## Total    1315   431.41 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# season, in.potato
subset.ps <- subset_samples(all.fung.1920.c.ps, in.potato == "Potato") # removing fall samples that are not in potato
subset.ps <- helltransform(subset.ps) # hellinger transformation
dm.subset <- phyloseq::distance(subset.ps, method = "bray") # bray-curtis distance matrix
adonis.season <- adonis2(formula = dm.subset ~ season, data = data.frame(subset.ps@sam_data)) # compare spring/summer seasons
adonis.season # 0.6% of variance between spring/summer
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = dm.subset ~ season, data = data.frame(subset.ps@sam_data))
    ##           Df SumOfSqs      R2      F Pr(>F)    
    ## season     1    1.887 0.00612 5.7266  0.001 ***
    ## Residual 930  306.510 0.99388                  
    ## Total    931  308.397 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# post.fumigation, for in.potato (unbalanced!)
# fix the WI fumigation fiasco:
subset.ps@sam_data[which(subset.ps@sam_data$state=="WI" & subset.ps@sam_data$year==20),"post.fumigation"] <- TRUE
subset.ps@sam_data[which(subset.ps@sam_data$state=="WI" & subset.ps@sam_data$year==19),"post.fumigation"] <- FALSE
adonis.fum <- adonis2(formula = dm.subset ~ post.fumigation, data = data.frame(subset.ps@sam_data)) 
adonis.fum # 3.8% of variance
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = dm.subset ~ post.fumigation, data = data.frame(subset.ps@sam_data))
    ##                  Df SumOfSqs      R2      F Pr(>F)    
    ## post.fumigation   1   13.473 0.04369 42.485  0.001 ***
    ## Residual        930  294.924 0.95631                  
    ## Total           931  308.397 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# potato in spring/summer vs rotation in fall (unbalanced!)
adonis.potato <- adonis2(formula = dm.fung1920.hel.bc ~ in.potato, data = data.frame(all.fung.1920.c.ps@sam_data))
adonis.potato # 0.6% of variance 
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = dm.fung1920.hel.bc ~ in.potato, data = data.frame(all.fung.1920.c.ps@sam_data))
    ##             Df SumOfSqs      R2      F Pr(>F)    
    ## in.potato    1     2.56 0.00579 7.8744  0.001 ***
    ## Residual  1353   439.49 0.99421                  
    ## Total     1354   442.05 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

In conclusion: state is by far the largest contributor of variance, then
it seems as though fumigation matters a little. Then year at just under
1%, and block the least at 0.4%. There’s probably a better way to do
this somewhere, but the main point is clear for now.

#### occupancy

``` r
asvc.fungi <- asvcumulate(all.fung.1920.c.ps) # get prevalence and mean abundances of each fungal ASV
summary(asvc.fungi %>% pull(pct_prevalence)) # nobody is core, most are in only a few samples
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0000  0.2214  0.8118  3.3723  2.8044 99.4096

``` r
# count the number of states in which an ASV was detected
num_states_per_asv_test <- num_states_per_asv(all.fung.1920.c.ps) # took like 8 minutes

# put the num_state count into the asvc dataframe
test_merge <- merge(asvc.fungi, data.frame(num_states_per_asv_test), by=0) 

# plot abundance occupancy
aostate <- ggplot(test_merge, aes(log10(pctreads), pct_prevalence, color=num_states_per_asv_test))+
  geom_point()+
  scale_color_continuous("Number of states detected")+
  scale_x_continuous("log mean abundance")+
  scale_y_continuous("percent occupancy")+
  ggtitle("Fungal occupancy / abundance")+theme_bw()
aostate
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
hist(test_merge %>% filter(pctreads > 0.1) %>% pull(pct_prevalence), breaks=25,
     main = "Occupancy of high-abundance ASVs", xlab = "Percent occupancy")
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->
Clear bimodal pattern: The ASVs on the left are really high-abundance,
but only in samples from a few states, while the ones on the right are
in pretty high abundance across the continent. 60% occupancy seems like
a natural cutoff.

``` r
ggplot(test_merge, aes(num_states_per_asv_test, log10(pctreads)))+
  geom_jitter()+
  scale_y_continuous("log10 percent abundance")+
  scale_x_continuous("Number of states ASV found in", breaks = c(1:9))+
  theme_bw()
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

#### bar plots of cosmopolitan and state/regionally specific taxa

``` r
# remove nonexistent ASVs that must have been from the blank samples
filter <- phyloseq::genefilter_sample(all.fung.1920.c.ps, filterfun_sample(function(x) x > 0)) # for each ASV get a true/false whether it has over 0 reads
all.fung.1920.c.ps <- prune_taxa(filter, all.fung.1920.c.ps) # remove the zero-count taxa from the ps object

# transform to percent abundance
fung1920.pa <- transform_sample_counts(all.fung.1920.c.ps, function(OTU) OTU/sum(OTU) * 100)

# try putting the asvc (as a matrix) into the tax table of the resulting object and graphing by "taxa" as number of states each ASV was found in
rownames(test_merge) <- test_merge$Row.names
state.count <- test_merge %>% select(num_states_per_asv_test) # does this need to be in the same order as taxa are in ps?
state.count$num_states_per_asv_test <- as.character(state.count$num_states_per_asv_test) # needs to be a character matrix for phyloseq tax table

# warning: stupid bullshit
# there's a glitch in phyloseq where you can't plot or subset by the leftmost column of a tax table, so I need to throw in another bogus one.
dummy <- rep("Eukarya", times=nrow(fung1920.pa@tax_table))
dummy <- as.matrix(dummy)
rownames(dummy) <- rownames(fung1920.pa@tax_table)
tmerge <- merge(dummy, as.matrix(state.count), by=0)
rownames(tmerge) <- tmerge$Row.names
tmerge <- tmerge %>% select(V1, num_states_per_asv_test)
# stupid bullshit over

# add in the stupid dataframe containing state-counts as our tax table
fung1920.pa.statecount <- phyloseq(tax_table(as.matrix(tmerge)),
                                   otu_table(fung1920.pa),
                                   sample_data(fung1920.pa))

# so we tax-glom before plotting to make it easier on our hardworking computer
fung1920.pa.statecount <- tax_glom(fung1920.pa.statecount, taxrank = "num_states_per_asv_test")

gg.bar.statecount <- plot_bar(fung1920.pa.statecount, fill="num_states_per_asv_test")+
   geom_bar(stat="identity", position="stack")+
   scale_fill_manual("ASVs present in\n this many states", values = natparks.pals("Acadia", 9))+
   facet_grid(~factor(state, levels = c("OR", "ID", "CO", "MN", "ND", "WI", "MI", "ME", "US"), labels = c("OR", "ID", "CO", "MN1", "MN2", "WI", "MI", "ME1", "ME2")), scales = "free", space = "free")+
   scale_y_continuous(limits = c(0,100))+
   ggtitle("")+
   theme(axis.text.x = element_blank())
gg.bar.statecount
```

    ## Warning: Removed 1 rows containing missing values (`geom_bar()`).
    ## Removed 1 rows containing missing values (`geom_bar()`).

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
# what is causing the variability in some states? (ME, ND, and WI in particular)
fung1920.pa.statecount.wi.20 <- subset_samples(fung1920.pa.statecount, state=="WI" & year == 20)
wi20state <- plot_bar(fung1920.pa.statecount.wi.20, fill="num_states_per_asv_test")+
   geom_bar(stat="identity", position="stack")+
   scale_fill_manual("ASVs present in\n this many states", values = natparks.pals("Acadia", 9))+
   facet_grid(~factor(season, levels = c("Spring", "Summer", "Fall")), scales = "free", space = "free")+
   scale_y_continuous(limits = c(0,100))+
   ggtitle("WI 2020 samples")+
   theme(axis.text.x = element_blank())
fung1920.pa.statecount.wi.19 <- subset_samples(fung1920.pa.statecount, state=="WI" & year == 19)
wi19state <- plot_bar(fung1920.pa.statecount.wi.19, fill="num_states_per_asv_test")+
   geom_bar(stat="identity", position="stack")+
   scale_fill_manual("ASVs present in\n this many states", values = natparks.pals("Acadia", 9))+
   facet_grid(~factor(season, levels = c("Spring", "Summer", "Fall")), scales = "free", space = "free")+
   scale_y_continuous(limits = c(0,100))+
   ggtitle("WI 2019 samples")+
   theme(axis.text.x = element_blank())
wi19state + wi20state + plot_layout(guides = "collect")
```

    ## Warning: Removed 1 rows containing missing values (`geom_bar()`).
    ## Removed 1 rows containing missing values (`geom_bar()`).

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

``` r
fung1920.pa.statecount.mi <- subset_samples(fung1920.pa.statecount, state=="MI" & year==20)
mi.state <- plot_bar(fung1920.pa.statecount.mi, fill="num_states_per_asv_test")+
   geom_bar(stat="identity", position="stack")+
   scale_fill_manual("ASVs present in\n this many states", values = natparks.pals("Acadia", 9))+
   facet_grid(~post.fumigation, scales = "free", space = "free")+
   scale_y_continuous(limits = c(0,100))+
   ggtitle("MI 2020 samples by fumigation")+
   theme(axis.text.x = element_blank())
mi.state
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-18-3.png)<!-- -->

#### rarefaction curves for WI diversity

``` r
wi.rarefy.df <- rarecurve(otu_table(subset_samples(all.fung.1920.c.ps, state=="WI")), step=50, cex=0.5, tidy = TRUE) # make a dataframe of rarefaction counts for WI samples
wi.rarefy.df$year <- as.numeric(substr(wi.rarefy.df$Site, 12, 13)) # assign the year
wi.rarefy.df[which(is.na(wi.rarefy.df$year)),"year"] <- 20
wi.rarefy.df$month <- as.numeric(gsub("_.*$","",sub("..............", "", wi.rarefy.df$Site))) # assign the month
table(wi.rarefy.df$month)
unique(wi.rarefy.df %>% filter(month==2) %>% pull(Site))
wi.rarefy.df[which(wi.rarefy.df$month==2),"month"] <- 11 # fix mislabeled month
wi.rarefy.df$postfum <- ifelse(wi.rarefy.df$year==19, "unfumigated",
                               ifelse(wi.rarefy.df$year==20 & wi.rarefy.df$month < 11, "fumigated", "unfumigated")) # fumigated category

# plot
wi.raref.19 <- ggplot(wi.rarefy.df %>% filter(year==19), aes(Sample, Species, color=as.factor(month)))+
  geom_line(size=0.07)+
  scale_color_discrete("month")+
  scale_y_continuous("ASV richness")+
  scale_x_continuous("read number")+
  ggtitle("WI 2019 ITS")+
  geom_vline(xintercept = 10000, linetype="dotted")+
  theme_bw()

wi.raref.20 <- ggplot(wi.rarefy.df %>% filter(year==20), aes(Sample, Species, color=as.factor(month)))+
  geom_line(size=0.07)+
  scale_color_discrete("month")+
  scale_y_continuous("ASV richness")+
  scale_x_continuous("read number")+
  ggtitle("WI 2020 ITS")+
  geom_vline(xintercept = 10000, linetype="dotted")+
  theme_bw()

wi.raref.19 + wi.raref.20  
```

This makes it seem like post-fumigation really did a number on the
diversity, but remember this is CONFOUNDED with field effects. The
two-year field could just have much lower diversity for other reasons…
but this is still very low diversity and I’m tempted to think fumigation
had to have made a difference.

#### Make an upset plot

``` r
all.states.asvs <- rownames(test_merge %>% arrange(rank) %>% filter(num_states_per_asv_test == 9)) # sequences present in each state

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
venn_obj.binary <- as.data.frame(venn_obj.binary) # creates a dataframe where rows are ASVs, columns are block_env categories, and entries are either 1 or 0
upset_order <- colnames(venn_obj.binary) # assigns the order of sets in the upset plot as the column names

# plot
state.upset <- upset(venn_obj.binary, 
      sets = rev(upset_order),
      mainbar.y.label = 'ASV count',
      sets.x.label = 'Observed ASVs by site',
      order.by = 'freq', nintersects = 30, # can change order.by to "degree", or adjust number of intersections to show 
      keep.order = T)
state.upset # 
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->
\#### Who are the cosmopolitan fungi/eukaryotes?

``` r
all.state.asvs.fung.1920.pa <- prune_taxa(all.states.asvs, fung1920.pa) # ps object with only the 9-state ASVs
length(unique(tax_table(all.state.asvs.fung.1920.pa)[,"Kingdom"]))
```

    ## [1] 4

``` r
# plot by kingdom
plot_bar(all.state.asvs.fung.1920.pa, fill="Kingdom")+
   geom_bar(stat="identity", position="stack")+
   scale_fill_manual("ASVs present in\n this many states", values = natparks.pals("Yellowstone", 4))+
   facet_grid(~state, scales = "free", space = "free")+
   scale_y_continuous(limits = c(0,100))+
   ggtitle("Cosmopolitan Eukaryotic ASVs")+
   theme(axis.text.x = element_blank())
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
not.fungi <- rownames(data.frame(tax_table(all.state.asvs.fung.1920.pa)) %>% filter(Kingdom !="Fungi"))
all.state.asvs.fung.1920.pa@tax_table[not.fungi,]
```

    ## Taxonomy Table:     [4 taxa by 7 taxonomic ranks]:
    ##         Kingdom         Phylum        Class             Order              
    ## ASV135  "Metazoa"       "Nematoda"    NA                NA                 
    ## ASV716  "Viridiplantae" "Chlorophyta" "Chlorophyceae"   "Chlamydomonadales"
    ## ASV1141 "Alveolata"     "Ciliophora"  "Phyllopharyngea" "Cyrtophorida"     
    ## ASV1317 "Alveolata"     "Ciliophora"  "Spirotrichea"    NA                 
    ##         Family              Genus               Species
    ## ASV135  NA                  NA                  NA     
    ## ASV716  "Chlorosarcinaceae" "Chlorosarcinopsis" "eremi"
    ## ASV1141 NA                  NA                  NA     
    ## ASV1317 NA                  NA                  NA

``` r
# only fungi
all.state.asvs.actuallyfungi.1920.pa <- prune_taxa(setdiff(all.states.asvs, not.fungi), fung1920.pa) # ps
plot_bar(all.state.asvs.actuallyfungi.1920.pa , fill="Phylum")+
   geom_bar(stat="identity", position="stack")+
   scale_fill_manual("Phylum", values = natparks.pals("GrandCanyon", 6))+
   facet_grid(~state, scales = "free", space = "free")+
   scale_y_continuous(limits = c(0,80))+
   ggtitle("Cosmopolitan Fungal ASVs")+
   theme(axis.text.x = element_blank())
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-21-2.png)<!-- -->

``` r
# by Class
plot_bar(all.state.asvs.actuallyfungi.1920.pa , fill="Class")+
   geom_bar(stat="identity", position="stack")+
   #scale_fill_manual("Phylum", values = natparks.pals("GrandCanyon", 6))+
   facet_grid(~state, scales = "free", space = "free")+
   scale_y_continuous(limits = c(0,80))+
   ggtitle("Cosmopolitan Fungal ASVs")+
   theme(axis.text.x = element_blank())
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-21-3.png)<!-- -->

``` r
all.state.asvs.fung.1920.pa@tax_table[rownames(test_merge[all.states.asvs,])[1:7],] # taxonomy of the most abundant cosmopolitan fungi
```

    ## Taxonomy Table:     [7 taxa by 7 taxonomic ranks]:
    ##        Kingdom Phylum              Class                Order                
    ## ASV578 "Fungi" "Ascomycota"        "Leotiomycetes"      "Thelebolales"       
    ## ASV77  "Fungi" "Ascomycota"        "Sordariomycetes"    "Glomerellales"      
    ## ASV17  "Fungi" "Mortierellomycota" "Mortierellomycetes" "Mortierellales"     
    ## ASV1   "Fungi" "Ascomycota"        "Sordariomycetes"    "Hypocreales"        
    ## ASV4   "Fungi" "Basidiomycota"     "Tremellomycetes"    "Cystofilobasidiales"
    ## ASV878 "Fungi" "Mortierellomycota" "Mortierellomycetes" "Mortierellales"     
    ## ASV516 "Fungi" "Mortierellomycota" "Mortierellomycetes" "Mortierellales"     
    ##        Family                 Genus              Species      
    ## ASV578 "Pseudeurotiaceae"     "Pseudogymnoascus" "pannorum"   
    ## ASV77  "Plectosphaerellaceae" "Gibellulopsis"    "piscis"     
    ## ASV17  "Mortierellaceae"      "Mortierella"      "elongata"   
    ## ASV1   "Nectriaceae"          "Gibberella"       "intricans"  
    ## ASV4   "Mrakiaceae"           "Tausonia"         "pullulans"  
    ## ASV878 "Mortierellaceae"      "Mortierella"      "minutissima"
    ## ASV516 "Mortierellaceae"      "Mortierella"      "hyalina"

Cosmopolitan eukaryotes consist largely of fungi, but there is one
Nematode, two members of Ciliophora (protists) and that green algae
Chlorosarcinopsis eremi.

When looking at Fungal ASVs only, the most dominant phylum is
Ascomycota. Western states except OR also have lots of Basidiomycota.
Maine, USDA, and North Dakota also have some Mortierellomycota.

It seems like little is known about the first two, P. pannorum is
closely related to the fungus that causes white-nose disease in bats. M.
elongata is apparently PGP for many diverse plants, [more info
here](https://edis.ifas.ufl.edu/publication/SS679). It is highly
abundant in soils, associates with plant roots, and harbors bacteria
inside its hyphae. The bacteria can modulate host metabolism… and a PGP
mechanism is currently unknown. So this would be a very good candidate
to examine cross-domain interactions: Is M. elongata associated with
tuber yield? If so, which bacteria associate with it, or if only in
certain states, which bacteria might be associated with it that could
drive a growth-association?

#### another attempt at plotting percent abundances of core eukaryotic ASVs at Phylum level

``` r
fung1920.pa # percent abundance ps object with all taxa
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 22035 taxa and 1355 samples ]
    ## sample_data() Sample Data:       [ 1355 samples by 28 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 22035 taxa by 7 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 22035 reference sequences ]

``` r
ITS.ASVs.asvc <- read.csv(file = "/Users/klas0061/Desktop/UMN/figures/from_server/ITS.ASVs.asvc") # info with state-counts
core.its.asvs <- ITS.ASVs.asvc %>% filter(num_states_per_asv_test == 9) %>% pull(Row.names) # finally the core its asvs
all.state.asvs.fung.1920.pa <- prune_taxa(core.its.asvs, fung1920.pa) # ps object with only the 9-state ASVs

gg.bar.core.its.asvs <- plot_bar(all.state.asvs.fung.1920.pa, fill="Phylum")+
   geom_bar(stat="identity", position="stack")+
   scale_fill_manual(values = natparks.pals("DeathValley", 10)[2:10])+
   facet_grid(~state, scales = "free", space = "free")+
   scale_y_continuous(limits = c(0,100))+
   ggtitle("Cosmopolitan Eukaryotic ASVs")+
   theme(axis.text.x = element_blank())+theme_bw()
gg.bar.core.its.asvs
```

![](05_serious_fungal_beginnings_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

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
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] UpSetR_1.4.0                NatParksPalettes_0.2.0     
    ##  [3] metagenomeSeq_1.40.0        RColorBrewer_1.1-3         
    ##  [5] glmnet_4.1-7                Matrix_1.5-3               
    ##  [7] limma_3.54.2                DESeq2_1.38.3              
    ##  [9] SummarizedExperiment_1.28.0 Biobase_2.58.0             
    ## [11] MatrixGenerics_1.10.0       matrixStats_0.63.0         
    ## [13] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
    ## [15] IRanges_2.32.0              S4Vectors_0.36.2           
    ## [17] BiocGenerics_0.44.0         vegan_2.6-4                
    ## [19] lattice_0.20-45             permute_0.9-7              
    ## [21] patchwork_1.1.2             phyloseq_1.42.0            
    ## [23] lubridate_1.9.2             forcats_1.0.0              
    ## [25] stringr_1.5.0               dplyr_1.1.1                
    ## [27] purrr_1.0.1                 readr_2.1.4                
    ## [29] tidyr_1.3.0                 tibble_3.2.1               
    ## [31] ggplot2_3.4.2               tidyverse_2.0.0            
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] colorspace_2.1-0       XVector_0.38.0         rstudioapi_0.14       
    ##  [4] farver_2.1.1           bit64_4.0.5            AnnotationDbi_1.60.2  
    ##  [7] fansi_1.0.4            codetools_0.2-19       splines_4.2.3         
    ## [10] cachem_1.0.7           geneplotter_1.76.0     knitr_1.42            
    ## [13] ade4_1.7-22            jsonlite_1.8.4         annotate_1.76.0       
    ## [16] cluster_2.1.4          png_0.1-8              compiler_4.2.3        
    ## [19] httr_1.4.5             fastmap_1.1.1          cli_3.6.1             
    ## [22] htmltools_0.5.5        tools_4.2.3            igraph_1.4.1          
    ## [25] gtable_0.3.3           glue_1.6.2             GenomeInfoDbData_1.2.9
    ## [28] reshape2_1.4.4         Rcpp_1.0.10            vctrs_0.6.1           
    ## [31] Biostrings_2.66.0      rhdf5filters_1.10.1    multtest_2.54.0       
    ## [34] ape_5.7-1              nlme_3.1-162           iterators_1.0.14      
    ## [37] xfun_0.38              timechange_0.2.0       lifecycle_1.0.3       
    ## [40] gtools_3.9.4           XML_3.99-0.14          zlibbioc_1.44.0       
    ## [43] MASS_7.3-58.2          scales_1.2.1           hms_1.1.3             
    ## [46] parallel_4.2.3         biomformat_1.26.0      rhdf5_2.42.0          
    ## [49] yaml_2.3.7             gridExtra_2.3          memoise_2.0.1         
    ## [52] stringi_1.7.12         RSQLite_2.3.1          highr_0.10            
    ## [55] foreach_1.5.2          caTools_1.18.2         BiocParallel_1.32.6   
    ## [58] shape_1.4.6            rlang_1.1.0            pkgconfig_2.0.3       
    ## [61] bitops_1.0-7           Wrench_1.16.0          evaluate_0.20         
    ## [64] Rhdf5lib_1.20.0        labeling_0.4.2         bit_4.0.5             
    ## [67] tidyselect_1.2.0       plyr_1.8.8             magrittr_2.0.3        
    ## [70] R6_2.5.1               gplots_3.1.3           generics_0.1.3        
    ## [73] DelayedArray_0.24.0    DBI_1.1.3              pillar_1.9.0          
    ## [76] withr_2.5.0            mgcv_1.8-42            survival_3.5-3        
    ## [79] KEGGREST_1.38.0        RCurl_1.98-1.12        crayon_1.5.2          
    ## [82] KernSmooth_2.23-20     utf8_1.2.3             tzdb_0.4.0            
    ## [85] rmarkdown_2.21         locfit_1.5-9.8         grid_4.2.3            
    ## [88] data.table_1.14.8      blob_1.2.4             digest_0.6.31         
    ## [91] xtable_1.8-4           munsell_0.5.0
