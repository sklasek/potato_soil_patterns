# dada2 processing of 16S amplicon sequence data for potato project
# scott klasek
# 5-17-2022

# load libraries
library(tidyverse)
library(dada2)
library(phyloseq)

# 1. Specify file paths and inspect read quality
fqpath <- getwd() # it's NOT set to take the fastq files from the previous directory, so this can be run in a directory called "sequence_processing" and output can be stored there
fnFs <- sort(list.files(fqpath, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(fqpath, pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)
filtFs <- file.path(fqpath, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(fqpath, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
# skipping for troubleshooting purposes: write out plots of quality score distributions:

# 2. Process reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(17,20), truncLen=c(280,205), maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
errF <- learnErrors(filtFs, multithread=TRUE, randomize = TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, randomize = TRUE)
pdf("error.plot.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE) 
seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))
seqtab.trim <- seqtab[,nchar(colnames(seqtab)) %in% 402:432]
table(nchar(getSequences(seqtab.trim)))

seqtab.nochim <- removeBimeraDenovo(seqtab.trim, method="consensus", multithread=TRUE, verbose=TRUE) 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
tdf <- as.data.frame(track)
rm(derepFs) # remove the large objects we don't need anymore
rm(derepRs)
rm(dadaFs)
rm(dadaRs)

# function to plot read track data
plot.track.reads <- function(tdf){ # tdf is a data frame of read counts with columns named input, filtered, denoisedF/R, merged, and nonchim
  tdf$name <- rownames(tdf)
  tdf$filtered_out <- tdf$input-tdf$filtered # reads that were filtered out
  tdf$noise_removed <- tdf$filtered-with(tdf, pmin(denoisedF, denoisedR)) # reads whose F or R was too noisy
  tdf$unmerged <- (tdf$filtered-tdf$noise_removed)-tdf$merged # unmerged reads 
  tdf$chimeras <- tdf$merged-tdf$nonchim # chimeras 
  tdfs <- tdf %>% select("name","filtered_out","noise_removed","unmerged","chimeras", "nonchim") # select only columns of importance
  tdfl <- pivot_longer(tdfs, !name, names_to = "step", values_to = "reads") # convert to long df
  tdfl$step <- factor(tdfl$step, levels = c("filtered_out","noise_removed","unmerged","chimeras", "nonchim"))
  gg <- ggplot(tdfl,aes(name,reads,fill=step))+ # plot
    geom_bar(stat="identity")+
    theme_bw()+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  return(gg)
}
pdf("track.reads.pdf")
plot.track.reads(tdf)
dev.off() 

# function to count read track data 
table.track.reads <- function(tdf){ # get counts of reads that passed each step and write them out in case they are needed in addition to the plot
  tdf$name <- rownames(tdf)
  tdf$filtered_out <- tdf$input-tdf$filtered # reads that were filtered out
  tdf$noise_removed <- tdf$filtered-with(tdf, pmin(denoisedF, denoisedR)) # reads whose F or R was too noisy
  tdf$unmerged <- (tdf$filtered-tdf$noise_removed)-tdf$merged # unmerged reads 
  tdf$chimeras <- tdf$merged-tdf$nonchim # chimeras 
  tdfs <- tdf %>% select("filtered_out","noise_removed","unmerged","chimeras", "nonchim") # select only columns of importance
  return(tdfs)
}
track.read.counts <- table.track.reads(tdf)
write_tsv(track.read.counts, file = "./track_reads_table.txt") # write the counts as a tsv

# 3. Assign taxonomy to reads
taxa <- assignTaxonomy(seqtab.nochim, "/home/kinkelll/klas0061/dbs/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/home/kinkelll/klas0061/dbs/silva_species_assignment_v138.1.fa.gz") # species is assigned only for an EXACT match 
taxa.print <- taxa
rownames(taxa.print) <- NULL

# 4. Get basic metadata from file names
metadata <- data.frame("state"=substring(sample.names, 1,2), # the first two characters in the name are the state
                       "objective"=as.numeric(substring(sample.names, 4,4)), # the fourth is the objective (1 or 2)
                       "rotation"=as.numeric(substring(sample.names, 6,6)), # the sixth is the rotation (2 or 3)
                       "plot"=as.numeric(substring(sample.names, 8,10)), # the 8th through 10th are the plot number, which contains treatment and replicate info
                       "year"=as.numeric(substring(sample.names, 12,13)), # the 12th and 13th are the year (last two digits)
                       "month"=gsub("_.*$","",sub("..............", "", sample.names)), # the month is the substring after the first 14 characters of the name, but before the first underscore thereafter
                       "code"=sub(".*_", "", sample.names)) # the code is the last substring after the final underscore
rownames(metadata) <- sample.names

# 5. Make phyloseq object
ps <- phyloseq(tax_table(taxa),
               sample_data(metadata),
               otu_table(seqtab.nochim, taxa_are_rows = FALSE))

# 6. Report and remove unwanted taxa
# report by domain, also include chloroplasts and mitochondria
levels(factor(tax_table(ps)[,1], exclude = NULL)) # domain-level taxonomic annotation
percentNAs <- sum(is.na(tax_table(ps)[,1]))/ntaxa(ps)*100 
taxa_breakdown <- c("Bacteria","Archaea","Eukaryotes","Unclassified","Chloroplasts","Mitochondria")
percent <- c(sum(tax_table(ps)[,1]=="Bacteria", na.rm=TRUE)/ntaxa(ps)*100,
             sum(tax_table(ps)[,1]=="Archaea", na.rm=TRUE)/ntaxa(ps)*100,
             sum(tax_table(ps)[,1]=="Eukaryota", na.rm=TRUE)/ntaxa(ps)*100,
             percentNAs,
             sum(tax_table(ps)[,4]=="Chloroplast", na.rm=TRUE)/ntaxa(ps)*100,
             sum(tax_table(ps)[,5]=="Mitochondria", na.rm=TRUE)/ntaxa(ps)*100)
taxa_annotation_breakdown <- data.frame(taxa_breakdown, percent)
paste("All domains and unclassified sequences add up to", sum(percent[1:3], percentNAs) ,"percent") # check that domains and NAs add up to 100
write_tsv(taxa_annotation_breakdown, file = "./taxa_annotation_breakdown.txt")

# remove unwanted taxonomies
paste("The phyloseq object contains", ntaxa(ps), "ASVs before removing Eukaryotes, Chloroplasts, and Mitochondria.")
ps.trim1 <- subset_taxa(ps, (Kingdom!="Eukaryota")) 
paste((ntaxa(ps) - ntaxa(ps.trim1)), "Eukaryotic and Unclassified sequences were removed.")
ps.trim2 <- subset_taxa(ps.trim1, (Order!="Chloroplast") | is.na(Order)) 
paste(ntaxa(ps.trim1) - ntaxa(ps.trim2), "Chloroplast sequences were removed.")
ps.trim3 <- subset_taxa(ps.trim2, (Family!="Mitochondria") | is.na(Family))
paste(ntaxa(ps.trim2) - ntaxa(ps.trim3), "Mitochondrial sequences were removed.")
paste("We are left with", ntaxa(ps.trim3), "ASVs.")

# 7. Display library sizes across samples
ps.trim3.df <- as.data.frame(sample_data(ps.trim3))
ps.trim3.df$LibrarySize <- sample_sums(ps.trim3) 
ps.trim3.df <- ps.trim3.df[order(ps.trim3.df$LibrarySize),]
libsizes <- ps.trim3.df[,"LibrarySize"]
write_tsv(libsizes, file = "./library_sizes.txt")

# 8. Write out phyloseq object
saveRDS(ps.trim3, "./processed.ps")

# 9. Report session information
sessionInfo()

