---
title: Tutorial
layout: home
---

## The Data

Prior to the workshop, we obtained a small, representative dataset from each BRC.  As a result, we have SSU sequencing datasets from corn, miscanthus, poplar, switchgrass, and sorghum.  These sequences are from different regions of the 16S rRNA gene, depending on the project source.  To check for this, I aligned these sequences to an E. coli 16S SSU gene reference.  I found that CABBI sequences represented region 533-668 (R1), and 783-644 (R2); CBI was 515-751 (R1) and 804-557 (R2); GLBRC was 514-798 (R1) and 804-517 (R2); and JBEI was 514-750 (R1) and 804-559 (R2).  These are all overlapping, which was helpful but different lengths.  I trimmed everything on the left side (but CABBI) by 19 bp, which would make the R1 and R2 similar starting points.  These starting points are “on average”, so we’re also “doing our best here”.  The goal is to get the sequences as close to simlar as we can.  For each dataset, I think trimmed, quality controlled, and merged reads.  I took only merged reads from ALL BRCs and then clustered them at 100% identity with CD-HIT using this [script](https://github.com/dsamoht/utility/tree/main?tab=readme-ov-file#collapse_asvpy).  This resulted in a merged file that we will use for this tutorial.

For this tutorial, we are assuming that your data has been processed through dada2 and then imported into phyloseq as a phyloseq object.  We are going to go through the pieces of phyloseq objects in the tutorial but you can read more about these details [here](https://joey711.github.io/phyloseq/import-data.html).

## “Core” Objectives

To explore a core, we’ve identified a few key “operations” that are useful.
* You need to be able to extract specific taxa from a phyloseq obeject
* You need to be able to look at the abundance distribution within a microbiome
* You need to be able to look at the prevalence distribution within a microbiome
* You need to be able to define a core and subset the taxa wihtin that core.
* You need to be able to compare that core to non-core data

## Breakout Group Exercises (To be brainstored as a group)

* How does the number of core taxa change based on different thresholds?  i.e., a plot of total number of core taxa based on increasing prevelance or abundance
* How much does the core change if defined by specific feedstocks?
* How much does the core change if looking at root vs not root?
* At what point do you not need any more samples to see a difference in core?

# Tutorial Coding Portion

## Get the Data
You will need to download the data onto your laptop.  This action requires that you KNOW where you download files when you interact with your browser -- it varies by your own personal settings.  Mine, for example, is into a "Downloads" folder.  

Navigate to this [website](https://github.com/germs-lab/brc-data) and click on the Code button to download a zip file of the dataset.  Unzip it!  Wee!  This directory now contains your files and you can make it your *working directory* for your tutorial.  For example, for me it would be `~/Downlaods/brc-data`.

## Prepare Personal Touches

First, set the working directory with *setwd* in R to your tutorial working directory.  You will also want to load your R libraries now.

```
setwd("~/Box Sync/sro-analysis-core/")
library(vegan)
library(phyloseq)
```

## Load Data

```
otu_csv <- read.csv("abundance-table-final.csv", header=TRUE)
rownames(otu_csv) <- otu_csv[,1]
otu_csv[,1] <- NULL
otu <- otu_table(as.matrix(otu_csv), taxa_are_rows = FALSE)
taxa_csv <- read.csv("tax-final.csv", header=TRUE)
rownames(taxa_csv) <- taxa_csv[,1]
taxa_csv[,1] <- NULL
tax <- tax_table(as.matrix(taxa_csv))
meta <- read.csv("meta.txt", header=TRUE, sep="\t")
rownames(meta) <- meta$X
meta$X <- NULL
meta_phy <- sample_data(meta)
phy <- phyloseq(tax, otu, meta_phy)
```


## Rename ASVs to be not crazy long 

```
dna <- Biostrings::DNAStringSet(taxa_names(phy))
names(dna) <- taxa_names(phy)
ps <- merge_phyloseq(phy, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
```

## Explore Data 

```
#Understanding data objects in phyloseq
ps
otu_table(ps)[1:5, 1:5]
tax_table(ps)[1:5]
sample_data(ps)[1:5]

#You can sort and subset by different things, like the total sum in a sample
myTaxa = names(sort(taxa_sums(ps), decreasing = TRUE)[1:10])
ex1 = prune_taxa(myTaxa, ps)
taxa_sums(ex1) #sum of all observations of each taxa
sample_sums(ex1) #sum of all taxa observe in each sample

#You can transform from total observations to relative abundances
pseq.rel <- microbiome::transform(ex1, "compositional")
taxa_sums(pseq.rel)

#You can make rarefaction curves
ps.rarefied <- rarefy_even_depth(ex1, rngseed=1, sample.size=15, replace=F)
otu.rarefied <- as.data.frame(otu_table(ps.rarefied))
otu.rarefied
otu.rarecurve = rarecurve(otu.rarefied)

#You can filter taxa by abundance
filtered_ps <- filter_taxa(ex1, function(x) mean(x) > 100, TRUE)

#You can add prevelance and abundance data to a dataframe
prevelance = apply(X = otu_table(ex1),
                     MARGIN = 2,
                     FUN = function(x){sum(x > 0)})
tot_samples <- dim(sample_data(ex1))[1]
dat = data.frame(Prevalence = prevelance, Prop_Prev = prevelance/tot_samples,
                          TotalAbundance = taxa_sums(ex1),
                          tax_table(ex1))

#You can explore distributions
library(ggplot2)
ggplot(data = dat, aes(x = Prevalence)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Prevalence",
       x = "Prevalence",
       y = "Frequency")

ggplot(data = dat, aes(x = Prop_Prev)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Prevalence",
       x = "Prevalence",
       y = "Frequency")

ggplot(data = dat, aes(x = TotalAbundance)) +
  geom_histogram(binwidth = 1000, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Prevalence",
       x = "Abundance",
       y = "Frequency")

ggplot(data = dat, aes(x = TotalAbundance, y = Prop_Prev)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Scatter Plot of Prevalence vs Total Abundance",
       x = "Prevalence",
       y = "Total Abundance")
```

## Core vs not core

```
myTaxa2 = names(sort(taxa_sums(ps), decreasing = TRUE)[1:100])
subset_ps = prune_taxa(myTaxa2, ps)
prev = apply(X = otu_table(subset_ps),
                     MARGIN = 2,
                     FUN = function(x){sum(x > 0)})
dim(sample_data(ex2))[1]
dat2 = data.frame(Prevalence = prev, Prop_Prev = prev/dim(sample_data(ex2))[1],
                           TotalAbundance = taxa_sums(ex2),
                           tax_table(ex2))

myTaxa #This is the select list to label as core
myTaxa2 #This is the NOT core
dat2$core <- ifelse(rownames(dat2) %in% myTaxa, "core", "other")
ggplot(data = dat2, aes(x = TotalAbundance, y = Prevalence, color=core)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Scatter Plot of Prevalence vs Total Abundance",
       x = "Prevalence",
       y = "Total Abundance")
```

## Getting the core as a phyloseq object
```
threshold_prev = 0.2
threshold_abund = 2000
prev_cutoff <- dat2[dat2$Prop_Prev > threshold_prev,] # Only the prevalance
prev_cutoff$Prop_Prev
abund_cutoff <- prev_cutoff[prev_cutoff$TotalAbundance > threshold_abund,]
myTaxa <- rownames(abund_cutoff)
length(myTaxa)
mycore = prune_taxa(myTaxa, ps)
```

## Subset by Metadata
```
ps_poplar <- subset_samples(ps, crop == "poplar")
ps.rarefied <- rarefy_even_depth(ex1, rngseed=1, sample.size=5, replace=F)
otu.rarefied <- as.data.frame(otu_table(ps.rarefied))
otu.rarecurve = rarecurve(otu.rarefied)
```

## Other Phyloseq Tutorials

[https://joey711.github.io/phyloseq/index.html](https://joey711.github.io/phyloseq/index.html)

[https://micca.readthedocs.io/en/latest/phyloseq.html](https://micca.readthedocs.io/en/latest/phyloseq.html)

[https://vaulot.github.io/tutorials/Phyloseq_tutorial.html](https://vaulot.github.io/tutorials/Phyloseq_tutorial.html)


