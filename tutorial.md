---
title: Tutorial
layout: home
---

## The Data

In class, Adina will introduce 16s rRNA amplicon sequencing.  There are several great review papers on the topic but to get you started here is [one](https://journals.asm.org/doi/10.1128/spectrum.00563-23).

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

The objective of this tutorial is to be able to identify a core with the criteria of ranked abundance and to evaluate its distribution of prevalence and abundance compared to 'non-core' taxa.

## Get the Data
You will need to download the data onto your laptop.  This action requires that you KNOW where you download files when you interact with your browser -- it varies by your own personal settings.  Mine, for example, is into a "Downloads" folder.  

Navigate to this [website](https://github.com/germs-lab/brc-data) and click on the Code button to download a zip file of the dataset.  Unzip it!  Wee!  This directory now contains your files and you can make it your *working directory* for your tutorial.  For example, for me it would be `~/Downloads/brc-data`.

## Prepare Personal Touches

First, set the working directory with *setwd* in R to your tutorial working directory.  You will also want to load your R libraries now.

```
setwd("~/Box Sync/sro-analysis-core/")
library(vegan)
library(phyloseq)
```

## Load Data
I have loaded all the data into a phyloseq object to save time.  The phyloseq object has 4 parts, and you can see that by executing this code.

```
ps <- readRDS('phyloseq-object.rds')
ps
otu_table(ps)[1:5, 1:5]  #Note that this code prints the first 5 rows and 5 columns of the matrix
tax_table(ps)[1:5]
sample_data(ps)[1:15]
```

### Task 1. 
With a partner, discuss the different contents of the phyloseq object.

## Explore Data and subsetting the Core
You can sort and subset by different things, like the total observations of an ASV.  Some of the functions you are using here in phyloseq include `taxa_sums` and `names`.  You are also using R functions like `sort`.  You can google many of these functions but here is the [manual for phyloseq](https://bioconductor.org/packages/release/bioc/manuals/phyloseq/man/phyloseq.pdf)

```
taxa_sums(ps)
sort(taxa_sums(ps), decreasing = TRUE)
names(sort(taxa_sums(ps), decreasing = TRUE)[1:10])
core_taxa_names <- names(sort(taxa_sums(ps), decreasing = TRUE)[1:10])
```

```
ps_core = prune_taxa(core_taxa_names, ps)
taxa_sums(ps_core) #sum of all observations of each taxa
sample_sums(ps_core) #sum of all taxa observe in each sample
```

### Task 2.
With a new partner, discuss the what the above code accomplished.  What did `sample_sums` do?

## Making another "core" phyloseq object
You now have a phyloseq object called `ps_core`, which is a type of core defined by the 10 most abundant taxa.  You can also define cores other ways.

#You can filter taxa by abundance - where the taxa have an average abundance among samples over 100 observations

```
ps_core_alt <- filter_taxa(ps, function(x) mean(x) > 100, TRUE)
```

## Let's put things in a table format that is readable for you - this will merge the abundance and taxa information

```
#You can add prevelance and abundance data to a dataframe
prevalance = apply(X = otu_table(ps_core),
                     MARGIN = 2,
                     FUN = function(x){sum(x > 0)})
tot_samples <- dim(sample_data(ps_core))[1] #this is the total number of samples
core_data = data.frame(Prevalence = prevalance, Prop_Prev = prevalance/tot_samples,
                          TotalAbundance = taxa_sums(ps_core),
                          tax_table(ps_core))

#You can explore distributions
library(ggplot2)
ggplot(data = core_data, aes(x = Prevalence)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Prevalence",
       x = "Prevalence",
       y = "Frequency")
```

### Task 3.
With a partner, make a histogram of total abundance -or- proportional prevalence.  Can you make a plot of abundance vs prevelance?  There is a hint below.

```
ggplot(data = dat, aes(x = ___, y = ___)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Scatter Plot of Prevalence vs Total Abundance",
       x = "Prevalence",
       y = "Total Abundance")
```

## Other Metrics to Identify a core
You can get a list of ASVs based on prevalence or abundance.

```
threshold_prev = 0.2
threshold_abund = 2000
prev_cutoff <- core_data[core_data$Prop_Prev > threshold_prev,] # Only the prevalance
prev_cutoff
abund_cutoff <- prev_cutoff[prev_cutoff$TotalAbundance > threshold_abund,]
taxa_to_filter <- rownames(abund_cutoff)
length(myTaxa)
mycore = prune_taxa(taxa_to_filter, ps)
```

### Task 4
Choose a way to define a core.  Can you create a phyloseq object with taxa that meet your prevalance and abundance thresholds?  You can either subset from teh core phyloseq object or from the larger phyloseq object.  [Hint:  `prune_taxa`]

## Adding a label of core or not core to your table

First, let's create two lists of taxa.  One will represent "core taxa" and the other "other taxa".

```
other_taxa = names(sort(taxa_sums(ps), decreasing = TRUE)[1:100])
other_ps = prune_taxa(other_taxa, ps)
prev = apply(X = otu_table(other_ps),
                     MARGIN = 2,
                     FUN = function(x){sum(x > 0)})
dim(sample_data(other_ps))[1]
data_table = data.frame(Prevalence = prev, Prop_Prev = prev/dim(sample_data(other_ps))[1],
                           TotalAbundance = taxa_sums(other_ps),
                           tax_table(other_ps))
threshold_prev = 0.2
threshold_abund = 2000
prev_cutoff <- core_data[core_data$Prop_Prev > threshold_prev,] # Only the prevalance
prev_cutoff
abund_cutoff <- prev_cutoff[prev_cutoff$TotalAbundance > threshold_abund,]
core_taxa <- rownames(abund_cutoff)
data_table$core_id <- ifelse(rownames(data_table) %in% core_taxa, "core", "other")
```

### Task 5 

Make a prevalence abundance plot of the core and non-core.

Hint:

```
ggplot(data = data_table, aes(x = ___, y = ___, color=core_id)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Scatter Plot of Prevalence vs Total Abundance",
       x = "Prevalence",
       y = "Total Abundance")
```       

## Phyloseq is now your playground

Now that you can create a phyloseq object of your core, you can do a lot of things with microbiome packages.  Phyloseq is a great place to start.

### Task 6

Make a bar chart of the of phylogenies in the core.  You can see this [tutorial](https://joey711.github.io/phyloseq/plot_bar-examples.html)


## Other Phyloseq Tutorials

[https://joey711.github.io/phyloseq/index.html](https://joey711.github.io/phyloseq/index.html)

[https://micca.readthedocs.io/en/latest/phyloseq.html](https://micca.readthedocs.io/en/latest/phyloseq.html)

[https://vaulot.github.io/tutorials/Phyloseq_tutorial.html](https://vaulot.github.io/tutorials/Phyloseq_tutorial.html)


