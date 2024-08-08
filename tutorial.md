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
