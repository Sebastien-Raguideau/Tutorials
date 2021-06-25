# 16S processing tutoriel
16S sequencing is done by selectively targeting some specific conserved regions of the 16s ribosomal RNA. This is done using primer designed to correspond exactly to the conserved regions so that all reads start exactly from the conserved region and span the variable region from which we can infer taxonomy. Primers/Amplicons are designed/selected so that Forward and Reverse reads starting from each side of the amplicon may meet in the midle and overlap. 

![alt tag](16s-gene.png)

### Plan 
1) Using dada2 
3) Phyloseq for diversity and graphs

## 1) DADA2
Dada2 is a denoising algorithm. It aims at distinguishing differences between variants and read sequencing substitution errors. It output ASV for amplicon sequence variants. ASV are higher resolutions than OTUs and aim at a variant level resolutions.

To start, let's create an output folder where we will do this analysis, AD_16S. 

<details><summary>Reveal commands</summary>
<p>

```
cd ~/Projects
mkdir AD_16S
cd AD_16S
```

</p>
</details>

Launch R and load library needed:

    library(dada2)

This library contains all the functions from dada2.

Let's define the path to our dataset as well as the path to our output folder

    path_data = "/home/ubuntu/Data/AD/16S/data"
    path_out = "/home/ubuntu/Projects/AD_16S"

### Read quality
Before starting any sort of analysis we need to be sure that the reads are of good quality. 
Dada2 allows us to do so easily 

    R1_files=sort(Sys.glob(paste(path_data,"*R1*",sep="/"), dirmark = FALSE))
    R2_files=gsub("R1","R2",R1_files)

**plotQualityProfile** to ... plot fastq files quality as a function of position alon the read

    pdf("quality.pdf")
    plotQualityProfile(R1_files)
    plotQualityProfile(R2_files)
    dev.off()

Use evince to look at the output. 

![alt tag](/figs/R1_qual_init.png)
![alt tag](/figs/R2_qual_init.png)

How is the quality of the reads? 
There is a pattern in quality, describe what can be seen.

We are going to filter quality using the  `filterAndTrim`function, first we need to create file path corresponding to clean/trimmed/filtered reads : 

    temp_filter_path = paste(path_out,"/temp",sep ="")
    dir.create(temp_filter_path)
    
    sample.name <- gsub("_R1.fastq","",basename(R1_files))
    Filtered_R1 <- file.path(temp_filter_path, paste0(sample.name, "_R1_Filtered.fastq"))
    Filtered_R2 <- file.path(temp_filter_path, paste0(sample.name, "_R2_Filtered.fastq"))

Bioinformatic practical : takes some time to read the documentation and try to devise by yourself how the function   filterAndTrim should be used. 

<details><summary> Answer</summary>
<p>
Documentation can be found there
https://letmegooglethat.com/?q=dada2+filterAndTrim&l=1
</p>
</details>
<details><summary> Answer part2</summary>
<p>
We definetely want to use the options :  trunclen, maxN, maxEE, truncQ, rm.phix, compress, verbose, multithread 
</p>
</details>

<details><summary> Answer part3</summary>
<p>

    Log_filtering =  filterAndTrim(R1_files,Filtered_R1,R2_files,Filtered_R2,truncLen=c(240,160),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, verbose=TRUE,multithread=TRUE)

</p>
</details>

Let's look at the quality of currents these filtered reads.

    pdf("quality_Filtered.pdf")
    plotQualityProfile(Filtered_R1)
    plotQualityProfile(Filtered_R2)
    dev.off()


### Learning the error rates
This is the most time consumming part of the pipeline. 


    Error_R1 <- learnErrors(Filtered_R1, multithread=TRUE)
    Error_R2 <- learnErrors(Filtered_R2, multithread=TRUE)

<details><summary> If it takes too long</summary>
<p>

    Error_R1 = readRDS(file  =  "/home/ubuntu/Prerun/AD16S/errors_R1.rds")
    Error_R2 = readRDS(file  =  "/home/ubuntu/Prerun/AD16S/errors_R2.rds")

</p>
</details>

Let's have a look at the error model that dada2 will be using : 

    pdf("Errors_learned.pdf")
    plotErrors(Error_R1, nominalQ=TRUE)
    plotErrors(Error_R2, nominalQ=TRUE)
    dev.off()

![alt tag](/figs/Error_learned.png)

Can you intuite what the red and black lines correspond to? 

### Dereplication

    Derep_R1 <- derepFastq(Filtered_R1, verbose=TRUE)
    Derep_R2 <- derepFastq(Filtered_R2, verbose=TRUE)

### Error correction

    dada_R1 <- dada(Derep_R1, err=Error_R1, multithread=TRUE,pool=TRUE)
    dada_R2 <- dada(Derep_R2, err=Error_R2, multithread=TRUE,pool=TRUE)

### Merge reads 

    mergers <- mergePairs(dada_R1, Derep_R1, dada_R2, Derep_R2, verbose=TRUE)


### Chimera

    Seqtab <- makeSequenceTable(mergers)
    Seqtab.nochim <- removeBimeraDenovo(Seqtab, method="consensus", verbose=TRUE)

### Summary

    getN <- function(x) sum(getUniques(x))
    track <- cbind(Log_filtering, sapply(dada_R1, getN), sapply(dada_R2, getN), sapply(mergers, getN), rowSums(Seqtab.nochim))
    colnames(track) <- c("input", "filtered", "denoised_R1", "denoised_R2", "merged", "nonchim")
    rownames(track) <- sample.name
    track

Why does the numbers or reads goes down as a result of merging? 
The number of chimera is in the unit why is does the counts diminution is bigger than that? 
How many different sequences do we end up with?

### Taxonomic annotation

    taxa <- assignTaxonomy(Seqtab.nochim, "~/Databases/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
    taxa <- addSpecies(taxa, "~/Databases/silva_species_assignment_v132.fa.gz")

#### Intermediary Results

    write.csv(Seqtab.nochim,paste(path_out,'sequence_table.csv',sep="/"))
    write.csv(taxa,paste(path_out,'taxonomy_table.csv',sep="/"))

# 2) Phyloseq
Phyloseq is a R library for handling taxonomy data, it contain multiple handy function for ploting diversity, taxonomic profile/diversity .... More documentation at https://joey711.github.io/phyloseq/
### Creating a phyloseq object
The first step needed is the most complicated one, creating the phyloseq object. 
For most application the phyloseq object need at least 
- an "otu_table" : a table of counts of otus/variants ... 
- a taxonomic table : corresponding taxonomic annotation for each otus/variant
- a metadata table : it allows to use the full power of the phyloseq and integrating/representing easily your data. 

First load the library and the metadata

     library(phyloseq)
     metadata <- read.table("/home/ubuntu/Data/AD/16S/metadata.tsv", sep="\t", header=TRUE, row.names=1)
     str(metadata)

Then let's create a phyoseq object

     ps=phyloseq(otu_table(Seqtab.nochim, taxa_are_rows=FALSE),tax_table(taxa),sample_data(metadata))


Bioinformatic pratical :  try to interpret the error message and solve the issue
<details><summary> Answer</summary>
<p>
Try first to find some information on the internet
<details><summary> Answer part2</summary>
<p>
The issue is related to samples names as defined inside Seqtab.nochim, check what the sample names are there and compare these to the sample names in metadata, you need to make it so they are identical.
<details><summary> Answer part3</summary>
<p>

    asv_table = otu_table(Seqtab.nochim, taxa_are_rows=FALSE)
    sample_names(asv_table) = gsub("_R1_Filtered.fastq","",sample_names(asv_table)) 
    ps=phyloseq(asv_table,tax_table(taxa),sample_data(metadata))

</p>
</details>
</p>
</details>
</p>
</details>

### Taxonomy profile
Let's look at phylum level taxonomy profile using the function **tax_glom**

    ps_phylum=tax_glom(ps, "Phylum")
    p1 = plot_bar(ps_phylum, fill="Phylum")

![alt tag](./figs/taxa.png)

This gives use all phylum, let's get only the 10 most abundant, using the function **prune_taxa**

    top10= names(sort(colSums(otu_table(ps_phylum)),decreasing=TRUE))[1:10]
    ps_phylumN =  transform_sample_counts(ps_phylum, function(variant) variant/sum(variant))
    ps_phylumN.top10=prune_taxa(top10, ps_phylumN)
    p2 = plot_bar(ps_phylumN.top10, x="week", fill="Phylum")


Save the plot : 

    pdf(paste(path_out,"phylum_profile.pdf",sep="/"))
    print(p1)
    print(p2)
    dev.off() 

### Richness plot
Using the function **plot_richness** allows to .... plot diverse richness measures

    pdf(paste(path_out,"Richness.pdf",sep="/"))
    plot_richness(ps,x="week",measures=c('Chao1','Simpson'),color="ch4...")
    dev.off()
    
![alt tag](./figs/diversity.png)

### NMDS plot
A Nmds plot is an ordination plot : a method to represent a high dimensional object in a 2 dimensional plane. We have 10 samples with ~ 500 coordinates and we want to represente that with only about 2 (X,Y). The 

    ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
    ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
    p1 = plot_ordination(ps.prop, ord.nmds.bray, title="Bray NMDS",label="week",color="ch4...")

What does rmse/resid means? 
Looking at this plot there seems to be an outlier. Let's remove it and redo this.

    ps_clean <- prune_samples((rownames(sample_data(ps)) !="AD7_W40"), ps)
    ps.prop <- transform_sample_counts(ps_clean, function(otu) otu/sum(otu))
    ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
    p2 = plot_ordination(ps.prop, ord.nmds.bray, title="Bray NMDS",label="week",color="ch4...")

![alt tag](./figs/NMDS.png)

Save the plots 

    pdf(paste(path_out,"NMDS.pdf",sep="/"))
    print(p1)
    print(p2)
    dev.off()
