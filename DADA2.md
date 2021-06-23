# 16S processing tutoriel
### Intro : 
- quick reminder of what the data : 16S region selection, paired reads, fastq files 
- Dataset is a time serie of 10 samples from AD, additional infos on that particular industrial AD, feed/....

### Plan 
1) Using dada2 
2) asv vs otu 
3) Phyloseq for diversity and graphs
4) Vegan with adonis/nmds

## 1) DADA2
In the terminal go to Projects and create the folder we will be using for this analysis, AD_16S 

<details><summary>Reveal commands</summary>
<p>

```
cd ~/Projects
mkdir AD_16S
```

</p>
</details>

Launch R and load library needed:

    library(dada2)

Let's define the path to our dataset as well as the path to our output folder

    path_data = "/home/ubuntu/Data/AD/16S/data"
    path_out = "/home/ubuntu/Projects/AD_16S"

### Read quality
Before starting any sort of analysis we need to be sure that the reads are of good quality. 
Dada2 allows us to do so easily 

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

    temp_filter_path = paste(out,"/temp",sep ="")
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
**Explain error rate estimation.**

    Error_R1 <- learnErrors(Filtered_R1, multithread=TRUE)
    Error_R2 <- learnErrors(Filtered_R2, multithread=TRUE)

<details><summary> If it takes too long</summary>
<p>

    readRDS(file  =  "/home/ubuntu/Prerun/AD16S/errors_R1.rds")
    readRDS(file  =  "/home/ubuntu/Prerun/AD16S/errors_R2.rds")

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
    rownames(track) <- sample.names
    track

### Taxonomic annotation

    taxa <- assignTaxonomy(Seqtab.nochim, "~/seb/Database/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
    taxa <- addSpecies(taxa, "~/seb/Database/silva_species_assignment_v132.fa.gz")

#### Results

    write.csv(Seqtab.nochim,paste(out,'sequence_table.csv',sep="/"))
    write.csv(taxa,paste(out,'taxonomy_table.csv',sep="/"))

