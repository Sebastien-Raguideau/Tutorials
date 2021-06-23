# 16S processing tutoriel
### Intro : 
- quick reminder of what the data : 16S region selection, paired reads, fastq files 
- Dataset is a time serie of 10 samples from AD, additional infos on that particular industrial AD, feed/....

### Plan 
1) Using dada2 
2) asv vs otu 
3) Phyloseq for diversity and graphs
4) Vegan with adonis/nmds

## DADA2
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

out =  filterAndTrim(R1_files,Filtered_R1,R2_files,Filtered_R2,truncLen=c(240,160),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, verbose=TRUE,multithread=TRUE)

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

readRDS(file  =  "/home/ubuntu/PreRun/AD16S/errors_R1.rds")
readRDS(file  =  "/home/ubuntu/PreRun/AD16S/errors_R2.rds")

</p>
</details>
