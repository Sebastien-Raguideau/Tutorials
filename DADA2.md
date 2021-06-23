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

