library(dada2)
#setwd("~/Documents/Bioinformatics_scripts/R_scripts/endoliths")
### ### ### ### ### 
### Read files  ### 
### ### ### ### ###  
cat("Reading files")
path <- "."
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

cat("Processing",length(sample.names),"samples:", sample.names)


#Inspect read quality profiles
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])

### ### ### ### ### ### 
### Filter and trim ### 
### ### ### ### ### ### 
cat("Filtering and trimming")
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
### set parameters  maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(290,210),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE) 
head(out)

### ### ### ### ### ### 
### Learn error rates ### 
### ### ### ### ### ### 

cat("Learning error rates")
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

### ### ### ### ### ### 
### Sample inference ### 
### ### ### ### ### ### 

#apply the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#Inspecting the returned dada-class object:
dadaFs[[1]]


### ### ### ### ### ### ### 
### Merge paired reads ### 
### ### ### ### ### ### ### 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


### ### ### ### ### ### ### 
### Contructuct ASV table ### 
### ### ### ### ### ### ### 

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


### ### ### ### ### ### ### 
### Remove chimeras ### 
### ### ### ### ### ### ### 

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

## calculate the  frequency of chimeric sequences 
sum(seqtab.nochim)/sum(seqtab)

### ### ### ### ### ### ### ## ### ###
### Track reads through the pipeline ### 
### ### ### ### ### ### ### ## ### ### 

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

### ### ### ### ### ###
### Assign taxonomy ### 
### ### ### ### ### ### 

taxa <- assignTaxonomy(seqtab.nochim, "/share/databases/SILVA_for_dada2/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
#taxa <- assignTaxonomy(seqtab.nochim, "~/Downloads/silva_nr_v132_train_set.fa", multithread=TRUE)
##optional: make species level assignments based on exact matching
taxa <- addSpecies(taxa, "/share/databases/SILVA_for_dada2/silva_species_assignment_v132.fa.gz")
#taxa <- addSpecies(taxa, "~/Downloads/silva_species_assignment_v132.fa")
#inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


### ### ### ### ### ###
### Export table  #### 
### ### ### ### ### ### 

samples.out <- rownames(seqtab.nochim)

asv.1=t(seqtab.nochim)
asv.1=cbind(asv.1, "sum"=rowSums(asv.1)) 
asv.2=merge(asv.1, taxa, by="row.names")
colnames(asv.2)[1]="Sequence"
asv.3=asv.2[,c(2:ncol(asv.2),1)]
asv.final=asv.3[order(-asv.3$sum),]
rownames(asv.final) = sprintf("ASV%04d", 1:nrow(asv.final))
write.table(asv.final, "coral_endoliths_ASV_table.txt",  quote = FALSE)
write.table(track, "coral_endoliths_ASV_stats.txt", quote = FALSE)
