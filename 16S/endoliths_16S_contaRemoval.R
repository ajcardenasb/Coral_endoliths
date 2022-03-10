###################################################################
#####Identifying and removing contaminant OTUs normalized data#####
###################################################################

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/16S/")

#read in OTU table generated in mothur
#otu = read.table("./Input_files/coral_endoliths_ASV_table.txt", header = T)
otu = read.table("./Input_files/coral_endoliths_pooling_ASV_table.txt", header = T)
names(otu)
tax=otu[,c(75:81)]
#Remove samples with low read number
otu.n=apply(otu[,c(1:73)], 2, as.numeric)
otu.o=otu.n[, colSums(otu.n) > 1000]
rownames(otu.o)=rownames(otu)
message(ncol(otu.o)," samples with > 1000 reads were retained out of ", ncol(otu.n), " total samples")

#Identify and removing contaminant ASVs raw data
otu.r=as.data.frame(sweep(otu.o,2,colSums(otu.o),`/`))
otu.r=transform(otu.r,  Sum = rowSums(otu.r[,2:ncol(otu.r)]))
names(otu.r)
otu.r=transform(otu.r,  SumNegs = rowSums(otu.r[,c(31,32)])) # define negative controls here
names(otu.r)
otu.r=transform(otu.r,  contaFactor=(otu.r$SumNegs/otu.r$Sum)*100)
rownames(otu.r)=rownames(otu)
Conta=subset(otu.r, otu.r$contaFactor > 10)
Conta$Family=otu$Family[match(rownames(Conta), rownames(otu))]
message("Number of total ASVs: ", nrow(otu))
message("Number of identified contaminant ASVs removed from the analysis: ", length(rownames(Conta)), "\n", Conta$Family[1],"\n", Conta$Family[2],"\n", Conta$Family[3],"\n", Conta$Family[4],"\n", Conta$Family[5])

# Export normalized and raw ASV tables
otu.noConta=subset(otu.o, !rownames(otu.o) %in% rownames(Conta))[,-c(31,32)]
colnames(otu.noConta)
otu.noConta.f=merge(otu.noConta, tax, by="row.names")
write.table(otu.noConta.f, "./outputs/ASVs_noContanoOut.raw.txt",  quote = FALSE, row.names=F, sep = "\t") #define sample range
message("Number of ASVs used in the analysis: ", length(rownames(otu.noConta.f)))
