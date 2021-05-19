
############################################################
##################### Alpha-diversity ######################
############################################################
library(vegan)
library(ggplot2)
library(plyr)
library(gridExtra)
library(dplyr)
library(GUniFrac)
setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/16S/ASVs/")

asv=read.table("./outputs/ASVs_noContanoOut.raw.txt", header = TRUE, row.names = 1)
map=read.table("./Input_files/skeleton16S-metadata.txt", header = TRUE, row.names = 1, sep ='\t')
map$Tissue=gsub("Green band", "Endolithic band", map$Tissue)
map$Tissue=gsub("White band", "Skeleton band", map$Tissue)

## aggregated to families
asv.tax.ag=aggregate(asv[, 1:65], by = list(asv[, 70]), FUN =  sum) #define sample range and group factor
fams=asv.tax.ag[,-1]
rownames(fams)=asv.tax.ag$Group.1

##rarefying
cnts=t(fams[, 1:65])
min(rowSums(cnts))
asv.rar=Rarefy(cnts, 1966)$otu.tab.rff
# asv.grp=merge(asv.rar, map, by="row.names")
# rownames(asv.grp)= asv.grp[,1]
# asv.grp=asv.grp[,-1]

###not rarefying
# cnts=t(asv[, 1:65])
# asv.grp=merge(cnts, map, by="row.names")
# rownames(asv.grp)= asv.grp[,1]
# asv.grp=asv.grp[,-1]


###plots
alpha=as.data.frame(t(estimateR(asv.rar,  smallsample = TRUE)))
alpha$Shannon=diversity(asv.rar, index = "shannon")#$shannon
alpha$treatment=map$Treatment[match(rownames(alpha),rownames(map))]
alpha$species=map$Species[match(rownames(alpha),rownames(map))]
alpha$tissue=map$Tissue[match(rownames(alpha),rownames(map))]
alpha$tissue=factor(alpha$tissue, levels = c("Tissue", "Endolithic band", "Skeleton band"))
alpha_cnts=subset(alpha, treatment == "Control")

# boxplots
pdf("./outputs/ASVs_skeleton16S_species.pdf", width=6.5,height=3, pointsize = 12)
ggplot(alpha_cnts, aes(x=species, y=S.ACE, fill=species)) + stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) +  scale_fill_manual(values=c("#8f2d56", "#218380"))  + facet_grid(~tissue) +  theme_classic() + labs( y= "Observed number of ASVs", x="") 
dev.off()
