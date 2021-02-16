
############################################################
##################### Alpha-diversity ######################
############################################################
library(vegan)
library(ggplot2)
library(plyr)
library(gridExtra)
library(ddply)
setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/16S/ASVs/")

asv=read.table("./outputs/ASVs_noContanoOut.raw.txt", header = TRUE, row.names = 1)
map=read.table("./Input_files/skeleton16S-metadata.txt", header = TRUE, row.names = 1, sep ='\t')
map$Tissue=gsub("Green band", "Endolithic band", map$Tissue)
map$Tissue=gsub("White band", "Skeleton band", map$Tissue)

###rarefying
# cnts=t(asv[, 1:65])
# min(rowSums(cnts))
# asv.rar=Rarefy(cnts, 4197)$otu.tab.rff
# asv.grp=merge(asv.rar, map, by="row.names")
# rownames(asv.grp)= asv.grp[,1]
# asv.grp=asv.grp[,-1]

###not rarefying
cnts=t(asv[, 1:65])
asv.grp=merge(cnts, map, by="row.names")
rownames(asv.grp)= asv.grp[,1]
asv.grp=asv.grp[,-1]


###plots
alpha=as.data.frame(t(estimateR(asv.grp[, 1:6491],  smallsample = TRUE)))
alpha$Shannon=diversity(asv.grp[, 1:6491], index = "shannon")#$shannon
alpha$treatment=map$Treatment[match(rownames(alpha),rownames(map))]
alpha$species=map$Species[match(rownames(alpha),rownames(map))]
alpha$tissue=map$Tissue[match(rownames(alpha),rownames(map))]
alpha$tissue=factor(alpha$tissue, levels = c("Tissue", "Endolithic band", "Skeleton band"))
alpha_cnts=subset(alpha, treatment == "Control")

# boxplots
pdf("./outputs/ASVs_skeleton16S_Shannon.pdf", width=6.5,height=3, pointsize = 12)
ggplot(alpha_cnts, aes(x=species, y=Shannon, fill=species)) + stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) +  scale_fill_manual(values=c("#8f2d56", "#218380"))  + facet_grid(~tissue) +  theme_classic() + labs( y= "Observed number of ASVs", x="") 
dev.off()




##################################################
##################### Stats ######################
##################################################


shapiro.test(alpha$S.ACE) # p-value > 0.05 implying we can assume the normality.
shapiro.test(alpha$S.obs) # p-value > 0.05 implying we can assume the normality.


### ACE
gon.t=subset(alpha, species == "Goniastrea" & tissue =="Tissue")
t.test(S.ACE~ treatment, gon.t )# p-value = 0.58
gon.g=subset(alpha, species == "Goniastrea" & tissue =="Green band")
t.test(S.ACE~ treatment,gon.g )# p-value = 0.18
gon.w=subset(alpha, species == "Goniastrea" & tissue =="White band")
t.test(S.ACE~ treatment,gon.w )#p-value = 0.59

por.t=subset(alpha, species == "Porites" & tissue =="Tissue")
t.test(S.ACE~ treatment, por.t )# p-value = 0.80
por.g=subset(alpha, species == "Porites" & tissue =="Green band")
t.test(S.ACE~ treatment,por.g )# p-value = 0.68
por.w=subset(alpha, species == "Porites" & tissue =="White band")
t.test(S.ACE~ treatment,por.w )#p-value = 0.72

#Observed

gon.t=subset(alpha, species == "Goniastrea" & tissue =="Tissue")
t.test(S.obs~ treatment, gon.t )# p-value = 0.97
gon.g=subset(alpha, species == "Goniastrea" & tissue =="Green band")
t.test(S.obs~ treatment,gon.g )# p-value = 0.11
gon.w=subset(alpha, species == "Goniastrea" & tissue =="White band")
t.test(S.obs~ treatment,gon.w )#p-value = 0.560

por.t=subset(alpha, species == "Porites" & tissue =="Tissue")
t.test(S.obs~ treatment, por.t )# p-value = 0.69
por.g=subset(alpha, species == "Porites" & tissue =="Green band")
t.test(S.obs~ treatment,por.g )# p-value = 0.56
por.w=subset(alpha, species == "Porites" & tissue =="White band")
t.test(S.obs~ treatment,por.w )#p-value = 0.47






############################################################
#####################  Family level  ######################
############################################################

asv=read.table("./outputs/ASVs_noContanoOut.raw.txt", header = TRUE, row.names = 1)
map=read.table("./Input_files/skeleton16S-metadata.txt", header = TRUE, row.names = 1, sep ='\t')
map$Tissue=gsub("Green band", "Endolithic band", map$Tissue)
map$Tissue=gsub("White band", "Skeleton band", map$Tissue)

fam=aggregate(asv[, 1:65], by = list(asv[, 70]), FUN =  sum) 
rownames(fam)= fam[,1]
fam=fam[,-1]
fam_t=t(fam)
write.table(fam, "outputs/family_counts_ASVs", quote =F, sep = "\t")
###plots
alpha=as.data.frame(t(estimateR(asv.grp[, 1:220],  smallsample = TRUE)))
alpha$Shannon=diversity(asv.grp[, 1:220], index = "shannon")#$shannon
alpha$treatment=map$Treatment[match(rownames(alpha),rownames(map))]
alpha$species=map$Species[match(rownames(alpha),rownames(map))]
alpha$tissue=map$Tissue[match(rownames(alpha),rownames(map))]
alpha$tissue=factor(alpha$tissue, levels = c("Tissue", "Endolithic band", "Skeleton band"))
alpha_cnts=subset(alpha, treatment == "Control")

# boxplots
pdf("./outputs/Alpha_diversity/ASVs_skeleton16S_ObservedAlpha.pdf", width=8,height=4, pointsize = 12)
ggplot(alpha_cnts, aes(x=species, y=Shannon, fill=species)) + stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) +  scale_fill_manual(values=c("#8f2d56", "#218380"))  + facet_grid(~tissue) +  theme_classic() + labs( y= "Observed number of ASVs", x="") 
dev.off()






