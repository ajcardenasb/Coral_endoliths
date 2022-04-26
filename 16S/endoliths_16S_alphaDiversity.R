
############################################################
##################### Alpha-diversity ######################
############################################################
library(vegan)
library(ggplot2)
library(plyr)
library(gridExtra)
library(dplyr)
library(GUniFrac)
library(lme4)
library(emmeans)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/16S/ASVs/")

asv=read.table("./outputs/ASVs_noContanoOut.raw.txt", header = TRUE, row.names = 1)
map=read.table("./Input_files/skeleton16S-metadata.txt", header = TRUE, row.names = 1, sep ='\t')
map$Tissue=gsub("Green band", "Endolithic band", map$Tissue)
map$Tissue=gsub("White band", "Skeleton", map$Tissue)

## aggregated to families
# asv.tax.ag=aggregate(asv[, 1:65], by = list(asv[, 70]), FUN =  sum) #define sample range and group factor
# fams=asv.tax.ag[,-1]
# rownames(fams)=asv.tax.ag$Group.1

##rarefying
cnts=t(asv[, 1:65])
min(rowSums(cnts))
asv.rar=Rarefy(cnts, min(rowSums(cnts)))$otu.tab.rff
#asv.grp=merge(asv.rar, map, by="row.names")
#rownames(asv.grp)= asv.grp[,1]
#asv.grp=asv.grp[,-1]

###not rarefying
# cnts=t(asv[, 1:65])
# asv.rar=cnts
# asv.grp=merge(cnts, map, by="row.names")
# rownames(asv.grp)= asv.grp[,1]
# asv.grp=asv.grp[,-1]


###plots
alpha=as.data.frame(t(estimateR(asv.rar)))
alpha$Shannon=diversity(asv.rar, index = "shannon")#$shannon
alpha$treatment=map$Treatment[match(rownames(alpha),rownames(map))]
alpha$species=map$Species[match(rownames(alpha),rownames(map))]
alpha$tissue=map$Tissue[match(rownames(alpha),rownames(map))]
alpha$genotype=map$Genotype[match(rownames(alpha),rownames(map))]
alpha$tissue=factor(alpha$tissue, levels = c("Tissue", "Endolithic band", "Skeleton"))
alpha_cnts=subset(alpha, treatment == "Control")

# boxplots
pdf("./outputs/ASVs_skeleton16S_alpha_families.pdf", width=6.5,height=3, pointsize = 12)
ggplot(alpha_cnts, aes(x=species, y=Shannon, fill=species)) + 
  stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) + 
  scale_fill_manual(values=c("#8f2d56", "#218380"))  + facet_grid(~tissue) + 
  theme_classic() + labs( y= "Shannon bacterial family diversity", x="") 
dev.off()


##################################################
##################### Stats ######################
##################################################

### full model = only controls
alpha_cnt=subset(alpha, treatment == "Control")
all_model=lmer(Shannon~species*tissue + (1 |genotype), data = alpha_cnt)
anova(all_model)
#comparing species
all_pairs = emmeans(all_model, pairwise ~ species|tissue, weights = "proportional", adjust="none")
rbind(all_pairs$contrasts, adjust="fdr")

#comparing tissues
all_pairs = emmeans(all_model, pairwise ~ tissue|species, weights = "proportional", adjust="none")
rbind(all_pairs$contrasts, adjust="fdr")
