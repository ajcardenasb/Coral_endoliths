library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/16S/ASVs/")

#############################################################
#####Taxonomic profiles of the 20 most abundant families#####
#############################################################

asv=read.table("./outputs/ASVs_noContanoOut.raw.txt", header = TRUE, row.names = 1)
map=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = T, row.names = 1, sep = "\t")
map$Tissue=gsub("Green band", "Endolithic band", map$Tissue)
map$Tissue=gsub("Green band", "Endolithic band", map$Tissue)

names(asv)
asv.tax.ag=aggregate(asv[, 1:65], by = list(asv[, 70]), FUN =  sum) #define sample range and group factor
#topFamilies=asv.tax.ag[order(rowSums(asv.tax.ag[, 2:ncol(asv.tax.ag)]),decreasing = TRUE),][1:20,1]
topFamilies=asv.tax.ag[order(rowSums(asv.tax.ag[, grepl("[PG]C[1-6][BS]", colnames(asv.tax.ag))]),decreasing = TRUE),][1:10,1] # top only in skeleton control samples 
fam.top=subset(asv.tax.ag, asv.tax.ag$Group.1 %in% topFamilies) 
fam.bot=subset(asv.tax.ag, !asv.tax.ag$Group.1 %in% topFamilies) 
fam.bot$Group.1=gsub(".*","zOthers", fam.bot$Group.1)
others=aggregate(fam.bot[, 2:ncol(fam.bot)], by = list(fam.bot[, 1]), FUN =  sum)
all.2 =rbind(fam.top, others)
all.l=melt(all.2, id.vars=c("Group.1"), variable.name = "Family", value.name = "Abundance")
colnames(all.l)=c("Family","Sample","Abundance")


## Add sample information
all.l$Tissue=map$Tissue[match(all.l$Sample, rownames(map))]
all.l$Treatment=map$Treatment[match(all.l$Sample, rownames(map))]
all.l$Species=map$Species[match(all.l$Sample, rownames(map))]
all.l$Genotype=map$Genotype[match(all.l$Sample, rownames(map))]
all.l$Tissue=factor(all.l$Tissue, levels = c("Tissue","Endolithic band", "Skeleton"))
all.l$Family=gsub("Family_XII", "Clostridiales_Family_XII", all.l$Family )
all.cnt=subset(all.l, Treatment = "Control")
final=all.cnt %>% group_by(Species, Family, Tissue) %>% summarise(Abundance=sum(Abundance))


## Plot
P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")
P10=c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A", "#66A61E", "#E6AB02" ,"#780116","#A6761D", "#2364aa", "#3da5d9", "#ababab")

pdf("./outputs/mean_barplots_Goniastrea_16S.pdf",  width = 7, height =3, pointsize = 12) 
ggplot() +geom_bar(aes(y = Abundance, x = Species, fill = Family), data = final, stat="identity", position = "fill") +  theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="Host colony") + scale_fill_manual(values=P10) + facet_grid(~Tissue) + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'bottom') +  theme_classic() + guides(fill=guide_legend(ncol=1))
dev.off() 

# pdf("./outputs/mean_barplots_Goniastrea_16S.pdf",  width = 7, height =7, pointsize = 12) 
# ggplot() +geom_bar(aes(y = Abundance, x = Species, fill = Family), data = final, stat="identity", position = "fill") +  theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="Host colony") + scale_fill_manual(values=P21) + facet_grid(~Tissue) + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'bottom') +  theme_classic() + guides(fill=guide_legend(ncol=1))
# dev.off() 

