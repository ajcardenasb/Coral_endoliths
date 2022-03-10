library(data.table)
library(ggplot2)
library(vegan)
library(VennDiagram)

############################################################
#####################  Venn diagrams  ######################
############################################################

#data in
map=read.table("./Input_files/skeleton16S-metadata.txt", header = TRUE, row.names = 1, sep ='\t')
map$Tissue=gsub("Green band", "Endolithic band", map$Tissue)
map$Tissue=gsub("White band", "Skeleton band", map$Tissue)

asvs=read.table("./outputs/family_counts_ASVs", header = TRUE, row.names = 1)
asvs_skel=asvs[,grepl( "[GP][CH][0-6][SB]", colnames(asvs))]
colnames(asvs_skel)=paste(colnames(asvs_skel), "16S", sep="_")

metag=read.table("../../metagenomes/taxonomy/Input_files/Families_counts_500",  sep = "\t" ,  row.names = 1, header = TRUE)
metadata_kai=read.table("../../metagenomes/taxonomy/Input_files/Endoliths_kaiju_taxonomy_500",  sep = "\t" , quote = "", header = TRUE)
metag_bac=subset(metag, rownames(metag) %in% subset(metadata_kai,  Superkingdom == "Bacteria" )$Family)
colnames(metag_bac)=paste(colnames(metag_bac), "MG", sep="_")

#calculate alpha diversity
alpha_asvs=as.data.frame(t(estimateR(t(asvs_skel),  smallsample = TRUE)))
alpha_asvs$Shannon=diversity(t(asvs_skel), index = "shannon")
alpha_metag=as.data.frame(t(estimateR(t(metag_bac),  smallsample = TRUE)))
alpha_metag$Shannon=diversity(t(metag_bac), index = "shannon")
alpha=rbind(alpha_asvs, alpha_metag)
alpha$Sample=gsub("_.*", "",rownames(alpha))
alpha$treatment=map$Treatment[match(alpha$Sample,rownames(map))]
alpha$species=map$Species[match(alpha$Sample,rownames(map))]
alpha$tissue=map$Tissue[match(alpha$Sample,rownames(map))]
alpha$tissue=factor(alpha$tissue, levels = c("Tissue", "Endolithic band", "Skeleton band"))
alpha$Library=ifelse(rownames(alpha) %like% "_16S", "16S", "Metagenomes")
alpha_cnts=subset(alpha, treatment == "Control")

# boxplots
pdf("./outputs/Alpha_diversity/16SvsMetaG_ObservedFams.pdf", width=7,height=3, pointsize = 12)
ggplot() +geom_bar(aes(y = S.obs, x = Sample, fill = Library), data = alpha_cnts, stat="identity", position = "dodge") +  scale_fill_manual(values=c("#8f2d56", "#218380"))  +  theme_classic() + labs( y= "Observed number of bacterial families", x="") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#Venn diagrams
asvs_skel$families=rownames(asvs_skel)
venn_asv=reshape2::melt(asvs_skel, id.vars=c("families"), variable.name = "Sample", value.name = "Abundance")
list_asvs_gon=subset(venn_asv,Abundance > 50 & Sample %like% "^G")$families
list_asvs_por=subset(venn_asv,Abundance > 50 & Sample %like% "^P")$families

metag_bac$families=rownames(metag_bac)
venn_metag=reshape2::melt(metag_bac, id.vars=c("families"), variable.name = "Sample", value.name = "Abundance")
list_metag_gon=subset(venn_metag,Abundance > 50 & Sample %like% "^G")$families
list_metag_por=subset(venn_metag,Abundance > 50 & Sample %like% "^P")$families

vennGon=venn.diagram(list(list_metag_gon, list_asvs_gon), NULL, category.names = c("MetaG", "16S"), rotation.degree = 90, fill = c("#8f2d56", "#D883A6"), lty="blank", cat.pos = c(180,0),  alpha = rep(0.8, 2), cat.fontfamily = "Helvetica",  fontfamily="Helvetica")
vennPor=venn.diagram(list(list_metag_por, list_asvs_por), NULL, category.names = c("MetaG", "16S"), rotation.degree = 90, fill = c("#218380", "#6CDAD7"), lty="blank", cat.pos = c(180,0),  alpha = rep(0.8, 2), cat.fontfamily = "Helvetica",  fontfamily="Helvetica")

pdf("./outputs/Observed_families_16SvsMetaG_vennDiagrams.pdf", width=4,height=6, pointsize = 20)
grid.arrange(grobTree(vennGon), grobTree(vennPor),  ncol = 1, nrow = 2)
dev.off()
