library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/taxonomy/")

map=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = TRUE, row.names = 1, sep ='\t')
kai_report=read.table("Input_files/Endoliths_kaiju_taxonomy_500", sep = "\t", header = T, fill = T, quote = "")
fam=read.table("Input_files/Families_counts_500", header = T, quote = "", sep = "\t")
dom=read.table("Input_files/Superkingdom_counts_500", header = T,  quote = "", sep = "\t")
mapp=read.table("../master_tables/mappedReads_500", header = T,  quote = "", sep = "\t")


### Superkingdom barplots
dom.l=reshape2::melt(dom, id.vars=c("Group.1"), variable.name = "Sample", value.name = "Abundance")
dom.l$Tissue=map$Tissue[match(dom.l$Sample, rownames(map))]
dom.l$Treatment=map$Treatment[match(dom.l$Sample, rownames(map))]
dom.l$Species=map$Species[match(dom.l$Sample, rownames(map))]
dom.l$Genotype=map$Genotype[match(dom.l$Sample, rownames(map))]

#only classified
col_tax=c("#f94144","#f8961e","#90be6d","#577590")
ggplot() +geom_bar(aes(y = Abundance, x = Species, fill = Group.1), data = subset(dom.l, Treatment == "Control" & !Group.1== "Short-contigs" &  !Group.1== "Unclassified"), stat="identity", position = "fill")  + labs( y= "Percentage of metagenomic ORFs", x="")  + facet_grid(~Tissue) +  theme_classic() + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'bottom', legend.key = element_blank(), strip.background = element_blank()) + guides(fill=guide_legend(ncol=5)) + scale_fill_manual(values=col_tax)

#with unclassified
col_tax=c("#f94144","#f8961e","#90be6d","#8d99ae","#6c757d","#577590")
ggplot() +geom_bar(aes(y = Abundance, x = Species, fill = Group.1), data = subset(dom.l, Treatment == "Control" ), stat="identity", position = "fill")  + labs( y= "Percentage of metagenomic ORFs", x="")  + facet_grid(~Tissue) +  theme_classic() + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'bottom', legend.key = element_blank(), strip.background = element_blank()) + guides(fill=guide_legend(ncol=5)) + scale_fill_manual(values=col_tax)


## with replicates
col_tax=c("#f94144","#f8961e","#90be6d","#8d99ae","#6c757d","#577590")
ggplot() +geom_bar(aes(y = Abundance, x = Genotype, fill = Group.1), data = subset(dom.l, Treatment == "Control" ), stat="identity", position = "fill")  + labs( y= "Percentage of metagenomic ORFs", x="")  + facet_grid(~Tissue) +  theme_classic() + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'bottom', legend.key = element_blank(), strip.background = element_blank(), legend.title = element_blank()) + guides(fill=guide_legend(ncol=6)) + scale_fill_manual(values=col_tax)



### Family barplots between species

# fixing family categories
# kai_report$Family=trimws(kai_report$Family, which = "left", whitespace = "[ \t\r\n]")
# kai_report$label=ifelse(kai_report$Family == " NA" | kai_report$Family == "NA" | is.na(kai_report$Family), as.character(paste("Unclassified", kai_report$Class, sep = " ")), as.character(kai_report$Family))
# kai_report$label=ifelse(kai_report$label == "Unclassified_NA", as.character(paste("Unclassified", kai_report$Phylum, sep = " ")), as.character(kai_report$label))
# kai_report$label=ifelse(kai_report$label == "Unclassified NA", as.character(paste("Unclassified", kai_report$Superkingdom, sep = " ")), as.character(kai_report$label))
# kai_report$label=ifelse(kai_report$label %in% c("3", "62", "89"), "Bacteria candidate division", as.character(kai_report$label))
# 

#as.character(wes_palette("Zissou1", 10, type = "continuous"))
#P10=c("#1f77b4","#ff7f0f","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#16bece","#D6D6D6")
#P20=c("#1f77b4","#ff7f0f","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#16bece", "#669900","#99cc33","#ccee66","#006699","#3399cc","#990066","#cc3399","#ff6600","#ff9900","#ffcc00","#D6D6D6")
#P10=rev(c("#ABB4C4","#2077b4","#aec7e8","#ff7f0e","#ffbb78","#2ca02c","#98df89","#17becf","#9edae5","#e377c2","#f7b6d2"))

P10=c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A", "#66A61E", "#E6AB02" ,"#780116","#A6761D", "#2364aa", "#3da5d9", "#ababab")

arc=subset(fam, fam$Group.1 %in% subset(kai_report, Superkingdom == "Archaea")$label) #98
topFamilies=arc[order(rowSums(arc[,2:ncol(arc)]),decreasing = TRUE),][c(1:10),1] #[c(3:12),1]
arc$Group.1=ifelse(arc$Group.1 %in% topFamilies, as.character(arc$Group.1), gsub(".*","zOthers", arc$Group.1))
arc.l=reshape2::melt(arc, id.vars=c("Group.1"), variable.name = "Sample", value.name = "Abundance")
arc.l$Tissue=map$Tissue[match(arc.l$Sample, rownames(map))]
arc.l$Treatment=map$Treatment[match(arc.l$Sample, rownames(map))]
arc.l$Species=map$Species[match(arc.l$Sample, rownames(map))]
arc.l$Genotype=map$Genotype[match(arc.l$Sample, rownames(map))]
arc.l=arc.l %>% group_by(Group.1, Treatment, Species, Tissue ) %>% summarise(Abundance = sum(Abundance))
arc_plot=ggplot() +geom_bar(aes(y = Abundance, x = Species, fill = Group.1), data = subset(arc.l, Treatment == "Control"), stat="identity", position = "fill")  + labs( y= "Percentage of archaeal ORFs", x="")  +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'none', legend.key = element_blank(), strip.background = element_blank(), legend.title = element_blank()) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=P10) +facet_grid(~Tissue, drop = T)
arc_leg=ggplot() +geom_bar(aes(y = Abundance, x = Species, fill = Group.1), data = subset(arc.l, Treatment == "Control"), stat="identity", position = "fill")  + labs( y= "Percentage of archaeal ORFs", x="")  +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'right', legend.key = element_blank(), strip.background = element_blank(), legend.title = element_blank()) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=P10) +facet_grid(~Tissue, drop = T)

bac=subset(fam, fam$Group.1 %in% subset(kai_report, Superkingdom == "Bacteria")$Family)#769 / 490
topFamilies=bac[order(rowSums(bac[,2:ncol(bac)]),decreasing = TRUE),][c(1:10),1]
bac$Group.1=ifelse(bac$Group.1 %in% topFamilies, as.character(bac$Group.1), gsub(".*","zOthers", bac$Group.1))
bac.l=reshape2::melt(bac, id.vars=c("Group.1"), variable.name = "Sample", value.name = "Abundance")
bac.l$Tissue=map$Tissue[match(bac.l$Sample, rownames(map))]
bac.l$Treatment=map$Treatment[match(bac.l$Sample, rownames(map))]
bac.l$Species=map$Species[match(bac.l$Sample, rownames(map))]
bac.l$Genotype=map$Genotype[match(bac.l$Sample, rownames(map))]
bac.l=bac.l %>% group_by(Group.1, Treatment, Species, Tissue ) %>% summarise(Abundance = sum(Abundance))

bac_plot=ggplot() +geom_bar(aes(y = Abundance, x = Species, fill = Group.1), data = subset(bac.l, Treatment == "Control"), stat="identity", position = "fill")  + labs( y= "Percentage of bacterial ORFs", x="")  +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'none', legend.key = element_blank(), strip.background = element_blank(), legend.title = element_blank()) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=P10) +facet_grid(~Tissue, drop = T)
bac_leg=ggplot() +geom_bar(aes(y = Abundance, x = Species, fill = Group.1), data = subset(bac.l, Treatment == "Control"), stat="identity", position = "fill")  + labs( y= "Percentage of bacterial ORFs", x="")  +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'right', legend.key = element_blank(), strip.background = element_blank(), legend.title = element_blank()) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=P10) +facet_grid(~Tissue, drop = T)

euk=subset(fam, fam$Group.1 %in% subset(kai_report, Superkingdom == "Eukaryota")$Family) #671
topFamilies=euk[order(rowSums(euk[,2:ncol(euk)]),decreasing = TRUE),][c(1:10),1]
euk$Group.1=ifelse(euk$Group.1 %in% topFamilies, as.character(euk$Group.1), gsub(".*","zOthers", euk$Group.1))
euk.l=reshape2::melt(euk, id.vars=c("Group.1"), variable.name = "Sample", value.name = "Abundance")
euk.l$Tissue=map$Tissue[match(euk.l$Sample, rownames(map))]
euk.l$Treatment=map$Treatment[match(euk.l$Sample, rownames(map))]
euk.l$Species=map$Species[match(euk.l$Sample, rownames(map))]
euk.l$Genotype=map$Genotype[match(euk.l$Sample, rownames(map))]
euk.l=euk.l %>% group_by(Group.1, Treatment, Species, Tissue ) %>% summarise(Abundance = sum(Abundance))
euk_plot=ggplot() +geom_bar(aes(y = Abundance, x = Species, fill = Group.1), data = subset(euk.l, Treatment == "Control"), stat="identity", position = "fill")  + labs( y= "Percentage of microbial eukaryotic ORFs", x="")  +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'none', legend.key = element_blank(), strip.background = element_blank(), legend.title = element_blank()) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=P10) +facet_grid(~Tissue, drop = T)
euk_leg=ggplot() +geom_bar(aes(y = Abundance, x = Species, fill = Group.1), data = subset(euk.l, Treatment == "Control"), stat="identity", position = "fill")  + labs( y= "Percentage of microbial eukaryotic ORFs", x="")  +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'right', legend.key = element_blank(), strip.background = element_blank(), legend.title = element_blank()) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=P10) +facet_grid(~Tissue, drop = T)

vir=subset(fam, fam$Group.1 %in% subset(kai_report, Superkingdom == "Viruses")$Family) #76
topFamilies=vir[order(rowSums(vir[,2:ncol(vir)]),decreasing = TRUE),][c(1:10),1]
vir$Group.1=ifelse(vir$Group.1 %in% topFamilies, as.character(vir$Group.1), gsub(".*","zOthers", vir$Group.1))
vir.l=reshape2::melt(vir, id.vars=c("Group.1"), variable.name = "Sample", value.name = "Abundance")
vir.l$Tissue=map$Tissue[match(vir.l$Sample, rownames(map))]
vir.l$Treatment=map$Treatment[match(vir.l$Sample, rownames(map))]
vir.l$Species=map$Species[match(vir.l$Sample, rownames(map))]
vir.l$Genotype=map$Genotype[match(vir.l$Sample, rownames(map))]
vir.l=vir.l %>% group_by(Group.1, Treatment, Species, Tissue ) %>% summarise(Abundance = sum(Abundance))
vir_plot=ggplot() +geom_bar(aes(y = Abundance, x = Species, fill = Group.1), data = subset(vir.l, Treatment == "Control"), stat="identity", position = "fill")  + labs( y= "Percentage of viral ORFs", x="")  +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'none', legend.key = element_blank(), strip.background = element_blank(), legend.title = element_blank()) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=P10) +facet_grid(~Tissue, drop = T)
vir_leg=ggplot() +geom_bar(aes(y = Abundance, x = Species, fill = Group.1), data = subset(vir.l, Treatment == "Control"), stat="identity", position = "fill")  + labs( y= "Percentage of viral ORFs", x="")  +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'right', legend.key = element_blank(), strip.background = element_blank(), legend.title = element_blank()) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=P10) +facet_grid(~Tissue, drop = T)

pdf("outputs/Families_barplots.pdf", width = 7, height = 7, pointsize = 12)
gridExtra::grid.arrange(arc_plot, bac_plot, euk_plot,  vir_plot, ncol=2, nrow=2)
dev.off()

pdf("outputs/Families_legends.pdf", width = 7, height = 7, pointsize = 12)
gridExtra::grid.arrange(arc_leg,bac_leg,euk_leg,vir_leg, ncol=2, nrow=2)
dev.off()


