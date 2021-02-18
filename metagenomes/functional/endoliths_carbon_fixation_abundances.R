setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/functional/")

library(compositions)
library(ggpubr)
library(tidyr)
library(data.table)

cnts=read.table("Input_files/KO_counts_500", header = T, sep = "\t", row.names = 1) 
cnts_clr=as.data.frame(apply(cnts,2,clr))
map=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = T, row.names = 1, sep = "\t")
kegg=read.table("~/Documents/Bioinformatics_scripts/KEGG_files/MostUpdated_hierarchy_ko00001", sep = "\t", header = T, fill = T)#, row.names = 5) 
kegg$module=trimws(kegg$module, which = "both", whitespace = "[ \t\r\n]")


#######################
### Carbon fixation ###
#######################

cnts$L3=kegg$L3[match(rownames(cnts), kegg$KO)]
l3=aggregate(cnts[1:48], by=list(cnts$L3), FUN=sum)
l3_clr=as.data.frame(apply(l3[2:49],2,clr))
rownames(l3_clr)=l3$Group.1
l3_t=as.data.frame(t(l3_clr))
l3_t$Species=map$Species[match(rownames(l3_t), rownames(map))]
l3_t$Treatment=map$Treatment[match(rownames(l3_t), rownames(map))]
l3_t$Group=paste(map$Species, map$Tissue, sep = "_")[match(rownames(l3_t), rownames(map))]
l3_cnt=subset(l3_t, Treatment == "Control")

l3_cnt$Group=factor(l3_cnt$Group, levels =c("Goniastrea_Endolithic band", "Porites_Endolithic band", "Goniastrea_Skeleton", "Porites_Skeleton"))

pdf("outputs/proka_cfix_metabolism_btwSpecies_1.pdf", width = 8, height = 5, pointsize = 12)
par(mfrow=c(1,2))
boxplot(l3_cnt$` Carbon fixation in photosynthetic organisms `~l3_cnt$Group, main = "Carbon fixation pathways in eukaryotes", col=c("#8f2d56","#218380"), ylab = "clr-trasformed counts", xlab = "")
boxplot(l3_cnt$` Carbon fixation pathways in prokaryotes ` ~l3_cnt$Group, main = "Carbon fixation pathways in prokaryotes", col=c("#8f2d56","#218380"), ylab = "clr-trasformed counts", xlab = "")
dev.off()

#### Carbon fixation pathways in prokaryotes

cnts$MD=kegg$MD[match(rownames(cnts), kegg$KO)]
md=aggregate(cnts[1:48], by=list(cnts$MD), FUN=sum)
md_clr=as.data.frame(apply(md[2:49],2,clr))
rownames(md_clr)=md$Group.1
md_t=as.data.frame(t(md_clr))
md_t$Species=map$Species[match(rownames(md_t), rownames(map))]
md_t$Treatment=map$Treatment[match(rownames(md_t), rownames(map))]
md_t$Group=paste(map$Species, map$Tissue, sep = "_")[match(rownames(md_t), rownames(map))]
md_cnt=subset(md_t, Treatment == "Control" )

md_cnt$Group=factor(md_cnt$Group, levels =c("Goniastrea_Endolithic band", "Porites_Endolithic band", "Goniastrea_Skeleton", "Porites_Skeleton"))

pdf("outputs/proka_cfix_metabolism_btwSpecies_2.pdf", width = 8, height = 3, pointsize = 12)
par(mfrow=c(1,3))
boxplot(md_cnt$M00173 ~md_cnt$Group, main = "Reductive citrate cycle\n(Arnon-Buchanan cycle)", col=c("#8f2d56","#218380"), ylab = "clr-trasformed counts", xlab = "")
boxplot(md_cnt$M00377 ~md_cnt$Group, main = "Reductive acetyl-CoA pathway\n(Wood-Ljungdahl pathway) ", col=c("#8f2d56","#218380"), ylab = "clr-trasformed counts", xlab = "")
boxplot(md_cnt$M00376 ~md_cnt$Group, main = "3-Hydroxypropionate\nbi-cycle", col=c("#8f2d56","#218380"), ylab = "clr-trasformed counts", xlab = "")
dev.off()


#### Carbon fixation pathways in prokaryotes stacked barplots
cnts=read.table("Input_files/KO_counts_500", header = T, sep = "\t", row.names = 1) 
map=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = T, row.names = 1, sep = "\t")
kegg=read.table("~/Documents/Bioinformatics_scripts/KEGG_files/MostUpdated_hierarchy_ko00001", sep = "\t", header = T, fill = T)#, row.names = 5) 
kegg$module=trimws(kegg$module, which = "both", whitespace = "[ \t\r\n]")
kegg$module=ifelse(kegg$L3 == " Carbon fixation in photosynthetic organisms ", "Calvin-Benson cycle", as.character(kegg$module))



cnts$MD=kegg$module[match(rownames(cnts), kegg$KO)]
md=aggregate(cnts[1:48], by=list(cnts$MD), FUN=sum)
car_fix_md=subset(md,Group.1 %like%  "Reductive acetyl-CoA pathway" | Group.1 %like%  "Calvin-Benson cycle" | Group.1 %like%  "Reductive citrate cycle" | Group.1 %like%  "3-Hydroxypropionate")
rownames(car_fix_md)=car_fix_md$Group.1
car_fix_l=reshape2::melt(car_fix_md, id.vars=c("Group.1"), variable.name = "Sample", value.name = "Abundance")

## normalize by number of KOs
kegg_cf=subset(kegg, module %in% rownames(car_fix_md))
length(unique((kegg %>% filter( module %like%  "Reductive citrate cycle"))$KO))
car_fix_l$Abundance=ifelse(car_fix_l$Group.1 %like% "Reductive acetyl-CoA pathway", car_fix_l$Abundance/6,ifelse(car_fix_l$Group.1 %like% "Calvin-Benson cycle", car_fix_l$Abundance/1,ifelse(car_fix_l$Group.1 %like% "3-Hydroxypropionate", car_fix_l$Abundance/7, ifelse(car_fix_l$Group.1 %like% "Reductive citrate cycle", car_fix_l$Abundance/19,as.numeric(as.character(car_fix_l$Abundance) )) )))

#add sample variables
car_fix_l$Species=map$Species[match(car_fix_l$Sample, rownames(map))]
car_fix_l$Treatment=map$Treatment[match(car_fix_l$Sample, rownames(map))]
car_fix_l$Group=paste(map$Species, map$Tissue, sep = "_")[match(car_fix_l$Sample, rownames(map))]
car_fix_cnt=subset(car_fix_l, Treatment == "Control" )
car_fix_cnt$Group=factor(car_fix_cnt$Group, levels =c("Goniastrea_Endolithic band", "Porites_Endolithic band", "Goniastrea_Skeleton", "Porites_Skeleton"))
car_fix_cnt$Sample=factor(car_fix_cnt$Sample, levels =c("GC1B", "GC2B", "GC3B", "GC4B", "GC5B", "GC6B","PC1B", "PC2B", "PC3B", "PC4B", "PC5B", "PC6B", "GC1S", "GC2S", "GC3S", "GC4S", "GC5S", "GC6S","PC1S", "PC2S", "PC3S", "PC4S", "PC5S", "PC6S"))
col_tax=c("#8f2d56","#d81159","#ffbc42","#006ba6")


pdf("outputs/proka_cfix_metabolism_btwSpecies_3.pdf", width = 8, height = 3, pointsize = 12)
ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Group.1), data = car_fix_cnt, stat="identity", position = "fill")  + labs( y= "Percentage of metagenomic ORFs", x="")   +  theme_classic() + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'bottom', legend.key = element_blank(), strip.background = element_blank()) + guides(fill=guide_legend(ncol=5)) + scale_fill_manual(values=col_tax)
dev.off()


