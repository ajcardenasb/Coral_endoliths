setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/taxonomy/")
library(tidyr)
library(ggplot2)
library(gridExtra)

ancom_band=read.table("outputs/ANCOM-BC_Families_species_band", sep = "\t", header = T)
ancom_skel=read.table("outputs/ANCOM-BC_Families_species_skel", sep = "\t", header = T)
cnts=read.table("Input_files/Families_counts_500", sep = "\t", header = T, row.names = 1)
cnts_clr=apply(cnts, 2, clr)

tax=rbind(ancom_band, ancom_skel)
tax=tidyr::separate(tax, Taxa, c("Superkingdom", "Phylum", "Class"), sep =";")

# replicated heatmaps
clr_cnt=cnts_clr[, grepl("[GP]C[0-9][BS]", colnames(cnts_clr))]
cnts_clr_sort=clr_cnt[, c(1,10,3,20,19,15,18,9,23,12,7,11,22,6,17,16,2,4, 24,5,14,21,8,13)]
clr_sig=subset(cnts_clr_sort, rownames(cnts_clr_sort) %in%  rownames(subset(ancom_band, Effect_size_W > 2)) |  rownames(cnts_clr_sort) %in%  rownames(subset(ancom_skel,  Effect_size_W > 2)) ) 
colour_taxa=data.frame(domain=tax$Superkingdom[match(rownames(clr_sig),rownames(tax))])
rownames(clr_sig)=paste(rownames(clr_sig), tax$Class[match(rownames(clr_sig),rownames(tax))], sep = " -" )
rownames(colour_taxa)=rownames(clr_sig)
#labels=paste(rownames(clr_sig), tax$Class[match(rownames(clr_sig),rownames(tax))], sep = " -" )

pdf("outputs/heatmaps_ANCOMBC_species_taxa.pdf", height = 12, width = 9)
pheatmap(clr_sig, color = colorRampPalette(c("#1f77b4", "white", "#f94144"))(50),  angle_col = "0",  scale = "row", cluster_row= T, cluster_col = FALSE, border_color = NA,  cellwidth = 5, cellheight =  5, fontsize_row = 6, fontsize_col= 8, legend = F,  annotation_row = colour_taxa, gaps_col = c(6,12,18), annotation_colors =  list(domain =c(Archaea="#f94144", Bacteria="#f8961e",Eukaryota="#90be6d",Viruses="#43aa8b") ))
dev.off()


## Volcano plots between treatments

ancom_gon_band=read.table("outputs/volcano_treatment_Gon_band", sep = "\t", header = T)
ancom_gon_skel=read.table("outputs/volcano_treatment_Gon_skel", sep = "\t", header = T)
ancom_por_band=read.table("outputs/volcano_treatment_Por_band", sep = "\t", header = T)
ancom_por_skel=read.table("outputs/volcano_treatment_Por_skel", sep = "\t", header = T)

tax_col=c( "#90be6d", "#000000")
a=ggplot(ancom_gon_band, aes(x = Beta, y = -log10(pval), color = Superkingdom)) + scale_colour_manual(values = tax_col) + ggtitle(label = "Goniastrea endolithic band", subtitle = "Enriched in control - Enriched in heat") + geom_point(size = 2.5, alpha = 1, na.rm = T) + theme_bw(base_size = 14) +  theme(legend.position = "right") + xlab(expression(log[2]("Effect size"))) + ylab(expression(-log[10]("p-val")))  + geom_hline(yintercept = 1.3, colour="#990000", linetype="dashed") + scale_y_continuous(trans = "log1p", limits = c(0,8))

tax_col=c("#f94144","#f8961e", "#90be6d","#000000","#000000", "#43aa8b")
b=ggplot(ancom_gon_skel, aes(x = Beta, y = -log10(pval), color = Superkingdom)) + scale_colour_manual(values = tax_col) + ggtitle(label = "Goniastrea skeleton", subtitle = "Enriched in control - Enriched in heat") + geom_point(size = 2.5, alpha = 1, na.rm = T) + theme_bw(base_size = 14) +  theme(legend.position = "right") + xlab(expression(log[2]("Effect size"))) + ylab(expression(-log[10]("p-val"))) +  geom_hline(yintercept = 1.3, colour="#990000", linetype="dashed") + scale_y_continuous(trans = "log1p", limits = c(0,8))

tax_col=c("#f8961e", "#90be6d","#000000", "#43aa8b")
c=ggplot(ancom_por_band, aes(x = Beta, y = -log10(pval), color = Superkingdom)) + scale_colour_manual(values = tax_col) + ggtitle(label = "Porites endolithic band", subtitle = "Enriched in control - Enriched in heat") + geom_point(size = 2.5, alpha = 1, na.rm = T) + theme_bw(base_size = 14) +  theme(legend.position = "right") + xlab(expression(log[2]("Effect size"))) + ylab(expression(-log[10]("p-val"))) +  geom_hline(yintercept = 1.3, colour="#990000", linetype="dashed") + scale_y_continuous(trans = "log1p", limits = c(0,8))

tax_col=c("#f8961e", "#90be6d","#000000", "#43aa8b")
d=ggplot(ancom_por_skel, aes(x = Beta, y = -log10(pval), color = Superkingdom)) + scale_colour_manual(values = tax_col) + ggtitle(label = "Porites skeleton", subtitle = "Enriched in control - Enriched in heat") + geom_point(size = 2.5, alpha = 1, na.rm = T) + theme_bw(base_size = 14) +  theme(legend.position = "right") + xlab(expression(log[2]("Effect size"))) + ylab(expression(-log[10]("p-val"))) +  geom_hline(yintercept = 1.3, colour="#990000", linetype="dashed") + scale_y_continuous(trans = "log1p", limits = c(0,8))

pdf("outputs/volcanos_ANCOMBC_treamtment_taxa.pdf", height = 7, width = 10, pointsize = 10)
grid.arrange(a,b,c,d,ncol=2)
dev.off()

