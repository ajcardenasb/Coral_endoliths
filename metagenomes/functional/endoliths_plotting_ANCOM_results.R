
library(pheatmap)
library(compositions)
library(data.table)
library(gridExtra)
setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/Functional/")


########################################################
##################### Volcano plots ####################
########################################################
#kos_res1=read.table("outputs/ANCOMBC_Allkos_betweenSpecies_band.txt", sep = "\t", header = T)
#kos_res1_metab=subset(kos_res1, L1 == " Metabolism")
kegg=read.table("~/Documents/Bioinformatics_scripts/KEGG_files/MostUpdated_hierarchy_ko00001", sep = "\t", header = T, fill = T)#kegg=read.table("~/scripts/hierarchy_ko00001", sep = "\t", header = T, fill = T)

pre_ancom_band=read.table("outputs/volcano_AllKos_species_band", sep = "\t", header = T)
pre_ancom_band$L1=kegg$L1[match(rownames(pre_ancom_band), kegg$KO)]
ancom_band=subset(pre_ancom_band, L1 %in% " Metabolism")
ancom_band$L2=ifelse(ancom_band$DA.SpeciesPorites=="FALSE", "no DA", as.character(ancom_band$L2))
ancom_band_s=subset(ancom_band, !L2 %in% c(" Not included in regular maps") & !is.na(L2))

pre_ancom_skel=read.table("outputs/volcano_AllKos_species_skel", sep = "\t", header = T)
pre_ancom_skel$L1=kegg$L1[match(rownames(pre_ancom_skel), kegg$KO)]
ancom_skel=subset(pre_ancom_skel, L1  %in% " Metabolism")
ancom_skel$L2=ifelse(ancom_skel$DA.SpeciesPorites=="FALSE", "no DA", as.character(ancom_skel$L2))
ancom_skel_s=subset(ancom_skel, !L2 %in% c(" Not included in regular maps") & !is.na(L2))
ancom_skel$L2=as.factor(ancom_skel$L2)

L2_col =c("#669900","#ccee66","#006699","#3399cc","#990066","#DA6CB6","#ff6600","#ff9900","#ffcc00","#fff275","#9AA0A8", "#000000")
#L2_col =c("#9AA0A8","#9AA0A8","#006699","#9AA0A8","#990066","#9AA0A8","#9AA0A8","#9AA0A8","#9AA0A8","#9AA0A8","#9AA0A8", "#000000")


a=ggplot(ancom_band_s, aes(x = Beta.SpeciesPorites, y = -log10(qval.SpeciesPorites), color = L2)) + scale_colour_manual(values = L2_col) + ggtitle(label = "Endolithic band", subtitle = "Enriched in Goniastrea - Enriched in Porites") + geom_point(size = 2.5, alpha = 1, na.rm = T) + theme_bw(base_size = 14) +  theme(legend.position = "right") + xlab(expression(log[2]("Effect size (W)"))) + ylab(expression(-log[10]("q-val"))) + geom_hline(yintercept = 1.30102, colour="#990000", linetype="dashed") +  geom_vline(xintercept = 1, colour="#990000", linetype="dashed") + geom_vline(xintercept = -1, colour="#990000", linetype="dashed") + scale_y_continuous(trans = "log1p", limits = c(0,8))
b=ggplot(ancom_skel_s, aes(x = Beta.SpeciesPorites, y = -log10(qval.SpeciesPorites), color = L2)) + scale_colour_manual(values = L2_col) + ggtitle(label = "Skeleton", subtitle = "Enriched in Goniastrea - Enriched in Porites") + geom_point(size = 2.5, alpha = 1, na.rm = T) + theme_bw(base_size = 14) +  theme(legend.position = "right") + xlab(expression(log[2]("Effect size (W)"))) + ylab(expression(-log[10]("q-val"))) + geom_hline(yintercept = 1.30102, colour="#990000", linetype="dashed") +  geom_vline(xintercept = 1, colour="#990000", linetype="dashed") + geom_vline(xintercept = -1, colour="#990000", linetype="dashed") + scale_y_continuous(trans = "log1p", limits = c(0,8))


pdf("outputs/volcanos_ANCOMBC_species_kos.pdf", height = 5, width = 20, pointsize = 10)
grid.arrange(a,b,ncol=2)
dev.off()

pdf("outputs/volcanos_ANCOMBC_legend.pdf", height = 5, width = 20, pointsize = 10)
ggplot(ancom_skel_s, aes(x = Beta.SpeciesPorites, y = -log10(pval.SpeciesPorites), color = L2)) + scale_colour_manual(values = L2_col) + ggtitle(label = "Skeleton", subtitle = "Enriched in Goniastrea - Enriched in Porites") + geom_point(size = 2.5, alpha = 1, na.rm = T) + theme_bw(base_size = 14) +  theme(legend.position = "bottom") + xlab(expression(log[2]("Effect size (W)"))) + ylab(expression(-log[10]("p-val"))) + geom_hline(yintercept = 1.30102, colour="#990000", linetype="dashed") + geom_hline(yintercept = 4.85, colour="#990000", linetype="dashed") + geom_vline(xintercept = 0.5, colour="#990000", linetype="dashed") + geom_vline(xintercept = -0.5, colour="#990000", linetype="dashed") + scale_y_continuous(trans = "log1p", limits = c(0,8))
dev.off()

#####################################################
##################### Heat maps  ####################
#####################################################


col_spe=c("#8f2d56","#218380")

kegg=read.table("~/Documents/Bioinformatics_scripts/KEGG_files/MostUpdated_hierarchy_ko00001", sep = "\t", header = T, fill = T)#, row.names = 5) 
kegg$module=trimws(kegg$module, which = "both", whitespace = "[ \t\r\n]")
map=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = T, sep = "\t")
cnts=read.table("Input_files/KO_counts_500", header = T, row.names = 1, sep = "\t")

cnts$L3=kegg$L3[match(rownames(cnts), kegg$KO)]
cnts$L3=ifelse(is.na(cnts$L3), "Others", as.character(cnts$L3))
cnts_gg=aggregate(cnts[,1:48], by=list(cnts$L3), FUN = sum)
rownames(cnts_gg)=cnts_gg[,1]
cnts_md=apply(cnts_gg[,-1],2,clr)

l3_res1=read.table("outputs/ANCOMBC_l3_betweenSpecies_band.txt", sep = "\t", header = T)
l3_res1_sub=subset(l3_res1, L1 %in% " Metabolism" )
# selection of samples band
map_control=subset(map, Treatment == "Control" )
spec_clr=as.data.frame(cnts_md[,which(colnames(cnts_md) %in% map_control$Sample)])
#spec_clr_sig=subset(spec_clr, rownames(spec_clr) %in% rownames(subset(md_res1, abs(W_Species) >2 )) )
spec_clr_sig=subset(spec_clr, rownames(spec_clr) %in% rownames(l3_res1_sub) )
spec_clr_sig$L2=kegg$L2[match(rownames(spec_clr_sig), kegg$L3)]

annotation_colors = list( L2= c(" Amino acid metabolism"="#669900"," Biosynthesis of other secondary metabolites"="#ccee66", " Carbohydrate metabolism"="#006699"," Energy metabolism" = "#3399cc"," Glycan biosynthesis and metabolism"="#990066", " Lipid metabolism"="#DA6CB6", " Metabolism of cofactors and vitamins"="#ff6600", " Metabolism of other amino acids"="#ff9900"," Metabolism of terpenoids and polyketides"="#ffcc00"," Nucleotide metabolism"= "#fff275", " Xenobiotics biodegradation and metabolism"="#9AA0A8"))
#annotation_colors = list( L2= c(" Amino acid metabolism"="#cb94b3"," Biosynthesis of other secondary metabolites"="#e9a2b4", " Carbohydrate metabolism"="#fdd9b7"," Energy metabolism" = "#fbfbc5"," Glycan biosynthesis and metabolism"="#c3e7bf", " Lipid metabolism"="#a1d3cc", " Metabolism of cofactors and vitamins"="#dbdac7", " Metabolism of terpenoids and polyketides"="#fce5c5"," Nucleotide metabolism"= "#fed8af", " Xenobiotics biodegradation and metabolism"="#febaa1"))



spec_clr_sig_contrast=subset(spec_clr_sig,  L2 %in% c(" Amino acid metabolism", " Carbohydrate metabolism"," Energy metabolism" , " Glycan biosynthesis and metabolism", " Lipid metabolism",  "Metabolism of terpenoids and polyketides", " Metabolism of other amino acids")  ) 
#spec_clr_sig_contrast=spec_clr_sig
color_rows=data.frame(L2=kegg$L2[match(rownames(spec_clr_sig_contrast), kegg$L3)])# %>% arrange(L2)
rownames(color_rows)=rownames(spec_clr_sig_contrast)
hp_spec=spec_clr_sig_contrast[,c(1,10,3,20,19,15,18,9,23,12,7,11)] #,22,6,17,16,2,4,24,5,14,21,8,13)]
pdf("outputs/ANCOM-BC_speciesL3_heatmap.pdf", width = 15, height = 5, pointsize = 12)
pheatmap(hp_spec, color = colorRampPalette(c("white","white", "black"))(50),  cellwidth = 10,  cellheight =  8, fontsize_row = 8, fontsize_col= 8,  legend = F, gaps_col = c(6,12),  cluster_rows = T,cluster_col = F, scale = "row", annotation_row = color_rows, annotation_colors = annotation_colors)
dev.off()

