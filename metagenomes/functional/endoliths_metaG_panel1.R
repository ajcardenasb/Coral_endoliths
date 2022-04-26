setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/functional/")
library(ggplot2)
library(phyloseq)
library(data.table)

###############################
##ordination between species ##
###############################


kos=read.table("Input_files/KO_counts", header = T, row.names = 1)

kegg=read.table("~/Documents/Bioinformatics_scripts/KEGG_files/hierarchy_ko00001", header = T, fill = T, sep = "\t")
map=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = TRUE, row.names = 1, sep ='\t')

metaKos=subset(kegg, L1 %in% " Metabolism")$KO
Ckos=read.table("~/Documents/Bioinformatics_scripts/KEGG_files/ko_C_metabolism", fill = T, sep = "\t")
Ckos$V2=gsub("ko:", "",Ckos$V2 )
Nkos=read.table("~/Documents/Bioinformatics_scripts/KEGG_files/ko_N_metabolism",  fill = T, sep = "\t")
Nkos$V2=gsub("ko:", "",Nkos$V2 )

otu.t= otu_table(kos, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax=data.frame(KOs=rownames(otu.t),Met=ifelse(rownames(otu.t) %in% metaKos, "yes", "no"), Cmet=ifelse(rownames(otu.t) %in% Ckos$V2, "yes", "no"), Nmet=ifelse(rownames(otu.t) %in% Nkos$V2, "yes", "no"))
rownames(tax)=rownames(otu.t)
tax.t= tax_table(as.matrix(tax))
phy.meta= phyloseq(otu.t, sam.t, tax.t)

P2=c("#8f2d56", "#218380")

phy.meta.t=microbiome::transform(phy.meta, transform = "clr", target = "OTU", shift = 0, scale = 1)

cnt=subset_samples(phy.meta.t, Treatment=="Control" )

### all metabolims KOs
met_kos=subset_taxa(cnt, Met== "yes" )
PCOA_met = ordinate(met_kos, method = "RDA", distance = "euclidean")
b=plot_ordination(met_kos,PCOA_met, color = "Species", shape = "Tissue") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P2) +  labs(  tag = "B", title = "Metabolism KOs") +  theme_classic() + theme(legend.position = "none") # + stat_ellipse(geom  = "polygon",  alpha = 0.1 )

### C metabolism KOs
c_kos=subset_taxa(cnt, Cmet == "yes" )
PCOA_c = ordinate(c_kos, method = "RDA", distance = "euclidean")
c=plot_ordination(c_kos,PCOA_c , color = "Species", shape = "Tissue") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P2) +  labs(  tag = "C", title = "Carbon metabolism KOs")  +  theme_classic() + theme(legend.position = "none")#+ theme(legend.key = element_blank()) #+ stat_ellipse(geom  = "polygon",  alpha = 0.1 )

### N metabolism KOs
n_kos=subset_taxa(cnt, Nmet == "yes" )
PCOA_n = ordinate(n_kos, method = "RDA", distance = "euclidean")
d=plot_ordination(n_kos,PCOA_n , color = "Species", shape = "Tissue") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P2) +  labs(  tag = "D", title = "Nitrogen metabolism KOs") +  theme_classic() + theme(legend.position = "none")#+ stat_ellipse(geom  = "polygon",  alpha = 0.1 )

pdf("outputs/ordination_panel_KOs.pdf", height = 3, width = 8, pointsize = 8)
grid.arrange(b,c,d, ncol=3, nrow=1)
dev.off()


########################################################
##################### Volcano plots ####################
########################################################

L2_col= c("#669900","#ccee66","#006699","#3399cc", "#990066", "#DA6CB6", "#ff6600","#ff9900", "#ffcc00",   "#fff275", "#9AA0A8", "#000000","#000000")
ancom_band=read.table("outputs/volcano_species_band", sep = "\t", header = T)
ancom_band$L1=kegg$L1[match(rownames(ancom_band), kegg$KO)]
ancom_band$L2=ifelse(ancom_band$DA.SpeciesPorites=="FALSE", "no DA", as.character(ancom_band$L2))
ancom_band_s=subset(ancom_band, L1 %in% " Metabolism" & !L2 == " Not included in regular maps")

ancom_skel=read.table("outputs/volcano_species_skel", sep = "\t", header = T)
ancom_skel$L1=kegg$L1[match(rownames(ancom_skel), kegg$KO)]
ancom_skel$L2=ifelse(ancom_skel$DA.SpeciesPorites=="FALSE", "no DA", as.character(ancom_skel$L2))
ancom_skel_s=subset(ancom_skel,  L1 %in% " Metabolism" & !L2 == " Not included in regular maps")


vol1=ggplot(ancom_band_s, aes(x = Beta.SpeciesPorites, y = -log10(qval.SpeciesPorites), color = L2)) + scale_colour_manual(values = L2_col) + ggtitle(label = "Endolithic band", subtitle = "Enriched in Goniastrea - Enriched in Porites") + geom_point(size = 2.5, alpha = 1, na.rm = T) + theme_bw(base_size = 14) +  theme(legend.position = "none") + xlab(expression(log[2]("Effect size (W)"))) + ylab(expression(-log[10]("p-val"))) + geom_hline(yintercept = 1.30102, colour="#990000", linetype="dashed")  + geom_vline(xintercept = 1, colour="#990000", linetype="dashed") + geom_vline(xintercept = -1, colour="#990000", linetype="dashed") + scale_y_continuous(trans = "log1p", limits = c(0,8))
vol2=ggplot(ancom_skel_s, aes(x = Beta.SpeciesPorites, y = -log10(qval.SpeciesPorites), color = L2)) + scale_colour_manual(values = L2_col) + ggtitle(label = "Skeleton", subtitle = "Enriched in Goniastrea - Enriched in Porites") + geom_point(size = 2.5, alpha = 1, na.rm = T) + theme_bw(base_size = 14) +  theme(legend.position = "none") + xlab(expression(log[2]("Effect size (W)"))) + ylab(expression(-log[10]("p-val"))) + geom_hline(yintercept = 1.30102, colour="#990000", linetype="dashed")  + geom_vline(xintercept = 0.5, colour="#990000", linetype="dashed") + geom_vline(xintercept = -0.5, colour="#990000", linetype="dashed") + scale_y_continuous(trans = "log1p", limits = c(0,8))

pdf("outputs/volcanos_ANCOMBC_species_kos.pdf", height = 4, width = 8, pointsize = 10)
grid.arrange(vol1,vol2,ncol=2)
dev.off()

pdf("outputs/volcanos_ANCOMBC_species_kos_greenbandOnly.pdf", height = 3.5, width = 4, pointsize = 10)
ggplot(ancom_band_s, aes(x = Beta.SpeciesPorites, y = -log10(qval.SpeciesPorites), color = L2)) + scale_colour_manual(values = L2_col) + ggtitle(label = "Endolithic band", subtitle = "Enriched in Goniastrea - Enriched in Porites") + geom_point(size = 2, alpha = 1, na.rm = T) + theme_bw(base_size = 14) +  theme(legend.position = "none") + xlab(expression(log[2]("Effect size (W)"))) + ylab(expression(-log[10]("p-val"))) + geom_hline(yintercept = 1.30102, colour="#990000", linetype="dashed")  + geom_vline(xintercept = 1, colour="#990000", linetype="dashed") + geom_vline(xintercept = -1, colour="#990000", linetype="dashed") + scale_y_continuous(trans = "log1p", limits = c(0,8))
dev.off()

# pdf("outputs/volcanos_ANCOMBC_legend.pdf", height = 5, width = 7.5, pointsize = 10)
# ggplot(ancom_band_s, aes(x = Beta.SpeciesPorites, y = -log10(pval.SpeciesPorites), color = L2)) + scale_colour_manual(values = L2_col) + ggtitle(label = "Skeleton", subtitle = "Enriched in Goniastrea - Enriched in Porites") + geom_point(size = 2.5, alpha = 1, na.rm = T) + theme_bw(base_size = 14) +  theme(legend.position = "right") + xlab(expression(log[2]("Effect size (W)"))) + ylab(expression(-log[10]("p-val"))) + geom_hline(yintercept = 1.30102, colour="#990000", linetype="dashed") + geom_hline(yintercept = 4.85, colour="#990000", linetype="dashed") + geom_vline(xintercept = 0.5, colour="#990000", linetype="dashed") + geom_vline(xintercept = -0.5, colour="#990000", linetype="dashed") + scale_y_continuous(trans = "log1p", limits = c(0,8))
# dev.off()

########################################################
##################### Nano SIMS ####################
########################################################

nan=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/nanoSIMS/input_files/nanoSIMS_data.txt", header = T, sep = "\t")
nan$Treatment=gsub("27", "Control", nan$Treatment)
nan$Treatment=gsub("33", "Heat", nan$Treatment)
Pspec=c("#8f2d56", "#218380")


### Species only
Species_C=subset(nan, compartment == "both"  & Treatment == "Control" & Experiment == "13C")
c_plot=ggplot(Species_C, aes(Species, average, color = Species, fill = Species)) + geom_col(data = Species_C, position = position_dodge(0.8), width = 1) + ylim(0, 0.05) + geom_errorbar( aes(ymin = average, ymax = average+se), data = Species_C , position = position_dodge(0.8)) + scale_color_manual(values = Pspec) + theme_classic() + scale_fill_manual(values = Pspec)+ theme(legend.position = "none")   + labs(y= "Atom% excess") 
Species_N=subset(nan, compartment == "both"  & Treatment == "Control" & Experiment == "15N")
n_plot=ggplot(Species_N, aes(Species, average, color = Species, fill = Species)) + geom_col(data = Species_N, position = position_dodge(0.8), width = 1) + ylim(0, 0.006) + geom_errorbar( aes(ymin = average, ymax = average+se), data = Species_N , position = position_dodge(0.8)) + scale_color_manual(values = Pspec) + theme_classic() + scale_fill_manual(values = Pspec)+ theme(legend.position = "none")   + labs(y= "Atom% excess") 
pdf("NanoSIMS_species.pdf", height = 2.5, width = 3, pointsize = 10)
gridExtra::grid.arrange(c_plot,n_plot,ncol=2)
dev.off()

#########################################################
############### Shannon per compartments ###############
#########################################################


map= read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = T, sep = "\t", row.names = 1)
kai=read.table("Input_files/KO_counts", header = T, row.names = 1, quote = "", sep = "\t")
kai.r=round(kai, digits = 0)#sweep(kai,2,colSums(kai),"/")

kegg=read.table("~/Documents/Bioinformatics_scripts/KEGG_files/hierarchy_ko00001", header = T, fill = T, sep = "\t")
metaKos=subset(kegg, L1 %in% " Metabolism")$KO
Ckos=read.table("~/Documents/Bioinformatics_scripts/KEGG_files/ko_C_metabolism", fill = T, sep = "\t")
Ckos$V2=gsub("ko:", "",Ckos$V2 )
Nkos=read.table("~/Documents/Bioinformatics_scripts/KEGG_files/ko_N_metabolism",  fill = T, sep = "\t")
Nkos$V2=gsub("ko:", "",Nkos$V2 )

meta_all=subset(kai.r, rownames(kai.r) %in% metaKos)
kai.a=as.data.frame(t(estimateR(t(meta_all))))
kai.a$Shannon=diversity(t(kai.r), index = "shannon")#$shannon
kai.a$Species=map$Species[match(rownames(kai.a),rownames(map))]
kai.a$Tissue=map$Tissue[match(rownames(kai.a),rownames(map))]
kai.a$Treatment=map$Treatment[match(rownames(kai.a),rownames(map))]
kai.a$Genotype=map$Genotype[match(rownames(kai.a),rownames(map))]
kai.s=subset(kai.a, Treatment == "Control")

pdf("outputs/alphaDiversity_species_KOs.pdf", width = 3, height = 2.5, pointsize = 12)
ggplot(data=kai.s, aes(x=Species, y=Shannon, fill=Species)) + scale_fill_manual(values = c("#8f2d56", "#218380"), name = "Species") + stat_boxplot(geom ='errorbar', width = 0.7, lwd=0.5)+  geom_boxplot(width=0.7, lwd=0.5, fatten=1)  +  theme_bw()   + labs(y = "KEGG orthologs Shannon diversity", x = "")  + theme( axis.line = element_line(colour = "grey20"), axis.ticks.length = unit(0.2 , "cm"), legend.position="none",  panel.spacing = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 0, vjust = 0.5, face = "bold"), axis.text.x = element_text(color = "grey20", size = 8, angle = 0, hjust = .5, vjust = .5, face = "bold") ) + facet_grid(~Tissue)
dev.off()

#####################################################
##################### Heat maps  ####################
#####################################################


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

map_control=subset(map, Treatment == "Control" )
spec_clr=as.data.frame(cnts_md[,which(colnames(cnts_md) %in% map_control$Sample)])
spec_clr_sig=subset(spec_clr, rownames(spec_clr) %in% rownames(l3_res1_sub) )
spec_clr_sig$L2=kegg$L2[match(rownames(spec_clr_sig), kegg$L3)]

annotation_colors = list( L2= c(" Amino acid metabolism"="#669900"," Biosynthesis of other secondary metabolites"="#ccee66", " Carbohydrate metabolism"="#006699"," Energy metabolism" = "#3399cc"," Glycan biosynthesis and metabolism"="#990066", " Lipid metabolism"="#DA6CB6", " Metabolism of cofactors and vitamins"="#ff6600", " Metabolism of other amino acids"="#ff9900"," Metabolism of terpenoids and polyketides"="#ffcc00"," Nucleotide metabolism"= "#fff275", " Xenobiotics biodegradation and metabolism"="#9AA0A8"))
spec_clr_sig_contrast=subset(spec_clr_sig,  L2 %in% c(" Amino acid metabolism", " Carbohydrate metabolism"," Energy metabolism" , " Glycan biosynthesis and metabolism", " Lipid metabolism",  "Metabolism of terpenoids and polyketides", " Metabolism of other amino acids")  ) 
color_rows=data.frame(L2=kegg$L2[match(rownames(spec_clr_sig_contrast), kegg$L3)])# %>% arrange(L2)
rownames(color_rows)=rownames(spec_clr_sig_contrast)
hp_spec=spec_clr_sig_contrast[,c(1,10,3,20,19,15,18,9,23,12,7,11)] #,22,6,17,16,2,4,24,5,14,21,8,13)]
test=as.matrix(hp_spec)
hp_spec_final=test[match(l3_res1_sub_sorted$L3, rownames(test)), ]
#pdf("outputs/ANCOM-BC_speciesL3_heatmap.pdf", width = 15, height = 5, pointsize = 12)
pheatmap(hp_spec_final, color = colorRampPalette(c("white","white", "black"))(50),  cellwidth = 9,  cellheight =  6, fontsize_row = 7, fontsize_col= 8,  legend = F, gaps_col = c(6,12),  cluster_rows = F,cluster_col = F, scale = "row", annotation_row = color_rows, annotation_colors = annotation_colors)
#dev.off()


### barplots with effect sizes
l3_res1_sub$L3=rownames(l3_res1_sub)
l3_res1_sub2=subset(l3_res1_sub, L3 %in% rownames(spec_clr_sig_contrast))
l3_res1_sub_sorted=l3_res1_sub2 %>% arrange(W_Species)
head(l3_res1_sub_sorted)
pdf("outputs/ANCOM-BC_speciesL3_mirror.pdf", width = 7, height = 7, pointsize = 12)
ggplot(data=l3_res1_sub2, aes(x=reorder(L3, W_Species), y=W_Species)) + geom_bar(stat="identity", position = "dodge") + coord_flip()   + guides(fill = guide_legend(reverse = TRUE)) + theme_classic() + labs(y="Effect size")
dev.off()


       
