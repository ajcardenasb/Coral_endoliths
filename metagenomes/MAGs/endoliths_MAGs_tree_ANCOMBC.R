library(ggtree)
library(compositions)
library(dplyr)
library(pheatmap)
library(ANCOMBC)
library(phyloseq)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/")

map=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/16S/ASVs/Input_files/skeleton16S-metadata.txt", header = TRUE, row.names = 1, sep ='\t')
#bin=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/Outputs/MAGs_abundance.txt")
bin=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/Outputs/MAGs_counts.txt")
meta.mag= read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/Input_files/tree_metadata.txt", header = T, row.names = 1)
tree = ggtree::read.tree("Input_files/Endolith_tree_edited.nwk")



#########################
#### Just phylo tree ####
#########################
Pphyla=c("#2077b4","#aec7e8","#ff7f0e", "#ffbb78", "#2ca02c", "#98df89","#17becf", "#9edae5", "#e377c2", "#f7b6d2", "#e33e33",   "#f1c453", "#efea5a", "#16db93" ,  "#B6B6B6")
meta.mag= read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/Input_files/tree_metadata.txt", header = T)
pdf("Outputs/MAGs_tree.pdf", height = 7, width = 7, pointsize = 10)
ggtree(tree,layout="radial") %<+% meta.mag +  geom_tippoint(aes(color=Taxa, size = 5)) + geom_treescale(x=0, y=0, width=0, color = "black") + scale_size_identity() + geom_tiplab(aes(angle=angle), size=2, align=F, linesize=.5, offset = 0.2) + scale_color_manual(values=Pphyla) 
dev.off()


##############################
### ANCOM between species ###
##############################
#bin=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/Outputs/MAGs_abundance.txt")
#tax.mag= read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/Input_files/metadata_MAGs_complete.txt", header = T, row.names = 1, sep = "\t")
bin=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/Outputs/MAGs_counts.txt")
bin.f=round(bin, digits = 0)
map=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = T, sep = "\t", row.names = 1)
#meta.mag= read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/Input_files/tree_metadata.txt", header = T, row.names = 1)
meta.mag= read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/Input_files/tree_metadata.txt", header = T, row.names = 1)
map$Species=as.factor(map$Species)
map$Treatment=as.factor(map$Treatment)

#cnts=read.table("Input_files/metabolic_Kos", header = T, row.names = 1, sep = "\t")

otu.t= otu_table(bin, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
phy_species= phyloseq(otu.t,  sam.t )#, tax.t)

######################
## between species ##
######################

phy_species_band=subset_samples(phy_species, Treatment  == "Control" & Tissue == "Endolithic band")
res1=ancombc(phyloseq=phy_species_band,formula="Species",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Species",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res1_df=data.frame(Beta=res1[["res"]][["beta"]], Beta_se=res1[["res"]][["se"]], W=res1[["res"]][["W"]],pval=res1[["res"]][["p_val"]], qval=res1[["res"]][["q_val"]],DA=res1[["res"]][["diff_abn"]])
res1_sig=subset(res1_df, DA.SpeciesPorites == "TRUE")[,c(1,2,6,8,10)]
colnames(res1_sig)=c("Beta_Intercept",	"Beta_Species",	"W_Species",	"pval_Species",	"qval_Species")
res1_sig$Diff_more_abundant=ifelse( res1_sig$W_Species < 0 , "Goniastrea", "Porites")
res1_sig$Phylum=meta.mag$Taxa[match(rownames(res1_sig),rownames(tax.mag))]
message("Number of DA KOs: ", nrow(res1_sig), "\nNumber of DA MAGs Goniastrea: ", nrow(subset(res1_sig, W_Species < 0 )), "\nNumber of DA MAGs enriched in Porites: ", nrow(subset(res1_sig, W_Species > 0 )))
write.table(res1_sig,  "outputs/ANCOMBC_mags_betweenSpecies_band.txt", sep = "\t", quote = F, row.names = T )

phy_species_skel=subset_samples(phy_species, Treatment  == "Control" & Tissue == "Skeleton")
res2=ancombc(phyloseq=phy_species_skel,formula="Species",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Species",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res2_df=data.frame(Beta=res2[["res"]][["beta"]], Beta_se=res2[["res"]][["se"]], W=res2[["res"]][["W"]],pval=res2[["res"]][["p_val"]], qval=res2[["res"]][["q_val"]],DA=res2[["res"]][["diff_abn"]])
res2_sig=subset(res2_df, DA.SpeciesPorites == "TRUE")[,c(1,2,6,8,10)]
colnames(res2_sig)=c("Beta_Intercept",	"Beta_Species",	"W_Species",	"pval_Species",	"qval_Species")
res2_sig$Diff_more_abundant=ifelse( res2_sig$W_Species < 0 , "Goniastrea", "Porites")
res2_sig$Phylum=meta.mag$Taxa[match(rownames(res2_sig),rownames(tax.mag))]
message("Number of DA KOs: ", nrow(res2_sig), "\nNumber of DA MAGs Goniastrea: ", nrow(subset(res2_sig, W_Species < 0 )), "\nNumber of DA MAGs enriched in Porites: ", nrow(subset(res2_sig, W_Species > 0 )))
write.table(res2_sig,  "outputs/ANCOMBC_mags_betweenSpecies_skel.txt", sep = "\t", quote = F, row.names = T )

######### temperatures

phy_gon_band=subset_samples(phy_species, Species  == "Goniastrea" & Tissue == "Endolithic band")
res3=ancombc(phyloseq=phy_gon_band,formula="Treatment",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Treatment",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res3_df=data.frame(Beta=res3[["res"]][["beta"]], Beta_se=res3[["res"]][["se"]], W=res3[["res"]][["W"]],pval=res3[["res"]][["p_val"]], qval=res3[["res"]][["q_val"]],DA=res3[["res"]][["diff_abn"]])
res3_sig=subset(res3_df, DA.TreatmentHeat == "TRUE")[,c(1,2,6,8,10)]
colnames(res3_sig)=c("Beta_Intercept",	"Beta_Treatment",	"W_Treatment",	"pval_Treatment",	"qval_Treatment")
res3_sig$Diff_more_abundant=ifelse( res3_sig$W_Treatment < 0 , "Control", "Heat")
res3_sig$Phylum=meta.mag$Taxa[match(rownames(res3_sig),rownames(tax.mag))]
message("Number of DA KOs: ", nrow(res3_sig), "\nNumber of DA MAGs enriched in Control: ", nrow(subset(res3_sig, W_Treatment < 0 )), "\nNumber of DA MAGs enriched in Heat: ", nrow(subset(res3_sig, W_Treatment > 0 )))
write.table(res3_sig,  "outputs/ANCOMBC_mags_betweenTreatment_gon_band.txt", sep = "\t", quote = F, row.names = T )

phy_gon_skel=subset_samples(phy_species, Species  == "Goniastrea" & Tissue == "Skeleton")
res4=ancombc(phyloseq=phy_gon_skel,formula="Treatment",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Treatment",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res4_df=data.frame(Beta=res4[["res"]][["beta"]], Beta_se=res4[["res"]][["se"]], W=res4[["res"]][["W"]],pval=res4[["res"]][["p_val"]], qval=res4[["res"]][["q_val"]],DA=res4[["res"]][["diff_abn"]])
res4_sig=subset(res4_df, DA.TreatmentHeat == "TRUE")[,c(1,2,6,8,10)]
colnames(res4_sig)=c("Beta_Intercept",	"Beta_Treatment",	"W_Treatment",	"pval_Treatment",	"qval_Treatment")
res4_sig$Diff_more_abundant=ifelse( res4_sig$W_Treatment < 0 , "Control", "Heat")
res4_sig$Phylum=meta.mag$Taxam[match(rownames(res4_sig),rownames(tax.mag))]
message("Number of DA KOs: ", nrow(res4_sig), "\nNumber of DA MAGs enriched in Control: ", nrow(subset(res4_sig, W_Treatment < 0 )), "\nNumber of DA MAGs enriched in Heat: ", nrow(subset(res4_sig, W_Treatment > 0 )))
write.table(res4_sig,  "outputs/ANCOMBC_mags_betweenTreatment_gon_skelt", sep = "\t", quote = F, row.names = T )

phy_por_band=subset_samples(phy_species, Species  == "Porites" & Tissue == "Endolithic band")
res5=ancombc(phyloseq=phy_por_band,formula="Treatment",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Treatment",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res5_df=data.frame(Beta=res5[["res"]][["beta"]], Beta_se=res5[["res"]][["se"]], W=res5[["res"]][["W"]],pval=res5[["res"]][["p_val"]], qval=res5[["res"]][["q_val"]],DA=res5[["res"]][["diff_abn"]])
res5_sig=subset(res5_df, DA.TreatmentHeat == "TRUE")[,c(1,2,6,8,10)]
colnames(res5_sig)=c("Beta_Intercept",	"Beta_Treatment",	"W_Treatment",	"pval_Treatment",	"qval_Treatment")
res5_sig$Diff_more_abundant=ifelse( res5_sig$W_Treatment < 0 , "Control", "Heat")
res5_sig$Phylum=meta.mag$Taxa[match(rownames(res5_sig),rownames(tax.mag))]
message("Number of DA KOs: ", nrow(res5_sig), "\nNumber of DA MAGs enriched in Control: ", nrow(subset(res5_sig, W_Treatment < 0 )), "\nNumber of DA MAGs enriched in Heat: ", nrow(subset(res5_sig, W_Treatment > 0 )))
write.table(res5_sig,  "outputs/ANCOMBC_mags_betweenTreatment_por_band.txt", sep = "\t", quote = F, row.names = T )

phy_por_skel=subset_samples(phy_species, Species  == "Porites" & Tissue == "Skeleton")
res6=ancombc(phyloseq=phy_por_skel,formula="Treatment",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Treatment",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F, alpha = 0.05,global = TRUE)
res6_df=data.frame(Beta=res6[["res"]][["beta"]], Beta_se=res6[["res"]][["se"]], W=res6[["res"]][["W"]],pval=res6[["res"]][["p_val"]], qval=res6[["res"]][["q_val"]],DA=res6[["res"]][["diff_abn"]])
res6_sig=subset(res6_df, DA.TreatmentHeat == "TRUE")[,c(1,2,6,8,10)]
colnames(res6_sig)=c("Beta_Intercept",	"Beta_Treatment",	"W_Treatment",	"pval_Treatment",	"qval_Treatment")
res6_sig$Diff_more_abundant=ifelse( res6_sig$W_Treatment < 0 , "Control", "Heat")
res6_sig$Phylum=meta.mag$Taxa[match(rownames(res6_sig),rownames(tax.mag))]
message("Number of DA KOs: ", nrow(res6_sig), "\nNumber of DA MAGs enriched in Control: ", nrow(subset(res6_sig, W_Treatment < 0 )), "\nNumber of DA MAGs enriched in Heat: ", nrow(subset(res6_sig, W_Treatment > 0 )))
write.table(res6_sig,  "outputs/ANCOMBC_mags_betweenTreatment_por_skelt", sep = "\t", quote = F, row.names = T )


meta.mag$DA_band=res1_sig$Diff_more_abundant[match(rownames(meta.mag),rownames(res1_sig) )]
meta.mag$DA_band=ifelse(is.na(meta.mag$DA_band), "None", as.character(meta.mag$DA_band))
message("Number of DA MAGs enriched in Goniastrea: ", nrow(subset(meta.mag, DA_band == "Goniastrea" )), "\nNumber of DA MAGs enriched in Porites: ", nrow(subset(meta.mag, DA_band == "Porites")) , "\nNumber of non-DA MAGs: ", nrow(subset(meta.mag, DA_band == "None")))

meta.mag$DA_skel=res2_sig$Diff_more_abundant[match(rownames(meta.mag),rownames(res2_sig) )]
meta.mag$DA_skel=ifelse(is.na(meta.mag$DA_skel), "None", as.character(meta.mag$DA_skel))
message("Number of DA MAGs enriched in Goniastrea: ", nrow(subset(meta.mag, DA_skel == "Goniastrea" )), "\nNumber of DA MAGs enriched in Porites: ", nrow(subset(meta.mag, DA_skel == "Porites")) , "\nNumber of non-DA MAGs: ", nrow(subset(meta.mag, DA_skel == "None")))

meta.mag$DA=meta.mag$DA_band
meta.mag$DA=ifelse(meta.mag$DA_band == "None", as.character(meta.mag$DA_skel), as.character(meta.mag$DA_band))
write.table(meta.mag, "Input_files/metadata_DA_MAGS.txt", quote = F, sep = "\t")

#################################################################
#### MAGs coverage % between species hierarchical clustering ####
#################################################################
#tax.mag= read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/Input_files/metadata_MAGs_complete.txt", header = T, row.names = 1, sep = "\t")
tax.mag= read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/Input_files/metadata_DA_MAGS.txt", header = T, row.names = 1, sep = "\t")
#meta.mag= read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/Input_files/tree_metadata.txt", header = T, row.names = 1)
map=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt",  sep = "\t", header = T, row.names = 1)
pro=read.table("Input_files/profile_endoliths", fill = T, sep = "\t", header = T, row.names = 1)
bin=pro[, grepl(pattern = "[SB]..mapped.reads", x = colnames(pro))][c(1:75),]
names(bin)=gsub("\\..*", "", names(bin))
#write.table(bin, "Outputs/MAGs_counts.txt")
bin_clr=apply(bin, 2, clr)
#rownames(bin_clr)=paste(rownames(bin_clr),tax.mag$Taxa,sep = " - ")


#both compartments

bin=bin_clr[,grepl("[PG]C[1-6][BS]", colnames(bin_clr))]
bin_sort=bin[,c(1,3,5,7,9,11,2,4,6,8,10,12,13,15,17,19,21,23,14,16,18,20,22,24)]
mag_colors=tax.mag[,c(1,4)]
mag_colors$DA=gsub("Goniastrea|Porites", "Enriched", mag_colors$DA)
rownames(mag_colors)=rownames(bin_sort)
pdf("Outputs/MAG_coverage_noTree_both.pdf",  width = 7, height =7, pointsize = 12)
pheatmap(bin_sort, color = colorRampPalette(c("#1f77b4", "white", "#f94144"))(50),  angle_col = "90",  scale = "row",  cluster_row= T, cluster_col = FALSE, border_color = NA,  cellwidth = 5, cellheight =  5, fontsize_row = 6, fontsize_col= 8, legend = F, gaps_col = c(6,12,18), annotation_row =mag_colors,  annotation_colors =  list(DA= c(Enriched = "black", None="white"), Taxa = c(Acidobacteriota="#2077b4", Bacteroidota="#aec7e8",Chloroflexota="#ff7f0e", Cyanobacteria="#ffbb78", Desulfobacterota="#2ca02c", Firmicutes="#98df89",Gemmatimonadota="#17becf",KSB1= "#9edae5", Myxococcota="#e377c2", Omnitrophota="#f7b6d2", Planctomycetota="#e33e33",  Proteobacteria= "#f1c453", SAR324 = "#efea5a", Spirochaetota = "#16db93" , Verrucomicrobiota= "#B6B6B6"), Species =c(Goniastrea="#B87A89",Porites=  "#66A39F")))
dev.off()
