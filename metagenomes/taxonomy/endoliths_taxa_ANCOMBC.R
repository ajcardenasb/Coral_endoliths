setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/taxonomy/")
library(phyloseq)
library(ANCOMBC)

map=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = T, row.names = 1, sep = "\t")
all=read.table("Input_files/Families_counts_500", header = T, sep = "\t", row.names = 1)
tax=read.table("Input_files/Endoliths_kaiju_taxonomy_500", sep = "\t", header = T, fill = T, quote = "")
tax.u=data.frame(Family=unique(tax$label)) 
tax.u$Domain=tax$Superkingdom[match(tax.u$Family, tax$label)]
rownames(tax.u)=tax.u$Family  

# ordination plot all families
all.fams.p= otu_table(all, taxa_are_rows=TRUE)
sam.p= sample_data(data.frame(map))
tax.p=tax_table(as.matrix(tax.u))
phy.all.fams= phyloseq(all.fams.p,  sam.p, tax.p)


## between species band 
C1=subset_samples(phy.all.fams, Treatment == "Control" & Tissue == "Endolithic band")
res1=ancombc(phyloseq=C1,formula="Species",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Species",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res1_df=data.frame(Beta=res1[["res"]][["beta"]], Beta_se=res1[["res"]][["se"]], W=res1[["res"]][["W"]],pval=res1[["res"]][["p_val"]], qval=res1[["res"]][["q_val"]],DA=res1[["res"]][["diff_abn"]])
res1_sig=subset(res1_df, SpeciesPorites.5 == "TRUE")[,3:5]
colnames(res1_sig)=c("Effect_size_W",	"pval",	"qval")
res1_sig$Diff_more_abundant=ifelse( res1_sig$Effect_size_W < 0 , "Goniastrea", "Porites")
res1_sig$Family=tax.u$Family[match(rownames(res1_sig),rownames(tax.u))]
res1_sig$Domain=tax.u$Domain[match(rownames(res1_sig),rownames(tax.u))]
message("Number of DA families: ", nrow(res1_sig), "\nNumber of DA families enriched in Goniastrea: ", nrow(subset(res1_sig, Diff_more_abundant == "Goniastrea" )), "\nNumber of DA families enriched in Porites: ", nrow(subset(res1_sig, Diff_more_abundant == "Porites" )))
#write.table(res1_sig, "outputs/ANCOM-BC_Families_species_band", quote = F,  sep = "\t")

skel_phy=subset_samples(phy.all.fams, Treatment == "Control" & Tissue == "Skeleton")
res2=ancombc(phyloseq=skel_phy,formula="Species",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Species",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res2_df=data.frame(Beta=res2[["res"]][["beta"]], Beta_se=res2[["res"]][["se"]], W=res2[["res"]][["W"]],pval=res2[["res"]][["p_val"]], qval=res2[["res"]][["q_val"]],DA=res2[["res"]][["diff_abn"]])
res2_sig=subset(res2_df, SpeciesPorites.5 == "TRUE")[,3:5]
colnames(res2_sig)=c("Effect_size_W",	"pval",	"qval")
res2_sig$Diff_more_abundant=ifelse( res2_sig$Effect_size_W < 0 , "Goniastrea", "Porites")
res2_sig$Family=tax.u$Family[match(rownames(res2_sig),rownames(tax.u))]
res2_sig$Domain=tax.u$Domain[match(rownames(res2_sig),rownames(tax.u))]
message("Number of DA families: ", nrow(res2_sig), "\nNumber of DA families enriched in Goniastrea: ", nrow(subset(res2_sig, Diff_more_abundant == "Goniastrea" )), "\nNumber of DA families enriched in Porites: ", nrow(subset(res2_sig, Diff_more_abundant == "Porites" )))
#write.table(res2_sig, "outputs/ANCOM-BC_Families_species_skel", quote = F, sep = "\t")


#########################
## between treatments  ##
#########################
gon_band_phy=subset_samples(phy.all.fams, Species == "Goniastrea" & Tissue == "Endolithic band")
res3=ancombc(phyloseq=gon_band_phy,formula="Treatment",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Treatment",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res3_df=data.frame(Beta=res3[["res"]][["beta"]], Beta_se=res3[["res"]][["se"]], W=res3[["res"]][["W"]],pval=res3[["res"]][["p_val"]], qval=res3[["res"]][["q_val"]],DA=res3[["res"]][["diff_abn"]])
res3_sig=subset(res3_df, TreatmentHeat.5 == "TRUE")[,3:5]
colnames(res3_sig)=c("Effect_size_W",	"pval",	"qval")
res3_sig$Diff_more_abundant=ifelse( res3_sig$Effect_size_W < 0 , "Control", "Heat")
res3_sig$Family=tax.u$Family[match(rownames(res3_sig),rownames(tax.u))]
res3_sig$Domain=tax.u$Domain[match(rownames(res3_sig),rownames(tax.u))]
message("Number of DA families: ", nrow(res3_sig), "\nNumber of DA families enriched in control: ", nrow(subset(res3_sig, Diff_more_abundant == "Control" )), "\nNumber of DA families enriched in heat: ", nrow(subset(res3_sig, Diff_more_abundant == "Heat" )))
#write.table(res3_sig, "outputs/ANCOM-BC_Families_treatment_Gon_band", quote = F,  sep = "\t")

gon_skel_phy=subset_samples(phy.all.fams, Species == "Goniastrea" & Tissue == "Skeleton")
res4=ancombc(phyloseq=gon_skel_phy,formula="Treatment",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Treatment",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res4_df=data.frame(Beta=res4[["res"]][["beta"]], Beta_se=res4[["res"]][["se"]], W=res4[["res"]][["W"]],pval=res4[["res"]][["p_val"]], qval=res4[["res"]][["q_val"]],DA=res4[["res"]][["diff_abn"]])
res4_sig=subset(res4_df, TreatmentHeat.5 == "TRUE")[,3:5]
colnames(res4_sig)=c("Effect_size_W",	"pval",	"qval")
res4_sig$Diff_more_abundant=ifelse( res4_sig$Effect_size_W < 0 , "Control", "Heat")
res4_sig$Family=tax.u$Family[match(rownames(res4_sig),rownames(tax.u))]
res4_sig$Domain=tax.u$Domain[match(rownames(res4_sig),rownames(tax.u))]
message("Number of DA families: ", nrow(res4_sig), "\nNumber of DA families enriched in control: ", nrow(subset(res4_sig, Diff_more_abundant == "Control" )), "\nNumber of DA families enriched in heat: ", nrow(subset(res4_sig, Diff_more_abundant == "Heat" )))
#write.table(res4_sig, "outputs/ANCOM-BC_Families_treatment_Gon_skel", quote = F,  sep = "\t")

por_band_phy=subset_samples(phy.all.fams, Species == "Porites" & Tissue == "Endolithic band")
res5=ancombc(phyloseq=por_band_phy,formula="Treatment",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Treatment",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res5_df=data.frame(Beta=res5[["res"]][["beta"]], Beta_se=res5[["res"]][["se"]], W=res5[["res"]][["W"]],pval=res5[["res"]][["p_val"]], qval=res5[["res"]][["q_val"]],DA=res5[["res"]][["diff_abn"]])
res5_sig=subset(res5_df, TreatmentHeat.5 == "TRUE")[,3:5]
colnames(res5_sig)=c("Effect_size_W",	"pval",	"qval")
res5_sig$Diff_more_abundant=ifelse( res5_sig$Effect_size_W < 0 , "Control", "Heat")
res5_sig$Family=tax.u$Family[match(rownames(res5_sig),rownames(tax.u))]
res5_sig$Domain=tax.u$Domain[match(rownames(res5_sig),rownames(tax.u))]
message("Number of DA families: ", nrow(res5_sig), "\nNumber of DA families enriched in control: ", nrow(subset(res5_sig, Diff_more_abundant == "Control" )), "\nNumber of DA families enriched in heat: ", nrow(subset(res5_sig, Diff_more_abundant == "Heat" )))
#write.table(res5_sig, "outputs/ANCOM-BC_Families_treatment_Por_band", quote = F,  sep = "\t")

por_skel_phy=subset_samples(phy.all.fams, Species == "Porites" & Tissue == "Skeleton")
res6=ancombc(phyloseq=por_skel_phy,formula="Treatment",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Treatment",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res6_df=data.frame(Beta=res6[["res"]][["beta"]], Beta_se=res6[["res"]][["se"]], W=res6[["res"]][["W"]],pval=res6[["res"]][["p_val"]], qval=res6[["res"]][["q_val"]],DA=res6[["res"]][["diff_abn"]])
res6_sig=subset(res6_df, TreatmentHeat.5 == "TRUE")[,3:5]
colnames(res6_sig)=c("Effect_size_W",	"pval",	"qval")
res6_sig$Diff_more_abundant=ifelse( res6_sig$Effect_size_W < 0 , "Control", "Heat")
res6_sig$Family=tax.u$Family[match(rownames(res6_sig),rownames(tax.u))]
res6_sig$Domain=tax.u$Domain[match(rownames(res6_sig),rownames(tax.u))]
message("Number of DA families: ", nrow(res6_sig), "\nNumber of DA families enriched in control: ", nrow(subset(res6_sig, Diff_more_abundant == "Control" )), "\nNumber of DA families enriched in heat: ", nrow(subset(res6_sig, Diff_more_abundant == "Heat" )))
#write.table(res6_sig, "outputs/ANCOM-BC_Families_treatment_Por_skel", quote = F,  sep = "\t")


#outputs for volcano plots

vol_gon_band=res3_df[,c(1,3,4,5)]
colnames(vol_gon_band)=c("Beta", "W", "pval",	"qval")
vol_gon_band$Superkingdom=tax.u$Domain[match(rownames(vol_gon_band),rownames(tax.u))]
vol_gon_band$Superkingdom=ifelse(vol_gon_band$qval > 0.05, "noDA", as.character(vol_gon_band$Superkingdom))
#write.table(vol_gon_band, "outputs/volcano_treatment_Gon_band", quote = F,  sep = "\t")

vol_gon_skel=res4_df[,c(1,3,4,5)]
colnames(vol_gon_skel)=c("Beta", "W", "pval",	"qval")
vol_gon_skel$Superkingdom=tax.u$Domain[match(rownames(vol_gon_skel),rownames(tax.u))]
vol_gon_skel$Superkingdom=ifelse(vol_gon_skel$qval > 0.05, "noDA", as.character(vol_gon_skel$Superkingdom))
#write.table(vol_gon_skel, "outputs/volcano_treatment_Gon_skel", quote = F,  sep = "\t")

vol_por_band=res5_df[,c(1,3,4,5)]
colnames(vol_por_band)=c("Beta", "W", "pval",	"qval")
vol_por_band$Superkingdom=tax.u$Domain[match(rownames(vol_por_band),rownames(tax.u))]
vol_por_band$Superkingdom=ifelse(vol_por_band$qval > 0.05, "noDA", as.character(vol_por_band$Superkingdom))
#write.table(vol_por_band, "outputs/volcano_treatment_Por_band", quote = F,  sep = "\t")

vol_por_skel=res6_df[,c(1,3,4,5)]
colnames(vol_por_skel)=c("Beta", "W", "pval",	"qval")
vol_por_skel$Superkingdom=tax.u$Domain[match(rownames(vol_por_skel),rownames(tax.u))]
vol_por_skel$Superkingdom=ifelse(vol_por_skel$qval > 0.05, "noDA", as.character(vol_por_skel$Superkingdom))
#write.table(vol_por_skel, "outputs/volcano_treatment_Por_skel", quote = F,  sep = "\t")

