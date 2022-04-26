library(ANCOMBC)
library(phyloseq)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/functional/")

map=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = T, sep = "\t", row.names = 1)
map$Species=as.factor(map$Species)
map$Species = relevel(map$Species, ref="Goniastrea")
map$Treatment=as.factor(map$Treatment)
kegg=read.table("~/Documents/Bioinformatics_scripts/KEGG_files/MostUpdated_hierarchy_ko00001", sep = "\t", header = T, fill = T)#, row.names = 5) 
kegg$module=trimws(kegg$module, which = "both", whitespace = "[ \t\r\n]")
cnts=read.table("Input_files/KO_counts_500", header = T, row.names = 1, sep = "\t")

cnts$MD=kegg$L3[match(rownames(cnts), kegg$KO)]
cnts$MD=ifelse(is.na(cnts$MD), "Others", as.character(cnts$MD))
cnts_md=aggregate(cnts[,1:48], by=list(cnts$MD), FUN = sum)
rownames(cnts_md)=cnts_md[,1]
cnts_md=cnts_md[,-1]

otu.t= otu_table(cnts_md, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
phy_species= phyloseq(otu.t,  sam.t )#, tax.t)


######################
## between species ##
######################

phy_species_band=subset_samples(phy_species,Treatment  == "Control" & Tissue == "Endolithic band")
res1=ancombc(phyloseq=phy_species_band,formula="Species",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Species",struc_zero =T,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res1_df=data.frame(Beta=res1[["res"]][["beta"]], Beta_se=res1[["res"]][["se"]], W=res1[["res"]][["W"]],pval=res1[["res"]][["p_val"]], qval=res1[["res"]][["q_val"]],DA=res1[["res"]][["diff_abn"]])
res1_sig=subset(res1_df, DA.SpeciesPorites == "TRUE")[,c(2,6,8,10)]
colnames(res1_sig)=c("Beta_Species",	"W_Species",	"pval_Species",	"qval_Species")
res1_sig$Diff_more_abundant=ifelse( res1_sig$W_Species < 0 , "Goniastrea", "Porites")
res1_sig$L2=kegg$L2[match(rownames(res1_sig), kegg$L3)] #res1_sig$L2=kegg$L2[match(rownames(res1_sig), kegg$L3)] 
res1_sig$L1=kegg$L1[match(rownames(res1_sig), kegg$L3)] #res1_sig$L2=kegg$L2[match(rownames(res1_sig), kegg$L3)] 
message("Number of DA  L3: ", nrow(res1_sig), "\nNumber of DA L3 enriched in Goniastrea: ", nrow(subset(res1_sig, W_Species < 0 )), "\nNumber of DA L3 enriched in Porites: ", nrow(subset(res1_sig, W_Species > 0 )))
message("Number of metabolic DA L3: ", nrow(subset(res1_sig, L1 == " Metabolism")), "\nNumber of DA fL3 enriched in Goniastrea: ", nrow(subset(res1_sig, L1 == " Metabolism" &  W_Species < 0 )), "\nNumber of DA L3 enriched in Porites: ", nrow(subset(res1_sig,  L1 == " Metabolism" & W_Species > 0 )))
write.table(res1_sig,  "outputs/ANCOMBC_L3_betweenSpecies_band.txt", sep = "\t", quote = F, row.names = T )

phy_species_skel=subset_samples(phy_species,Treatment  == "Control" & Tissue == "Skeleton")
res2=ancombc(phyloseq=phy_species_skel,formula="Species",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Species",struc_zero =T,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res2_df=data.frame(Beta=res2[["res"]][["beta"]], Beta_se=res2[["res"]][["se"]], W=res2[["res"]][["W"]],pval=res2[["res"]][["p_val"]], qval=res2[["res"]][["q_val"]],DA=res2[["res"]][["diff_abn"]])
res2_sig=subset(res2_df, DA.SpeciesPorites == "TRUE")[,c(1,2,6,8,10)]
colnames(res2_sig)=c("Beta_Intercept",	"Beta_Species",	"W_Species",	"pval_Species",	"qval_Species")
res2_sig$Diff_more_abundant=ifelse( res2_sig$W_Species < 0 , "Goniastrea", "Porites")
res2_sig$L2=kegg$L2[match(rownames(res2_sig), kegg$L3)]
res2_sig$L1=kegg$L1[match(rownames(res2_sig), kegg$L3)]
message("Number of DA  L3: ", nrow(res2_sig), "\nNumber of DA L3 enriched in Goniastrea: ", nrow(subset(res2_sig, W_Species < 0 )), "\nNumber of DA L3 enriched in Porites: ", nrow(subset(res2_sig, W_Species > 0 )))
message("Number of metabolic DA L3: ", nrow(subset(res2_sig, L1 == " Metabolism")), "\nNumber of DA fL3 enriched in Goniastrea: ", nrow(subset(res2_sig, L1 == " Metabolism" &  W_Species < 0 )), "\nNumber of DA L3 enriched in Porites: ", nrow(subset(res2_sig,  L1 == " Metabolism" & W_Species > 0 )))
write.table(res2_sig,  "outputs/ANCOMBC_L3_betweenSpecies_skel.txt", sep = "\t", quote = F, row.names = T )


#########################
## between treatments  ##
#########################

## Gon band
phy_treat_gon_band=subset_samples(phy_species, Species  == "Goniastrea" & Tissue == "Endolithic band")
res3=ancombc(phyloseq=phy_treat_gon_band,formula="Treatment",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Treatment",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res3_df=data.frame(Beta=res3[["res"]][["beta"]], Beta_se=res3[["res"]][["se"]], W=res3[["res"]][["W"]],pval=res3[["res"]][["p_val"]], qval=res3[["res"]][["q_val"]],DA=res3[["res"]][["diff_abn"]])
res3_sig=subset(res3_df, DA.TreatmentHeat == "TRUE")[,c(1,2,6,8,10)]
colnames(res3_sig)=c("Beta_Intercept",	"Beta_Treatment",	"W_Treatment",	"pval_Treatment",	"qval_Treatment")
res3_sig$Diff_more_abundant=ifelse( res3_sig$W_Treatment < 0 , "Control", "Heat")
res3_sig$L2=kegg$L2[match(rownames(res3_sig), kegg$L3)]
res3_sig$L1=kegg$L1[match(rownames(res3_sig), kegg$L3)]
message("Number of DA KOs: ", nrow(res3_sig), "\nNumber of DA families enriched in Control: ", nrow(subset(res3_sig, W_Treatment < 0 )), "\nNumber of DA families enriched in Heat: ", nrow(subset(res3_sig, W_Treatment > 0 )))
message("Number of metabolic DA L3: ", nrow(subset(res3_sig, L1 == " Metabolism")), "\nNumber of DA L3 enriched in Goniastrea: ", nrow(subset(res3_sig, L1 == " Metabolism" &  W_Treatment < 0 )), "\nNumber of DA L3 enriched in Porites: ", nrow(subset(res3_sig,  L1 == " Metabolism" & W_Treatment > 0 )))
write.table(res3_sig,  "outputs/ANCOMBC_L3_betweenTreatments_gon_band.txt", sep = "\t", quote = F, row.names = T )


## Gon skel
phy_treat_gon_skel=subset_samples(phy_species, Species  == "Goniastrea" & Tissue == "Skeleton")
res4=ancombc(phyloseq=phy_treat_gon_skel,formula="Treatment",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Treatment",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res4_df=data.frame(Beta=res4[["res"]][["beta"]], Beta_se=res4[["res"]][["se"]], W=res4[["res"]][["W"]],pval=res4[["res"]][["p_val"]], qval=res4[["res"]][["q_val"]],DA=res4[["res"]][["diff_abn"]])
res4_sig=subset(res4_df, DA.TreatmentHeat == "TRUE")[,c(1,2,6,8,10)]
colnames(res4_sig)=c("Beta_Intercept",	"Beta_Treatment",	"W_Treatment",	"pval_Treatment",	"qval_Treatment")
res4_sig$Diff_more_abundant=ifelse( res4_sig$W_Treatment < 0 , "Control", "Heat")
res4_sig$L2=kegg$L2[match(rownames(res4_sig), kegg$L3)]
res4_sig$L1=kegg$L1[match(rownames(res4_sig), kegg$L3)]
message("Number of DA KOs: ", nrow(res4_sig), "\nNumber of DA families enriched in Control: ", nrow(subset(res4_sig, W_Treatment < 0 )), "\nNumber of DA families enriched in Heat: ", nrow(subset(res4_sig, W_Treatment > 0 )))
message("Number of metabolic DA L3: ", nrow(subset(res4_sig, L1 == " Metabolism")), "\nNumber of DA L3 enriched in Goniastrea: ", nrow(subset(res4_sig, L1 == " Metabolism" &  W_Treatment < 0 )), "\nNumber of DA L3 enriched in Porites: ", nrow(subset(res4_sig,  L1 == " Metabolism" & W_Treatment > 0 )))
write.table(res4_sig,  "outputs/ANCOMBC_L3_betweenTreatments_gon_skel.txt", sep = "\t", quote = F, row.names = T )

## Por band
phy_treat_por_band=subset_samples(phy_species, Species  == "Porites" & Tissue == "Endolithic band")
res5=ancombc(phyloseq=phy_treat_por_band,formula="Treatment",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Treatment",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res5_df=data.frame(Beta=res5[["res"]][["beta"]], Beta_se=res5[["res"]][["se"]], W=res5[["res"]][["W"]],pval=res5[["res"]][["p_val"]], qval=res5[["res"]][["q_val"]],DA=res5[["res"]][["diff_abn"]])
res5_sig=subset(res5_df, DA.TreatmentHeat == "TRUE")[,c(1,2,6,8,10)]
colnames(res5_sig)=c("Beta_Intercept",	"Beta_Treatment",	"W_Treatment",	"pval_Treatment",	"qval_Treatment")
res5_sig$Diff_more_abundant=ifelse( res5_sig$W_Treatment < 0 , "Control", "Heat")
res5_sig$L2=kegg$L2[match(rownames(res5_sig), kegg$L3)]
res5_sig$L1=kegg$L1[match(rownames(res5_sig), kegg$L3)]
message("Number of DA KOs: ", nrow(res5_sig), "\nNumber of DA families enriched in Control: ", nrow(subset(res5_sig, W_Treatment < 0 )), "\nNumber of DA families enriched in Heat: ", nrow(subset(res5_sig, W_Treatment > 0 )))
message("Number of metabolic DA L3: ", nrow(subset(res5_sig, L1 == " Metabolism")), "\nNumber of DA L3 enriched in Goniastrea: ", nrow(subset(res5_sig, L1 == " Metabolism" &  W_Treatment < 0 )), "\nNumber of DA L3 enriched in Porites: ", nrow(subset(res5_sig,  L1 == " Metabolism" & W_Treatment > 0 )))
write.table(res5_sig,  "outputs/ANCOMBC_L3_betweenTreatments_por_band.txt", sep = "\t", quote = F, row.names = T )


## Por skel
phy_treat_gon_skel=subset_samples(phy_species, Species  == "Porites" & Tissue == "Skeleton")
res6=ancombc(phyloseq=phy_treat_gon_skel,formula="Treatment",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Treatment",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res6_df=data.frame(Beta=res6[["res"]][["beta"]], Beta_se=res6[["res"]][["se"]], W=res6[["res"]][["W"]],pval=res6[["res"]][["p_val"]], qval=res6[["res"]][["q_val"]],DA=res6[["res"]][["diff_abn"]])
res6_sig=subset(res6_df, DA.TreatmentHeat == "TRUE")[,c(1,2,6,8,10)]
colnames(res6_sig)=c("Beta_Intercept",	"Beta_Treatment",	"W_Treatment",	"pval_Treatment",	"qval_Treatment")
res6_sig$Diff_more_abundant=ifelse( res6_sig$W_Treatment < 0 , "Control", "Heat")
res6_sig$L2=kegg$L2[match(rownames(res6_sig), kegg$L3)]
res6_sig$L1=kegg$L1[match(rownames(res6_sig), kegg$L3)]
message("Number of DA KOs: ", nrow(res6_sig), "\nNumber of DA families enriched in Control: ", nrow(subset(res6_sig, W_Treatment < 0 )), "\nNumber of DA families enriched in Heat: ", nrow(subset(res6_sig, W_Treatment > 0 )))
message("Number of metabolic DA L3: ", nrow(subset(res6_sig, L1 == " Metabolism")), "\nNumber of DA L3 enriched in Goniastrea: ", nrow(subset(res6_sig, L1 == " Metabolism" &  W_Treatment < 0 )), "\nNumber of DA L3 enriched in Porites: ", nrow(subset(res6_sig,  L1 == " Metabolism" & W_Treatment > 0 )))
write.table(res6_sig,  "outputs/ANCOMBC_L3_betweenTreatments_por_skel.txt", sep = "\t", quote = F, row.names = T )

