#library(selbal)
setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/functional/")
source("~/Documents/Bioinformatics_scripts/R_scripts/orphan/remove_rare.R")
library("selbal")

ancom_band=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/functional/outputs/ANCOMBC_L3_betweenSpecies_band.txt", header = T, sep = "\t", row.names = 1)
ancom_skel=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/functional/outputs/ANCOMBC_L3_betweenSpecies_skel.txt", header = T, sep = "\t", row.names = 1)

met=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = T, sep = "\t", row.names = 1)
kegg=read.table("~/Documents/Bioinformatics_scripts/KEGG_files/MostUpdated_hierarchy_ko00001", sep = "\t", header = T, fill = T)#, row.names = 5) 
kegg$module=trimws(kegg$module, which = "both", whitespace = "[ \t\r\n]")
meta_l3=unique(subset(kegg, L1 == " Metabolism" & L2 %in% c(" Energy metabolism" ," Carbohydrate metabolism"))$L3)
#meta_l3=unique(subset(kegg, L1 == " Metabolism" )$L3)
cnts=read.table("Input_files/KO_counts_500", header = T, row.names = 1, sep = "\t")

cnts$MD=kegg$L3[match(rownames(cnts), kegg$KO)]
cnts$MD=ifelse(is.na(cnts$MD), "Others", as.character(cnts$MD))
cnts_md=aggregate(cnts[,1:48], by=list(cnts$MD), FUN = sum)
rownames(cnts_md)=cnts_md[,1]
cnts_md=cnts_md[,-1]
#cnts_met_l3_pre=subset(cnts_md, rownames(cnts_md) %in% meta_l3)
#cnts_met_l3=as.data.frame(t(remove_rare(table=cnts_met_l3_pre, cutoff_pro=0.1))) removing processes that are not in at least 10%
cnts_met_l3=as.data.frame(t(subset(cnts_md, rownames(cnts_md) %in% rownames(ancom_band) & rownames(cnts_md)  %in% meta_l3)))
n_col=ncol(cnts_met_l3)
cnts_met_l3$Species=met$Species[match(rownames(cnts_met_l3),rownames(met))]
cnts_met_l3$Tissue=met$Tissue[match(rownames(cnts_met_l3),rownames(met))]
cnts_met_l3$Treatment=met$Treatment[match(rownames(cnts_met_l3),rownames(met))]
cnts_met_l3$Genotype=met$Genotype[match(rownames(cnts_met_l3),rownames(met))]

band_cnt=subset(cnts_met_l3, Treatment == "Control" & Tissue == "Endolithic band") # band
band_cnt=subset(cnts_met_l3, Treatment == "Control" & Tissue == "Skeleton") # skel
band_cnt$Species=as.factor(band_cnt$Species)
sebal_Res = selbal.cv(as.matrix(band_cnt[,1:n_col]), band_cnt$Species, zero.rep = "one")
sebal_Res$accuracy.nvar # the smallest number of variables with an accuracy within the minimum accuray plus one standard error.
sebal_Res$var.barplot # the frequency of the variables selected in some step of the CV process
grid.draw(sebal_Res$global.plot)
plot.tab(sebal_Res$cv.tab) #  summary of the CV-procedure.variables included either in the Global balance (second column) or the three most frequent balances in the CV process (last three columns). The first column provides the percentage of times each variable has appeared in CV and, the last row points the proportion of times the most repeated balances have appeared.
summary(sebal_Res$cv.accuracy) #test prediction or classification accuracy (AUC in this case) for each cv iteration (length equal to n.fold*n.iter)
sebal_Res$global.balance#he variables included in the Global Balance. The first column, Taxa, provides the names of the selected taxa while variable Group specifies wether the taxon is included in the numerator (NUM) or the denominator (DEN) of the balance.
sebal_Res$glm# result of a regression model. Balance selection is implemented through a regression model where the balance itself and additional covariates are the explanatory variables, and the variable of interest the response. The user has access to the regression model in order to obtain and analyze any of its characteristics.
summary(sebal_Res$glm)
sebal_Res$opt.nvar# the last element represents the optinal number of variables considered in the balance as the optimal.

#plots

#species
band_cnt=subset(cnts_met_l3, Treatment == "Control" & Tissue == "Endolithic band") # band
skel_cnt=subset(cnts_met_l3, Treatment == "Control" & Tissue == "Skeleton") # skel
band_cnt$Species=as.factor(band_cnt$Species)
skel_cnt$Species=as.factor(skel_cnt$Species)
sebal_band = selbal.cv(as.matrix(band_cnt[,1:n_col]), band_cnt$Species, zero.rep = "one")
sebal_skel = selbal.cv(as.matrix(skel_cnt[,1:n_col]), skel_cnt$Species, zero.rep = "one")

pdf("./outputs/selbal_output_band.pdf",  width = 12, height = 7, pointsize = 10) 
grid.draw(sebal_band$global.plot)
dev.off()

pdf("./outputs/selbal_output_skel.pdf",  width = 12, height = 7, pointsize = 10) 
skel=grid.draw(sebal_skel$global.plot)
dev.off()
