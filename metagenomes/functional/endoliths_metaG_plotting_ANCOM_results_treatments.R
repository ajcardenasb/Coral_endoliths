setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/functional/")

library(compositions)
library(ggplot2)

map=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = T, sep = "\t")
cnts=read.table("Input_files/metabolic_KO_counts_500", header = T, sep = "\t", row.names = 1) # for KOs
cnts=round(cnts, digits = 0)

#######################
## between treatments #
#######################

##comparison 1: goniastrea endolithic band
cnts_gon_band=apply(cnts[,grepl("G[CH][1-6]B", colnames(cnts))],2,clr)
colnames(cnts_gon_band)=gsub("[0-9]", "", colnames(cnts_gon_band))
means_gon_band=as.data.frame(t(apply(cnts_gon_band, 1, function(x) tapply(x, colnames(cnts_gon_band), mean))))
means_gon_band$enrichment=ifelse(means_gon_band$GCB > means_gon_band$GHB, "Control", "Heat")
gon_band=read.table("outputs/ANCOM_results/metabolims_KO_counts_ANCOM_betweenTretat_Gon_band.txt", sep = "\t", header = T)
gon_band_sig=subset(gon_band, detected_0.6 == "TRUE")
gon_band_sig$Enrichment=means_gon_band$enrichment[match(gon_band_sig$taxa_id, rownames(means_gon_band))]
message("Number of enriched features in control: ", nrow(subset(gon_band_sig, Enrichment == "Control")), "\nNumber of enriched features in Heat: ", nrow(subset(gon_band_sig, Enrichment == "Heat")))
#write.table(gon_band_sig, "outputs/ANCOM_results/Final_ANCOM_betweenTreatment_Gonband.txt",  sep = "\t", quote = F, row.names = F )

##comparison 2: goniastrea skeleton
cnts_gon_skel=apply(cnts[,grepl("G[CH][1-6]S", colnames(cnts))],2,clr)
colnames(cnts_gon_skel)=gsub("[0-9]", "", colnames(cnts_gon_skel))
means_gon_skel=as.data.frame(t(apply(cnts_gon_skel, 1, function(x) tapply(x, colnames(cnts_gon_skel), mean))))
means_gon_skel$enrichment=ifelse(means_gon_skel$GCS > means_gon_skel$GHS, "Control", "Heat")
gon_skel=read.table("outputs/ANCOM_results/metabolims_KO_counts_ANCOM_betweenTretat_Gon_skel.txt", sep = "\t", header = T)
gon_skel_sig=subset(gon_skel, detected_0.6 == "TRUE")
gon_skel_sig$Enrichment=means_gon_skel$enrichment[match(gon_skel_sig$taxa_id, rownames(means_gon_skel))]
message("Number of enriched features in control: ", nrow(subset(gon_skel_sig, Enrichment == "Control")), "\nNumber of enriched features in Heat: ", nrow(subset(gon_skel_sig, Enrichment == "Heat")))
#write.table(gon_skel_sig, "outputs/ANCOM_results/Final_ANCOM_betweenTreatment_Gonskel.txt",  sep = "\t", quote = F, row.names = F )

##comparison 3: Porites band
cnts_por_band=apply(cnts[,grepl("P[CH][1-6]B", colnames(cnts))],2,clr)
colnames(cnts_por_band)=gsub("[0-9]", "", colnames(cnts_por_band))
means_por_band=as.data.frame(t(apply(cnts_por_band, 1, function(x) tapply(x, colnames(cnts_por_band), mean))))
means_por_band$enrichment=ifelse(means_por_band$PCB > means_por_band$PHB, "Control", "Heat")
por_band=read.table("outputs/ANCOM_results/metabolims_KO_counts_ANCOM_betweenTretat_Por_band.txt", sep = "\t", header = T)
por_band_sig=subset(por_band, detected_0.6 == "TRUE")
por_band_sig$Enrichment=means_por_band$enrichment[match(por_band_sig$taxa_id, rownames(means_por_band))]
message("Number of enriched features in control: ", nrow(subset(por_band_sig, Enrichment == "Control")), "\nNumber of enriched features in Heat: ", nrow(subset(por_band_sig, Enrichment == "Heat")))
#write.table(por_band_sig, "outputs/ANCOM_results/Final_ANCOM_betweenTreatment_Gonskel.txt",  sep = "\t", quote = F, row.names = F )

##comparison 3: Porites skel
cnts_por_skel=apply(cnts[,grepl("P[CH][1-6]S", colnames(cnts))],2,clr)
colnames(cnts_por_skel)=gsub("[0-9]", "", colnames(cnts_por_skel))
means_por_skel=as.data.frame(t(apply(cnts_por_skel, 1, function(x) tapply(x, colnames(cnts_por_skel), mean))))
means_por_skel$enrichment=ifelse(means_por_skel$PCS > means_por_skel$PHS, "Control", "Heat")
por_skel=read.table("outputs/ANCOM_results/metabolims_KO_counts_ANCOM_betweenTretat_Por_skel.txt", sep = "\t", header = T)
por_skel_sig=subset(por_skel, detected_0.6 == "TRUE")
por_skel_sig$Enrichment=means_por_skel$enrichment[match(por_skel_sig$taxa_id, rownames(means_por_skel))]
message("Number of enriched features in control: ", nrow(subset(por_skel_sig, Enrichment == "Control")), "\nNumber of enriched features in Heat: ", nrow(subset(por_skel_sig, Enrichment == "Heat")))
#write.table(por_skel_sig, "outputs/ANCOM_results/Final_ANCOM_betweenTreatment_Gonskel.txt",  sep = "\t", quote = F, row.names = F )



#######################################################
### Between species plot 1 : absolute abundance L2 ###
#######################################################
kegg=read.table("~/Documents/Bioinformatics_scripts/KEGG_files/hierarchy_ko00001", sep = "\t", header = T, fill = T, quote = "")

por_skel_sig$comparison="Porites skeleton"
por_band_sig$comparison="Porites band"
gon_skel_sig$comparison="Goniastrea skeleton"
gon_band_sig$comparison="Goniastrea band"
all=rbind(por_skel_sig,por_band_sig,gon_skel_sig,gon_band_sig)
all$L2=kegg$L2[match(all$taxa_id, kegg$KO)]
all$L3=kegg$L3[match(all$taxa_id, kegg$KO)]

summ_l2= all %>% group_by(Enrichment, L2, comparison) %>%  tally()
ggplot(summ_l2, aes(x = L2, y=n, fill=Enrichment)) + geom_bar( stat="identity", alpha =0.8, position = "dodge") + scale_fill_manual(values = c("#247BA0", "#ffe066"))   + labs(x="KEGG L2 category",y="Number of DA KOs") + theme_classic() + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1))+ coord_flip() + facet_grid(~comparison)

top_categories = all %>% filter(Enrichment == "Heat") %>% group_by(Enrichment, L2, ) %>%  tally()
top_categories2 = all %>% group_by(Enrichment, L3, comparison) %>%  tally()

tmp =all %>% filter(Enrichment == "Heat") %>% group_by( L2, L3, comparison) %>%  tally()
out2=reshape2::dcast(tmp, L2+L3~comparison,  value.var= "n",  fun.aggregate = sum)

tmp =all %>% filter(Enrichment == "Heat") %>% group_by(L2, comparison ) %>%  tally()
out1=reshape2::dcast(tmp, L2~comparison,  value.var= "n",  fun.aggregate = sum)
write.table(out1, "outputs/L2_betwTreatments_DAKOssummary.txt", sep = "\t", quote = F, row.names = F)

#####################################################################
### Between treatments plot 2 : absolute abundance L2 procesess ###
#####################################################################

l3_carbo=all %>% filter(L2 == " Carbohydrate metabolism") %>% group_by(Enrichment, L3, comparison) %>%  tally()
a=ggplot(l3_carbo, aes(x = L3, y=n, fill=Enrichment)) + geom_bar( stat="identity", alpha =0.8, position = "dodge") + scale_fill_manual(values = c("#247BA0", "#ffe066"))   + labs(x="Carbohydrate metabolism",y="Number of DA KOs") + theme_classic() + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1))+ coord_flip() + facet_grid(~comparison)

l3_terpe=all %>% filter(L2 == " Metabolism of terpenoids and polyketides") %>% group_by(Enrichment, L3, comparison) %>%  tally()
b=ggplot(l3_terpe, aes(x = L3, y=n, fill=Enrichment)) + geom_bar( stat="identity", alpha =0.8, position = "dodge") + scale_fill_manual(values = c("#247BA0", "#ffe066"))   + labs(x="Terpenoids and polyketides metabolism",y="Number of DA KOs") + theme_classic() + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1))+ coord_flip() + facet_grid(~comparison)

l3_amino=all %>% filter(L2 == " Amino acid metabolism") %>% group_by(Enrichment, L3, comparison) %>%  tally()
c=ggplot(l3_amino, aes(x = L3, y=n, fill=Enrichment)) + geom_bar( stat="identity", alpha =0.8, position = "dodge") + scale_fill_manual(values = c("#247BA0", "#ffe066"))   + labs(x="Amino acid metabolism",y="Number of DA KOs") + theme_classic() + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1))+ coord_flip() + facet_grid(~comparison)

l3_vita=all %>% filter(L2 == " Metabolism of cofactors and vitamins") %>% group_by(Enrichment, L3, comparison) %>%  tally()
d=ggplot(l3_vita, aes(x = L3, y=n, fill=Enrichment)) + geom_bar( stat="identity", alpha =0.8, position = "dodge") + scale_fill_manual(values = c("#247BA0", "#ffe066"))   + labs(x="Metabolism of cofactors and vitamins",y="Number of DA KOs") + theme_classic() + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1))+ coord_flip() + facet_grid(~comparison)

l3_ener=all %>% filter(L2 == " Energy metabolism") %>% group_by(Enrichment, L3, comparison) %>%  tally()
e=ggplot(l3_ener, aes(x = L3, y=n, fill=Enrichment)) + geom_bar( stat="identity", alpha =0.8, position = "dodge") + scale_fill_manual(values = c("#247BA0", "#ffe066"))   + labs(x="Energy metabolism",y="Number of DA KOs") + theme_classic() + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1))+ coord_flip() + facet_grid(~comparison)

l3_metab=all %>% filter(L2 == " Biosynthesis of other secondary metabolites") %>% group_by(Enrichment, L3, comparison) %>%  tally()
f=ggplot(l3_metab, aes(x = L3, y=n, fill=Enrichment)) + geom_bar( stat="identity", alpha =0.8, position = "dodge") + scale_fill_manual(values = c("#247BA0", "#ffe066"))   + labs(x="Biosynthesis of other secondary metabolites",y="Number of DA KOs") + theme_classic() + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1))+ coord_flip() + facet_grid(~comparison)

gridExtra::grid.arrange(a,b,c,d,e,f, ncol=3,nrow=2)


##################################################
### Between treatments plot 2 : distribution  ###
#################################################
top_categories2 = all %>% group_by(Enrichment, L3, L2,comparison) %>%  tally()
heat=subset(top_categories2, Enrichment == "Heat")
heat$Species=ifelse(heat$comparison %like% "Goniastrea", "Goniastrea", "Porites")
heat$Tissue=ifelse(heat$comparison %like% "band", "Endolithic band", "Skeleton")
heat$label=paste(heat$L2, heat$L3 ,sep = " | ")
top20= heat %>% group_by(comparison)  %>% slice_max(n, n = 10, with_ties = FALSE)


heat.m=reshape2::dcast(heat, label ~ comparison, value.var= "n", fun.aggregate = sum)
rownames(heat.m)=heat.m$label
heat.m=heat.m[,-1]

library(pheatmap)
pheatmap(heat.m, color = colorRampPalette(c("white","red", "red"))(50),  cellwidth = 7, cellheight =  5, fontsize_row = 6, fontsize_col= 8, legend = F, gaps_col = c(2),  cluster_rows = F,cluster_col = F)


