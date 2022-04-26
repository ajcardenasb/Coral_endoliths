library(ggplot2)
library(phyloseq)
setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/functional/")
########################################################
##################### Ordination ####################
########################################################

kos=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/functional/Input_files/KO_counts_500", header = T, row.names = 1)

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

Ptreat=c("#247BA0", "#ffe066")

phy.meta.t=microbiome::transform(phy.meta, transform = "clr", target = "OTU", shift = 0, scale = 1)

gon_band=subset_samples(phy.meta.t, Species== "Goniastrea") #&  Tissue == "Endolithic band")
gon_band_met =subset_taxa(gon_band, Met== "yes" )
ord_gon_band = ordinate(gon_band_met, method = "RDA", distance = "euclidean")
ordi_a=plot_ordination(gon_band_met,ord_gon_band, color = "Treatment", shape = "Tissue") + geom_point(size = 3, alpha = 1) + theme_bw() + scale_colour_manual(values=Ptreat) +  labs(  title = "Goniastrea endolithic band") +  theme_classic() + theme(legend.position = "none") # + stat_ellipse(geom  = "polygon",  alpha = 0.1 )

por_band=subset_samples(phy.meta.t, Species== "Porites") #&  Tissue == "Endolithic band")
por_band_met =subset_taxa(por_band, Met== "yes" )
ord_por_band = ordinate(por_band_met, method = "RDA", distance = "euclidean")
ordi_b=plot_ordination(por_band_met,ord_por_band, color = "Treatment", shape = "Tissue") + geom_point(size = 3, alpha = 1) + theme_bw() + scale_colour_manual(values=Ptreat) +  labs(  title = "Porites endolithic band") +  theme_classic() + theme(legend.position = "none") # + stat_ellipse(geom  = "polygon",  alpha = 0.1 )

#pdf("outputs/KEGG_ordination_treatment.pdf", width = 7, height = 7, pointsize = 12)
#grid.arrange(a,b,ncol=2)
#dev.off()
## 4 compartments
gon_band=subset_samples(phy.meta.t, Species== "Goniastrea" &  Tissue == "Endolithic band")
gon_band_met =subset_taxa(gon_band, Met== "yes" )
ord_gon_band = ordinate(gon_band_met, method = "RDA", distance = "euclidean")
ordi_a=plot_ordination(gon_band_met,ord_gon_band, color = "Treatment") + geom_point(size = 3, alpha = 1) + theme_bw() + scale_colour_manual(values=Ptreat) +  labs(  title = "") +  theme_classic() + theme(legend.position = "none") # + stat_ellipse(geom  = "polygon",  alpha = 0.1 )

gon_skel=subset_samples(phy.meta.t, Species== "Goniastrea" &  Tissue == "Skeleton")
gon_skel_met =subset_taxa(gon_skel, Met== "yes" )
ord_gon_skel = ordinate(gon_skel_met, method = "RDA", distance = "euclidean")
ordi_b=plot_ordination(gon_skel_met,ord_gon_skel, color = "Treatment") + geom_point(size = 3, alpha = 1) + theme_bw() + scale_colour_manual(values=Ptreat) +  labs(  title = "") +  theme_classic() + theme(legend.position = "none") # + stat_ellipse(geom  = "polygon",  alpha = 0.1 )

por_band=subset_samples(phy.meta.t, Species== "Porites" &  Tissue == "Endolithic band")
por_band_met =subset_taxa(por_band, Met== "yes" )
ord_por_band = ordinate(por_band_met, method = "RDA", distance = "euclidean")
ordi_c=plot_ordination(por_band_met,ord_por_band, color = "Treatment") + geom_point(size = 3, alpha = 1) + theme_bw() + scale_colour_manual(values=Ptreat) +  labs(  title = "") +  theme_classic() + theme(legend.position = "none") # + stat_ellipse(geom  = "polygon",  alpha = 0.1 )

por_skel=subset_samples(phy.meta.t, Species== "Porites" &  Tissue == "Skeleton")
por_skel_met =subset_taxa(por_skel, Met== "yes" )
ord_por_skel = ordinate(por_skel_met, method = "RDA", distance = "euclidean")
ordi_d=plot_ordination(por_skel_met,ord_por_skel, color = "Treatment") + geom_point(size = 3, alpha = 1) + theme_bw() + scale_colour_manual(values=Ptreat) +  labs(  title = "") +  theme_classic() + theme(legend.position = "none") # + stat_ellipse(geom  = "polygon",  alpha = 0.1 )


########################################################
##################### Nano SIMS ####################
########################################################

nan=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/nanoSIMS/input_files/nanoSIMS_data.txt", header = T, sep = "\t")
nan$Treatment=gsub("27", "Control", nan$Treatment)
nan$Treatment=gsub("33", "Heat", nan$Treatment)
Ptreat=c("#247BA0", "#ffe066")

### C assimiation between treatments
Carb=subset(nan, compartment == "both" & Experiment== "13C")
nan_a=ggplot(Carb, aes(Treatment, average, color = Treatment, fill =Treatment)) + geom_col(data = Carb, position = position_dodge(0.8), width = 1) +  geom_errorbar( aes(ymin = average, ymax = average+se), data = Carb , position = position_dodge(0.8)) + scale_color_manual(values = Ptreat)  + theme_classic() + scale_fill_manual(values = Ptreat)+ theme(legend.position = "none") + facet_grid(~Species)  + labs(y= "13C assimilation (Atom% excess)")  + ylim(0, 0.15)

### C assimiation between treatments
Nitro=subset(nan, compartment == "both" & Experiment== "15N")
nan_b=ggplot(Nitro, aes(Treatment, average, color = Treatment, fill =Treatment)) + geom_col(data = Nitro, position = position_dodge(0.8), width = 1) +  geom_errorbar( aes(ymin = average, ymax = average+se), data = Nitro , position = position_dodge(0.8)) + scale_color_manual(values = Ptreat)  + theme_classic() + scale_fill_manual(values = Ptreat)+ theme(legend.position = "none") + facet_grid(~Species)   + labs(y= "15N assimilation (Atom% excess)") + ylim(0, 0.3)

#grid.arrange(nan_a,nan_b,ncol=2)
### 
nan_b=ggplot(Nitro, aes(Treatment, average, color = Treatment, fill =Treatment)) + geom_col(data = Nitro, position = position_dodge(0.8), width = 1) +  geom_errorbar( aes(ymin = average, ymax = average+se), data = Nitro , position = position_dodge(0.8)) + scale_color_manual(values = Ptreat)  + theme_classic() + scale_fill_manual(values = Ptreat)+ theme(legend.position = "none") + facet_grid(~Species)   + labs(y= "15N assimilation (Atom% excess)") + ylim(0, 0.002)
pdf("outputs/nano_FuncPanel2.pdf", height = 3, width = 4, pointsize = 10)
grid.arrange(nan_a,nan_b,ncol=2)
dev.off()
#######################################################
#################### Volcano plots ####################
#######################################################

#L2_col =c("#669900","#ccee66","#006699","#3399cc","#990066","#DA6CB6","#ff6600","#ff9900","#ffcc00","#fff275", "#000000", "#000000")
Ptreat=c("#247BA0", "#ffe066", "#000000")
vol_gon_band=read.table("outputs/volcano_gon_band", sep = "\t", header = T)
vol_gon_band$color=ifelse(vol_gon_band$DA.TreatmentHeat=="FALSE", "no DA", ifelse(vol_gon_band$W.TreatmentHeat  > 0, "Heat", "Control"))
vol_gon_band_s=subset(vol_gon_band, L1 %in% " Metabolism")# & !L2 %in% c(" Not included in regular maps") & !is.na(L2))
vol1=ggplot(vol_gon_band_s, aes(x = Beta.TreatmentHeat, y = -log10(qval.TreatmentHeat+0.001), color = color)) + scale_colour_manual(values = Ptreat) + ggtitle(label = "", subtitle = "") + geom_point(size = 2.5, alpha = 1, na.rm = T) + theme_bw(base_size = 14) +  theme(legend.position = "none") + xlab(expression(log[2]("Effect size (W)"))) + ylab(expression(-log[10]("q-val"))) + geom_hline(yintercept = 1.30102, colour="#990000", linetype="dashed")  + geom_vline(xintercept = 1, colour="#990000", linetype="dashed") + geom_vline(xintercept = -1, colour="#990000", linetype="dashed") + scale_y_continuous(trans = "log1p", limits = c(0,8))

vol_gon_skel=read.table("outputs/volcano_gon_skel", sep = "\t", header = T)
vol_gon_skel$color=ifelse(vol_gon_skel$DA.TreatmentHeat=="FALSE", "no DA", ifelse(vol_gon_skel$W.TreatmentHeat  > 0, "Heat", "Control"))
vol_gon_skel_s=subset(vol_gon_skel, L1 %in% " Metabolism" & !L2 %in% c(" Not included in regular maps") & !is.na(L2))
vol2=ggplot(vol_gon_skel_s, aes(x = Beta.TreatmentHeat, y = -log10(qval.TreatmentHeat), color = color)) + scale_colour_manual(values = Ptreat) + ggtitle(label = "", subtitle = "") + geom_point(size = 2.5, alpha = 1, na.rm = T) + theme_bw(base_size = 14) +  theme(legend.position = "none") + xlab(expression(log[2]("Effect size (W)"))) + ylab(expression(-log[10]("q-val"))) + geom_hline(yintercept = 1.30102, colour="#990000", linetype="dashed") + geom_vline(xintercept = 1, colour="#990000", linetype="dashed") + geom_vline(xintercept = -1, colour="#990000", linetype="dashed") + scale_y_continuous(trans = "log1p", limits = c(0,8))

vol_por_band=read.table("outputs/volcano_por_band", sep = "\t", header = T)
vol_por_band$color=ifelse(vol_por_band$DA.TreatmentHeat=="FALSE", "no DA", ifelse(vol_por_band$W.TreatmentHeat  > 0, "Heat", "Control"))
vol_por_band_s=subset(vol_por_band, L1 %in% " Metabolism" & !L2 %in% c(" Not included in regular maps") & !is.na(L2))
vol3=ggplot(vol_por_band_s, aes(x = Beta.TreatmentHeat, y = -log10(qval.TreatmentHeat), color =color)) + scale_colour_manual(values = Ptreat) + ggtitle(label = "", subtitle = "") + geom_point(size = 2.5, alpha = 1, na.rm = T) + theme_bw(base_size = 14) +  theme(legend.position = "none") + xlab(expression(log[2]("Effect size (W)"))) + ylab(expression(-log[10]("q-val"))) + geom_hline(yintercept = 1.30102, colour="#990000", linetype="dashed")  + geom_vline(xintercept = 1, colour="#990000", linetype="dashed") + geom_vline(xintercept = -1, colour="#990000", linetype="dashed") + scale_y_continuous(trans = "log1p", limits = c(0,8))

vol_por_skel=read.table("outputs/volcano_por_skel", sep = "\t", header = T)
vol_por_skel$color=ifelse(vol_por_skel$DA.TreatmentHeat=="FALSE", "no DA", ifelse(vol_por_skel$W.TreatmentHeat  > 0, "Heat", "Control"))
vol_por_skel_s=subset(vol_por_skel, L1 %in% " Metabolism" & !L2 %in% c(" Not included in regular maps") & !is.na(L2))
vol4=ggplot(vol_por_skel_s, aes(x = Beta.TreatmentHeat, y = -log10(qval.TreatmentHeat+0.001), color = color)) + scale_colour_manual(values = Ptreat) + ggtitle(label = "", subtitle = "") + geom_point(size = 2.5, alpha = 1, na.rm = T) + theme_bw(base_size = 14) +  theme(legend.position = "none") + xlab(expression(log[2]("Effect size (W)"))) + ylab(expression(-log[10]("q-val"))) + geom_hline(yintercept = 1.30102, colour="#990000", linetype="dashed") + geom_vline(xintercept = 1, colour="#990000", linetype="dashed") + geom_vline(xintercept = -1, colour="#990000", linetype="dashed") + scale_y_continuous(trans = "log1p", limits = c(0,8))


#pdf("outputs/volcanos_ANCOMBC_species_kos.pdf", height = 5, width = 20, pointsize = 10)
grid.arrange(vol1,vol2,vol3,vol4,ncol=2)
#dev.off()

# pdf("outputs/volcanos_ANCOMBC_legend.pdf", height = 5, width = 20, pointsize = 10)
# ggplot(vol_por_skel_s, aes(x = Beta.TreatmentHeat, y = -log10(qval.TreatmentHeat+0.001), color = L2)) + scale_colour_manual(values = L2_col) + ggtitle(label = "Porites skeleton", subtitle = "Enriched in Control - Enriched in Heat") + geom_point(size = 2.5, alpha = 1, na.rm = T) + theme_bw(base_size = 14) +  theme(legend.position = "right") + xlab(expression(log[2]("Effect size (W)"))) + ylab(expression(-log[10]("p-val"))) + geom_hline(yintercept = 1.30102, colour="#990000", linetype="dashed") + geom_vline(xintercept = 1, colour="#990000", linetype="dashed") + geom_vline(xintercept = -1, colour="#990000", linetype="dashed") + scale_y_continuous(trans = "log1p", limits = c(0,8))
# dev.off()

##################################################
#################### Barplots ####################
##################################################
anc_gon_band=read.table("outputs/ANCOMBC_kos_betweenTreatments_gon_band.txt", sep = "\t", header = T)
anc_gon_band_s=subset(anc_gon_band, L1 %in% " Metabolism")
anc_por_band=read.table("outputs/ANCOMBC_kos_betweenTreatments_por_band.txt", sep = "\t", header = T)
anc_por_band_s=subset(anc_por_band, L1 %in% " Metabolism")

anc_gon_band_filtered = anc_gon_band_s %>% group_by(Diff_more_abundant, L2, L3) %>% filter(Diff_more_abundant == "Heat") %>% filter( L2 %in% c(" Amino acid metabolism", " Carbohydrate metabolism", " Energy metabolism", " Glycan biosynthesis and metabolism")) %>%   dplyr::summarise(Effect_size = max(W_Treatment))
anc_por_band_filtered = anc_por_band_s %>% group_by(Diff_more_abundant, L2, L3) %>% filter(Diff_more_abundant == "Heat") %>% filter( L2 %in% c(" Amino acid metabolism", " Carbohydrate metabolism", " Energy metabolism", " Glycan biosynthesis and metabolism")) %>%   dplyr::summarise(Effect_size = max(W_Treatment))

L2_col1 =c("#669900","#006699","#3399cc","#990066","#DA6CB6","#ff9900","#ffcc00","#fff275", "#000000", "#000000")
bar_a=ggplot(anc_gon_band_filtered, aes(x=reorder(L3,Effect_size), y=Effect_size, label = L3)) + geom_text(aes(y = 0) ) + geom_bar(stat = "identity", aes(fill = L2)) + coord_flip()  + scale_x_discrete(breaks = NULL) + theme_bw() + scale_fill_manual(values =L2_col1) + ylim(-7,7) + labs(x="", title = "Goniastrea") + theme_classic() + theme(axis.text.y = element_text(size = 1))          
L2_col2 =c("#669900","#006699","#3399cc","#990066","#ffcc00","#fff275", "#000000", "#000000")
bar_b=ggplot(anc_por_band_filtered, aes(x=reorder(L3,Effect_size), y=Effect_size, label = L3)) + geom_text(aes(y = 0) ) + geom_bar(stat = "identity", aes(fill = L2)) + coord_flip()  + scale_x_discrete(breaks = NULL) + theme_bw() + scale_fill_manual(values =L2_col2) + ylim(-7,7) + labs(x="", title = "Porites") + theme_classic() + theme(axis.text.y = element_text(size = 1))          

grid.arrange(bar_a,bar_b,ncol=2)



##############
### Plots ####
##############

# pdf("outputs/ordi_FuncPanel2.pdf", height = 2, width = 3, pointsize = 10)
# grid.arrange(ordi_a,ordi_b,ncol=1)
# dev.off()

pdf("outputs/ordi_FuncPanel2.pdf", height = 4, width = 4, pointsize = 10)
grid.arrange(ordi_a,ordi_b,ordi_c,ordi_d,ncol=2)
dev.off()

pdf("outputs/nano_FuncPanel2.pdf", height = 3, width = 4, pointsize = 10)
grid.arrange(nan_a,nan_b,ncol=2)
dev.off()

pdf("outputs/volca_FuncPanel2.pdf", height = 7, width = 4.5, pointsize = 10)
grid.arrange(vol1,vol2,vol3,vol4,ncol=2)
dev.off()

pdf("outputs/bar_FuncPanel2.pdf", height = 6, width = 4.5, pointsize = 10)
grid.arrange(bar_a,bar_b,ncol=1)
dev.off()


##############
### donuts ####
##############

data_A=subset(nan, Species == "Goniastrea" & Treatment == "Control" & Experiment == "13C" & !compartment == "both" )
data_A$fraction = data_A$average / sum(data_A$average) # Compute percentages
data_A$ymax = cumsum(data_A$fraction) # Compute the cumulative percentages (top of each rectangle)
data_A$ymin = c(0, head(data_A$ymax, n=-1)) # Compute the bottom of each rectangle
data_A$labelPosition <- (data_A$ymax + data_A$ymin) / 2
data_A$label = paste0(data_A$compartment, "\n value: ", data_A$fraction)
dona_1=ggplot(data_A, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=compartment)) +  geom_rect() + geom_text( x=2, aes(y=labelPosition, label=label, color=compartment), size=6) + scale_fill_brewer(palette=3) +  scale_color_brewer(palette=3) + coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")

data_E=subset(nan, Species == "Goniastrea" & Treatment == "Heat" & Experiment == "13C" & !compartment == "both" )
data_E$fraction = data_E$average / sum(data_E$average) # Compute percentages
data_E$ymax = cumsum(data_E$fraction) # Compute the cumulative percentages (top of each rectangle)
data_E$ymin = c(0, head(data_E$ymax, n=-1)) # Compute the bottom of each rectangle
data_E$labelPosition <- (data_E$ymax + data_E$ymin) / 2
data_E$label = paste0(data_E$compartment, "\n value: ", data_E$fraction)
dona_2=ggplot(data_E, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=compartment)) +  geom_rect() + geom_text( x=2, aes(y=labelPosition, label=label, color=compartment), size=6) + scale_fill_brewer(palette=3) +  scale_color_brewer(palette=3) + coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")


data_B=subset(nan, Species == "Porites" & Treatment == "Control" & Experiment == "13C" & !compartment == "both" )
data_B$fraction = data_B$average / sum(data_B$average) # Compute percentages
data_B$ymax = cumsum(data_B$fraction) # Compute the cumulative percentages (top of each rectangle)
data_B$ymin = c(0, head(data_B$ymax, n=-1)) # Compute the bottom of each rectangle
data_B$labelPosition <- (data_B$ymax + data_B$ymin) / 2
data_B$label = paste0(data_B$compartment, "\n value: ", data_B$fraction)
dona_3=ggplot(data_B, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=compartment)) +  geom_rect() + geom_text( x=2, aes(y=labelPosition, label=label, color=compartment), size=6) + scale_fill_brewer(palette=3) +  scale_color_brewer(palette=3) + coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")


data_G=subset(nan, Species == "Porites" & Treatment == "Heat" & Experiment == "13C" & !compartment == "both" )
data_G$fraction = data_G$average / sum(data_G$average) # Compute percentages
data_G$ymax = cumsum(data_G$fraction) # Compute the cumulative percentages (top of each rectangle)
data_G$ymin = c(0, head(data_G$ymax, n=-1)) # Compute the bottom of each rectangle
data_G$labelPosition <- (data_G$ymax + data_G$ymin) / 2
data_G$label = paste0(data_G$compartment, "\n value: ", data_G$fraction)
dona_4=ggplot(data_G, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=compartment)) +  geom_rect() + geom_text( x=2, aes(y=labelPosition, label=label, color=compartment), size=6) + scale_fill_brewer(palette=3) +  scale_color_brewer(palette=3) + coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")

data_C=subset(nan, Species == "Goniastrea" & Treatment == "Control" & Experiment == "15N" & !compartment == "both" )
data_C$fraction = data_C$average / sum(data_C$average) # Compute percentages
data_C$ymax = cumsum(data_C$fraction) # Compute the cumulative percentages (top of each rectangle)
data_C$ymin = c(0, head(data_C$ymax, n=-1)) # Compute the bottom of each rectangle
data_C$labelPosition <- (data_C$ymax + data_C$ymin) / 2
data_C$label = paste0(data_C$compartment, "\n value: ", data_C$fraction)
dona_5=ggplot(data_C, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=compartment)) +  geom_rect() + geom_text( x=2, aes(y=labelPosition, label=label, color=compartment), size=6) + scale_fill_brewer(palette=3) +  scale_color_brewer(palette=3) + coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")

data_F=subset(nan, Species == "Goniastrea" & Treatment == "Heat" & Experiment == "15N" & !compartment == "both" )
data_F$fraction = data_F$average / sum(data_F$average) # Compute percentages
data_F$ymax = cumsum(data_F$fraction) # Compute the cumulative percentages (top of each rectangle)
data_F$ymin = c(0, head(data_F$ymax, n=-1)) # Compute the bottom of each rectangle
data_F$labelPosition <- (data_F$ymax + data_F$ymin) / 2
data_F$label = paste0(data_F$compartment, "\n value: ", data_F$fraction)
dona_6=ggplot(data_F, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=compartment)) +  geom_rect() + geom_text( x=2, aes(y=labelPosition, label=label, color=compartment), size=6) + scale_fill_brewer(palette=3) +  scale_color_brewer(palette=3) + coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")

data_D=subset(nan, Species == "Porites" & Treatment == "Control" & Experiment == "15N" & !compartment == "both" )
data_D$fraction = data_D$average / sum(data_D$average) # Compute percentages
data_D$ymax = cumsum(data_D$fraction) # Compute the cumulative percentages (top of each rectangle)
data_D$ymin = c(0, head(data_D$ymax, n=-1)) # Compute the bottom of each rectangle
data_D$labelPosition <- (data_D$ymax + data_D$ymin) / 2
data_D$label = paste0(data_D$compartment, "\n value: ", data_D$fraction)
dona_7=ggplot(data_D, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=compartment)) +  geom_rect() + geom_text( x=2, aes(y=labelPosition, label=label, color=compartment), size=6) + scale_fill_brewer(palette=3) +  scale_color_brewer(palette=3) + coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")


data_H=subset(nan, Species == "Porites" & Treatment == "Heat" & Experiment == "15N" & !compartment == "both" )
data_H$fraction = data_H$average / sum(data_H$average) # Compute percentages
data_H$ymax = cumsum(data_H$fraction) # Compute the cumulative percentages (top of each rectangle)
data_H$ymin = c(0, head(data_H$ymax, n=-1)) # Compute the bottom of each rectangle
data_H$labelPosition <- (data_H$ymax + data_H$ymin) / 2
data_H$label = paste0(data_H$compartment, "\n value: ", data_H$fraction)
dona_8=ggplot(data_H, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=compartment)) +  geom_rect() + geom_text( x=2, aes(y=labelPosition, label=label, color=compartment), size=6) + scale_fill_brewer(palette=3) +  scale_color_brewer(palette=3) + coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")



pdf("outputs/NanoSIMS_donut_treatment.pdf", height = 8, width = 3, pointsize = 10)
plot=gridExtra::grid.arrange(dona_1,dona_2,dona_3,dona_4,dona_5,dona_6,dona_7,dona_8,ncol=2)
dev.off()

