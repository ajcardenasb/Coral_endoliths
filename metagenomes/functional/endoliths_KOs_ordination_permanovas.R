setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/functional/")
library(ggplot2)
library(phyloseq)

kos=read.table("Input_files/metabolic_KO_counts_500", header = T, row.names = 1)
map=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = TRUE, row.names = 1, sep ='\t')

Ckos=read.table("~/Documents/Bioinformatics_scripts/KEGG_files/ko_C_metabolism", fill = T, sep = "\t")
Ckos$V2=gsub("ko:", "",Ckos$V2 )
Nkos=read.table("~/Documents/Bioinformatics_scripts/KEGG_files/ko_N_metabolism",  fill = T, sep = "\t")
Nkos$V2=gsub("ko:", "",Nkos$V2 )

otu.t= otu_table(kos, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax=data.frame(KOs=rownames(otu.t), Cmet=ifelse(rownames(otu.t) %in% Ckos$V2, "yes", "no"), Nmet=ifelse(rownames(otu.t) %in% Nkos$V2, "yes", "no"))
rownames(tax)=rownames(otu.t)
tax.t= tax_table(as.matrix(tax))
phy.meta= phyloseq(otu.t, sam.t, tax.t)

P2=c("#8f2d56", "#218380")

phy.meta.t=microbiome::transform(phy.meta, transform = "clr", target = "OTU", shift = 0, scale = 1)

###############################
##ordination between species ##
###############################

cnt=subset_samples(phy.meta.t, Treatment=="Control" )


### metab.  KOs
PCOA_cnt = ordinate(cnt, method = "RDA", distance = "euclidean")
a=plot_ordination(cnt,PCOA_cnt, color = "Species", shape = "Tissue") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P2) +  labs(  tag = "A", title = "All KEGG KOs") +  theme_classic() + theme(legend.position = "none") # + stat_ellipse(geom  = "polygon",  alpha = 0.1 )

### C metabolism KOs
c_kos=subset_taxa(cnt, Cmet == "yes" )
PCOA_c = ordinate(c_kos, method = "RDA", distance = "euclidean")
c=plot_ordination(c_kos,PCOA_c , color = "Species", shape = "Tissue") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P2) +  labs(  tag = "C", title = "Carbon metabolism KOs")  +  theme_classic() + theme(legend.position = "none")#+ theme(legend.key = element_blank()) #+ stat_ellipse(geom  = "polygon",  alpha = 0.1 )

### C metabolism KOs
n_kos=subset_taxa(cnt, Nmet == "yes" )
PCOA_n = ordinate(n_kos, method = "RDA", distance = "euclidean")
d=plot_ordination(n_kos,PCOA_n , color = "Species", shape = "Tissue") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P2) +  labs(  tag = "D", title = "Nitrogen metabolism KOs") +  theme_classic() + theme(legend.position = "none")#+ stat_ellipse(geom  = "polygon",  alpha = 0.1 )


#pdf("outputs/KEGG_ordination_Species.pdf", width = 7, height = 7, pointsize = 12)
grid.arrange(a,c,d, ncol=2)
#dev.off()



##################################
##ordination between treatments ##
##################################

P2=c("#247BA0", "#ffe066")

gon_band=subset_samples(phy.meta.t, Species== "Goniastrea" &  Tissue == "Endolithic band")
ord_gon_band = ordinate(gon_band, method = "RDA", distance = "euclidean")
a=plot_ordination(gon_band,ord_gon_band, color = "Treatment") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P2) +  labs(  tag = "A", title = "Goniastrea endolithic band") +  theme_classic() + theme(legend.position = "none") # + stat_ellipse(geom  = "polygon",  alpha = 0.1 )

gon_skel=subset_samples(phy.meta.t, Species== "Goniastrea" &  Tissue == "Skeleton")
ord_gon_skel = ordinate(gon_skel, method = "RDA", distance = "euclidean")
b=plot_ordination(gon_skel,ord_gon_skel, color = "Treatment") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P2) +  labs(  tag = "B", title = "Goniastrea endolithic skeleton") +  theme_classic() + theme(legend.position = "none") # + stat_ellipse(geom  = "polygon",  alpha = 0.1 )

por_band=subset_samples(phy.meta.t, Species== "Porites" &  Tissue == "Endolithic band")
ord_por_band = ordinate(por_band, method = "RDA", distance = "euclidean")
c=plot_ordination(por_band,ord_por_band, color = "Treatment") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P2) +  labs(  tag = "C", title = "Porites endolithic band") +  theme_classic() + theme(legend.position = "none") # + stat_ellipse(geom  = "polygon",  alpha = 0.1 )

por_skel=subset_samples(phy.meta.t, Species== "Porites" &  Tissue == "Skeleton")
ord_por_skel = ordinate(por_skel, method = "RDA", distance = "euclidean")
d=plot_ordination(por_skel,ord_por_skel, color = "Treatment") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P2) +  labs(  tag = "D", title = "Porites endolithic skeleton") +  theme_classic() + theme(legend.position = "none") # + stat_ellipse(geom  = "polygon",  alpha = 0.1 )

#pdf("outputs/KEGG_ordination_treatment.pdf", width = 7, height = 7, pointsize = 12)
grid.arrange(a,b,c,d, ncol=2)
#dev.off()


#######################################################
############# Permanovas full models ##################
#######################################################

met_kos_df=as.data.frame(t(otu_table(phy.meta.t)))
met_kos_df$species=map$Species[match(rownames(met_kos_df), rownames(map))]
met_kos_df$tissue=map$Tissue[match(rownames(met_kos_df), rownames(map))]
met_kos_df$treatment=map$Treatment[match(rownames(met_kos_df), rownames(map))]
met_kos_df$genotype=map$Genotype[match(rownames(met_kos_df), rownames(map))]
met_kos_df$group=paste(met_kos_df$species, met_kos_df$tissue,met_kos_df$treatment, sep = "_")
ado=adonis(met_kos_df[,1:3975]~met_kos_df$species*met_kos_df$tissue*met_kos_df$treatment, method = "euclidean")
temp=as.data.frame(ado$aov.tab)

#species and compartments
met_controls=subset(met_kos_df, treatment == "Control")
met_controls$group2=paste(met_controls$species, met_controls$tissue, sep = "_")
pairwise.adonis(met_controls[,1:3974],met_controls$group2, sim.method = "euclidean")

#treatments
met_gon=subset(met_kos_df,  species == "Goniastrea")
pairwise.adonis(met_gon[,1:3975],met_gon$group, sim.method = "euclidean")
met_por=subset(met_kos_df,  species == "Porites")
pairwise.adonis(met_por[,1:3975],met_por$group, sim.method = "euclidean")

