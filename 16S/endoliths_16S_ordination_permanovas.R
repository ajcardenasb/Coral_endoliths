############################################################
##################### phyloseq #############################
############################################################
library(phyloseq)
library(ggplot2)
library(plyr)
library(gridExtra)
setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/16S/ASVs/")

asv=read.table("./outputs/ASVs_noContanoOut.raw.txt", header = TRUE, row.names = 1)
map=read.table("./Input_files/skeleton16S-metadata.txt", header = TRUE, row.names = 1, sep ='\t')
map$Tissue=gsub("Green band", "Endolithic band", map$Tissue)
map$Tissue=gsub("Green band", "Endolithic band", map$Tissue)

otu.t= otu_table(asv[, 1:65], taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t= tax_table(as.matrix(asv[, 66:ncol(asv)]))

phy.all= phyloseq(otu.t, tax.t,  sam.t)
P3=c("#e36600ff", "#008000ff", "#ccccccff")

phy.t=microbiome::transform(phy.all, transform = "log10", target = "OTU", shift = 0, scale = 1)

gon=subset_samples(phy.t, Species=="Goniastrea" & Treatment == "Control")# & !Genotype =="G3")
PCOA_gon = ordinate(gon, method = "RDA", distance = "euclidean")

por=subset_samples(phy.t, Species=="Porites"& Treatment == "Control")
PCOA_por = ordinate(por, method = "RDA", distance = "euclidean")

o.g=plot_ordination(phy.t,PCOA_gon, color = "Tissue") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P3) + ggtitle("Goniastrea") +  theme_classic()
o.p=plot_ordination(phy.t,PCOA_por, color = "Tissue") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P3) + ggtitle("Porites") +  theme_classic()

pdf("./outputs/skeleton16S_ordination.pdf", width=8,height=3, pointsize = 12)
grid.arrange(o.g,o.p, ncol=2, nrow=1)
dev.off()


#####################################################
######## Stats on community composition ############
####################################################
library(vegan)
library(pairwiseAdonis)
asv.n=as.data.frame(t(sweep(asv[, 1:65],2,colSums(asv[, 1:65]),`/`)))
asv.n$tissue=map$Tissue[match(rownames(asv.n), rownames(map))]
asv.n$species=map$Species[match(rownames(asv.n), rownames(map))]
asv.n$treatment=map$Treatment[match(rownames(asv.n), rownames(map))]
asv.n$genotype=map$Genotype[match(rownames(asv.n), rownames(map))]

##overal model
asv_adonis=adonis(asv.n[,1:6491]~ asv.n$species * asv.n$tissue * asv.n$treatment ) # 0.6666667
asv_adonis_df=as.data.frame(asv_adonis[["aov.tab"]])
write.table(asv_adonis_df, "outputs/endoliths_overal_adonis_ASVs.txt", sep = "\t", row.names = T, quote = F)

#comparing species in tissue pools
asv.n$groups=paste(asv.n$species, asv.n$tissue, sep = "-")
controls=subset(asv.n,  treatment == "Control") 
species_pairWadonis_df=pairwise.adonis(controls[,1:6491], controls$groups,  sim.method = "bray", p.adjust.m = "fdr", perm = 999) # 0.6666667
write.table(species_pairWadonis_df, "outputs/endoliths_pairwiseAdonis_speceis_ASVs.txt", sep = "\t", row.names = F, quote = F)


#comparing treatment in tissue pools

tiss.G=subset(asv.n, species == "Goniastrea"  & treatment == "Control" & tissue =="Tissue") 
adonis(tiss.G[,1:6491]~ tiss.G$species) # 0.6666667

gree.G=subset(asv.n, species == "Goniastrea"  & treatment == "Control" & tissue =="Green band") 
adonis(gree.G[,1:6491]~ gree.G$treatment) # 0.675 

whit.G=subset(asv.n, species == "Goniastrea"  & treatment == "Control" & tissue =="White band") 
adonis(whit.G[,1:6491]~whit.G$treatment) # 0.903

tiss.P=subset(asv.n, species == "Porites"  &   tissue =="Tissue" & treatment == "Control") 
adonis(tiss.P[,1:6491]~ tiss.P$treatment ) # 0.756 

gree.P=subset(asv.n, species == "Porites"  & tissue =="Green band" & treatment == "Control") 
adonis(gree.P[,1:6495]~ gree.P$treatment) # 0.945  

whit.P=subset(asv.n, species == "Porites"  &  tissue =="White band" & treatment == "Control") 
adonis(whit.P[,1:6495]~ whit.P$treatment ) # 0.952


#comparing tissues from controls

por=subset(asv.n, species == "Porites" & treatment =="Control")
pairwiseAdonis::pairwise.adonis(por[,1:6491], por$tissue,  sim.method = 'bray', p.adjust.m ='fdr',perm=999) # 0.6666667

gon=subset(asv.n, species == "Goniastrea" & treatment =="Control")
pairwiseAdonis::pairwise.adonis(gon[,1:6491], gon$tissue,  sim.method = 'bray', p.adjust.m ='fdr',perm=999) # 0.6666667

