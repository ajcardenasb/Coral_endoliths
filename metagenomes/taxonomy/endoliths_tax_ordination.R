setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/taxonomy/")
library(phyloseq)

map=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = T, row.names = 1, sep = "\t")
all=read.table("Input_files/Families_counts", header = T, sep = "\t", row.names = 1)
tax=read.table("Input_files/Endoliths_kaiju_taxonomy", sep = "\t", header = T, fill = T, quote = "")
tax.u=data.frame(Family=unique(tax$label)) 
tax.u$Domain=tax$Superkingdom[match(tax.u$Family, tax$label)]
rownames(tax.u)=tax.u$Family  

# ordination plot all families
all.fams.p= otu_table(all, taxa_are_rows=TRUE)
sam.p= sample_data(data.frame(map))
tax.p=tax_table(as.matrix(tax.u))
phy.all.fams= phyloseq(all.fams.p,  sam.p, tax.p)
phy.all.fams.t=microbiome::transform(phy.all.fams, transform = "clr", target = "OTU", shift = 0, scale = 1)

##species
P2=c("#8f2d56", "#218380")
phy.all.fams.s=subset_samples(phy.all.fams.t, Treatment == "Control")
PCOA_all.famst = ordinate(phy.all.fams.s, method = "RDA", distance = "euclidean")
pdf("outputs/ordination_allFamilies_species.pdf", width = 7, height = 7, pointsize = 12)
plot_ordination(phy.all.fams.s,PCOA_all.famst, color = "Species", shape = "Tissue") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P2) + ggtitle("") +  theme_classic() + labs(tag = "A", title = "All families") + theme(legend.position = "none") 
dev.off()

##treatments
treat_col=c("#247BA0", "#ffe066")
phy_gon=subset_samples(phy.all.fams.t, Species == "Goniastrea")
phy_por=subset_samples(phy.all.fams.t, Species == "Porites")

PCOA_phy_gon = ordinate(phy_gon, method = "RDA", distance = "euclidean")
PCOA_phy_por = ordinate(phy_por, method = "RDA", distance = "euclidean")

gon_plot=plot_ordination(phy_gon, PCOA_phy_gon, color = "Treatment", shape = "Tissue") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=treat_col) + ggtitle("") +  theme_classic() + facet_grid(~Species)  + labs(tag = "A", title = "Goniastrea") + theme(legend.position = "none") 
por_plot=plot_ordination(phy_por, PCOA_phy_por, color = "Treatment", shape = "Tissue") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=treat_col) + ggtitle("") +  theme_classic() + facet_grid(~Species)  + labs(tag = "B", title = "Porites") + theme(legend.position = "none") 

pdf("outputs/ordination_allFamilies_treatments.pdf", width = 7, height = 4, pointsize = 12)
gridExtra::grid.arrange(gon_plot,por_plot,ncol=2,nrow=1)
dev.off()


### per compartment

arch.p=subset_taxa(phy.all.fams.t, Domain == "Archaea")
PCOA_arch = ordinate(arch.p, method = "RDA", distance = "euclidean")
arch_plot=plot_ordination(arch.p,PCOA_arch, color = "Species", shape = "Tissue") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P2) +  theme_classic() + labs(tag = "B", title = "Archaeal families") + theme(legend.position = "none") 

bac.p=subset_taxa(phy.all.fams.t, Domain == "Bacteria")
PCOA_bac = ordinate(bac.p, method = "RDA", distance = "euclidean")
bac_plot=plot_ordination(bac.p,PCOA_bac, color = "Species", shape = "Tissue") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P2) +  theme_classic() + labs(tag = "C", title = "Bacterial families")  + theme(legend.position = "none") 

euk.p=subset_taxa(phy.all.fams.t, Domain == "Eukaryota")
PCOA_euk = ordinate(euk.p, method = "RDA", distance = "euclidean")
euk_plot=plot_ordination(euk.p,PCOA_euk, color = "Species", shape = "Tissue") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P2) +  theme_classic() + labs(tag = "D", title = "Microbial eukaryotic families") + theme(legend.position = "none") 

vir.p=subset_taxa(phy.all.fams.t, Domain == "Viruses")
PCOA_vir = ordinate(vir.p, method = "RDA", distance = "euclidean")
vir_plot=plot_ordination(vir.p,PCOA_vir, color = "Species", shape = "Tissue") + geom_point(size = 4, alpha = 1) + theme_bw() + scale_colour_manual(values=P2) +  theme_classic() + labs(tag = "E", title = "Viral families") + theme(legend.position = "none") 

pdf("outputs/ordination_allFamilies_perDomain.pdf", width = 7, height = 7, pointsize = 12)
gridExtra::grid.arrange(arch_plot, bac_plot, euk_plot,  vir_plot, ncol=2, nrow=2)
dev.off()

#############################################
############# Permanovas ####################
##############################################

all_fams_df=as.data.frame(t(otu_table(phy.all.fams.t)))
all_fams_df$species=map$Species[match(rownames(all_fams_df), rownames(map))]
all_fams_df$tissue=map$Tissue[match(rownames(all_fams_df), rownames(map))]
all_fams_df$treatment=map$Treatment[match(rownames(all_fams_df), rownames(map))]
all_fams_df$genotype=map$Genotype[match(rownames(all_fams_df), rownames(map))]
all_fams_df$group1=paste(all_fams_df$species, all_fams_df$tissue,all_fams_df$treatment, sep = "_")
all_fams_df$group2=paste(all_fams_df$species, all_fams_df$tissue, sep = "_")

ado=adonis(all_fams_df[,1:1644]~all_fams_df$species*all_fams_df$tissue*all_fams_df$treatment, method = "euclidean", by = )
as.data.frame(ado$aov.tab)

#tissues
all_controls=subset(all_fams_df, treatment == "Control")
pairwise.adonis(all_controls[,1:1644],all_controls$group2, sim.method = "euclidean")

#treatments
all_gon=subset(all_fams_df,  species == "Goniastrea")
pairwise.adonis(all_gon[,1:1644],all_gon$group1, sim.method = "euclidean")
all_por=subset(all_fams_df,  species == "Porites")
pairwise.adonis(all_por[,1:1644],all_por$group1, sim.method = "euclidean")


#Archaea
arc_fams_df=as.data.frame(t(otu_table(arch.p)))
arc_fams_df$species=map$Species[match(rownames(arc_fams_df), rownames(map))]
arc_fams_df$tissue=map$Tissue[match(rownames(arc_fams_df), rownames(map))]
arc_fams_df$treatment=map$Treatment[match(rownames(arc_fams_df), rownames(map))]
arc_fams_df$genotype=map$Genotype[match(rownames(arc_fams_df), rownames(map))]
arc_fams_df$group1=paste(arc_fams_df$species, arc_fams_df$tissue,arc_fams_df$treatment, sep = "_")
arc_fams_df$group2=paste(arc_fams_df$species, arc_fams_df$tissue, sep = "_")
ado=adonis(arc_fams_df[,1:98]~arc_fams_df$species*arc_fams_df$tissue*arc_fams_df$treatment, method = "euclidean", by = )
as.data.frame(ado$aov.tab)

#tissues
arc_controls=subset(arc_fams_df, treatment == "Control")
pairwise.adonis(arc_controls[,1:98],arc_controls$group2, sim.method = "euclidean")

#treatments
arc_gon=subset(arc_fams_df,  species == "Goniastrea")
pairwise.adonis(arc_gon[,1:98],arc_gon$group1, sim.method = "euclidean")
arc_por=subset(arc_fams_df,  species == "Porites")
pairwise.adonis(arc_por[,1:98],arc_por$group1, sim.method = "euclidean")


#Bacteria
bac_fams_df=as.data.frame(t(otu_table(bac.p)))
bac_fams_df$species=map$Species[match(rownames(bac_fams_df), rownames(map))]
bac_fams_df$tissue=map$Tissue[match(rownames(bac_fams_df), rownames(map))]
bac_fams_df$treatment=map$Treatment[match(rownames(bac_fams_df), rownames(map))]
bac_fams_df$genotype=map$Genotype[match(rownames(bac_fams_df), rownames(map))]
bac_fams_df$group1=paste(bac_fams_df$species, bac_fams_df$tissue,bac_fams_df$treatment, sep = "_")
bac_fams_df$group2=paste(bac_fams_df$species, bac_fams_df$tissue, sep = "_")
ado=adonis(bac_fams_df[,1:806]~bac_fams_df$species*bac_fams_df$tissue*bac_fams_df$treatment, method = "euclidean", by = )
as.data.frame(ado$aov.tab)

#tissues
bac_controls=subset(bac_fams_df, treatment == "Control")
pairwise.adonis(bac_controls[,1:806],bac_controls$group2, sim.method = "euclidean")

#treatments
bac_gon=subset(bac_fams_df,  species == "Goniastrea")
pairwise.adonis(bac_gon[,1:806],bac_gon$group1, sim.method = "euclidean")
bac_por=subset(bac_fams_df,  species == "Porites")
pairwise.adonis(bac_por[,1:806],bac_por$group1, sim.method = "euclidean")


#eukaryotes
euk_fams_df=as.data.frame(t(otu_table(euk.p)))
euk_fams_df$species=map$Species[match(rownames(euk_fams_df), rownames(map))]
euk_fams_df$tissue=map$Tissue[match(rownames(euk_fams_df), rownames(map))]
euk_fams_df$treatment=map$Treatment[match(rownames(euk_fams_df), rownames(map))]
euk_fams_df$genotype=map$Genotype[match(rownames(euk_fams_df), rownames(map))]
euk_fams_df$group1=paste(euk_fams_df$species, euk_fams_df$tissue,euk_fams_df$treatment, sep = "_")
euk_fams_df$group2=paste(euk_fams_df$species, euk_fams_df$tissue, sep = "_")
ado=adonis(euk_fams_df[,1:639]~euk_fams_df$species*euk_fams_df$tissue*euk_fams_df$treatment, method = "euclidean", by = )
as.data.frame(ado$aov.tab)

#tissues
euk_controls=subset(euk_fams_df, treatment == "Control")
pairwise.adonis(euk_controls[,1:639],euk_controls$group2, sim.method = "euclidean")

#treatments
euk_gon=subset(euk_fams_df,  species == "Goniastrea")
pairwise.adonis(euk_gon[,1:639],euk_gon$group1, sim.method = "euclidean")
euk_por=subset(euk_fams_df,  species == "Porites")
pairwise.adonis(euk_por[,1:639],euk_por$group1, sim.method = "euclidean")


#viruses
vir_fams_df=as.data.frame(t(otu_table(vir.p)))
vir_fams_df$species=map$Species[match(rownames(vir_fams_df), rownames(map))]
vir_fams_df$tissue=map$Tissue[match(rownames(vir_fams_df), rownames(map))]
vir_fams_df$treatment=map$Treatment[match(rownames(vir_fams_df), rownames(map))]
vir_fams_df$genotype=map$Genotype[match(rownames(vir_fams_df), rownames(map))]
vir_fams_df$group1=paste(vir_fams_df$species, vir_fams_df$tissue,vir_fams_df$treatment, sep = "_")
vir_fams_df$group2=paste(vir_fams_df$species, vir_fams_df$tissue, sep = "_")
ado=adonis(vir_fams_df[,1:100]~vir_fams_df$species*vir_fams_df$tissue*vir_fams_df$treatment, method = "euclidean", by = )
as.data.frame(ado$aov.tab)


#tissues
vir_controls=subset(vir_fams_df, treatment == "Control")
pairwise.adonis(vir_controls[,1:100],vir_controls$group2, sim.method = "euclidean")

#treatments
vir_gon=subset(vir_fams_df,  species == "Goniastrea")
pairwise.adonis(vir_gon[,1:100],vir_gon$group1, sim.method = "euclidean")
vir_por=subset(vir_fams_df,  species == "Porites")
pairwise.adonis(vir_por[,1:100],vir_por$group1, sim.method = "euclidean")

