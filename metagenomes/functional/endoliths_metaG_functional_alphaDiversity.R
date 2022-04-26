library(ggplot2)
library(vegan)
library(ggpubr)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/functional/")
############################################################
##################### alpha div ############################
############################################################
map= read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = T, sep = "\t", row.names = 1)
kai=read.table("Input_files/metabolic_KO_counts_500", header = T, row.names = 1, quote = "", sep = "\t")
kai.r=round(kai, digits = 0)#sweep(kai,2,colSums(kai),"/")
kai.r[ kai.r < 20] <- 0
##rarefying
cnts=t(kai.r[, 1:48])
min_val=min(rowSums(cnts)) #686534
kai.rar=Rarefy(cnts, min_val)$otu.tab.rff

P2=c("#8f2d56", "#218380")

kai.a=as.data.frame(t(estimateR(kai.rar)))
kai.a$Shannon=diversity(kai.rar, index = "shannon")
kai.a$Species=map$Species[match(rownames(kai.a),rownames(map))]
kai.a$Tissue=map$Tissue[match(rownames(kai.a),rownames(map))]
kai.a$Treatment=map$Treatment[match(rownames(kai.a),rownames(map))]
kai.a$Genotype=map$Genotype[match(rownames(kai.a),rownames(map))]
kai.s=subset(kai.a, Treatment == "Control")


all_model=lmer(Shannon~Treatment*Species*Tissue + (1 |Genotype), data = kai.a)
anova(all_model) 
all_pairs = emmeans(all_model, pairwise ~ Treatment|Species|Tissue, weights = "proportional", adjust="none")
rbind(all_pairs$contrasts, adjust="fdr")

### full model = only controls
all_model=lmer(Shannon~Species*Tissue + (1 |Genotype), data = kai.s)

#pairwise comaring Species
all_pairs = emmeans(all_model, pairwise ~ Species|Tissue, weights = "proportional", adjust="none")
rbind(all_pairs$contrasts, adjust="fdr")

#pairwise comaring Tissues
all_pairs = emmeans(all_model, pairwise ~ Tissue|Species, weights = "proportional", adjust="none")
rbind(all_pairs$contrasts, adjust="fdr")

pdf("outputs/alphaDiversity_KOs.pdf", width = 5, height = 5, pointsize = 12)
shannon_plot=ggplot(data=kai.s, aes(x=Species, y=Shannon, fill=Species)) + scale_fill_manual(values = c ("#8f2d56", "#218380"), name = "Species") + stat_boxplot(geom ='errorbar', width = 0.7, lwd=0.5)+  geom_boxplot(width=0.7, lwd=0.5, fatten=1) + facet_grid(~Tissue, space = "free", scales = "free") +  theme_bw()   + labs(y = "KEGG orthologs Shannon diversity", x = "")  + theme( axis.line = element_line(colour = "grey20"), axis.ticks.length = unit(0.2 , "cm"), legend.position="none",  panel.spacing = unit(1, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 0, vjust = 0.5, face = "bold"), axis.text.x = element_text(color = "grey20", size = 8, angle = 0, hjust = .5, vjust = .5, face = "bold") ) #+ theme_classic2()
dev.off()
