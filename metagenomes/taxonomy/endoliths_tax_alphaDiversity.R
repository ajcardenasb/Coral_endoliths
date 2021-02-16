library(ggplot2)
library(vegan)
library(ggpubr)
library(GUniFrac)
library(lme4)
library(emmeans)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/taxonomy/")
############################################################
##################### alpha div ############################
############################################################
map= read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = T, sep = "\t", row.names = 1)
#kai=read.table("Input_files/Families_counts", header = T, row.names = 1, quote = "", sep = "\t")
kai=read.table("Input_files/Families_counts_500", header = T, row.names = 1, quote = "", sep = "\t")

###rarefying
cnts=t(round(kai[, 1:48], digits = 0))
min(rowSums(cnts)) # determine sample with lowest counts
kai.r=t(Rarefy(cnts, 698219)$otu.tab.rff)

#kai.r=kai
P2=c("#8f2d56", "#218380")
kai.a=as.data.frame(t(estimateR(t(kai.r))))
kai.a$Shannon=diversity(t(kai.r), index = "shannon")
kai.a$Species=map$Species[match(rownames(kai.a),rownames(map))]
kai.a$Tissue=map$Tissue[match(rownames(kai.a),rownames(map))]
kai.a$Treatment=map$Treatment[match(rownames(kai.a),rownames(map))]
kai.a$Genotype=map$Genotype[match(rownames(kai.a),rownames(map))]
kai.s=subset(kai.a, Treatment == "Control")

######################
##### Observed #######
######################

### full model + temperature
all_model=lmer(S.chao1~Species*Tissue*Treatment* + (1 |Genotype), data = kai.a)
anova(all_model)

#pairwise comaring  treatment per species and tissues
all_pairs = emmeans(all_model, pairwise ~ Treatment|Species|Tissue, weights = "proportional", adjust="none")
rbind(all_pairs$contrasts, adjust="fdr")

### full model = only controls
all_model=lmer(S.chao1~Species*Tissue + (1 |Genotype), data = kai.s)
anova(all_model)

#pairwise comaring Species
all_pairs = emmeans(all_model, pairwise ~ Species|Tissue, weights = "proportional", adjust="none")
rbind(all_pairs$contrasts, adjust="fdr")

#pairwise comaring Tissues
all_pairs = emmeans(all_model, pairwise ~ Tissue|Species, weights = "proportional", adjust="none")
rbind(all_pairs$contrasts, adjust="fdr")


observed_plot=ggplot(data=kai.s, aes(x=Species, y=S.obs, fill=Species)) + scale_fill_manual(values = c ("#8f2d56", "#218380"), name = "Species") + stat_boxplot(geom ='errorbar', width = 0.7, lwd=0.5)+  geom_boxplot(width=0.7, lwd=0.5, fatten=1) + facet_grid(~Tissue, space = "free", scales = "free") +  theme_bw()   + labs(y = "Observed richness", x = "")  + theme( axis.line = element_line(colour = "grey20"), axis.ticks.length = unit(0.2 , "cm"), legend.position="none",  panel.spacing = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 0, vjust = 0.5, face = "bold"), axis.text.x = element_text(color = "grey20", size = 8, angle = 0, hjust = .5, vjust = .5, face = "bold") )# + theme_classic2()
#observed_plot=ggplot(data=kai.a, aes(x=Treatment, y=S.obs, fill=Treatment)) + scale_fill_manual(values = c ("#247BA0", "#ffe066"), name = "Species") + stat_boxplot(geom ='errorbar', width = 0.7, lwd=0.5)+  geom_boxplot(width=0.7, lwd=0.5, fatten=1) + facet_grid(~Species+Tissue, space = "free", scales = "free") +  theme_bw()   + labs(y = "Observed richness", x = "")  + theme( axis.line = element_line(colour = "grey20"), axis.ticks.length = unit(0.2 , "cm"), legend.position="none",  panel.spacing = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 0, vjust = 0.5, face = "bold"), axis.text.x = element_text(color = "grey20", size = 8, angle = 0, hjust = .5, vjust = .5, face = "bold") )# + theme_classic2()


#####################
##### Shannon #######
##################### 

all_model=lmer(Shannon~Treatment*Species*Tissue + (1 |Genotype), data = kai.a)
anova(all_model) 


### full model = only controls
all_model=lmer(Shannon~Species*Tissue + (1 |Genotype), data = kai.s)

#pairwise comaring Species
all_pairs = emmeans(all_model, pairwise ~ Species|Tissue, weights = "proportional", adjust="none")
rbind(all_pairs$contrasts, adjust="fdr")

#pairwise comaring Tissues
all_pairs = emmeans(all_model, pairwise ~ Tissue|Species, weights = "proportional", adjust="none")
rbind(all_pairs$contrasts, adjust="fdr")

shannon_plot=ggplot(data=kai.s, aes(x=Species, y=Shannon, fill=Species)) + scale_fill_manual(values = c ("#8f2d56", "#218380"), name = "Species") + stat_boxplot(geom ='errorbar', width = 0.7, lwd=0.5)+  geom_boxplot(width=0.7, lwd=0.5, fatten=1) + facet_grid(~Tissue, space = "free", scales = "free") +  theme_bw()   + labs(y = "Shannon diversity", x = "")  + theme( axis.line = element_line(colour = "grey20"), axis.ticks.length = unit(0.2 , "cm"), legend.position="none",  panel.spacing = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 0, vjust = 0.5, face = "bold"), axis.text.x = element_text(color = "grey20", size = 8, angle = 0, hjust = .5, vjust = .5, face = "bold") ) #+ theme_classic2()
#shannon_plot=ggplot(data=kai.s, aes(x=Species, y=Shannon, fill=Species)) + scale_fill_manual(values = c ("#8f2d56", "#218380"), name = "Species") + stat_boxplot(geom ='errorbar', width = 0.7, lwd=0.5)+  geom_boxplot(width=0.7, lwd=0.5, fatten=1) + facet_grid(~Tissue, space = "free", scales = "free") +  theme_bw()   + labs(y = "Shannon diversity", x = "")  + theme( axis.line = element_line(colour = "grey20"), axis.ticks.length = unit(0.2 , "cm"), legend.position="none",  panel.spacing = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 0, vjust = 0.5, face = "bold"), axis.text.x = element_text(color = "grey20", size = 8, angle = 0, hjust = .5, vjust = .5, face = "bold") ) #+ theme_classic2()


pdf("outputs/alphaDiversity_taxa.pdf", width = 5, height = 5, pointsize = 12)
gridExtra::grid.arrange(observed_plot, shannon_plot, ncol=1)
dev.off()

# pdf("outputs/alphaDiversity_taxa2.pdf", width = 5, height = 5, pointsize = 12)
# gridExtra::grid.arrange(observed_plot, shannon_plot, ncol=1)
# dev.off()