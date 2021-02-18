
setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/")

####################
#### Binning % #####
####################

library(tidyverse)
library(ggtree)

options(scipen = 5)

map=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = TRUE, row.names = 1, sep ='\t')
pro=read.table("Input_files/profile_endoliths", fill = T, sep = "\t", header = T, row.names = 1)
bact=as.data.frame(t(read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/taxonomy/Input_files/Superkingdom_counts_500", header = T, row.names = 1)))
#bact2=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/taxonomy/outputs/NumberOfBacterialreads.txt", header = T,   sep = " ")
com=pro[, grepl(pattern = "community", x = colnames(pro))]
names(com)=gsub("\\..*", "", names(com))

### reads mapped to bins table = counts for ANCOM
# mapped=pro[, grepl(pattern = "[SB]..mapped.reads", x = colnames(pro))][c(1:75),]
# names(mapped)=gsub("\\..*", "", names(mapped))
# write.table(mapped, "~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/Outputs/MAGs_counts.txt", sep = "\t", quote = F)

bact$bin.per=colSums(com[c(1:75),])
bact$unb.per=t(com[76,])
bact$percentageBac=(bact$Bacteria*100)/rowSums(bact[,c(1:5)])
bact$Unbinned.final=bact$unb.per-(bact$percentageBac-bact$bin.per)
bact$Binned.final=bact$bin.per
bact$UnbinBac.final=bact$percentageBac-bact$bin.per
bact$Sample=rownames(bact)
bact[bact < 0] <- 0

long=reshape2::melt(bact[,c(9:12)], id.vars=c("Sample"), variable.name = "Class", value.name = "Abundance")
long$Species=map$Species[match(long$Sample, rownames(map))]
long$Tissue=map$Tissue[match(long$Sample, rownames(map))]
long$Treatment=map$Treatment[match(long$Sample, rownames(map))]
long$Class=gsub("Unbinned.final", "Unbinned" ,long$Class)
long$Class=gsub("Binned.final", "Bacterial binned" ,long$Class)
long$Class=gsub("UnbinBac.final", "Bacterial unbinned" ,long$Class)
long_sub=subset(long, Treatment == "Control")
long_fixed=long_sub %>% group_by(Species, Class, Tissue) %>% summarise(Abundance=sum(Abundance))
long_sub %>% group_by(Species, Class, Tissue) %>% summarise(Abundance=sum(Abundance)) %>% filter(Tissue == "Endolithic band" & Species == "Goniastrea"  )  %>%  mutate(percentage=Abundance/sum(Abundance)) 

P21=c("#001233", "#023e7d", "#979dac")
pdf("Outputs/binning_stats.pdf",  width = 7, height =3, pointsize = 12) 
ggplot() +geom_bar(aes(y = Abundance, x = Species, fill = Class), data = long_fixed, stat="identity", position = "fill") +  theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of metagenomic sequences", x="") + scale_fill_manual(values=P21) + facet_grid(~Tissue) + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'bottom') +  theme_classic() +facet_grid(~Tissue)
dev.off() 

# perfecntage of bined fraction
100-mean(as.numeric(com[76,]))

## numbers##
# long2=long
# long2$file=gsub("[0-9]", "", long2$file)
# long2 %>% group_by(file, Class) %>% summarize(mean = mean(Abundance))
# stats=aggregate(Abundance ~ file+Class, long2 , mean)
# stats2=dcast(stats, file~Class, value.var= "Abundance", fun.aggregate = mean)
# stats2$summary=paste("(", stats2$Binned.final,"/", (stats2$Binned.final+stats2$UnbinBac.final), ")")
# rowSums(stats2[,2:4])

######################################################################
#### MAGs coverage by phylum with respect to the binned fraction ####
######################################################################

bin=pro[, grepl(pattern = "binned", x = colnames(pro))][c(1:75),]
names(bin)=gsub("\\..*", "", names(bin))
bin$mag=rownames(bin)
bin.l=melt(bin, id.vars=c("mag"), variable.name = "Sample", value.name = "Abundance")
meta.mag= read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/Input_files/tree_metadata.txt", header = T)
bin.l$Phylum=meta.mag$Taxa[match(bin.l$mag, meta.mag$MAG)]
bin.l$Species=map$Species[match(bin.l$Sample, rownames(map))]
bin.l$Tissue=map$Tissue[match(bin.l$Sample, rownames(map))]
bin.l$Treatment=map$Treatment[match(bin.l$Sample, rownames(map))]
bin_fixed= bin.l %>% filter(Treatment == "Control") %>%group_by(Species, Phylum, Tissue ) %>% summarise(Abundance=sum(Abundance))


Pphyla=c("#2077b4","#aec7e8","#ff7f0e", "#ffbb78", "#2ca02c", "#98df89","#17becf", "#9edae5", "#e377c2", "#f7b6d2", "#e33e33",   "#f1c453", "#efea5a", "#16db93" ,  "#B6B6B6")
pdf(file = "./Outputs/MAGs_coverage_binned.pdf", width = 4, height =3, pointsize = 12) 
ggplot() +geom_bar(aes(y = Abundance, x = Species, fill = Phylum), data = bin_fixed, stat="identity", position = "fill") +  theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "MAG relative abundance", x="") + scale_fill_manual(values=Pphyla) +  theme_classic() + facet_grid(~Tissue) + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'none') +facet_grid(~Tissue)
dev.off()

