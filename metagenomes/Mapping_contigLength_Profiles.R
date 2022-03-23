
### @Zygote:~/projects/endoliths/metaGs/community_analyses/master_table

cont_len=read.table("endoliths_contig_length", sep = "\t")#read.table("endoliths_contig_length", sep = "\t")
cnts=read.table("~/projects/endoliths/metaGs/community_analyses/functional_profiles/mapping2faa/endoliths_allContigs_counts_table", header =T, row.names = 1)

#100
cnts$group=ifelse(rownames(cnts) %in% subset(cont_len, V2 >= 100)$V1, "Long", "Short")
gg_100=aggregate(cnts[,1:48], by = list(cnts$group), FUN=sum)
write.table(gg_100, "mappedReads_100", sep = "\t", quote = F, row.names = F)

#150
cnts$group=ifelse(rownames(cnts) %in% subset(cont_len, V2 >= 150)$V1, "Long", "Short")
gg_150=aggregate(cnts[,1:48], by = list(cnts$group), FUN=sum)
write.table(gg_150, "mappedReads_150", sep = "\t", quote = F, row.names = F)

#200
cnts$group=ifelse(rownames(cnts) %in% subset(cont_len, V2 >= 200)$V1, "Long", "Short")
gg_200=aggregate(cnts[,1:48], by = list(cnts$group), FUN=sum)
write.table(gg_200, "mappedReads_200", sep = "\t", quote = F, row.names = F)

#250
cnts$group=ifelse(rownames(cnts) %in% subset(cont_len, V2 >= 250)$V1, "Long", "Short")
gg_250=aggregate(cnts[,1:48], by = list(cnts$group), FUN=sum)
write.table(gg_250, "mappedReads_250", sep = "\t", quote = F, row.names = F)

#500
cnts$group=ifelse(rownames(cnts) %in% subset(cont_len, V2 >= 500)$V1, "Long", "Short")
gg_500=aggregate(cnts[,1:48], by = list(cnts$group), FUN=sum)
write.table(gg_500, "mappedReads_500", sep = "\t", quote = F, row.names = F)

cnts_250=subset(cnts, rownames(cnts) %in% subset(cont_len, V2 >= 250)[,1])
write.table(cnts_250[,-49], "Endoliths_final_contig_count_table_over250bp", row.names = T, quote = F)

cnts_500=subset(cnts, rownames(cnts) %in% subset(cont_len, V2 >= 500)[,1])
write.table(cnts_500[,-49], "Endoliths_final_contig_count_table_over500bp", row.names = T, quote = F)



## export lists
write.table(subset(cont_len, V2 >= 500)[,1], "Contigs_over_500_list", row.names = F)
write.table(subset(cont_len, V2 >= 200)[,1], "Contigs_over_200_list", row.names = F)
write.table(subset(cont_len, V2 >= 250)[,1], "Contigs_over_250_list", row.names = F)
write.table(subset(cont_len, V2 >= 100)[,1], "Contigs_over_100_list", row.names = F)
write.table(subset(cont_len, V2 >= 150)[,1], "Contigs_over_150_list", row.names = F)



#local
setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/master_tables/") 
P2=c("#A8C69F", "#99A1A6")
c100=read.table("mappedReads_100",  sep = "\t" ,  header =T)
c100_l=reshape2::melt(c100, id.vars=c("Group.1"), variable.name = "Sample", value.name = "Abundance")
c150=read.table("mappedReads_150",  sep = "\t" ,  header =T)
c150_l=reshape2::melt(c150, id.vars=c("Group.1"), variable.name = "Sample", value.name = "Abundance")
c200=read.table("mappedReads_200",  sep = "\t" ,  header =T)
c200_l=reshape2::melt(c200, id.vars=c("Group.1"), variable.name = "Sample", value.name = "Abundance")
c250=read.table("mappedReads_250",  sep = "\t" ,  header =T)
c250_l=reshape2::melt(c250, id.vars=c("Group.1"), variable.name = "Sample", value.name = "Abundance")
c500=read.table("mappedReads_500",  sep = "\t" ,  header =T)
c500_l=reshape2::melt(c500, id.vars=c("Group.1"), variable.name = "Sample", value.name = "Abundance")

a=ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Group.1), data = c100_l, stat="identity", position = "fill")  + labs( title = "Cut-off = 100bp",y= "Relative abundance in the metagenome", x="")  +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'right', legend.key = element_blank(), strip.background = element_blank(), legend.title = element_blank()) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=P2) 
b=ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Group.1), data = c150_l, stat="identity", position = "fill")  + labs( title = "Cut-off = 150bp",y= "Relative abundance in the metagenome", x="")  +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'right', legend.key = element_blank(), strip.background = element_blank(), legend.title = element_blank()) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=P2) 
c=ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Group.1), data = c200_l, stat="identity", position = "fill")  + labs( title = "Cut-off = 200bp",y= "Relative abundance in the metagenome", x="")  +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'right', legend.key = element_blank(), strip.background = element_blank(), legend.title = element_blank()) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=P2) 
d=ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Group.1), data = c250_l, stat="identity", position = "fill")  + labs( title = "Cut-off = 250bp",y= "Relative abundance in the metagenome", x="")  +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'right', legend.key = element_blank(), strip.background = element_blank(), legend.title = element_blank()) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=P2) 
e=ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Group.1), data = c500_l, stat="identity", position = "fill")  + labs( title = "Cut-off = 500bp",y= "Relative abundance in the metagenome", x="")  +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'right', legend.key = element_blank(), strip.background = element_blank(), legend.title = element_blank()) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=P2) 

gridExtra::grid.arrange(a,b,c,d,e, ncol=3)


#control band per species
c500_r=sweep(c500[2:49],2,colSums(c500[2:49]),"/")
c500_r=cbind(Group.1=c500[,1], c500_r )
c500_l=reshape2::melt(c500_r, id.vars=c("Group.1"), variable.name = "Sample", value.name = "Abundance")
c500_s=subset(c500_l, Sample %like% "[GP]C[0-6]B" )
c500_s$Sample=factor(c500_s$Sample, levels=c("GC1B", "GC2B", "GC3B", "GC4B", "GC5B", "GC6B", "PC1B", "PC2B", "PC3B", "PC4B", "PC5B", "PC6B"))
ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Group.1), data = c500_s, stat="identity", position = "fill")  + labs( title = "Cut-off = 500bp",y= "Relative abundance in the metagenome", x="")  +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'right', legend.key = element_blank(), strip.background = element_blank(), legend.title = element_blank()) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=P2) 
c500_s$Species=ifelse(c500_s$Sample %like% "^G", "Goniastrea", "Porites")
c500_s2=subset(c500_s, Group.1 == "Long")
boxplot(c500_s2$Abundance~c500_s2$Species,  main = "500 bp", ylim=c(0,1))


#control band per species
c250_r=sweep(c250[2:49],2,colSums(c250[2:49]),"/")
c250_r=cbind(Group.1=c250[,1], c250_r )
c250_l=reshape2::melt(c250_r, id.vars=c("Group.1"), variable.name = "Sample", value.name = "Abundance")
c250_s=subset(c250_l, Sample %like% "[GP]C[0-6]B" )
c250_s$Sample=factor(c250_s$Sample, levels=c("GC1B", "GC2B", "GC3B", "GC4B", "GC5B", "GC6B", "PC1B", "PC2B", "PC3B", "PC4B", "PC5B", "PC6B"))
ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Group.1), data = c250_s, stat="identity", position = "fill")  + labs( title = "Cut-off = 500bp",y= "Relative abundance in the metagenome", x="")  +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'right', legend.key = element_blank(), strip.background = element_blank(), legend.title = element_blank()) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=P2) 
c250_s$Species=ifelse(c250_s$Sample %like% "^G", "Goniastrea", "Porites")
c250_s2=subset(c250_s, Group.1 == "Long")
boxplot(c250_s2$Abundance~c250_s2$Species, main = "250 bp",ylim=c(0,1))

#control band per species
c100_r=sweep(c100[2:49],2,colSums(c100[2:49]),"/")
c100_r=cbind(Group.1=c100[,1], c100_r )
c100_l=reshape2::melt(c100_r, id.vars=c("Group.1"), variable.name = "Sample", value.name = "Abundance")
c100_s=subset(c100_l, Sample %like% "[GP]C[0-6]B" )
c100_s$Sample=factor(c100_s$Sample, levels=c("GC1B", "GC2B", "GC3B", "GC4B", "GC5B", "GC6B", "PC1B", "PC2B", "PC3B", "PC4B", "PC5B", "PC6B"))
ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Group.1), data = c100_s, stat="identity", position = "fill")  + labs( title = "Cut-off = 500bp",y= "Relative abundance in the metagenome", x="")  +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'right', legend.key = element_blank(), strip.background = element_blank(), legend.title = element_blank()) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=P2) 
c100_s$Species=ifelse(c100_s$Sample %like% "^G", "Goniastrea", "Porites")
c100_s2=subset(c100_s, Group.1 == "Long")
boxplot(c100_s2$Abundance~c100_s2$Species, main = "100 bp",ylim=c(0,1))

par(mfrow=c(1,3))
boxplot(c100_s2$Abundance~c100_s2$Species, main = "100 bp",ylim=c(0,1), ylab= "Relative abundance of contigs included", xlab="")
boxplot(c250_s2$Abundance~c250_s2$Species, main = "250 bp",ylim=c(0,1), ylab= "Relative abundance of contigs included", xlab="")
boxplot(c500_s2$Abundance~c500_s2$Species,  main = "500 bp", ylim=c(0,1), ylab= "Relative abundance of contigs included", xlab="")
