setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/")
library(pheatmap)
library(reshape2)
library(dplyr)
library(gplots)

map=read.table("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/endoliths_metadata.txt", header = TRUE, row.names = 1, sep ='\t')
#meta.mag= read.table("Input_files/tree_metadata.txt", header = T, row.names = 1)
pht=read.table("Outputs/MAGs_complementary_annotations_prokka", sep = "\t", header = T)
tax.mag= read.table("Input_files/metadata_DA_MAGS.txt", header = T, row.names = 1, sep = "\t")
meta.mag_sorted=tax.mag[ order(tax.mag$DA_band, tax.mag$DA_skel), ] #sorting by enrichment
#meta.mag_sorted=tax.mag[ order(tax.mag$Taxa), ] #sorting by taxa

################################################
#### MAGs metabolic %1 ####
################################################

class=read.table("Input_files/MAGs_class.txt", fill = T, sep = "\t", header = T, row.names = 1)
class$Labels=paste(rownames(class), class$Class, sep = " - ")
MAG_labels=t(class[-1])
met=t(read.table("Input_files/Metabolic_FinalTable", fill = T, sep = "\t", header = T, row.names = 1))
rownames(met)=gsub(".R_input.txt", "", rownames(met))
met.s=as.data.frame(met[,order(colnames(met))]) #sorting
met.s.f=met.s[,-c(18:22)] #filtering processes out
colnames(met.s.f)=gsub("-S-0", "", colnames(met.s.f))
colnames(met.s.f)=gsub(".*:", "", colnames(met.s.f))
column_colors=meta.mag_sorted[,-4]
met_phy=merge(pht,met.s.f, by.x ="MAG", by.y="row.names")
rownames(met_phy)=met_phy[,1]
met_phy=met_phy[,-c(1:3)]

metabolic_final=met_phy[match(rownames(meta.mag_sorted), rownames(met_phy)), ]
#pdf(file =  "./outputs/MAG_metabolic_potential.pdf",  width = 12, height = 12, pointsize = 12) 
pheatmap(metabolic_final, labels_row= MAG_labels, color = colorRampPalette(c("white", "black"))(50),  cellwidth = 8, cellheight =  5, fontsize_row = 6, fontsize_col= 8, legend = F,   cluster_rows = T,cluster_col = F, annotation_row   = column_colors, annotation_colors =  list(Taxa = c(Acidobacteriota="#2077b4", Bacteroidota="#aec7e8",Chloroflexota="#ff7f0e", Cyanobacteria="#ffbb78", Desulfobacterota="#2ca02c", Firmicutes="#98df89",Gemmatimonadota="#17becf",KSB1= "#9edae5", Myxococcota="#e377c2", Omnitrophota="#f7b6d2", Planctomycetota="#e33e33", Proteobacteria= "#f1c453", SAR324 = "#efea5a", Spirochaetota = "#16db93" , Verrucomicrobiota= "#B6B6B6"), DA_band =c(Goniastrea="#8f2d56",Porites=  "#218380", None = "white"),  DA_skel =c(Goniastrea="#8f2d56",Porites=  "#218380", None = "white")))
#dev.off()

#extract data
long_metabolic=melt(metabolic_final)%>% filter(value == "1" )
long_metabolic_complete=merge(long_metabolic,tax.mag, by.x ="Var2", by.y="row.names")

process="Carbon fixation"
long_metabolic_complete %>% filter(Var1 == process & DA_band == "Goniastrea") %>% arrange(Taxa)
long_metabolic_complete %>% filter(Var1 == process & DA_band == "Porites") %>% arrange(Taxa)
long_metabolic_complete %>% filter(Var1 == process & DA_skel == "Goniastrea") %>% arrange(Taxa)
long_metabolic_complete %>% filter(Var1 == process & DA_skel == "Porites") %>% arrange(Taxa)

#measure redundancy
redundancy=data.frame(Process=unique(long_metabolic_complete$Var1))
gon_redu=long_metabolic_complete %>% filter(DA_band == "Goniastrea" |  DA_skel == "Goniastrea") %>% group_by(Var1) %>% tally() 
por_redu=long_metabolic_complete %>% filter(DA_band == "Porites" |  DA_skel == "Porites") %>% group_by(Var1) %>% tally() 
redundancy$Goniastrea=gon_redu$n[match(redundancy$Process, gon_redu$Var1)]
redundancy$Porites=por_redu$n[match(redundancy$Process, por_redu$Var1)]
redundancy[is.na(redundancy)]<- 0
redundancy$ratio=redundancy$Goniastrea/redundancy$Porites
redundancy_sub=subset(redundancy, ratio > 1 & !ratio =="Inf")
mean(redundancy_sub$ratio)
# redundancy$Goniastrea=(redundancy$Goniastrea*100)/15
# redundancy$Porites=(redundancy$Porites*100)/15
# redundancy_l=melt(redundancy)
# boxplot(redundancy_l$value~redundancy_l$variable)

##############################
#### metabolic potential 2 ####
#############################
pro=read.table("Input_files/profile_endoliths", fill = T, sep = "\t", header = T, row.names = 1)
bin=pro[, grepl(pattern = "binned", x = colnames(pro))][c(1:75),]
names(bin)=gsub("\\..*", "", names(bin))
names(bin)=gsub("[0-9]", "", names(bin))
bin.mean=t(apply(bin, 1, function(x) tapply(x, colnames(bin), mean)))
bin.mean_filt=bin.mean[,c(1,5,2,6)]

met_phy2=met_phy[,7:ncol(met_phy)]
met_phy2$`Anoxygenic photosynthesis`=met_phy$RCI+met_phy$RCII
met_phy2$`Photosynthesis`=met_phy$PS1+met_phy$PS2
pr.ab=lapply(met_phy2, function(x) ifelse(x >= 1, "present", "absent")) %>% bind_rows()

all=cbind(bin.mean_filt, pr.ab)
combined=lapply(all[, 5:29], function(x) aggregate(all[, 1:4], by = list(x), FUN =  sum)[2,]) %>% bind_rows()
pathways=names(all[, 5:29])
combined=combined[,-1]
rownames(combined)=pathways
combined.f=combined[complete.cases(combined), ]

#pdf(file =  "./outputs/MAG_metabolism_abundances.pdf",  width = 7, height = 7, pointsize = 12) 
pheatmap(combined.f, color = colorRampPalette(c("white", "gray", "black"))(50),  angle_col = "90",  scale = "row",  cluster_row= T, cluster_col = FALSE, border_color = NA,  cellwidth = 8, cellheight =  6, fontsize_row = 6, fontsize_col= 8, legend = T, gaps_col = c(2))
#dev.off()
