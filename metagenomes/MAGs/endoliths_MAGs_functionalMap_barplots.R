library(lemon)
library(reshape2)
library(ggplot2)
setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/")

#Analyze enrichM data to exclude pathways < 70% completeness
enr=read.table("Outputs/module_completeness.tsv",header =T, sep = "\t", quote = "")
map=read.table("Input_files/metadata_DA_MAGS.txt",header =T, sep = "\t", quote = "")
enr$DA_band=map$DA_band[match(enr$Genome_name, rownames(map))]
enr$DA_skel=map$DA_skel[match(enr$Genome_name, rownames(map))]
enr$Taxa=map$Taxa[match(enr$Genome_name, rownames(map))]

kos=read.table("Input_files/endolith_MAGs_KOs",header = F, row.names = 1)
kos$genome=gsub("_.*", "",rownames(kos))
kos$DA_band=map$DA_band[match(kos$genome, rownames(map))]
kos$DA_skel=map$DA_skel[match(kos$genome, rownames(map))]
kos$Taxa=map$Taxa[match(kos$genome, rownames(map))]

#Pphyla=c("#2077b4","#aec7e8","#ff7f0e", "#ffbb78", "#2ca02c", "#98df89","#17becf", "#9edae5", "#e377c2", "#f7b6d2", "#e33e33",   "#f1c453", "#efea5a", "#16db93" ,  "#B6B6B6")

#### from metabolic

met=t(read.table("Input_files/Metabolic_FinalTable", fill = T, sep = "\t", header = T, row.names = 1))
rownames(met)=gsub(".R_input.txt", "", rownames(met))
met.s=as.data.frame(met[,order(colnames(met))]) #sorting
met.s.f=met.s[,-c(18:22)] #filtering processes out
colnames(met.s.f)=gsub("-S-0", "", colnames(met.s.f))
colnames(met.s.f)=gsub(".*:", "", colnames(met.s.f))
met.s.f$MAG=rownames(met.s.f)
met_long=reshape2::melt(met.s.f)
met_merg=merge(met_long,map, by.x="MAG" ,by.y="row.names")

#### all from metabolic
s_modules=c("Sulfide oxidation", "Sulfite oxidation" ,"Sulfur oxidation", "Sulfate reduction" , "Sulfite reduction", "Sulfur reeduction", "Sulfite oxidation", "Thiosulfate disproportionation","Thiosulfate oxidation" )
n_modules=c( "Nitrogen fixation" ,"Nitrite reduction", "Anammox" , "Ammonia oxidation", "Nitric oxide reduction", "Nitrite oxidation", "Nitrous oxide reduction","Nitrate reduction", "Nitrite ammonification" )
c_modules=c( "Fermentation" ,"Ethanol oxidation","Methanogenesis","Acetate oxidation","Methanotrophy", "Hydrogen generation","Hydrogen oxidation")
phot_modules=c( "M00598" ,"M00597", "M00163" , "M00161")
Cfix_modules=c("M00173", "M00165", "M00377" ,"M00376")

band_photo=subset(enr, Module_id %in% phot_modules & Percent_steps_found >= 70)[,c(1,3,6,7,9)] %>% group_by(DA_band, Taxa, Module_name) %>%filter(!DA_band == "None") %>% tally() 
band_cfix=subset(enr, Module_id %in% Cfix_modules & Percent_steps_found >= 70)[,c(1,3,6,7,9)] %>% group_by(DA_band, Taxa, Module_name) %>%filter(!DA_band == "None") %>% tally()
band_cfix_all=subset(enr, Module_id %in% Cfix_modules & Percent_steps_found >= 70)[,c(1,3,6,7,9)] %>% group_by(DA_band, Taxa) %>%filter(!DA_band == "None") %>% tally()
band_cfix_all$Module_name="Carbon fixation"
band_all=subset(met_merg, variable %in% c_modules | variable %in% n_modules | variable %in% s_modules) %>% filter(!DA_band == "None" ) %>% group_by(DA_band, Taxa, variable) %>% summarise(n=sum(value))
colnames(band_all)[3]="Module_name"
band_final=rbind(band_all, band_cfix, band_photo, band_cfix_all) %>% filter(!Module_name %in% c("Sulfite reduction" , "Anammox" ,"Ethanol oxidation" ,  "Sulfite reduction",  "Ammonia oxidation" , "Nitrite oxidation" , "Methanogenesis", "Photosystem II" ))
colnames(band_final)=c("Species","Taxa","Module_name", "n")
band_final$Compartment="Band"

skel_photo=subset(enr, Module_id %in% phot_modules & Percent_steps_found >= 70)[,c(1,3,6,8,9)] %>% group_by(DA_skel, Taxa, Module_name) %>%filter(!DA_skel == "None") %>% tally() 
skel_cfix=subset(enr, Module_id %in% Cfix_modules & Percent_steps_found >= 70)[,c(1,3,6,8,9)] %>% group_by(DA_skel, Taxa, Module_name) %>%filter(!DA_skel == "None") %>% tally() 
skel_cfix_all=subset(enr, Module_id %in% Cfix_modules & Percent_steps_found >= 70)[,c(1,3,6,8,9)] %>% group_by(DA_skel, Taxa) %>%filter(!DA_skel == "None") %>% tally() 
skel_cfix_all$Module_name="Carbon fixation"
skel_all=subset(met_merg, variable %in% c_modules | variable %in% n_modules | variable %in% s_modules) %>% filter(!DA_skel == "None" ) %>% group_by(DA_skel, Taxa, variable) %>% summarise(n=sum(value)) 
colnames(skel_all)[3]="Module_name"
skel_final=rbind(skel_all, skel_cfix, skel_photo, skel_cfix_all) %>% filter(!Module_name %in% c("Sulfite reduction" , "Anammox" ,"Ethanol oxidation" ,  "Sulfite reduction",  "Ammonia oxidation" , "Nitrite oxidation" , "Methanogenesis"))
colnames(skel_final)=c("Species","Taxa","Module_name", "n")
skel_final$Compartment="Skeleton"

final=rbind(band_final, skel_final)
Pphyla=c("#aec7e8","#ff7f0e", "#ffbb78", "#2ca02c", "#98df89","#17becf", "#9edae5", "#e377c2",  "#e33e33",   "#f1c453", "#efea5a", "#16db93" ,  "#B6B6B6")

pdf("Outputs/barplots_MAGmap.pdf", height = 3.5, width = 20, pointsize = 10)
ggplot() +geom_bar(aes(y = n, x = Species, fill = Taxa), data = final, stat="identity", position = "stack") + 
  geom_col(width = 0.8, position = position_dodge2(width = 0.8)) + labs( y= "", x="") +
  scale_fill_manual(values=Pphyla) +  facet_wrap(Compartment~ Module_name, scales='free_y', nrow = 2 ) + theme_classic() +
  theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'none', axis.text.x=element_blank()) + 
  guides(fill=guide_legend(ncol=1)) + scale_y_continuous( expand = c(0, 0)) 
dev.off()



#measure redundancy
final2=final  %>% group_by(Species, Module_name ) %>% summarise(n=sum(n)) 
final_contst=reshape2::dcast(final2, Module_name ~Species, value.var ="n")
final_contst[is.na(final_contst)]<- 0
message("percentage of processes more redundat in Goniastrea (number): ", nrow(subset(final_contst, Goniastrea > Porites))/nrow(final_contst)*100, "(", nrow(subset(final_contst, Goniastrea > Porites)),")")
final_contst_sub=subset(final_contst, Goniastrea > Porites)
final_contst_sub$ratio=final_contst_sub$Goniastrea/final_contst_sub$Porites
final_contst_sub$ratio=ifelse(final_contst_sub$ratio == "Inf", as.numeric(final_contst_sub$Goniastrea), as.numeric(final_contst_sub$ratio))
mean(final_contst_sub$ratio)

#band
final2=final %>% filter(Compartment == "Band") %>% group_by(Species, Module_name ) %>% summarise(n=sum(n)) 
final_contst=reshape2::dcast(final2, Module_name ~Species, value.var ="n")
final_contst[is.na(final_contst)]<- 0
message("percentage of processes more redundat in Goniastrea (number): ", nrow(subset(final_contst, Goniastrea > Porites))/nrow(final_contst)*100, "(", nrow(subset(final_contst, Goniastrea > Porites)),")")
final_contst_sub=subset(final_contst, Goniastrea > Porites)
final_contst_sub$ratio=final_contst_sub$Goniastrea/final_contst_sub$Porites
final_contst_sub$ratio=ifelse(final_contst_sub$ratio == "Inf", as.numeric(final_contst_sub$Goniastrea), as.numeric(final_contst_sub$ratio))
mean(final_contst_sub$ratio)

#skeleton
final2=final %>% filter(Compartment == "Skeleton") %>% group_by(Species, Module_name ) %>% summarise(n=sum(n)) 
final_contst=reshape2::dcast(final2, Module_name ~Species, value.var ="n")
final_contst[is.na(final_contst)]<- 0
message("percentage of processes more redundat in Porites (number): ", nrow(subset(final_contst, Porites > Goniastrea))/nrow(final_contst)*100, "(", nrow(subset(final_contst,  Porites > Goniastrea)),")")
final_contst_sub=subset(final_contst, Porites > Goniastrea)
final_contst_sub$ratio=final_contst_sub$Porites/final_contst_sub$Goniastrea
final_contst_sub$ratio=ifelse(final_contst_sub$ratio == "Inf", as.numeric(final_contst_sub$Porites), as.numeric(final_contst_sub$ratio))
mean(final_contst_sub$ratio)
