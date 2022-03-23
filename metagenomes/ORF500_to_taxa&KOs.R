#@Zygote:~/projects/endoliths/metaGs/community_analyses/master_table
library(tidyr)
cnts=read.table("Endoliths_final_contig_count_table_over500bp", header =T, row.names = 1)
cnts=round(cnts, digits = 0)
kai=read.table("endoliths_allContigs_KaijuOut_mem_withNames", sep = "\t", quote = "", fill = T)
kai=separate(kai, V4, sep = ";", into = c('Superkingdom','Phylum','Class','Order','Family', 'genus', 'species'), remove = FALSE)
kegg=read.table("endoliths_fullKEGG_annotations_allContigs", sep = "\t")


#fix with unclassified all fields
unique(kai$Superkingdom)
kai$Superkingdom=sub("^$", "Unclassified",kai$Superkingdom)
kai$Superkingdom=ifelse(kai$Superkingdom == "NA" | is.na(kai$Superkingdom) , "Unclassified", as.character(kai$Superkingdom))
unique(kai$Superkingdom)

kai$Family=trimws(kai$Family, which = "both" , whitespace = "[ \t\r\n]")
kai$Order=trimws(kai$Order, which = "both", whitespace = "[ \t\r\n]")
kai$Phylum=trimws(kai$Phylum, which = "both", whitespace = "[ \t\r\n]")
kai$Class=trimws(kai$Class, which = "both", whitespace = "[ \t\r\n]")

kai$label=ifelse(kai$Family == "NA" | is.na(kai$Family), as.character(paste("Unclassified", kai$Order, sep = "_")), as.character(kai$Family))
kai$label=ifelse(kai$label == "Unclassified_NA", as.character(paste("Unclassified", kai$Class, sep = "_")), as.character(kai$label))
kai$label=ifelse(kai$label == "Unclassified_NA", as.character(paste("Unclassified", kai$Phylum, sep = "_")), as.character(kai$label))
kai$label=ifelse(kai$label == "Unclassified_NA", as.character(paste("Unclassified", kai$Superkingdom, sep = "_")), as.character(kai$label))
unique(kai$label)


cnts$Family=kai$label[match(rownames(cnts), kai$V2)]
cnts$Family=ifelse(is.na(cnts$Family), "Unclassified", as.character(cnts$Family))
cnts$Superkingdom=kai$Superkingdom[match(rownames(cnts), kai$V2)]
cnts$Superkingdom=ifelse(is.na(cnts$Superkingdom), "Unclassified", as.character(cnts$Superkingdom))

### ORF to tax and kos
#Domain
all_dom=aggregate(cnts[,1:48], by = list(cnts$Superkingdom), FUN=sum)
write.table(all_dom, "Superkingdom_counts_500", sep = "\t", quote = F, row.names = F)

#Families
all_fams=aggregate(cnts[,1:48], by = list(cnts$Family), FUN=sum)
write.table(all_fams, "Families_counts_500", sep = "\t", quote = F, row.names = F)

#Kos
cnts$KO=kegg$V2[match(rownames(cnts), kegg$V1)]
cnts$KO=ifelse(is.na(cnts$KO), "Unclassified", as.character(cnts$KO))
all_KOs=aggregate(cnts[,1:48], by = list(cnts$KO), FUN=sum)
write.table(all_KOs, "KO_counts_500", sep = "\t", quote = F, row.names = F)

tpm=read.table("~/projects/endoliths/metaGs/community_analyses/functional_profiles/mapping2faa/endoliths_allContigs_TPM_table", header =T, row.names = 1)
tpm=round(tpm, digits = 0)
tpm$KO=kegg$V2[match(rownames(tpm), kegg$V1)]
tpm$KO=ifelse(is.na(tpm$KO), "Unclassified", as.character(tpm$KO))
tpm_KOs=aggregate(tpm[,1:48], by = list(tpm$KO), FUN=sum)
write.table(tpm_KOs, "KO_tpm_500", sep = "\t", quote = F, row.names = F)


#metabolism Kos
kegg2=read.table("~/scripts/hierarchy_ko00001", header= T, sep="\t", fill =T)
pre_meta_kos=all_KOs
pre_meta_kos$Group.1=ifelse(pre_meta_kos$Group.1 %in% subset(kegg2, L1 == " Metabolism")$KO, as.character(pre_meta_kos$Group.1), "Others")
meta_KOs=aggregate(pre_meta_kos[,2:49], by = list(pre_meta_kos$Group.1), FUN=sum)
write.table(meta_KOs, "metabolic_KO_counts_500", sep = "\t", quote = F, row.names = F)



#export kaiju taxonomy
final_kaiju=unique(kai[,-c(1,2,3)])
write.table(final_kaiju,  "Endoliths_kaiju_taxonomy_500", sep = "\t", quote = F, row.names = F )

### QC
sum(all_dom$GC1B)
sum(all_fams$GC1B)
sum(all_KOs$GC1B)
sum(meta_KOs$GC1B)





