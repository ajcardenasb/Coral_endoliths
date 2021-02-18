

## prepare for enrichM
setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/")
kos=read.table("Input_files/endolith_MAGs_KOs",header = F, row.names = 1)
kos$genome=gsub("_.*", "",rownames(kos))
kos_matrix=reshape2::dcast(kos, V2~genome)
write.table(kos_matrix,"Outputs/input_for_enrichM", quote = F, row.names = F, sep = "\t")

#Analyze enrichM data to explide pathways < 80% completeness
enr=read.table("Outputs/module_completeness.tsv",header =T, sep = "\t", quote = "")
map=read.table("Input_files/metadata_DA_MAGS.txt",header =T, sep = "\t", quote = "")
enr$DA_band=map$DA_band[match(enr$Genome_name, rownames(map))]
enr$DA_skel=map$DA_skel[match(enr$Genome_name, rownames(map))]
enr$Taxa=map$Taxa[match(enr$Genome_name, rownames(map))]

## carbon fixation
subset(enr, Module_id == "M00173" & Percent_steps_found >= 70)[,c(1,3,6,7,9)]%>% arrange(DA_band, Taxa)#Reductive citrate cycle (Arnon-Buchanan cycle) 
subset(enr, Module_id == "M00165" & Percent_steps_found >= 70)[,c(1,3,6,7,9)] %>% arrange(DA_band, Taxa)#Calvin cycle
subset(enr, Module_id == "M00377" & Percent_steps_found >= 70)[,c(1,3,6,7,9)] %>% arrange(DA_band, Taxa)#Wood-Ljungdahl pathway
subset(enr, Module_id == "M00376" & Percent_steps_found >= 70)[,c(1,3,6,7,9)] %>% arrange(DA_band, Taxa)#3-hydroxypropionate bicycle

subset(enr, Module_id == "M00173" & Percent_steps_found >= 70)[,c(1,3,6,8,9)]%>% arrange(DA_skel, Taxa)#Reductive citrate cycle (Arnon-Buchanan cycle) 
subset(enr, Module_id == "M00165" & Percent_steps_found >= 70)[,c(1,3,6,8,9)] %>% arrange(DA_skel, Taxa)#Calvin cycle
subset(enr, Module_id == "M00377" & Percent_steps_found >= 70)[,c(1,3,6,8,9)] %>% arrange(DA_skel, Taxa)#Wood-Ljungdahl pathway
subset(enr, Module_id == "M00376" & Percent_steps_found >= 70)[,c(1,3,6,8,9)] %>% arrange(DA_skel, Taxa)#3-hydroxypropionate bicycle

## photosynthesis
subset(enr, Module_id == "M00598" & Percent_steps_found >= 70)[,c(1,3,6,7,9)]%>% arrange(DA_band, Taxa)#Anoxygenic photosystem I
subset(enr, Module_id == "M00597" & Percent_steps_found >= 70)[,c(1,3,6,7,9)] %>% arrange(DA_band, Taxa)#Anoxygenic photosystem II
subset(enr, Module_id == "M00163" & Percent_steps_found >= 70)[,c(1,3,6,7,9)] %>% arrange(DA_band, Taxa)#	Photosystem I
subset(enr, Module_id == "M00161" & Percent_steps_found >= 70)[,c(1,3,6,7,9)] %>% arrange(DA_band, Taxa)#Photosystem II

subset(enr, Module_id == "M00598" & Percent_steps_found >= 70)[,c(1,3,6,8,9)]%>% arrange(DA_skel, Taxa)#Anoxygenic photosystem I
subset(enr, Module_id == "M00597" & Percent_steps_found >= 70)[,c(1,3,6,8,9)] %>% arrange(DA_skel, Taxa)#Anoxygenic photosystem II
subset(enr, Module_id == "M00163" & Percent_steps_found >= 70)[,c(1,3,6,8,9)] %>% arrange(DA_skel, Taxa)#	Photosystem I
subset(enr, Module_id == "M00161" & Percent_steps_found >= 70)[,c(1,3,6,8,9)] %>% arrange(DA_skel, Taxa)#Photosystem II

## metanogenesis
subset(enr, Module_id == "M00567" & Percent_steps_found >= 70)[,c(1,3,6,7,9)]%>% arrange(DA_band, Taxa)#CO2 to methane
subset(enr, Module_id == "M00356" & Percent_steps_found >= 70)[,c(1,3,6,7,9)]%>% arrange(DA_band, Taxa)#methanol to methane
subset(enr, Module_id == "M00357" & Percent_steps_found >= 70)[,c(1,3,6,7,9)]%>% arrange(DA_band, Taxa)#acetate to methane 
subset(enr, Module_id == "M00358" & Percent_steps_found >= 70)[,c(1,3,6,7,9)]%>% arrange(DA_band, Taxa)#CoM biosynthesis, serves as a methyl group carrier in key reactions within the pathway of methane formation from C1 precursors

subset(enr, Module_id == "M00567" & Percent_steps_found >= 70)[,c(1,3,6,8,9)]%>% arrange(DA_skel, Taxa)#CO2 to methane
subset(enr, Module_id == "M00356" & Percent_steps_found >= 70)[,c(1,3,6,8,9)] %>% arrange(DA_skel, Taxa)#methanol to methane
subset(enr, Module_id == "M00357" & Percent_steps_found >= 70)[,c(1,3,6,8,9)] %>% arrange(DA_skel, Taxa)#	acetate to methane 
subset(enr, Module_id == "M00358" & Percent_steps_found >= 70)[,c(1,3,6,8,9)] %>% arrange(DA_skel, Taxa)#	CoM biosynthesis, serves as a methyl group carrier in key reactions within the pathway of methane formation from C1 precursors

## metanothrophy

subset(enr, Module_id == "M00346" & Percent_steps_found >= 70)[,c(1,3,6,7,9)]%>% arrange(DA_band, Taxa)#serine pathway
subset(enr, Module_id == "M00345" & Percent_steps_found >= 70)[,c(1,3,6,7,9)]%>% arrange(DA_band, Taxa)#ribulose monophosphate pathway 
subset(enr, Module_id == "M00344" & Percent_steps_found >= 70)[,c(1,3,6,7,9)]%>% arrange(DA_band, Taxa)#xylulose monophosphate pathway
subset(enr, Module_id == "M00563" & Percent_steps_found >= 70)[,c(1,3,6,7,9)]%>% arrange(DA_band, Taxa)#trimethylamine metabolism

subset(enr, Module_id == "M00346" & Percent_steps_found >= 70)[,c(1,3,6,8,9)]%>% arrange(DA_skel, Taxa)#serine pathway
subset(enr, Module_id == "M00345" & Percent_steps_found >= 70)[,c(1,3,6,8,9)] %>% arrange(DA_skel, Taxa)#ribulose monophosphate pathway 
subset(enr, Module_id == "M00344" & Percent_steps_found >= 70)[,c(1,3,6,8,9)] %>% arrange(DA_skel, Taxa)#	xylulose monophosphate pathway
subset(enr, Module_id == "M00563" & Percent_steps_found >= 70)[,c(1,3,6,8,9)] %>% arrange(DA_skel, Taxa)#	trimethylamine metabolism

