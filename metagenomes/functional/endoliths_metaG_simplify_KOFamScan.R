setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/functional/")

# Rscript simplifying_KofamScan_output.R KEGG_table output
library(data.table)
args=commandArgs(TRUE)
kegg=read.table(args[1], sep = "\t") #"Input_files/test"
kegg.f=subset(kegg, V6 < 0.001)[,c(2:3,6)]
smallesEval=setDT(kegg.f)[, .SD[which.min(V6)], V2][,c(1:2)]
write.table(smallesEval, args[2], row.names = F, col.names = F, quote = F, sep = "\t") 

