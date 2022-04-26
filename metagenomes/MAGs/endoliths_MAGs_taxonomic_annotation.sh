
#########################
### Phylogenomic tree ###
#########################

# Phylogenetic tree using UBCG pipeline on complete and near complete bins @BinnerMachine:/home/cardena/checkm-metaG-Skeleton/tree
#ran in BinnerMachine, zygote had problems with this particular jar


#Converting fasta files to bcg files using prodigal and hmmsearch tools.
cat 75MAGs | while read bin ; do java -jar UBCG.jar extract -i ./endoliths/${bin}.fasta -bcg_dir ./endoliths_ubcg_out -label ${bin}  ; done
#Generating multiple alignments from bcg files
java -jar UBCG.jar align -bcg_dir ./endoliths_ubcg_out  -out_dir ./endoliths_ubcg_out2 -t 32 -raxml -prefix Phylogenomic_endoliths1
#transferring to zygote
scp -r endoliths_ubcg_out* cardena@134.34.126.43:/home/cardena/projects/endoliths/metaGs/mags/phylo_tree/ubcg_out/

The final tree marked with GSI was written to './endoliths_ubcg_out2/Phylogenomic_endoliths1/Phylogenomic_endoliths1.UBCG_gsi(92).codon.50.label.nwk'


############################
### Taxonomic annotation ###
############################
## GTDB
export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate gtdbtk
export GTDBTK_DATA_PATH=/home/cardena/databases/GTDB/release95/

 gtdbtk classify_wf --genome_dir . --out_dir ./output_GTDB --cpus 64 -x fa


#########################
### Metabolic profile ###
#########################

# 1. Metabolic profiling using METABOLIC pipeline on complete and near complete bins @BinnerMachine:/home/cardena/checkm-metaG-Skeleton/metabolic
#Run Prokka first to have final faa files

export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate prokka
for i in *fa ; do mv ${i} ${i//metabat/Endolith} ; done
cat 75MAGs | while read line ; do  prokka --prefix ${line} --force --outdir ${line} --locustag ${line} --cpus 132 --rawproduct ~/projects/endoliths/metaGs/mags/final_bins/${line}.fasta ; done

#Run METABOLIC, https://github.com/AnantharamanLab/METABOLIC/blob/master/README.md
##Perl libraries werent properly installed but still got the output!
export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate metabolic
perl /home/cardena/software/METABOLIC/METABOLIC-C.pl -in-gn . -t 56 -o metabolic_75 -m /home/cardena/software/METABOLIC/
cd /metabolic_75/R_input/
python ~/middleChild/scripts/MergingMultipleTables.py . .
scp cardena@134.34.126.43:/home/cardena/projects/endoliths/metaGs/mags/metabolic/metabolic_75/R_input/FinalTable .
