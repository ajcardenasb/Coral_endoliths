
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

#2. photosyntesis proteins were annotated from Prokka
#grep "Photosystem I" ./*/*.tsv | grep -v "Photosystem II" | grep -v "Ycf3" | cut -f1 -d ":" | sort | uniq # 26 MAGs have any protein. none have psaA, just cyanobacteria
#grep "Photosystem II" ./*/*.tsv | grep -v "Ycf3" | cut -f1 -d ":" | sort | uniq # 3 MAGs have any protein
grep "bch"  ./*/*.tsv | cut -f1 -d ":" | sort | uniq # 45 MAGs have any protein related to bacteriochlorophyll synthesis
grep "puf[LM]"  ./*/*.tsv | cut -f1 -d ":" | sort | uniq # 2 MAGs  have any protein related to photosynthetic unit forming call_genes
grep "psc"  ./*/*.tsv | cut -f1 -d ":" | sort | uniq #Photosystem1-like gene in chlorobi (pscC gene)
grep "psa"  ./*/*.tsv | cut -f1 -d ":" | sort | uniq # 30 PS1
grep "psb"  ./*/*.tsv | cut -f1 -d ":" | sort | uniq # 3 PS2
grep "crt[EBIYZX]"  ./*/*.tsv | cut -f1 -d ":" | sort | uniq # 43 MAGs encode for genes associated to carotenoids

# 2. using KOfamscan + metabolic to add methanogens
export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate  KOfamscan
cat list | while read line ;  do /home/cardena/software/kofam_scan-1.3.0/exec_annotation -f detail-tsv  -o ${line}_KOfamScan  ~/projects/endoliths/metaGs/mags/prokka/${line}/${line}.faa ; done
for i in *_KOfamScan; do Rscript ~/scripts/simplifying_KofamScan_output.R $i ${i}_simplified ; done

for i in *_simplified ; do grep -f methanogenesis_KOs $i  > ${i//_KOfamScan_simplified/_methanogenesistKOs} ; done
cat  *_methanogenesistKOs > MAGs_methanogenesis_results
scp cardena@134.34.126.43:~/projects/endoliths/metaGs/mags/KOfamscan/MAGs_methanogenesis_results ~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/Input_files/

# 3. from hmms
export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate prodigal # this includes hmmer version 3
for i in ~/projects/endoliths/metaGs/mags/prokka/*/*faa ; do hmmsearch  --noali --cut_ga --cpu 64 --tblout $(basename $i | sed 's/.faa/_PF10643_results/g')  PF10643.hmm $i ; done
for i in ~/projects/endoliths/metaGs/mags/prokka/*/*faa ; do hmmsearch  --noali --cut_ga --cpu 64 --tblout $(basename $i | sed 's/.faa/_TIGR01157_results/g')  TIGR01157.hmm $i ; done
for i in ~/projects/endoliths/metaGs/mags/prokka/*/*faa ; do hmmsearch  --noali --cut_ga --cpu 64 --tblout $(basename $i | sed 's/.faa/_TIGR01115_results/g')  TIGR01115.hmm $i ; done
for i in ~/projects/endoliths/metaGs/mags/prokka/*/*faa ; do hmmsearch  --noali --cut_ga --cpu 64 --tblout $(basename $i | sed 's/.faa/_TIGR01335_results/g')  TIGR01335.hmm $i ; done
for i in ~/projects/endoliths/metaGs/mags/prokka/*/*faa ; do hmmsearch  --noali --cut_ga --cpu 64 --tblout $(basename $i | sed 's/.faa/_TIGR01336_results/g')  TIGR01336.hmm $i ; done
for i in ~/projects/endoliths/metaGs/mags/prokka/*/*faa ; do hmmsearch  --noali --cut_ga --cpu 64 --tblout $(basename $i | sed 's/.faa/_PF00223_results/g')  PF00223.hmm $i ; done
for i in ~/projects/endoliths/metaGs/mags/prokka/*/*faa ; do hmmsearch  --noali --cut_ga --cpu 64 --tblout $(basename $i | sed 's/.faa/_TIGR01151_results/g')  TIGR01151.hmm $i ; done
for i in ~/projects/endoliths/metaGs/mags/prokka/*/*faa ; do hmmsearch  --noali --cut_ga --cpu 64 --tblout $(basename $i | sed 's/.faa/_TIGR03039_results/g')  TIGR03039.hmm $i ; done
for i in ~/projects/endoliths/metaGs/mags/prokka/*/*faa ; do hmmsearch  --noali --cut_ga --cpu 64 --tblout $(basename $i | sed 's/.faa/_PF00421_results/g')  PF00421.hmm $i ; done

cat *_results  | grep Endolith | grep -v "#"


##################################
### HMMs for energy metabolism ### ~/projects/endoliths/metaGs/mags/energy_metabolism
##################################
grep "^P" all_hmm > pfam_hmm
grep "^K" all_hmm > kegg_hmm
grep "^T" all_hmm > tigr_hmm

cat kegg_hmm | while read line ; do cp ~/databases/KOFam/profiles/${line}.hmm . ; done
cat pfam_hmm | while read line ; do wget https://pfam.xfam.org/family/${line}/hmm -O ${line}.hmm ; done
cat tigr_hmm | while read line ; do cp ~/databases/HMM_models/TIGRFAM/hmms/${line}.HMM ${line}.hmm ; done

cat *.hmm > energy_metabolism.hmm
for i in ~/projects/endoliths/metaGs/mags/prokka/*/*faa ; do hmmsearch  --noali --cut_ga --cpu 64 --tblout $(basename $i | sed 's/.faa/_energy_metabolism_results/g')  energy_metabolism.hmm $i ; done
cat *_energy_metabolism_results  | grep Endolith | grep -v "#" > All_MAGs_energy_metabolism_results
scp cardena@134.34.126.43:~/projects/endoliths/metaGs/mags/energy_metabolism/All_MAGs_energy_metabolism_results ~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/MAGs/Input_files/
