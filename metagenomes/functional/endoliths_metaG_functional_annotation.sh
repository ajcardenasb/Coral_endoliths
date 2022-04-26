
#cardena@134.34.126.43:/home/cardena/projects/endoliths/metaGs/community_analyses/functional_profiles/

#########################
#### ORF prediction #####
#########################
export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate prodigal # this includes hmmer version 3
prodigal -a endoliths.assembly.faa -d endoliths.assembly.fna  -i ../endoliths.assembly.fa -o endoliths.assembly.gff -p meta

######################################
#### quantification using salmon ##### @~/projects/endoliths/metaGs/community_analyses/functional_profiles/mapping2faa
######################################

## all contigs
salmon index -t ../prodigal/endoliths.assembly.fna -i endoliths_assembly_all_salmonDB -p 56
for fastq in ~/projects/endoliths/metaGs/trimmed/*_1P.fastq ; do  salmon quant -i endoliths_assembly_all_salmonDB  -l A -1 ${fastq} -2 ${fastq//1P/2P}  -p 56 -o $(basename $fastq |  sed 's/_1P.fastq/_mapping_salmon/') --gcBias  ; done
for i in ./*_salmon/quant.sf ; do cut -f1,5 $i > ./counts/$(dirname "$i" | sed 's/.\///g' | sed 's/_mapping_salmon//'); done
python /home/cardena/scripts/MergingMultipleTables.py . .
for i in ./*_salmon/quant.sf ; do cut -f1,4 $i > ./TPM/$(dirname "$i" | sed 's/.\///g' | sed 's/_mapping_salmon//') ; done
python /home/cardena/scripts/MergingMultipleTables.py . .


#################
### KOfamscan ### @~/projects/endoliths/metaGs/community_analyses/functional_profiles/KOfamSca
#################
export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate  KOfamscan
## need to set parameters on /home/cardena/software/kofam_scan-1.3.0/config.yml
### the profile is defined on the profile file as:   and he ko_list as : ko_list_map01120
/home/cardena/software/kofam_scan-1.3.0/exec_annotation -f detail-tsv  -o endoliths_fullKEGG_KOfamScan_allContigs  /home/cardena/projects/endoliths/metaGs/community_analyses/functional_profiles/prodigal/endoliths.assembly.faa

export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate R.v4.0.2
### simplify KEGG annotations by choosing the minimum E-val pero contig
Rscript ~/scripts/simplifying_KofamScan_output.R endoliths_fullKEGG_KOfamScan_allContigs endoliths_fullKEGG_annotations_allContigs


#######################################################
### create master table before binning ORFs into KOS ## @~/projects/endoliths/metaGs/community_analyses/master_table
#######################################################
export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate bioawk

bioawk -c fasta '{ print $name, length($seq) }' < ~/projects/endoliths/metaGs/community_analyses/functional_profiles/prodigal/endoliths.assembly.fna > endolits_contig_length


### bin ORFs into KOs
Rscript  ~/scripts/ORFs2modules.R ~/projects/endoliths/metaGs/community_analyses/functional_profiles/mapping2faa/endoliths_allContigs_counts_table ~/projects/endoliths/metaGs/community_analyses/functional_profiles/mapping2faa/endoliths_allContigs_TPM_table endoliths_fullKEGG_annotations_allContigs
## get all KOs for C and N metabolism
wget  http://rest.kegg.jp/link/ko/map01120  -O map01200
cat map01200 | cut -f2 | sed 's/ko://' > KOs_map01200
wget http://rest.kegg.jp/link/ko/map00910 -O map00910
cat map00910 | cut -f2 | sed 's/ko://' > KOs_map00910
### bin ORFs into  C and N KOs
Rscript  ~/scripts/ORFs2CNKOs.R ~/projects/endoliths/metaGs/community_analyses/functional_profiles/mapping2faa/endoliths_allContigs_counts_table endoliths_fullKEGG_annotations_allContigs
### bin ORFs into kegg pathways
Rscript ~/scripts/ORFs2pathways.R  ~/projects/endoliths/metaGs/community_analyses/functional_profiles/mapping2faa/endoliths_allContigs_counts_table endoliths_fullKEGG_annotations_allContigs
### bin ORFs into metabolims KOs
Rscript ~/scripts/ORFs2metabolKOs.R ~/projects/endoliths/metaGs/community_analyses/functional_profiles/mapping2faa/endoliths_allContigs_counts_table endoliths_fullKEGG_annotations_allContigs ~/projects/endoliths/metaGs/community_analyses/functional_profiles/mapping2faa/endoliths_allContigs_counts_table



####################
### EnrichM ### @/home/cardena/projects/roseobacters/genomes/enrichM
#################

export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate enrichm
export ENRICHM_DB=/home/cardena/databases/enrichm_db/enrichm_database_v10/
#run rscript endoliths_prepare_for_enrichM
#scp input_for_enrichM cardena@134.34.126.43:/home/cardena/projects/endoliths/metaGs/mags/enrichM/

enrichm classify --genome_and_annotation_matrix input_for_enrichM
#scp cardena@134.34.126.43:~/projects/endoliths/metaGs/community_analyses/functional_profiles/enrichM/2021-02-03_15-42-enrichm_classify_output/module_completeness.tsv


############################################
### Photosystems, methanogens and RuBisCO ### ~/projects/endoliths/metaGs/community_analyses/functional_profiles/PFAM_hmms/
############################################

export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate prodigal # this includes hmmer version 3

for i in *hmm ; do hmmsearch  --noali --cut_ga --cpu 64 --tblout ${i//.hmm/_results} $i ~/projects/endoliths/metaGs/community_analyses/functional_profiles/prodigal/endoliths.assembly.faa ; done
### running wiht updated families
hmmsearch  --noali --cut_ga --cpu 64 --tblout  PF10657_results PF10657.hmm ~/projects/endoliths/metaGs/community_analyses/functional_profiles/prodigal/endoliths.assembly.faa
hmmsearch  --noali --cut_ga --cpu 64 --tblout  PF10643_results PF10643.hmm ~/projects/endoliths/metaGs/community_analyses/functional_profiles/prodigal/endoliths.assembly.faa

## cholobi pscA
grep "O07091" ../blastp/endoliths_assembly_sprot.outfmt6 | cut -f1 > pscA_uniprot

#Heliobacteria and Firmicutes pshA
grep "Q3IC07\|Q3ID16" ../blastp/endoliths_assembly_sprot.outfmt6 | cut -f1 > pshA_uniprot

  cat *_results  | grep "k1"  > endoliths_photosystems_contigs


cat PF00223_results PF00796_results PF07465_results PF05479_results PF03244_results PF02605_results PF02507_results PF02427_results PF01701_results PF01241_results > P1_Results
cut -f1 -d " " P1_Results | grep "k1" | sort | uniq > P1_contigs

cat PF07123_results PF04725_results PF02468_results PF18240_results PF13326_results PF06596_results PF06514_results PF06298_results PF05969_results PF05151_results PF02533_results PF02532_results PF02419_results PF01789_results PF01788_results PF01405_results PF00737_results PF00421_results > P2_Results
cut -f1 -d " " P2_Results | grep "k1" | sort | uniq > P2_contigs

cat PF00016_results PF00101_results PF02788_results > rubisco_Results
cut -f1 -d " " rubisco_Results | grep "k1" | sort | uniq > rubisco_contigs


##methanogens
export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate prodigal # this includes hmmer version 3
cat tigr_accessions | while read line ; do hmmsearch  --noali --cut_ga --cpu 64 --tblout ${line}_results /home/cardena/databases/HMM_models/TIGRFAM/hmms/${line}.HMM ~/projects/endoliths/metaGs/community_analyses/functional_profiles/prodigal/endoliths.assembly.faa ; done
cat *results | grep "k141" > endoliths_methanogenesis_contigs
scp cardena@134.34.126.43:/home/cardena/projects/endoliths/metaGs/community_analyses/functional_profiles/hmms_photo_rubisco/methanogenesis/endoliths_methanogenesis_contigs ~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/functional/Input_files/

## Carbon Fixers
export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate prodigal # this includes hmmer version 3
cat rubisco_Accessions | while read line ; do hmmsearch  --noali --cut_ga --cpu 64 --tblout ${line//.hmm/}_results ${line} ~/projects/endoliths/metaGs/community_analyses/functional_profiles/prodigal/endoliths.assembly.faa ; done
cat *results | grep "k141" > endoliths_rubisco_contigs


############################################
### Carbohydrate and amino acid metabolism ### ~/projects/endoliths/metaGs/community_analyses/functional_profiles/
############################################

## see categories
cut -f1 FOAM-onto_rel1.tsv | sort | uniq

#using FOAM structure cardena@zygote:~/databases/HMM_models/FOAM
grep "09_Carbohydrate" FOAM-onto_rel1.tsv | cut -f1,2,5 > KOs_Carbohydrate_Active_enzyme_CAZy
grep "04_Utililization" FOAM-onto_rel1.tsv | cut -f1,2,5 > KOs_Sugar_metabolism
grep "12_Transporters" FOAM-onto_rel1.tsv | cut -f1,2,3,4,5 > KOs_transporters
grep "06_Amino" FOAM-onto_rel1.tsv | cut -f1,2,5 > KOs_aminoacids

scp cardena@134.34.126.43:~/databases/HMM_models/FOAM/KOs_* ~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/metagenomes/functional/Input_files/
