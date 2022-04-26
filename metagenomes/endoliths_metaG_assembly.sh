#########################
### Trimming adaptors ###
#########################

for i in *_R1_001.fastq.gz ; do java -jar /home/cardena/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 16 $i ${i//_R1_/_R2_} -baseout ./trimmed/$(basename $i | sed 's/_R1_001.fastq.gz/.fastq/'g) ILLUMINACLIP:/home/cardena/software/Trimmomatic-0.38/adapters/Truseq_V.3_edited.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50; done
for i in GC[1234][BS]*1P.fastq*_R1_001.fastq.gz ; do java -jar /home/cardena/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 16 $i ${i//_R1_/_R2_} -baseout ./trimmed/$(basename $i | sed 's/_R1_001.fastq.gz/.fastq/'g) ILLUMINACLIP:/home/cardena/software/Trimmomatic-0.38/adapters/Truseq_V.3_edited.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50; done


#renaming files
for i in *P.fastq ; do mv $i $(basename $i | sed 's/M_19_1[12][0-9][0-9]_//' | sed 's/_[A-Z][0][0-9]_L00[1-7]//') ; done
for i in *.fastq.gz ; do cp $i $(basename $i | sed 's/M_19_1[12][0-9][0-9]_//' | sed 's/_[A-Z][0][0-9]_L00[1-7]//' | sed 's/_001//' | sed 's/_R1/_1P/' | sed 's/_R2/_2P/')  ; done

#################
### Assembly ###(middlechild)
#################

/home/cardena/software/megahit/megahit -t 16 -m 0.9 -1 GC1B_1P.fastq,GC1S_1P.fastq,GC2B_1P.fastq,GC2S_1P.fastq,GC3B_1P.fastq,GC3S_1P.fastq,GC4B_1P.fastq,GC4S_1P.fastq,GC5B_1P.fastq,GC5S_1P.fastq,GC6B_1P.fastq,GC6S_1P.fastq,GH1B_1P.fastq,GH1S_1P.fastq,GH2B_1P.fastq,GH2S_1P.fastq,GH3B_1P.fastq,GH3S_1P.fastq,GH4B_1P.fastq,GH4S_1P.fastq,GH5B_1P.fastq,GH5S_1P.fastq,GH6B_1P.fastq,GH6S_1P.fastq,PC1B_1P.fastq,PC1S_1P.fastq,PC2B_1P.fastq,PC2S_1P.fastq,PC3B_1P.fastq,PC3S_1P.fastq,PC4B_1P.fastq,PC4S_1P.fastq,PC5B_1P.fastq,PC5S_1P.fastq,PC6B_1P.fastq,PC6S_1P.fastq,PH1B_1P.fastq,PH1S_1P.fastq,PH2B_1P.fastq,PH2S_1P.fastq,PH3B_1P.fastq,PH3S_1P.fastq,PH4B_1P.fastq,PH4S_1P.fastq,PH5B_1P.fastq,PH5S_1P.fastq,PH6B_1P.fastq,PH6S_1P.fastq -2 GC1B_2P.fastq,GC1S_2P.fastq,GC2B_2P.fastq,GC2S_2P.fastq,GC3B_2P.fastq,GC3S_2P.fastq,GC4B_2P.fastq,GC4S_2P.fastq,GC5B_2P.fastq,GC5S_2P.fastq,GC6B_2P.fastq,GC6S_2P.fastq,GH1B_2P.fastq,GH1S_2P.fastq,GH2B_2P.fastq,GH2S_2P.fastq,GH3B_2P.fastq,GH3S_2P.fastq,GH4B_2P.fastq,GH4S_2P.fastq,GH5B_2P.fastq,GH5S_2P.fastq,GH6B_2P.fastq,GH6S_2P.fastq,PC1B_2P.fastq,PC1S_2P.fastq,PC2B_2P.fastq,PC2S_2P.fastq,PC3B_2P.fastq,PC3S_2P.fastq,PC4B_2P.fastq,PC4S_2P.fastq,PC5B_2P.fastq,PC5S_2P.fastq,PC6B_2P.fastq,PC6S_2P.fastq,PH1B_2P.fastq,PH1S_2P.fastq,PH2B_2P.fastq,PH2S_2P.fastq,PH3B_2P.fastq,PH3S_2P.fastq,PH4B_2P.fastq,PH4S_2P.fastq,PH5B_2P.fastq,PH5S_2P.fastq,PH6B_2P.fastq,PH6S_2P.fastq -o ../assemblies/megahit
#Before mapping remove spaces and special characters from headers

### SPAdes
export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate  spades

cat /home/cardena/projects/endoliths/metaGs/trimmed/*1P.fastq > all_1P.fastq
cat /home/cardena/projects/endoliths/metaGs/trimmed/*2P.fastq > all_2P.fastq
spades.py -1 all_1P.fastq  -2 all_2P.fastq -k 127 -t 64 -m 900 -o spades_assembly --meta

#Spades assembly failed error code 255
#check if they have the same pairs
head -1051 FRIDec16A_C10_decon_R1.fastq >test_R1.fastq
head -1076 FRIDec16A_C10_decon_R2.fastq >test_R2.fastq

/home/cardena/software/fastq-pair/build/fastq_pair all_1P.fastq all_2P.fastq

awk '{s++}END{print s/4}' all_1P.fastq
awk '{s++}END{print s/4}' all_1P.fastq

#https://github.com/rrwick/Unicycler/issues/152
### #IDBA
#IDBA-UD requires paired-end reads stored in single FastA file and a pair of reads is in consecutive two lines.
export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate idba

for i in /home/cardena/projects/endoliths/metaGs/trimmed/*1P.fastq; do fq2fa --merge $i ${i//1P.fastq/2P.fastq} $(basename $i | sed 's/1P.fastq/_IDBA_input.fasta/g')  ; done
for i in /home/cardena/projects/endoliths/metaGs/trimmed/GC[1234][BS]_merged.fastq.extendedFrags.fastq ; do fq2fa  $i  $(basename $i | sed 's/merged.fastq.extendedFrags.fastq/_IDBA_input.fasta/g')  ; done

cat *_IDBA_input.fasta > IDBA_input.fasta
idba_ud -r IDBA_input.fasta --mink 80 --maxk 140 --num_threads 64 -o assembly_idba_ut

#m_threads 64 -o assembly_idba_ut
#number of threads 64
#reads 2,190,124,857
#long reads 0
#extra reads 0
#read_length 151
#kmer 80

#Quast
export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate quast
/home/cardena/software/quast-5.0.2/quast.py ./megahit/final.contigs.fa ./idba_ut/scaffold.fa ./spades/scaffolds.fasta -o ./quast


###############
### Mapping ###(middlechild)
###############

bowtie2-build ~/Skeleton_diversity/metaGs/assemblies/megahit/final.contigs.fa skeletonMetaG_DB
for i in ~/Skeleton_diversity/metaGs/trimmed/*1P.fastq; do bowtie2 -x skeletonMetaG_DB  -1 $i -2 ${i//1P.fastq/2P.fastq} -S $(basename $i | sed 's/_1P.fastq/.sam/') -p 16 ; done
for i in *.sam ; do samtools view -bS $i  | samtools sort > ${i//sam/bam} ; done
for i in *.bam ; do samtools index  $i ; done
