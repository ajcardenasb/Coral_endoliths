
###############
### Binning ###(middlechild)
###############

/home/cardena/software/metabat/jgi_summarize_bam_contig_depths --outputDepth skeletonMetaG_contig_depths  *.bam
/home/cardena/software/metabat/metabat -i ~/Skeleton_diversity/metaGs/assemblies/megahit/final.contigs.fa  -a  skeletonMetaG_contig_depths -o metabat -t 16

#@binnermachine. CheckM is still missing GLIBC_2.23 in Middlechild.
checkm lineage_wf  --tab_table --pplacer_threads 16 -t 16 -x fa -f checkM.txt ./ ./checkM
# 580 cluster were formed.

###############
### Refining ###
###############

# Refining using https://github.com/snayfach/MAGpurify/tree/master/magpurify
export PATH=$PATH:/home/cardena/software/MAGpurify/
export MAGPURIFYDB=/home/cardena/databases/MAGpurify/MAGpurify-db-v1.0/

run_qc.py phylo-markers metabat.100.fa metabat.100_refining
run_qc.py clade-markers metabat.100.fa metabat.100_refining
run_qc.py tetra-freq metabat.100.fa metabat.100_refining
run_qc.py gc-content metabat.100.fa metabat.100_refining
run_qc.py known-contam metabat.100.fa metabat.100_refining
run_qc.py clean-bin metabat.100.fa metabat.100_refining

checkm lineage_wf  --tab_table --pplacer_threads 16 -t 16 -x fna -f checkM.txt ./ ./checkM

# refining for all 232 “refining bins”
cat analysis | while read line ; do for bin in *.fa ; do run_qc.py  $line $bin ${bin//.fa/_refining} -t 32 ; done ; done
ls -1 *.fa > refining_bins
mkdir checkm
cat refining_bins | while read bins ; do cp ./${bins//.fa/}_refining/cleaned_bin.fna ./checkm/${bins} ; done
checkm lineage_wf  --tab_table --pplacer_threads 32 -t 32 -x fa -f checkM.txt ./ ./checkM
scp cardena@10.254.145.13:/home/cardena/checkm-metaG-Skeleton/refining/second_refining/checkM.txt .
cut -f 1,2,12,13 checkem_ref5.txt | tr " " "_"> checkm_ref5_ed.txt
 mkdir ../sixth_refining/
cat need_more_refining | while read bins ; do cp ${bins}.fa ../../sixth_refining/${bins}.fa ; done
#### after first refining 9 new bins were processed as final bins and 223 were refined. After the second refining, none was treated as final fin but 66 were discarded as have completeness < 75, 157 were refined again.

### refining in refineM @MiddeChild:/home/cardena/Skeleton_diversity/metaGs/refineM
Trying it now on a subset (~/checkm-metaG-Skeleton/refining/sixth_refininf/subset/). When I tried for all 48 bins it ran for 4 days and then crashed, dont know exactly why
#You’ll need the stats output from checkm, so first, one more run
checkm lineage_wf  --tab_table --pplacer_threads 32 -t 32 -x fa -f checkM.txt ./ ./checkM
#Remove outliers based on genomic signatures
## calculate scaffold coverage with refineM. Remember to call bam files instead of the folder
export PATH=/home/cardena/miniconda2/bin:$PATH
source activate refinem
refinem scaffold_stats -c 32 ../../../final.contigs.fa . ./checkM/lineage.ms ../../../bams/*bam -x fa
refinem outliers <stats_output_dir>/scaffold_stats.tsv <outlier_output_dir>
refinem filter_bins <bin_dir> <outlier_output_dir>/outliers.tsv <filtered_output_dir>

#Remove outliers based on taxonomic assignments
refinem call_genes -c 40 . refinem_gene-output -x fa
refinem taxon_profile -c 40 refinem_gene-output <stats_output_dir>/scaffold_stats.tsv <reference_db> <reference_taxonomy> <taxon_profile_output_dir>


### refining in MAGpy/MAGwash ( cardena@rsse:~/checkm-metaG-Skeleton/refining/sixth_refininf/subset)
#https://github.com/WatsonLab/MAGpy/blob/master/install.md
#export PATH=/home/cardena/anaconda3/bin:$PATH
#conda activate snakemake
#This seems to be connected to the pipeline below :(  trying to get that ready then

##Plotting and other tools for MAG exploration
#export PATH=/home/cardena/anaconda2/bin:$PATH
#source activate magpy_install
#Only missing last step https://github.com/WatsonLab/MAGpy/blob/master/install.md


###################
### Final stats ###
###################

export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate checkm

### Final 75 bins
checkm lineage_wf  --tab_table --pplacer_threads 56 -t 56 -x fasta -f checkM.txt . checkM_out

## calculate bin coverage using checkM, https://github.com/Ecogenomics/CheckM/issues/129
export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate checkm
checkm coverage ~/projects/endoliths/metaGs/mags/final_bins/ coverage_endoliths /home/cardena/projects/endoliths/metaGs/binning/*.bam -x fasta -t 128

#2. Calculate profiles
checkm profile --tab_table -f profile_endoliths  coverage_endoliths
#3. Calculate extended stats
checkm qa -o 2 -t 56 -f extended_statistics_endoliths  --coverage_file coverage_endoliths ../final_bins/checkM_out/lineage.ms ../final_bins/checkM_out/
