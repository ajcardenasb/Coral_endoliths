# Repository - Coral endoliths

This repository contains the scripts used to analyze data and create figures for the manuscript "Greater functional diversity and redundancy of coral endolithic microbiomes align with lower coral bleaching susceptibility"

Raw sequencing data are deposited in the NCBI Sequence Read Archive (SRA) under BioProject [PRJNA757245](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA757245)

## Workflow

### Metabarcoding - 16S rRNA data analysis
1. Amplicon Sequence variance (ASV) were inferred using [dada2](https://github.com/benjjneb/dada2) using the script `/16S/endoliths_16S_dada2.R`
2. Quality checks (i.e., removal of putatively contaminant ASVs and removal of samples with < 1000 reads) were done using the script `/16S/endoliths_16S_contaRemoval.R`
3. Bar plots of most abundant bacterial taxa were created using the script`/16S/endoliths_16S_barplots.R`
4. Ordination plots and PERMANOVAs were done using [Vegan](https://github.com/vegandevs/vegan) with the script `/16S/endoliths_16S_ordination_permanovas.R`
5. Alpha diversity estimates and statistical comparisons were done using [Vegan](https://github.com/vegandevs/vegan) following the script `/16S/endoliths_16S_alphaDiversity.R`
6. Differentially abundance analysis was done using [ANCOMBC](https://github.com/FrederickHuangLin/ANCOMBC) following the script `/16S/endoliths_16S_ANCOMBC.R`
7. Venn diagrams were plotted to compare metagenomic with metabarcoding results using the script `/16S/endoliths_16S_vs_metaG_vennDiagram.R`


### Metagenomics

#### Assembly, ORF prediction, gene quantification and annotation.
1. Assembly, ORF prediction and estimation of gene counts were done using the script `endoliths_metaG_assembly.sh`
2. Assess contig length cutoffs to keep for downstream analysis 'ORFs_master_to_taxa_KOs' and `Mapping_contigLength_Profiles`
3. Create KO and KEGG module count matrix removing contigs < 250 bp 'ORF250_to_taxa&KOs.R'

#### community-based metagenomic taxonomic profiles
5. ORF-based taxonomic annotation was done using [kaiju](https://github.com/bioinformatics-centre/kaiju) as described in 'metaG_community_taxonomy.sh'
6. Relative abundance of top families per superkingdom were plotted using  'Barplots_topFamilies.R'
7. Ordination plots of all families were done using 'All_families_ordination_plots.R'
8. Differentially abundant families were determined by ANCOM-BC using the scripts '/metagenomes/taxonomy/ANCOMBC_tax.R'
9. Plots to represent differentially abundant famies were done using the scripts ''

#### community-based metagenomic functional profiles
10. ORF-based functional annotation was done using KOFamScan as described in 'metaG_community_funProfiles' and 'simplifying_KofamScan_output'
11. Ordination plots of KEGG orthologs were done using 'endoliths_KEGG_ordination'
12.  Differentially abundant KEGG orthologs and modules were determined by ANCOM-BC using the scripts ''
13. Plots to represent differentially abundant KOs and modules were done using the scripts '' 'endoliths_plotting_ANCOM_KOS.R
14. 'Photosystems_RuBisCO_TPMs.R' -> boxplots of normalized abundance of contigs annotated as Photosystems I and II and RuBisCO

### genomic-centric approaches (MAGs)
15. binning done using 'metaG_binning.sh'
15. 'endoliths_binning_stats.R' Plots relative abundance of MAGs by phylum in relation to the entire community and also in relation to the binned fraction
16. 'ggtree_endoliths.R' and 'metaG_phyloTree_functions.sh' plots the phylogenomic inference
18. 'endoliths_MAGs_ANCOMBC_betwSpec.R' identifies DA MAGs between species and treatments using ANCOM-BC and plots heatmaps comparing MAGs abundance between species
19. 'MAG_photo_and_auto.R' identifies photosynthesis and carbon fixation genes based on KO annotations
20. 'endoliths_MAGs_plotting_metabolic_potential.R' plots heatmaps with metabolic potential of MAGs
