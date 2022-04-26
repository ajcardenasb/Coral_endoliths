# Repository - Coral endoliths

This repository contains the scripts used to analyze data and create figures for the manuscript "Greater functional diversity and redundancy of coral endolithic microbiomes align with lower coral bleaching susceptibility"

Raw sequencing data are deposited in the NCBI Sequence Read Archive (SRA) under BioProject [PRJNA757245](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA757245)

## Workflow

### Metabaracoding - 16S rRNA data analysis
1. Amplicon Sequence variance (ASV) were inferred using [dada2](https://github.com/benjjneb/dada2) using the script `/16S/endoliths_16S_dada2.R`
2. Quality checks (i.e., removal of putatively contaminant ASVs and removal of samples with < 1000 reads) were done using the script `/16S/endoliths_16S_contaRemoval.R`
3. Bar plots of most abundant bacterial taxa were created using the script`/16S/endoliths_16S_barplots.R`
4. Ordination plots and PERMANOVAs were done using [Vegan](https://github.com/vegandevs/vegan) with the script `/16S/endoliths_16S_ordination_permanovas.R`
5. Alpha diversity estimates and statistical comparisons were done using [Vegan](https://github.com/vegandevs/vegan) following the script `/16S/endoliths_16S_alphaDiversity.R`
6. Differentially abundance analysis was done using [ANCOMBC](https://github.com/FrederickHuangLin/ANCOMBC) following the script `/16S/endoliths_16S_ANCOMBC.R`
7. Venn diagrams were plotted to compare metagenomic with metabarcoding results using the script `/16S/endoliths_16S_vs_metaG_vennDiagram.R`


### Metagenomics

#### Assembly, ORF prediction, gene quantification and annotation.
1. Assembly, ORF prediction and estimation of gene counts were done using the script `/metagenomes/endoliths_metaG_assembly.sh`
2. Assess contig length cutoffs to keep for downstream analysis 'ORFs_master_to_taxa_KOs' and `/metagenomes/Mapping_contigLength_Profiles.R`
3. Create KO and KEGG module count matrix removing contigs < 500 bp `/metagenomes/ORF500_to_taxa&KOs.R`

#### community-based metagenomic taxonomic profiles
5. ORF-based taxonomic annotation was done using [kaiju](https://github.com/bioinformatics-centre/kaiju) as described in `/metagenomes/taxonomy/endoliths_metaG_taxonomic_annotation.sh`
6. Relative abundance of top families per superkingdom were plotted using  'Barplots_topFamilies.R'
7. Ordination plots of all families were done using 'All_families_ordination_plots.R'
8. Differentially abundant families were determined by ANCOM-BC using the scripts '/metagenomes/taxonomy/ANCOMBC_tax.R'
9. Plots to represent differentially abundant famies were done using the scripts ''

#### community-based metagenomic functional profiles
10. ORF-based functional annotation was done using KOFamScan as described in `/metagenomes/functional/endoliths_metaG_functional_annotation.sh` and `/metagenomes/functional/endoliths_metaG_simplify_KOFamScan.R`
11. Ordination plots and PERMANOVAs of KEGG orthologs were done using the script `metagenomes/functional/endoliths_metaG_functional_ordination_permanovas.R`
12. KEGG orthologs alpha diversity was calculated usign the script `/metagenomes/functional/endoliths_metaG_functional_alphaDiversity.R`
13. Differentially abundant KEGG orthologs and modules were determined by ANCOM-BC using the scripts `/metagenomes/functional/endoliths_metaG_ANCOMBC_L3.R` `/metagenomes/functional/endoliths_metaG_ANCOMBC_KOs.R`
14. Selbal analysis to calculate biomarker features for each treatment/species was done using the script `/metagenomes/functional/endoliths_metaG_selbal.R`
15. Plots to represent differentially abundant KOs and modules were done using the scripts `/metagenomes/functional/endoliths_metaG_panel1.R`, `/metagenomes/functional/endoliths_metaG_panel2.R` `/metagenomes/functional/endoliths_metaG_plotting_ANCOM_results_treatments.R`
16. The script `/metagenomes/functional/endoliths_metaG_Cfixation.R` was used to plot carbon fixation pathways across species

### genomic-centric approaches (MAGs)
15. Metagenomic binning was done using the script `/metagenomes/MAGs/endoliths_MAGs_binning.sh`
16. MAG taxonomic annotation and phylogenomic tree were done using the script `/metagenomes/MAGs/endoliths_MAGs_taxonomic_annotation.sh`
17. 'endoliths_binning_stats.R' Plots relative abundance of MAGs by phylum in relation to the entire community and also in relation to the binned fraction
18. 'ggtree_endoliths.R' and 'metaG_phyloTree_functions.sh' plots the phylogenomic inference
19. 'endoliths_MAGs_ANCOMBC_betwSpec.R' identifies DA MAGs between species and treatments using ANCOM-BC and plots heatmaps comparing MAGs abundance between species
20. 'MAG_photo_and_auto.R' identifies photosynthesis and carbon fixation genes based on KO annotations
21. 'endoliths_MAGs_plotting_metabolic_potential.R' plots heatmaps with metabolic potential of MAGs
