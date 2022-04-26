#######################################
### Taxonomic profiles on all contigs ###
#######################################

/home/cardena/software/kaiju/bin/kaiju -a mem -z 64 -t ~/databases/kaiju_nr_euk_2020-05-25/nodes.dmp -f ~/databases/kaiju_nr_euk_2020-05-25/kaiju_db_nr_euk.fmi -i /home/cardena/projects/endoliths/metaGs/community_analyses/functional_profiles/prodigal/endoliths.assembly.fna  -o endoliths_allContigs_KaijuOut_mem
/home/cardena/software/kaiju/bin/kaiju2table -t ~/databases/kaiju_nr_euk_2020-05-25/nodes.dmp -n  ~/databases/kaiju_nr_euk_2020-05-25/names.dmp  -o endoliths_allContigs_KaijuOut_mem_full_report -r family -l superkingdom,phylum,class,order,family,genus,species  endoliths_allContigs_KaijuOut_mem -e

#needed for master table
/home/cardena/software/kaiju/bin/kaiju-addTaxonNames -t ~/databases/kaiju_nr_euk_2020-05-25/nodes.dmp -n  ~/databases/kaiju_nr_euk_2020-05-25/names.dmp -i endoliths_allContigs_KaijuOut_mem -o endoliths_allContigs_KaijuOut_mem_withNames -r superkingdom,phylum,class,order,family,genus,species


/home/cardena/software/kaiju/bin/kaiju-addTaxonNames -t ~/databases/kaiju_nr_euk_2020-05-25/nodes.dmp -n  ~/databases/kaiju_nr_euk_2020-05-25/names.dmp -i endoliths_allContigs_KaijuOut_mem -o endoliths_allContigs_KaijuOut_mem_withNames_full -p
/home/cardena/software/kaiju/bin/kaiju-addTaxonNames -t ~/databases/kaiju_nr_euk_2020-05-25/nodes.dmp -n  ~/databases/kaiju_nr_euk_2020-05-25/names.dmp -i endoliths_allContigs_KaijuOut_mem -o endoliths_allContigs_KaijuOut_mem_withNames_full2 -r superkingdom,kingdom,phylum,class,order,family,genus,species
