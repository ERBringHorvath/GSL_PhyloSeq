#!/bin/zsh

source ~/.zshrc

conda activate anvio

##Repeat workflow for each genomic analysis
#reformat FASTA files for anvio
for f in *.fna
do
anvi-script-reformat-fasta $f --simplify-names -o ${f%.fna}.renamed.fa -r ${f%.fna}.contig.names.txt
done

#generate contigs databases
for f in *.renamed.fa
do
anvi-gen-contigs-database --num-threads 8 -f $f -o ${f%.renamed.fa}.contigs.db --skip-mindful-splitting –
  project-name "$f"
done

#populate databases
for f in *.contigs.db
do
anvi-run-ncbi-cogs -T 6 -c $f
anvi-run-hmms -T 6 -c $f
anvi-run-scg-taxonomy -T 6 -c $f
done

#gen genomes storage; replace FILE with appropriate name
anvi-gen-genomes-storage -e FILE.txt -o FILE_bins-GENOMES.db

#generate pangenome, excluding partial gene calls for cleaner downstream analysis; replace FILE with appropriate name
anvi-pan-genome -g FILE_bins-GENOMES.db --project-name “NAME” --num-threads 6 --mcl-inflation 6 
    --exclude-partial-gene-calls

#extract all singleton sequences; replace FILE with appropriate name
anvi-get-sequences-for-gene-clusters -g FILE_bins-GENOMES.db -p FILE/FILE-PAN.db 
  --min-num-genomes-gene-cluster-occurs 1 --max-num-genomes-gene-cluster-occurs 1 --min-functional-homogeneity-index 0.75 -o All_Genes.fasta

#move into new directory; split master multi-FASTA file into single files
cat All_Genes.fasta | grep '>' | split -d -l 1 -a 5

##remove all files that aren’t GSL17-111 or GSL17-113
for f in *
do
grep -l ‘S_’ $f | xargs rm -f #All non-GSL genomes start with 'S_'
done

