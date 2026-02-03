#!/usr/bin/env bash

# Download SILVA database files (release 138.2)
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/taxonomy/tax_slv_ssu_138.2.txt.gz
gunzip tax_slv_ssu_138.2.txt.gz

wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.2.txt.gz
gunzip taxmap_slv_ssu_ref_nr_138.2.txt.gz

wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/taxonomy/tax_slv_ssu_138.2.tre.gz
gunzip tax_slv_ssu_138.2.tre.gz

# Transform to qza
# Import the Taxonomy Rank file:
qiime tools import --type 'FeatureData[SILVATaxonomy]' --input-path tax_slv_ssu_138.2.txt --output-path taxranks-silva-138.2-ssu-nr99.qza
# Import the Taxonomy Mapping file:
qiime tools import --type 'FeatureData[SILVATaxidMap]' --input-path taxmap_slv_ssu_ref_nr_138.2.txt --output-path taxmap-silva-138.2-ssu-nr99.qza
# Import the Taxonomy Hierarchy Tree file:
qiime tools import  --type 'Phylogeny[Rooted]'  --input-path tax_slv_ssu_138.2.tre  --output-path taxtree-silva-138.2-nr99.qza
# Import the sequence file:
qiime tools import --type 'FeatureData[RNASequence]' --input-path SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta --output-path silva-138.2-ssu-nr99-rna-seqs.qza
# Need to reverse transcribe the 16S RNA sequences
qiime rescript reverse-transcribe --i-rna-sequences silva-138.2-ssu-nr99-rna-seqs.qza --o-dna-sequences silva-138.2-ssu-nr99-seqs.qza 

# Prepare the silva taxonomy:
qiime rescript parse-silva-taxonomy --i-taxonomy-tree taxtree-silva-138.2-nr99.qza --i-taxonomy-map taxmap-silva-138.2-ssu-nr99.qza --i-taxonomy-ranks taxranks-silva-138.2-ssu-nr99.qza --o-taxonomy silva-138.2-ssu-nr99-tax.qza --p-include-species-labels
# remove sequences that contain 5 or more ambiguous bases:
qiime rescript cull-seqs --i-sequences silva-138.2-ssu-nr99-seqs.qza --o-clean-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza
# cut sequences Archaea (16S) >= 900 bp, Bacteria (16S) >= 1200 bp, and any Eukaryota (18S) >= 1400 bp (for quality)
qiime rescript filter-seqs-length-by-taxon --i-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza --i-taxonomy silva-138.2-ssu-nr99-tax.qza --p-labels Archaea Bacteria Eukaryota --p-min-lens 900 1200 1400 --o-filtered-seqs silva-138.2-ssu-nr99-seqs-filt.qza --o-discarded-seqs silva-138.2-ssu-nr99-seqs-discard.qza
# Dereplicating in uniq mode
qiime rescript dereplicate --i-sequences silva-138.2-ssu-nr99-seqs-filt.qza --i-taxa silva-138.2-ssu-nr99-tax.qza --p-mode uniq --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-derep-uniq.qza --o-dereplicated-taxa silva-138.2-ssu-nr99-tax-derep-uniq.qza
# Train naive bayes classifier:
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads silva-138.2-ssu-nr99-seqs-derep-uniq.qza --i-reference-taxonomy silva-138.2-ssu-nr99-tax-derep-uniq.qza --o-classifier silva-138.2-ssu-nr99-classifier.qza
