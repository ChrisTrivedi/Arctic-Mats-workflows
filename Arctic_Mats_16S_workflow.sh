#!/bin/bash

## 2018.02.21
## QIIME/VSearch workflow of Arctic Mat samples from 2014 & 2017
## Concatenating together two data sets and reclustering.

## Running from /Users/Chris/GoogleDrive/Working_Qiime/20180221_Arctic_Mats

## First I need to do Split Libraries for the NanoSeq run and then I can join the two SL files together.

# Validate mapping file
# To test for duplicate barcodes.
validate_mapping_file.py -o Mapping_Validate -m 20171018_Arctic_18_Parada_Master_Map.txt

# Join paired ends - about 3 mins
time join_paired_ends.py -f Raw/BS-L1_S1_L001_Run4_R1.fastq -r Raw/BS-L1_S1_L001_Run4_R2.fastq -o seq/Joined/ -j 50 -v

# Use PEAR to join paired ends - 4 min - almost 100 MB more joined
#pear -f Raw/BS-L1_S1_L001_Run4_R1.fastq -r #Raw/BS-L1_S1_L001_Run4_R2.fastq -o seq/Joined/PEAR -p 0.001 -v #50 -m 450 -n 250 -y 500m -j 4

# Extract barcodes using the mapping file - 1 min
time extract_barcodes.py -f seq/Joined/fastqjoin.join.fastq -m 20171018_Arctic_18_Parada_Master_Map.txt -l 12 -o seq/Prepped -a -v

# Split libraries (separate out samples from fastq) - 2 mins
time split_libraries_fastq.py --store_demultiplexed_fastq --phred_quality_threshold 20 -i seq/Prepped/reads.fastq -b seq/Prepped/barcodes.fastq -m 20171018_Arctic_18_Parada_Master_Map.txt --barcode_type 12 -o seq/SlOut/ -v

## Concatenante Split Libraries files from the BFP14_16_17 set and the NanoSeq 2018 run.

# Concatenate SL files
# seqs.fastq
cat seq/SlOut_14_16_17/seqs_14_16_17.fastq seq/SlOut/seqs.fastq > seq/SlOut_All_Arctic/seqs_Mats.fastq

# seqs.fna
cat seq/SlOut_14_16_17/seqs_14_16_17.fna seq/SlOut/seqs.fna > seq/SlOut_All_Arctic/seqs_Mats.fna

mkdir seq/VsearchOut

# Get quality stats - 30 sec
vsearch -fastq_stats seq/SlOut_All_Arctic/seqs_Mats.fastq --log seq/VsearchOut/seqs_stats.log

# Remove low quality reads (trimming not required for paired-end data) - 15 sec
time vsearch -fastx_filter seq/SlOut_All_Arctic/seqs_Mats.fna -fastaout seq/VsearchOut/seqs_filtered.fasta --fastq_maxee 0.5 --threads 4

# Dereplicate seqs - 12 sec
time vsearch -derep_fulllength seq/VsearchOut/seqs_filtered.fasta --output seq/VsearchOut/seqs_filtered_derep.fasta --sizeout --minuniquesize 2 --threads 4

# Reference chimera check - 4 mins
time vsearch -uchime_ref seq/VsearchOut/seqs_filtered_derep.fasta --db /macqiime/DB/gold.fasta --strand plus --nonchimeras seq/VsearchOut/seqs_filtered_derep_nochimeras.fasta --threads 4

# Cluster OTUs @ 97% - 1 min
time vsearch -cluster_fast seq/VsearchOut/seqs_filtered_derep_nochimeras.fasta --centroids seq/VsearchOut/seqs_filtered_derep_nochimeras_repset.fasta --sizein --xsize --relabel OTU_ --id 0.97 --threads 4

# Make an otus folder
mkdir otus/

# Copy this file to use as RepSet at a later time
cp seq/VsearchOut/seqs_filtered_derep_nochimeras_repset.fasta otus/RepSet.fna

# Map the original quality filtered reads back to OTUs - 11 mins
time vsearch -usearch_global seq/VsearchOut/seqs_filtered.fasta --db seq/VsearchOut/seqs_filtered_derep_nochimeras_repset.fasta --strand plus --id 0.97 -uc  seq/VsearchOut/OTU_map.uc --threads 4

# Modify OTU table for input into QIIME - 30 sec
python /macqiime/bin/uc2otutab.py seq/VsearchOut/OTU_map.uc > seq/VsearchOut/seqs_filtered_derep_nochimeras_repset_OTU-table.txt

# Convert to HDF5 biom type - instant
biom convert --table-type="OTU table" -i seq/VsearchOut/seqs_filtered_derep_nochimeras_repset_OTU-table.txt --to-hdf5 -o otus/Arctic_mats.biom

# Summarize BIOM table to check general stats - instant
biom summarize-table -i otus/Arctic_mats.biom -o otus/Arctic_mats_Summary.txt

# QIIME - Assign taxonomy - 1 min
time assign_taxonomy.py -t /macqiime/SILVA/taxonomy/taxonomy_all/97/taxonomy_7_levels.txt -r /macqiime/SILVA/rep_set/rep_set_all/97/97_otus.fasta -i otus/RepSet.fna -o otus/TaxonomyOut/ -v

# Add taxonomy to BIOM table - instant
echo “Adding Metadata”
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp otus/TaxonomyOut/RepSet_tax_assignments.txt -i otus/Arctic_mats.biom -o otus/Arctic_mats_otuTable.biom

# Computing Summaries of Taxa
echo "Computing Summaries"
biom summarize-table -i otus/Arctic_mats_otuTable.biom -o otus/Arctic_mats_otuTable_Summary.txt

# Align seqs (default: pynast - can also do in parallel) - 4 mins
time align_seqs.py -i otus/RepSet.fna -t /macqiime/SILVA/core_alignment/core_alignment_SILVA123.fasta -o otus/RepSet_Aligned/ -v

# Filter alignment to certain quality - 12 mins
time filter_alignment.py -i otus/RepSet_Aligned/RepSet_aligned.fasta -o otus/RepSet_Aligned/ -e 0.01

# Make tree file (.tre) for downstream use - 4 mins
time make_phylogeny.py -i otus/RepSet_Aligned/RepSet_aligned_pfiltered.fasta -o otus/RepSet_Aligned/rep_set.tre -l otus/RepSet_Aligned/tree_log.txt

# Filter OTU table and convert for analysis in R
filter_samples_from_otu_table.py -i otus/Arctic_mats_otuTable.biom -o otus/Arctic_mats_otuTable_filtered.biom --sample_id_fp Arctic_mat_IDs.txt

# Sort OTU table by sample type (sort Mapping file first & then use SortOrder)
sort_otu_table.py -i otus/Arctic_mats_otuTable_filtered.biom -o otus/Arctic_mats_otuTable_filtered_sorted.biom -m 20180221_Arctic_mats_Mapping_sorted.txt -s SortOrder

# Summarize BIOM table to check general stats - instant
biom summarize-table -i otus/Arctic_mats_otuTable_filtered_sorted.biom -o otus/Arctic_mats_sorted_Summary.txt

# Summarize Taxa/plots on this new OTU table and then we will filter for our diversity analyses (i.e. heatmap, BetaDiversity, etc.) - 5 min
time summarize_taxa_through_plots.py -i otus/Arctic_mats_otuTable_filtered_sorted.biom -o Analyses/Arctic_mats_TaxaSummary -v







