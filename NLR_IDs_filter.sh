#!/bin/env bash

set -xeuo pipefail # this is "bash safe mode", which catches many errors. use it!

#This script is used to filter out the NLR-Integrated domains from the PfamScan output <- It is general purpose script that can extract IDs for every species!
#Things to fix
#USERS should substitute the species basenanme in line # 78-81, and species name as prefix of major NLR clads in line 109-112
#The user must provide the right path to the inputs and output directories of their analysis. Putting everything in one location is prefferred. .

# USAGE: bash script "fasta file of genome"

#Use fasta file of genome as first input
genome=$1

#Reuse the genome as a variable

echo ${genome}

genomebase=$(basename ${genome} .fa)

echo ${genomebase}

#################################################################################################################
#THE FILLOWING COMMAND LINES ARE USED TO POOLOUT THE PROTEIN SEQUENCES OF NLR-INTEGRATED DOMAINS FOR EACH SPECIES
##################################################################################################################

# First, Sarris et al., 2016 Pfam-parser scripts (https://github.com/krasileva-group/plant_rgenes) were used as follows to parse the PfamScan outputs 
# to extracts all domains for each proteins and removes redundant nested hits with larger e-values

perl ~/analysis/nlr_annotation_scripts/K-parse_Pfam_domains_v3.1.pl \
       --pfam ${genomebase}.protein.fa_pfamscan.txt --evalue 0.001 \
       --output ${genomebase}.protein.fa_pfamscan_parsed.verbose --verbose T

# Parsing the output of PfamScan output parser using the script "K-parse_Pfam_domains_NLR-fusions-v2.4.1.pl"
# The script generates:
# Summary of number of NLRs and NLR-IDs
# Summary of integrated domains with species list for each domain
# Abundace list of IDs counted once for each family
# Contingency table per ID domain for each species
# Reference: Sarris et al., 2016

mkdir -p ${genomebase}_pfam-parser
mv ${genomebase}.protein.fa_pfamscan_parsed.verbose ./${genomebase}_pfam-parser/

perl ~/analysis/nlr_annotation_scripts/K-parse_Pfam_domains_NLR-fusions-v2.4.2.pl \
        --indir ./${genomebase}_pfam-parser/ --evalue 0.001 \
        -o ./${genomebase}_pfam-parser/ \
        -d db_descriptions.txt

# Create a directory to store integrted domain filters for each species and NLR-archtectures
mkdir -p ${genomebase}_NLR_IDs
                                                                                                                                                                                                                   # Remove comment lines and the first default blank line that come after grepping the comment lines
grep -v '#' ../${genomebase}.protein.fa_pfamscan.txt | awk '/^$/ && !f{f=1;next}1' \
        > ./${genomebase}_NLR_IDs/${genomebase}_augustusprotein_pfamscan.txt

# Filter out the NLR PfamScan outputs based on the filtered uniq NLR-integred domains from Pfam-parser script
awk 'NR==FNR{for (i=1;i<=NF;i++) a[$i];next} FNR==1 || ($7 in a)' \
       ../${genomebase}_pfam-parser/*_wordcloud_*.txt ./${genomebase}_NLR_IDs/${genomebase}_augustusprotein_pfamscan.txt \
        | sort | uniq > ./${genomebase}_NLR_IDs/${genomebase}_protein_nlrid_pfamscan.txt

# Extract uniq gene-id of NLR-IDs of either from NLR_annotation of a species that help to pool
# out their protein sequences
cat ./${genomebase}_NLR_IDs/${genomebase}_protein_nlrid_pfamscan.txt | cut -d ' ' -f1 \
        | sort | uniq > ./${genomebase}_NLR_IDs/${genomebase}_protein_nlrid_geneid.uniq.txt

# Pool out the protein sequences of IDs using the extracted corresponding uniq gene-id
seqtk subseq ../${genomebase}.NB-ARC_hmmsearch_perseqhit_proteinseq.fa \
	./${genomebase}_NLR_IDs/${genomebase}_protein_nlrid_geneid.uniq.txt \
        > ../../NLR_ID_proteins/${genomebase}_NLR_uniqid_protein.fa

# Filter out uniq NLR-integrated domains for each species with their corresponding Pfam-domain IDs
awk '{print $1 " "$6 " " $7}' ./${genomebase}_NLR_IDs/${genomebase}_protein_nlrid_pfamscan.txt \
        > ./${genomebase}_NLR_IDs/${genomebase}_protein_nlrid_pfamid_multoccurrence.txt

# THE FOLLOWING SCRIPT IS TO DETERMINE THE NUMBER OF NLR-ARCITECTURES WITH AND WITHOUT INTEGRATED DOMAINS                                                                                                          ################################################################################################################

# Filter out the gene id of each amino acid sequence, ignore the prefix (base/species) name and save in separate file to be used to extract linked NLR-IDs
cat ${genomebase}.NBLRR.aa | grep '^>' | sed 's/>Egrandis_101_V1_//g' > ./${genomebase}_NLR_IDs/${genomebase}_NBLRR_geneid.txt
cat ${genomebase}.NBTIRs.aa | grep '^>' | sed 's/>Egrandis_101_V1_//g' > ./${genomebase}_NLR_IDs/${genomebase}_NBTIR_geneid.txt
cat ${genomebase}.NBCoils.aa | grep '^>' | sed 's/>Egrandis_101_V1_//g' > ./${genomebase}_NLR_IDs/${genomebase}_NBCoil_geneid.txt
cat ${genomebase}.NBRNLs.aa | grep '^>' | sed 's/>Egrandis_101_V1_//g' > ./${genomebase}_NLR_IDs/${genomebase}_NBRNL_geneid.txt

# Navigate the the directory
cd ./${genomebase}_NLR_IDs

# Extract Pfam output with gene-id of specific NLR architectures with column 1,6 and 7 (gene id, pfam id, and integrated domains)
awk 'NR==FNR{for (i=1;i<=NF;i++) a[$i];next} FNR==1 || ($1 in a)' ${genomebase}_NBLRR_geneid.txt ${genomebase}_augustusprotein_pfamscan.txt \
	| sort | uniq | awk '{print $1, $6, $7}' > ${genomebase}_NBLRRid_geneid_pfid.txt
awk 'NR==FNR{for (i=1;i<=NF;i++) a[$i];next} FNR==1 || ($1 in a)' ${genomebase}_NBTIR_geneid.txt ${genomebase}_augustusprotein_pfamscan.txt \
	| sort | uniq | awk '{print $1, $6, $7}' > ${genomebase}_NBTIRid_geneid_pfid.txt
awk 'NR==FNR{for (i=1;i<=NF;i++) a[$i];next} FNR==1 || ($1 in a)' ${genomebase}_NBCoil_geneid.txt ${genomebase}_augustusprotein_pfamscan.txt \
	| sort | uniq | awk '{print $1, $6, $7}' > ${genomebase}_NBCOILid_geneid_pfid.txt
awk 'NR==FNR{for (i=1;i<=NF;i++) a[$i];next} FNR==1 || ($1 in a)' ${genomebase}_NBRNL_geneid.txt ${genomebase}_augustusprotein_pfamscan.txt \
	| sort | uniq | awk '{print $1, $6, $7}' > ${genomebase}_NBRNLid_geneid_pfid.txt

# Read all lines from file1 into the array arr[], and then check for each line in file2 if it already exists within the array (i.e. file1).
# The lines that are found will be printed in the order in which they appear in file2. Note that the comparison in arr uses the entire line
# from file2 as index to the array, so it will only report exact matches on entire lines.
awk 'NR==FNR{arr[$0];next} $0 in arr' ${genomebase}_NBLRRid_geneid_pfid.txt ${genomebase}_protein_nlrid_pfamid_multoccurrence.txt \
	> ${genomebase}_NBLRR_id_number.txt
awk 'NR==FNR{arr[$0];next} $0 in arr' ${genomebase}_NBTIRid_geneid_pfid.txt ${genomebase}_protein_nlrid_pfamid_multoccurrence.txt \
	> ${genomebase}_NBTIR_id_number.txt
awk 'NR==FNR{arr[$0];next} $0 in arr' ${genomebase}_NBCOILid_geneid_pfid.txt ${genomebase}_protein_nlrid_pfamid_multoccurrence.txt \
	> ${genomebase}_NBCOIL_id_number.txt
awk 'NR==FNR{arr[$0];next} $0 in arr' ${genomebase}_NBRNLid_geneid_pfid.txt ${genomebase}_protein_nlrid_pfamid_multoccurrence.txt \
	> ${genomebase}_NBRNL_id_number.txt

# Extracting columns containing integrated domains identified in eadh NLR archtectures 
cat ${genomebase}_NBLRR_id_number.txt | cut -d ' ' -f3 | sed -e $'1i\\\nEgrandis_NBLRRid' > ../../../NLRs_with_IDs/${genomebase}_NBLRR_multiids.txt
cat ${genomebase}_NBTIR_id_number.txt | cut -d ' ' -f3 | sed -e $'1i\\\nEgrandis_NBTIRid' > ../../../NLRs_with_IDs/${genomebase}_NBTIR_multiids.txt 
cat ${genomebase}_NBCOIL_id_number.txt | cut -d ' ' -f3 | sed -e $'1i\\\nEgrandis_NBCOILid' > ../../../NLRs_with_IDs/${genomebase}_NBCOIL_multiids.txt 
cat ${genomebase}_NBRNL_id_number.txt | cut -d ' ' -f3 | sed -e $'1i\\\nEgrandis_NBRNLid' > ../../../NLRs_with_IDs/${genomebase}_NBRNL_multiids.txt

# Detection for NLR gene clusters and head-to-head pairs. 
# Two adjacent genes separated by a short intergenic distance and oriented in divergent (−+) transcription configuration
# This intergenic distance is called promoter region and shared between the two genes
# The two adjacent head-to-head clustered genes are expected tobe regulated by this bidirectional promoter

# Filtering and sorting the input file by chromosome and start position

cat ../${genomebase}.NBLRR.gtf | awk '$3=="gene"{print $0}' | sort -k1,1 -k4,4n \
	> ${genomebase}.NBLRR_filtered_sortedgene.gtf

# Generating intergenic distance between head-to-head clustered genes
bedtools closest -a ${genomebase}.NBLRR_filtered_sortedgene.gtf \
	-b ${genomebase}.NBLRR_filtered_sortedgene.gtf -nonamecheck -io -id -S -D ref \
	> ${genomebase}.NBLRRgenes_h2hdistance.txt

# Filter out only important columns for gene clusters
 cat ${genomebase}.NBLRRgenes_h2hdistance.txt | cut -f1,9,10,18,19 > ${genomebase}.Distance_between_NBLRRgenes.txt

# Filter out those head-to-head NLR gene pairs with less than or equal to 15000 kb intergenic/promoter distance
awk '$5 < -1 && $5 >= -15000{print}' > ../../../Clustered_NLR_genes/${genomebase}.LT15kb_h2hdistance.txt \
	${genomebase}.Distance_between_NBLRRgenes.txt 

# Filter out those head-to-head NLR gene pairs with less than or equal to 200000 kb intergenic/promoter distance
awk '$5 < -1 && $5 >= -200000{print}' > ../../../Clustered_NLR_genes/${genomebase}.LT200kb_h2hdistance.txt \
	 ${genomebase}.Distance_between_NBLRRgenes.txt


