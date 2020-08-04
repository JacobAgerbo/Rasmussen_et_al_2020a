#!/bin/bash
#Written by Lara Puetz, April 5, 2017
#lara.c.puetz@gmail.com
#modified by Sarah Mak 27 Oct 2017

#This script is for the analyses of sequence data for the ITS2 amplificon libraries amplifying. The library was sequenced with the Illumina MiSeq 250PE chemistry. All samples where dual indexed with the Illumina Nextera indexing kits
#PCR triplicates were pooled prior to indexing (i.e. the triplicates for a given sample have the same index)


#----------
#Overview
#----------

#1. Download MiSeq sequencing files from the sequencing centre on to the HPC hologenomic servers
#2. Initial FastQC on zipped files
#3. Merge R1 and R2 reads with stagger (without assigning the sequence length)
#4. FastQC on merged reads (all the reads and also on each of them individually)
#5. Check the length distribution of the merged sequences
#6. Do quality filtering for each merged sample
#7. Remove primer and adaptors from the filtered sequences
#8. Move intermediate files into appropriate folders
#9. Check the length distribution again on the filtered sequences
#10. Dereplication of sequences using full-length matching
#11. Sort by size and remove singletons
#12. OTU clustering (deNovo OTU clustering)



#------------------------------------------------------------------------------------------------
#1. Transfer MiSeq sequencing files from the sequencing centre on to the HPC hologenomic servers
#------------------------------------------------------------------------------------------------
module load vsearch/v2.1.2  python/v2.7.12  cutadapt/v1.11   usearch/v9.0.2132  fastx-toolkit/v0.0.13  seqtk/v1.2  blast+/v2.5.0 #in case you don't these in your hpc server
cd data/
mkdir VF
mkdir VF_ITS2
cd VF_ITS2
mkdir 1-rawSequences
mv *.fastq.gz 1-rawSequences/

#------------------------------------------------------------------------------------------------------------
#2. Initial FastQC on zipped files
#------------------------------------------------------------------------------------------------------------
#Look at the FastQC files generatated from the seq centre and look at the summary data on each one of our samples (quality of data and how many reads)

#create a new directory called FastQC_initial and move all fastqc outputs from raw sequences into this folder
module load java/v1.8.0_131  fastqc/v0.11.5
mkdir FastQC_initial
screen
fastqc *.fastq.gz ../FastQC_initial/

#look at some of the html files to get an idea of the quality of the data generated.
#the quality is poor at the beginning of read 1... not so surprising because of the primer region
#also the quality in general is not the greatest BUT look at the data after merging
#how many reads does each file have.... you can look at the first 10 lines of each file to get the answer
cd FastQC_initial/
unzip '*fastqc.zip'
cd ../
mkdir read_counts
head -n10 FastQC_initial/*_fastqc/fastqc_data.txt > read_counts/1-raw_read_count.txt

#remove the unzipped folders for now to save on storage.... will no longer use
cd FastQC_initial/
rm -r *fastqc/
cd ../1-rawSequences/

#--------------------------------------------------------------------------------------------------------------------------------------------
#3. Merge R1 and R2 reads
#Merge read pairs with a 400 bp minimum length for the merged reads using usearch
#--------------------------------------------------------------------------------------------------------------------------------------------
#for usearch help see http://www.drive5.com/usearch/manual/
# and specifically for "fastq_mergepairs" http://www.drive5.com/usearch/manual/cmd_fastq_mergepairs.html

##Starting unzip the fasta files
cd 1-rawSequences/
gunzip *.fastq.gz

rename "pos-" "pos" *

#create a sample list and then the variable sample list within the terminal
find *.fastq > temp
sed 's/_S[0-9]*_L001_R[0-9]_001.fastq//g' temp > temp2
uniq temp2 > sample_list.txt
rm -f temp*

sample_list=$(cat sample_list.txt)
echo $sample_list #to make sure variable was created
wc -l sample_list.txt # count the number of lines... 98 samples

#use stagger script
#in my bin folder: /groups/hologenomics/mak/bin/vsearch_merge_ITS_tagger.sh
xsbatch -c 1 --mem-per-cpu 8000 -- bash /groups/hologenomics/mak/bin/vsearch_merge_ITS_tagger.sh sample_list.txt

#create one file with all of the merged sequences
cat *_merged.fq > merged_all.fq


mv 1-rawSequences

#move all necesarry files into the parent directory and continue from there
cd ..
mv 1-rawSequences/*merged.fq .
mv 1-rawSequences/sample_list.txt .

#--------------------------------------------------------------------------------------------------------------------------------------------
#4. FastQC on merged reads (all of them and also on each individual sample)
#--------------------------------------------------------------------------------------------------------------------------------------------
#look at the quality of the merged sequences
mkdir FastQC_merged
xsbatch -c 1 -- bash /groups/hologenomics/mak/bin/fastqc_merged.sh

#create a new directory called FastQC_merged_all and move all fastqc outputs from raw sequences into this folder
cd FastQC_merged

firefox merged_all_fastqc.html &
#sequence quality is not the greatest... but we still have some filtereing steps to do below

#look at the quality scores for all of the merged data but also individual samples. Did the quality improve after merging?

#how many reads does each file have.... you can look at the first 10 lines of each file to get the answer
#will also see the answer for the merged reads fastqc as well
unzip '*.zip'
head -n10 *_fastqc/fastqc_data.txt > 2-read_count_merged_data.txt
rm -r *_fastqc/
cd ..
mkdir read_counts
mv FastQC_merged/2-read_count_merged_data.txt read_counts/

#--------------------------------------------------------------------------------------------------------------------------------------------
#4. FastQC on merged reads (all of them and also on each individual sample)
#--------------------------------------------------------------------------------------------------------------------------------------------
#look at the quality of the merged sequences
mkdir FastQC_merged
xsbatch -c 1 -- bash ~/bin/fastqc_merged.sh

#create a new directory called FastQC_merged_all and move all fastqc outputs from raw sequences into this folder
cd FastQC_merged

#### to create a new directory for results generated by tagger script
unzip '*.zip'
head -n10 *_fastqc/fastqc_data.txt > 2-read_count_merged_data.txt
rm -r *_fastqc/
mkdir read_counts_tagger
mv 2-read_count_merged_data.txt read_counts_tagger/


cd ..
mkdir merged_tagger
mv *merged.gq merged_tagger/

#-------------------------------------------------------------------------------------------------------------------------------------------
#5. Check the length distribution of the merged sequences (not necessary but a nice visualization of the length distribution of sequences before more filtering steps to come)
#--------------------------------------------------------------------------------------------------------------------------------------------

#Check the length distribution
mkdir length_distributions

#xsbatch -c 1 --time 30:00 -- seqtk seq -a merged_all.fq > length_distributions/merged_all.fa
seqtk seq -a merged_all.fq > length_distributions/merged_all.fa

# altering sequence file (removing the illumima quality score data from the file)
vsearch -sortbylength length_distributions/merged_all.fa --output length_distributions/merged_all_sorted.fasta -sizeout
#Convert to tab format
fasta_formatter -t -i length_distributions/merged_all_sorted.fasta -o length_distributions/merged_all_sorted.tab
cd length_distributions/
mkdir 1-length_dist_merged/
mv * 1-length_dist_merged/
cd 1-length_dist_merged/
awk '{print NR,length($2)}' merged_all_sorted.tab > merged_all_for_plot

######## HERE! ######

gnuplot
set terminal png size 600,300 enhanced font 'Verdana,10'
set output "merged_all_combine_lgth_dist.png"
plot "merged_all_combine_for_plot" using 1:2 with points
exit



#--------------------------------------------------------------------------------------------------------------------------------------------
#6. Do quality filtering for each merged sample
# see here for description and help http://www.drive5.com/usearch/manual/cmd_fastq_filter.html
#ee=exp error; Small E means high quality, large E means low quality; 0.5 more stringent than 1
#--------------------------------------------------------------------------------------------------------------------------------------------
for b in $sample_list
do
vsearch -fastq_filter "$b"_merged.fq -fastaout "$b"_filtered.fa -fastq_maxee 1 --threads 16
done

#--------------------------------------------------------------------------------------------------------------------------------------------
#7. Remove primer and adaptors from the sequence
# http://www.drive5.com/usearch/manual/cmd_search_pcr.html
#--------------------------------------------------------------------------------------------------------------------------------------------

#Strip PCR primers #change from Tue's file.... because we need to reverse compliment the reverse primer
echo '>ITS7
GTGAGTCATCGAATCTTTG
>ITS4R
GCATATCAATAAGCGGAGGA' > ITS2_primers.fa


#look for D2 primers
echo '>U1F
GTGAAATTGTTGAAAGGGAA
>NL4R
CCGTCTTGAAACACGGACC' > D2_primers.fa

#remove primer and adaptors from the sequence
for a in $sample_list
do
usearch32 -search_pcr "$a"_filtered.fa -db ITS2_primers.fa -maxdiffs 3 -pcr_strip_primers -strand plus -ampout "$a"_primstrip.fa -pcrout "$a"_hits.txt -minamp 400 -maxamp 2000
done

## try -minamp 200 (29Oct2017) & it works :) <- need to adjust for D2
screen
for a in $sample_list
do
usearch32 -search_pcr "$a"_filtered.fa -db D2_primers.fa -maxdiffs 3 -pcr_strip_primers -strand plus -ampout "$a"_primstrip.fa -pcrout "$a"_hits.txt -minamp 200 -maxamp 2000
done


#And then, labelling, put a sample label in the FASTA headers (last step before OTU clustering!!)
rm -f *.fasta
screen
sample_list=$(cat sample_list.txt) #whenever enter screen, need to recall the sample list
echo $sample_list

for c in $sample_list
do
sed "-es/^>\(.*\) />\1;barcodelabel=$c /" < "$c"_primstrip.fa > "$c".fasta
done


cat *.fasta > reads.fa
grep -c ">" reads.fa
##Removing any sequences below 200 bp using tagger script: 9656589



#--------------------------------------------------------------------------------------------------------------------------------------------
#8. Move intermediate files into appropriate folders
#--------------------------------------------------------------------------------------------------------------------------------------------

#create directory for all individually merged samples and gzip for future use
mkdir 2-merged_unfiltered
mv *_merged.fq 2-merged_unfiltered/
cd 2-merged_unfiltered/
ls -1 | wc -l #to see how many files are in this directory (98)
gzip *
cd ..

#make directory for filtered data
mkdir 3-filtered_seq
mv *_filtered.fa 3-filtered_seq/
cd 3-filtered_seq
ls -1 | wc -l #to see how many files are in this directory (98)
cd ..

#move all files created from primer removal
mkdir 4-primer_strip_output
mv *_primstrip.fa 4-primer_strip_output/
mv *_hits.txt 4-primer_strip_output/
cd 4-primer_strip_output/
gzip *
cd ..

#make directory for merged, filtered, primer removal and labelled seq files
mkdir 5-merged_filtered_primstrip_labelled/
mv *.fasta 5-merged_filtered_primstrip_labelled/
cd 5-merged_filtered_primstrip_labelled/
gzip *
cd ..

#--------------------------------------------------------------------------------------------------------------------------------------------
#9. Check the length distribution again on the filtered sequences
#--------------------------------------------------------------------------------------------------------------------------------------------
vsearch --sortbylength reads.fa --output length_distributions/reads_sorted2.fasta -sizeout
cd length_distributions/
mkdir 2.1-length_dist_filtered_reads-400/
mv reads_sorted2.fasta 2.1-length_dist_filtered_reads-400/
#Convert to tab format
cd 2.1-length_dist_filtered_reads-400/
fasta_formatter -t -i reads_sorted2.fasta -o reads_sorted.tab
awk '{print NR,length($2)}' reads_sorted.tab > reads_sorted_for_plot

gnuplot
set terminal png size 600,300 enhanced font 'Verdana,10'
set output "reads_sorted_lgth_dist.png"
plot "reads_sorted_for_plot" using 1:2 with points
exit

cd ../../

#--------------------------------------------------------------------------------------------------------------------------------------------
#10. Dereplication of sequences using full-length matching
# creates a file that contains each unique sequence once and sorts them by decreasing cluster size.
# http://www.drive5.com/usearch/manual/cmd_derep_fulllength.html
#derep.fa file contains each unique seq once with an identifier of how often it occurs in the dataset
#--------------------------------------------------------------------------------------------------------------------------------------------
mkdir pipeline

#Dereplication with vsearch
vsearch -derep_fulllength reads.fa  --output pipeline/derepV.fa -sizeout


grep -c ">" pipeline/derepV.fa #341638 unique sequences
#341845

#--------------------------------------------------------------------------------------------------------------------------------------------
#11. Sort by size and remove singletons
#--------------------------------------------------------------------------------------------------------------------------------------------
#Sort by size and remove singletons usearch
usearch32 -sortbysize pipeline/derepV.fa -fastaout pipeline/sorted.fa -minsize 2

# to compare the number of reads still contained within each file:
grep -c ">" pipeline/derepV.fa	# --> 341,638 (total # of unique reads including singletons)
grep -c ">" pipeline/sorted.fa	# --> 106,390 (total # of unique reads excluding singletons)
# the large reduction in number of seqs shows us that there were a lot of singletons (~68.9% of the sequences were singletons)


#--------------------------------------------------------------------------------------------------------------------------------------------
#12. OTU clustering (deNovo OTU clustering) ### is this the way we want to do it? what happens when we do the open reference OTU picking?
#cluster_otus command http://www.drive5.com/usearch/manual/cmd_cluster_otus.html
#The cluster_otus command performs OTU clustering using the UPARSE-OTU algorithm
#--------------------------------------------------------------------------------------------------------------------------------------------
usearch32 -cluster_otus pipeline/sorted.fa -otus pipeline/otus1.fa
# NOTE: look at the output in terminal to tell you how many OTUs and chimeras and make note of it somewhere:   936 OTUs, 23058 chimeras
#default clustering % is 97%


#--------------------------------------------------------------------------------------------------------------------------------------------
#14. label the OTUs as "OTU1" "OTU2".... a nice script to give numbers to the otus
#
#--------------------------------------------------------------------------------------------------------------------------------------------
grep -c ">" pipeline/otus1.fa #516 otus
#Need to first allow for any programs in the drive5_py folder to run from anywhere/any directory within the virtual machine
#
python /groups/hologenomics/mak/bin/drive5_py/fasta_number.py pipeline/otus1.fa  OTU_ > pipeline/otusn.fasta

#check file
cd pipeline
head otusn.fasta
tail otusn.fasta	#check that the number of OTUs is 516 (see results above)
cd ..

#--------------------------------------------------------------------------------------------------------------------------------------------
#15. #Map reads back to OTUs (including singletons)
#--------------------------------------------------------------------------------------------------------------------------------------------
#use 97% id
usearch32 -usearch_global reads.fa -db pipeline/otusn.fasta -strand plus -id 0.97 -uc pipeline/readmap.uc
#use 99% id
usearch32 -usearch_global reads.fa -db pipeline/otusn.fasta -strand plus -id 0.99 -uc pipeline/readmap_99.uc
#97.8% matched

#use "vsearch" if the file is too big for usearch32, 99% id
vsearch --usearch_global reads.fa --db pipeline/otusn_combineall.fasta --strand plus --id 0.99 --uc pipeline/readmap_combine99.uc
#Matching query sequences: 12182492 of 12472893 (97.67%)


#--------------------------------------------------------------------------------------------------------------------------------------------
#16. Create OTU table
#--------------------------------------------------------------------------------------------------------------------------------------------

python /groups/hologenomics/mak/bin/drive5_py/uc2otutab.py  pipeline/readmap_99.uc > pipeline/otu_table.txt


#--------------------------------------------------------------------------------------------------------------------------------------------
#17. Convert the otu table into BIOM format for assigning taxonomy in the next step and take a look at the summary output to see how many number of counts each sample has... this will help inform what cutoffs to use for the next filtering step
#--------------------------------------------------------------------------------------------------------------------------------------------
biom convert -i pipeline/otu_table.txt -o pipeline/otu_table.biom --table-type="OTU table" --to-json

#basic stats
biom summarize-table -i pipeline/otu_table.biom -o pipeline/otu_table_summary.txt
#look at table
cat pipeline/otu_table_summary.txt

#--------------------------------------------------------------------------------------------------------------------------------------------
#18. Assign Taxonomy
#use blast to align sequences using UNITE reference data base on the input file we created above "labeled OTUs fasta file as "OTU1" "OUT2""
#NOTE: UNITE has three sets of QIIME files released, corresponding to the SHs resulting from clustering at the 97% and 99% threshold levels. The third set of files is the result of a dynamic use of clustering thresholds, such that some SHs are delimited at the 97% level, some at the 97.5% level, some at the 98% level, and so on; these choices were made manually by experts of those particular lineages of fungi. The syntax is the same throughout the three sets of files.
#For the trial, I tried all three releases
#--------------------------------------------------------------------------------------------------------------------------------------------
module unload blast+/v2.5.0
module load blast/v2.2.26
module load qiime/vMod

#Try 99% threshold level for ITS2 database
assign_taxonomy.py -i table_16S.txt.blast.txt -m blast -o test


#### to try to produce 3 highest blast results
####
mkdir moreTaxa

#try in 99% threshold level
assign_taxonomy.py -i otusn.fasta -m blast -N 10 -n 1 -r /groups/hologenomics/mak/database/D2_database/D2_new/sequence_gi_filtered.fasta -t /groups/hologenomics/mak/database/D2_database/D2_new/fasta_lineage.txt -o 1-otus_test_OTU2_newdb_Top50
more otusn_tax_assignments.txt


#--------------------------------------------------------------------------------------------------------------------------------------------
#19a. add assignmentsto OTU tables using the output information from the UNITE dynamic clustering database
#--------------------------------------------------------------------------------------------------------------------------------------------

### add assignments to OTU tables

## using taxonomic assignment with dynamic threshold
biom add-metadata -i pipeline/otu_table.biom -o pipeline/2-otus_table_UNITEdyn.biom --observation-metadata-fp pipeline/1-otus_test_OTU2_newdb_Top50/otusn_tax_assignments.txt --observation-header OTUId,taxonomy,confidence,ref --sc-separated taxonomy
ll -htr

#--------------------------------------------------------------------------------------------------------------------------------------------
#19b. add assignmentsto OTU tables using the output information from the UNITE dynamic clustering database
#--------------------------------------------------------------------------------------------------------------------------------------------


#use biom convert to get this file into a "usable" tab delimited format
biom convert -i 2-otus_table_UNITEdyn.biom -o 2-otus_table_UNITEdyn.tab --to-tsv --header-key taxonomy


mkdir 2-otus_table_UNITEdyn
mv 2-otus_table_UNITEdyn.* 2-otus_table_UNITEdyn/

#--------------------------------------------------------------------------------------------------------------------------------------------
#19c. add assignmentsto OTU tables using the output information from the UNITE dynamic clustering database
#--------------------------------------------------------------------------------------------------------------------------------------------

#also get a summary of your new table... should be the same details as above
biom summarize-table -i 2-otus_table_UNITEdyn/2-otus_table_UNITEdyn.biom -o 2-otus_table_UNITEdyn_summary.txt


#--------------------------------------------------------------------------------------------------------------------------------------------
#24. Try LULU script to remove error OTUs (optional)
#check details with Tobias GitHub: https://github.com/tobiasgf/lulu
#--------------------------------------------------------------------------------------------------------------------------------------------
## to create the match list for LULU from OTUs table without annotation using VSEARCH
vsearch --usearch_global pipeline/otusn.fasta --db pipeline/otusn.fasta --self --id .84 --iddef 1 --userout LULU/match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10

## move the stuffs to R and continue


#--------------------------------------------------------------------------------------------------------------------------------------------
#20. Let's start with having a look at the biom table using summarize_taxa_through_plots.py -h
#--------------------------------------------------------------------------------------------------------------------------------------------
summarize_taxa_through_plots.py -h
summarize_taxa_through_plots.py -i 2-otus_table_UNITEdyn.biom -o 3-summary

nohup firefox summary/taxa_summary_plots/bar_charts.html &

#--------------------------------------------------------------------------------------------------------------------------------------------
#21. validate mapping file (not run yet)
#--------------------------------------------------------------------------------------------------------------------------------------------

#run the following program to check my file in the ITS2_sorted directory


validate_mapping_file.py -m mappingSTUcombine.txt -o validate_mapping_file_output_STUcombine

#--------------------------------------------------------------------------------------------------------------------------------------------
#22. Run core diversity analyses (output e.g. barcharts, rarefactions, betadiversity indices) (NOT run Yet)
#--------------------------------------------------------------------------------------------------------------------------------------------

##run core diversity analyses
core_diversity_analyses.py -i otus_table_UNITEdyn.biom -o core_diversity -m mappingfileAnts.txt -e 500 --nonphylogenetic_diversity --suppress_group_significance

#--------------------------------------------------------------------------------------------------------------------------------------------
#23. Filter OTU table (Remove the plant OTUs from the analysis and also all OTUS that had no taxonomic assignment)
#--------------------------------------------------------------------------------------------------------------------------------------------
#Remove the plant OTUs from the analysis and also all OTUS that had no taxonomic assignment
filter_taxa_from_otu_table.py -i 2-otus_table_UNITEdyn/2-otus_table_UNITEdyn.biom -o 2-otu_table_filter_tax.biom -n "No blast hit"
### "k__Plantae"
biom summarize-table -i 2-otu_table_filter_tax.biom -o 2-otu_table_filter_tax_summary.txt



#removed samples from table with very low counts (<1000) NOTE: this step removes all of the negative controls as well
filter_samples_from_otu_table.py -i 2-otu_table_filter_tax.biom -o 3-otu_table_filt_taxDepth.biom -n 1000
biom summarize-table -i 3-otu_table_filt_taxDepth.biom -o 3-otu_table_filt_taxDepth_summary.txt

#use biom convert to get this file into a "usable" tab delimited format
biom convert -i 3-otu_table_filt_taxDepth.biom -o 3-otu_table_filt_taxDepth.tab --to-tsv --header-key taxonomy


## just retain "o__Saccharomycetales" in OTUs table
filter_taxa_from_otu_table.py -i 2-otus_table_UNITEdyn/2-otus_table_UNITEdyn.biom -o 2-otus_table_o_Sac_only.biom -p o__Saccharomycetales
biom summarize-table -i 2-otus_table_o_Sac_only.biom -o 2-otu_table_o_Sac_only_summary.txt
biom convert -i 2-otus_table_o_Sac_only.biom -o 2-otus_table_o_Sac_only.tab --to-tsv --header-key taxonomy

## just retain "c__Saccharomycetes" in OTUs table
filter_taxa_from_otu_table.py -i 2-otus_table_UNITEdyn/2-otus_table_UNITEdyn.biom -o 2-otus_table_class_Sac_only.biom -p c__Saccharomycetes
biom convert -i 2-otus_table_class_Sac_only.biom -o 2-otus_table_class_Sac_only.tab --to-tsv --header-key taxonomy
