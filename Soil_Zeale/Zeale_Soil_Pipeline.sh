##########################################################################################
#                                         Initial QC                                     #
##########################################################################################

### Create FASTQC output
FASTQ_DIR=''
WORK_DIR=''

module load java fastqc/v0.11.8a
fastqc -o $WORK_DIR/ $FASTQ_DIR/Zeale_S1.collapsed

### Create FASTQC output from all samples
cd $FASTQ_DIR/
find *_fastqc.zip > temp
sed 's/_fastqc.zip//g' temp > temp2
uniq temp2 > sample_list.txt
rm -f temp*
sample_list=$(cat sample_list.txt)
cd $WORK_DIR
for a in $sample_list
  do
  cd "$a"_fastqc
    total_seqs=`cat fastqc_data.txt | grep 'Total Sequences' | cut -f 2`
    gc_percent=`cat fastqc_data.txt | grep '%GC' | cut -f 2`
    seq_length=`cat fastqc_data.txt | grep 'Sequence length' | cut -f 2`
    seq_qual=`cat fastqc_data.txt | awk '/>>Per base sequence quality/,/>>END_MODULE/' | tail -n +3 | head -n -1 | awk '{total+=$2} END {print total/NR}'`
    n_count=`cat fastqc_data.txt | awk '/>>Per base N content/,/>>END_MODULE/' | tail -n +3 | head -n -1 | awk '{total+=$2} END {print total/NR}'`
		Dedub_level=`cat fastqc_data.txt | grep '#Total Deduplicated Percentage' | cut -f 2`
		echo -e "File Name:\t"$a"\nNumber of Sequences:\t${total_seqs}\nGC%:\t${gc_percent}\nSequence Length:\t${seq_length}\nAverage per base sequence quality:\t${seq_qual}\nRemaining data (%) after deduplication:\t${Dedub_level}\nN%\t${n_count}" > ../"$a"_short.txt
    cd ..
    rm -r "$a"_fastqc
done
rm sample_list.txt
cat *.txt > FastQC_PostTrim_report.txt  # Short report of fastQC output
rm *_short.txt

##########################################################################################
#                           Trimming with AdapterRemoval                                 #
##########################################################################################

module load AdapterRemoval
FASTA_DIR='/groups/hologenomics/jagerbo/data/02-Metabarcoding/Insects/Zeale/data'
WORK_DIR='/groups/hologenomics/jagerbo/data/02-Metabarcoding/V2'

cd $FASTA_DIR/
find *fastq.gz > temp
sed 's/_L001_R[1-2]_001.fastq.gz//g' temp > temp2
sed 's/TOG-VYJV-//g' temp2 > temp3
uniq temp3 > sample_list.txt
rm -f temp*
sample_list=$(cat sample_list.txt)
cd $WORK_DIR

for a in $sample_list
  do
    AdapterRemoval --file1 $FASTA_DIR/TOG-VYJV-"$a"_L001_R1_001.fastq.gz --file2 $FASTA_DIR/TOG-VYJV-"$a"_L001_R2_001.fastq.gz \
    --threads 28 \
    --mm 0.05 \
    --minlength 30 \
    --shift 5 \
    --trimns \
    --trimqualities \
    --qualitybase 33 \
    --minquality 28 \
    --basename $WORK_DIR/"$a" \
    --collapse
done

##########################################################################################
#                                   Sort files with Begum                                #
##########################################################################################

module load python/v2.7.12 Begum/v1.0

for i in 'Soil'
  do
  Begum sort -l pool_${i}.txt -p primers_zeale.txt -t tags_zeale.txt -s PSinfo_${i}.txt -d sorted -o sorted_${i}
done


##########################################################################################
#                                    Filter with Begum                                   #
##########################################################################################

module load python/v2.7.12 Begum/v1.0
for i in 'Soil'
  do
  Begum filter -i sorted_${i} -s ../PSinfo_${i}.txt -p 1 -l 2 -d ../filtered_3_of_3/ -o ${i}
done

##########################################################################################
#                                    Cluster with SumaClust                              #
##########################################################################################

module load sumaclust/v1.0.20

for i in 'Soil'
  do
    python /groups/hologenomics/software/DAMe/v0.9/bin/convertToUSearch.py -i ${i}.fna -lmin 100 -lmax 300
    sumaclust -e FilteredReads.forsumaclust.fna  -F OTUs_${i}.fna -B OTUs_${i}.biom -t 0.97
    python /groups/hologenomics/software/DAMe/v0.9/bin/tabulateSumaclust.py -i OTUs_${i}.fna -o table_${i}.txt -blast
done


##########################################################################################
#                                    BLAST with blatsn                                   #
##########################################################################################

module load vsearch/v2.1.2  python/v2.7.12  cutadapt/v1.11   usearch/v9.0.2132  fastx-toolkit/v0.0.13  seqtk/v1.2 blast/v2.2.26 qiime/vMod

for i in 'Soil'
  do
    biom summarize-table -i OTUs_${i}.biom -o OTUs_${i}_summary.txt
done

for i in Soil
  do
    blastn -query table_${i}.txt.blast.txt -db nt -num_alignments 1 -qcov_hsp_perc 80 -perc_identity 84 > ${i}_hits.blast
done

### LULU
for i in 'Soil'
  do
    makeblastdb -in table_${i}.txt.blast.txt -parse_seqids -dbtype nucl
    blastn -db table_${i}.txt.blast.txt -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query table_${i}.txt.blast.txt
done

##########################################################################################
#                                    Tax Assignment                                      #
##########################################################################################
# Upload table_Soil.txt.blast.txt to gbif database (https://www.gbif.org/tools/sequence-id)
# Tax file will be curated with LULU match list and merged in R
