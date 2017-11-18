#!/bin/sh

#These are passed as inputs to this script:
LIB=$1 # <libdir>
MAF=$2 #tumor maf file that contains point mutation information
BAM=$3 #BAM file
REFo=$4 #original reference seqeunce
NAME=$5 #sample id
xbase=${BAM##*/}
STUB=${xbase%.*}



echo "LIB: $LIB"
echo "MAF: $MAF"
echo "BAM: $BAM"
echo "REFo: $REFo"
echo "NAME: $NAME"
echo "STUB: $STUB"

#WORK_DIR=$(echo pwd) #working directory

#check to see if maflite or maf file:
#maf lite will have columns: contig, start_position, end_position
#maf will have columns: Chromosome Start_position End_position

##########################################################
#Generate Intervals
##########################################################

echo "generate intervals from maf or maflite"
title=`cat $MAF | head -n100 | grep chr | sed 's/\t/\n/g' | grep -P ^chr`

if [ ! -z "$title" ];
then
echo "this is a maflite"
grep -v "#" $MAF | gcol chr start end | sed -e 's/\s\+/:/'| sed -e 's/\s\+/-/' | sed '/^#/ d' | sed -e '1d' | sed 's/^M:/MT:/' > snp_mutations.intervals
if [ "$?" -ne 0 ]; then echo "command failed"; exit 1; fi

fi

titlemaf=`cat $MAF | head -n100 | grep position | sed 's/\t/\n/g' | grep -P ^Chromosome`
if [ ! -z "$titlemaf" ];
then
echo "this is a maf file"
grep -v "#" $MAF | gcol Chromosome Start_position End_position | sed -e 's/\s\+/:/'| sed -e 's/\s\+/-/' | sed '/^#/ d' | sed -e '1d' | sed -e 's/^M:/MT:/' > snp_mutations.intervals
if [ "$?" -ne 0 ]; then echo "command failed"; exit 1; fi

fi

#print out intervals
sed 's/:/'"$(printf '\011')"'/g' snp_mutations.intervals | sed 's/-/'"$(printf '\011')"'/g' > snp_mutations.intervals.bed

#######################################################
#Extracting Paired Reads from BAM
#######################################################

echo "extracting paired reads from bam"

samtools view -L snp_mutations.intervals.bed $BAM | cut -f1 > IDs_all.txt
if [ "$?" -ne 0 ]; then echo "command extracting read IDs failed"; exit 1; fi

java -Xmx7g -jar $LIB/FilterSamReads.jar I=$BAM O=tmp_bam.bam READ_LIST_FILE=IDs_all.txt FILTER=includeReadList WRITE_READS_FILES=false VALIDATION_STRINGENCY=LENIENT
if [ "$?" -ne 0 ]; then echo "command generating BAM based on read list failed"; exit 1; fi


echo "convert bam to fastq"
samtools view -H tmp_bam.bam | sed '$d' - > tmp_header_T.sam
if [ "$?" -ne 0 ]; then echo "command failed"; exit 1; fi

samtools view tmp_bam.bam | awk '$2 < 2040 { print }' > tmp0_T.sam
if [ "$?" -ne 0 ]; then echo "command failed"; exit 1; fi

cat tmp_header_T.sam tmp0_T.sam > tmp_filteredbamT.sam

java -Xmx7g -jar $LIB/SamToFastq.jar I=tmp_filteredbamT.sam F=${NAME}_tmp_sequence_1.fastq F2=${NAME}_tmp_sequence_2.fastq VALIDATION_STRINGENCY=LENIENT
if [ "$?" -ne 0 ]; then echo "command failed"; exit 1; fi


##########################################################
#Writing Files to Annotations:
##########################################################

echo "writing out annotation file for upload"
echo -e "$PWD/${NAME}_tmp_sequence_1.fastq\n$PWD/${NAME}_tmp_sequence_2.fastq"   >> ${NAME}.rna_reads_fastq_list.list
echo "done"
