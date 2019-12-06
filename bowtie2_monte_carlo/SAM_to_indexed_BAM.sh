#!/bin/bash
MINARGS=1

INFILE=$1
OUTFILE=$2
FASTQ_FILE=$3
SAMPLE_NAME=$4

TEMP_NAME="temp"
SAM_OPTIONS="-b -S"
PIPELINE="$HOME/bin/pipeline"

function die_unless () 
{
  if [ ! -f $1 ]
  then  
    echo -e >&2 "Cannot find file '$1'"
    move_on_to_next_file 
    exit ${2:-1}
#  else 
#    echo -e >&2 "Checked that file '$1' exists"
  fi 
}

function add_read_group ()
{
echo "add_read_group( $1, $2, $3 )" 
current_file=$1
outfile=`echo $current_file | sed 's/.bam/_rg.bam/'`
file_name=$2
fastq_file=${3-"${file_name}_R1_001.fastq"}

echo -e "current_file=$current_file\n outfile=$outfile\n file_name=$filename\n fastq_file=$fastq_file"

#Sample ID
SAMPLE_NAME=`echo "$file_name}" | cut -d "_" -f 1` 

#Sequencing Number
SEQ_NUMBER=`echo "${file_name}" | cut -d "_" -f 2`

#Find first sequence name in fastq file
EXAMPLE_SEQ_NAME=`head $fastq_file -n 1`
echo -e "EXAMPLE_SEQ_NAME:" $EXAMPLE_SEQ_NAME "\n"

#Name of machine, can be had from fastq
MACHINE_NAME=`echo "${EXAMPLE_SEQ_NAME}" | cut -d ":" -f 1 | tr -d "@"`

#Lane number, can be had from fastq
LANE=`echo "${EXAMPLE_SEQ_NAME}" | cut -d ":" -f 4`

PLATFORM="ILLUMINA"

LIBRARY="LIB1"

RG_OPTIONS="RGID=${MACHINE_NAME}.LANE${LANE}.${SEQ_NUMBER} RGLB=${LIBRARY} RGPL=${PLATFORM} RGPU=${MACHINE_NAME}.LANE${LANE}.${SAMPLE_NAME} SM=${SAMPLE_NAME}"

java -jar $PIPELINE/third_party_programs/picard/build/libs/picard.jar AddOrReplaceReadGroups \
I=${current_file} \
O=$OUTFILE \
$RG_OPTIONS 

}


#######################################################################
###########################################
### Convert genome SAMfile to BAM format and sort. 
if ([ ! -f $(echo $OUTFILE) ] )
then
echo -e "\nconverting genome match to BAM, double single end, $current_file"
die_unless $INFILE

samtools view $SAM_OPTIONS \
-o $TEMP_NAME.bam \
$INFILE

samtools sort "$TEMP_NAME.bam" \
-o ${TEMP_NAME}_sorted.bam 

rm $TEMP_NAME.bam

samtools index ${TEMP_NAME}_sorted.bam

add_read_group ${TEMP_NAME}_sorted.bam $SAMPLE_NAME $FASTQ_FILE 

samtools index $OUTFILE

rm ${TEMP_NAME}_sorted.bam

else
echo $current_file already exists, skipping.
fi
