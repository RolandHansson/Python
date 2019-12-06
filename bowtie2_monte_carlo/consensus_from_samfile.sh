
############################
### Make consensus sequence.
INFILE=$1
REFERENCE=$2
OUTFILE=$3
PIPELINE=$HOME/bin/pipeline

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


######################################################
####################################
### MAIN
if [ ! -f $(echo $OUTFILE) ]
then
echo -e "\nmaking consensus sequence"
die_unless $INFILE
die_unless $REFERENCE

mkdir scaffolds
touch $OUTFILE
echo -e "starting...\n\n"

#Get SAM files of exons from genome match

cat $REFERENCE | grep -o "Ht.*" | \
while read line; do 
name=`echo $line | sed "s/\s//"`
if [ `samtools view $INFILE $name -c` -eq 0 ]; then
echo -en "\e[1A"; echo -e "\e[0K\rskipping $name"
continue 
fi

echo -e "\e[1A" "\e[0K\rFinding reads from <${name}>";\
scaffold_length=`cat $REFERENCE | grep $name -A1 | tail -n1 | tr -d "\n" | wc -c`

{
samtools mpileup -uf $REFERENCE $INFILE -r $name -d 100000 | bcftools call -c --ploidy 1 > scaffolds/${name}_consensus.vcf;\
} > /dev/null 2> /dev/null #Silence this script, so it doesn't report every time it's run.
echo did $FOLDER/scaffolds/${name}_consensus.fastq; \

scaffold_sequence=`cat scaffolds/${name}_consensus.vcf | perl $PIPELINE/third_party_programs/vcfutils.pl vcf2fq | seqtk seq -A | tail -n1 | tr -d "\n"` 
length_difference=$[scaffold_length-`echo $scaffold_sequence | wc -c`] 
echo "length=$scaffold_length, length_difference=$length_difference"

if [ $length_difference -gt 0 ] ; then # samtools mpileup mislays the end of the sequence unless it has reads. 
  echo "padding sequence with extra 'n'"
  extra_n=`printf '%*s' "$length_difference" | tr ' ' "n"`
  scaffold_sequence=`echo "${scaffold_sequence}$extra_n"`
fi

echo -e ">$name\n$scaffold_sequence" >> $OUTFILE

echo -en "\e[1A"; echo "created $OUTFILE" 
done
echo -en "\e[1A"; echo -e "\e[0K\rDone";\
else
echo $OUTFILE already exists, skipping.
fi
