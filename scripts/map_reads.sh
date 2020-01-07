#!/bin/bash
# Map fastq files to reference genome. Requires (1) ref genome, (2) mapping algorithm, (3) full path to read file 1, (4) full path to read file 2 (if exists), (5) prefix for output files (if paired end reads).

# Read from command line: ref genome, fastq 1, fastq 2.
ref=${1} # 1st input is full path to reference genome
mapper=$2 # 2nd input is mapping algorithm 
p1=$3  # 3th input is full path to read 1
p2=$4  # 4th input is full path to read 2 (if paired-end)
BAMS_DIR=$5

# Set up environment.
source config.txt

# Move to vars_dir if variable is set, otherwise set VARS_DIR to current working directory. 
echo $BAMS_DIR
if [ -z "$BAMS_DIR" ]; then
  echo 'no output directory specified'
  BAMS_DIR=$(pwd)
  else 
  echo $BAMS_DIR specified
  cd $BAMS_DIR
fi

# Define prefix.
prefix=$(basename ${p1%%.*})

# Create temp directory specific for files (so no overwriting of other temp files).
TMP_DIR=${TMP_DIR}${prefix}_${mapper}/
mkdir -p ${TMP_DIR}

# Set names of intermediate files.
ref_index=${ref%.*}
ref_prefix=$(basename $ref_index)

# remove suffix
sam=${prefix}_${mapper}_${ref_prefix}'.sam' 
bam=${prefix}_${mapper}_${ref_prefix}'.bam' 
rgbam=${prefix}_${mapper}_${ref_prefix}'.rg.bam' 
sortbam=${prefix}_${mapper}_${ref_prefix}'.rg.sorted.bam'
rmdupbam=${prefix}_${mapper}_${ref_prefix}'.rmdup.bam'

# help message if no inputs are entered (-z checks if 1st arugment string is null).
if [ -z ${1} ]
then
  echo "This script will use Bowtie 2 to align two paired-end fastq files"
  echo "$(basename $0) <refGenome> <mapper> <dataset> <readFile1> <readFile2>"
  echo "First input is the Mtb ref genome"
  echo "Second input is the read mapper"
  echo "Third input is the dataset"
  echo "Forth input is location of file containing read 1"
  echo "Fifth input is location of file containing read 2 (if paired-end)" 
  exit
fi

# if no p2 given 
if [ -z ${p2} ]
then
  echo "read pair not specified"
fi

# get machine id_lane from SRA read identifier (old Illumina fastq format)
# if gzipped fastq
if [[ ${p1} == *.gz ]]; then
  seqid=$(zcat ${p1} | head -n 1)
else
# if not gzipped
  seqid=$(cat ${p1} | head -n 1)
fi
seqid="$(echo $seqid | cut -d' ' -f1)"
seqid="$(echo $seqid | cut -d':' -f3)"

id_lane=${seqid:-readgroup1} # default is "readgroup1"

#### MAPPING ####
# bowtie2 mapping
if [ $mapper == 'bowtie' ] || [ $mapper == 'bowtie2' ] ; then

  # if no indexing, index reference genome
  if [ ! -f ${ref%.*}".1.bt2" ] ; then
  echo "bowtie2 indexing $ref" >&2
  ${BOWTIE2_BUILD} ${ref} ${ref_index}
  fi 

  # map
  #if paired-end reads
  echo "mapping with bowtie2" >&2
  if [ ! -z ${p2} ]; then 
  ${BOWTIE2} --threads 4 -X 1100 -x ${ref_index} -1 ${p1} -2 ${p2} -S ${sam}
  # -x basename of reference index 
  # --un-gz gzips sam output
  # p is number of threads
  # end-to-end mapping is the default mode
  # -X 2000 increase maximum fragment length from default (500) to allow longer fragments (paired-end only) -X 2000 **used for roetzer data ***
  # ART simulations X = 1100 (mean fragment length = 650bp + 3 x 150-bp stdev)
  
  # if single-end reads
  elif [ -z ${p2} ]; then
    echo "single reads"
  ${BOWTIE2} --time --threads 4 -X 1100 -x ${ref_index} -U ${p1} -S ${sam}
  # -U for unpaired reads
  fi
  # Error handling
  if [ "$?" != "0" ]; then
    echo "[Error]" $LINENO "bowtie mapping failed ${p1}!" 1>&2
    exit 1
  fi
fi

# bwa mapping
if [ $mapper == 'bwa' ]; then
  # if no indexing, index reference genome
  if [ ! -f ${ref}".sa" ] ; then
  echo "bwa indexing $ref" >&2
  ${BWA} index ${ref} 
  fi

  # map
  echo "mapping with bwa" >&2
  # if paired-end reads
  if [ ! -z ${p2} ]; then 
    ${BWA} mem -t 4 ${ref} ${p1} ${p2} > ${sam}
  #  -t no. of threads.
  # if single-end reads
  elif [ -z ${p2} ]; then
      echo "single reads"
    ${BWA} mem -t 4 ${ref} ${p1} >  ${sam}
  fi
  # Error handling
  if [ "$?" != "0" ]; then
    echo "[Error]" $LINENO "bwa mapping failed ${p1}!" 1>&2
    exit 1
  fi
fi

# smalt mapping
if [ $mapper == 'smalt' ]; then
  # if no indexing, index reference genome
  # builds a hash index for the TB genome. Words of 14 base pair length are sampled at every 8th position in the genome. 
  if [ ! -f ${ref_index}".sma" ] ; then
  echo "smalt indexing $ref" >&2
  ${SMALT} index ${ref_index} ${ref}
  #leave wordlength at default (k=13) and skipstep=wordlength as default
  fi
  
  # map
  echo "mapping with smalt" >&2
  # if paired-end reads
   if [ ! -z ${p2} ]; then 
    ${SMALT} map -n 8 -i 1100 -o ${sam} ${ref_index} ${p1} ${p2} 
    # for ART simulations: -i max insert size: 1100bp = mean insert size (650-bp) + 3 x 150bp (they actually mean fragment length)
   elif [  -z ${p2} ]; then 
     echo "single reads"
     ${SMALT} map -n 4 -i 1100 -o ${sam} ${ref_index} ${p1}  
  fi
  # Error handling
  if [ "$?" != "0" ]; then
    echo "[Error]" $LINENO "smalt mapping failed ${p1}!" 1>&2
    exit 1
  fi
fi

#### POST-PROCESSING ####

# Convert sam to bam 
${SAMBAMBA} view -t 7 -S -h ${sam} -f bam -o ${bam}
#-S auto-detects input format, -h includes header, -o directs output

# add/replace read groups for post-processing with GATK
${PICARD} AddOrReplaceReadGroups \
  INPUT=${bam} \
  OUTPUT=${rgbam} \
  RGID=${id_lane} \
  RGLB=library1 \
  RGPU=${id_lane} \
  RGPL="illumina" \
  RGSM=${prefix}

# Sort the BAM 
${SAMBAMBA} sort ${rgbam}

# Index BAM
${SAMBAMBA} index -t4 ${sortbam}

# Remove duplicates. ## here.
${SAMBAMBA} markdup -r -p -t4  ${sortbam} ${rmdupbam} --tmpdir=${TMP_DIR} 

# Error handling
if [ "$?" != "0" ]; then
    echo "[Error]" $LINENO "Remove dups failed ${p1}!" 1>&2
    exit 1
fi

# Remove intermediate files.
rm ${sam} ${bam} ${rgbam} ${sortbam} ${sortbam}.bai

#### PRINT OUTPUT ####
echo "Done====" >&2
