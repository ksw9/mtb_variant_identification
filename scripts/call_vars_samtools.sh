#!/bin/bash
# Variant calling with Samtools. Requires: (1) reference genome, (2) bam file, (3) ploidy, (4) optionally, output directory.

# Read input arguments
ref=$1 
bam=$2 
ploidy=$3
VARS_DIR=$4 

# Set up environment.
source config.txt

# Move to vars_dir if variable is set, otherwise set VARS_DIR to current working directory. 
echo $VARS_DIR
if [ -z "$VARS_DIR" ]; then
  echo 'no output directory specified'
  VARS_DIR=$(pwd)
  else 
  echo $VARS_DIR specified
  cd $VARS_DIR
fi

# List of BAMs that need to be combined. (Only if list is not provided).
if [ -z ${bam} ]; then 
  ls ${BAMS_DIR}*${mapper}.rmdup.bam > ${mapper}_bam_list.txt
  bam=${mapper}_bam_list.txt
fi

## singly genotype with samtools.
echo 'singly calling'
prefix=${bam/.*}
prefix=$(basename $prefix)
    
# Call variants singly. Adds additional annotations for downstream filtering. 
${BCFTOOLS} mpileup -Ou -f ${ref} --threads 8 -a INFO/ADF,INFO/ADR,FORMAT/AD,DP,SP ${bam}  | ${BCFTOOLS} call --threads 8 --ploidy ${ploidy} -m -O z --gvcf 0 -o ${prefix}_samtools.vcf.gz 
	
# bgzip file and index
tabix -p vcf ${prefix}_samtools.vcf.gz
	
fi

if [ "$?" != "0" ]; then
	echo "[Error]" $LINENO "samtools variant calling failed!" 1>&2
	exit 1
fi