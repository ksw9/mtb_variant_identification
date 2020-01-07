#!/bin/bash
# Variant calling with Pilon. Requires: (1) reference genome, (2) bam file, (3) optionally, output directory.

# Read input arguments. 
ref=$1 
bam=$2
VARS_DIR=$3

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

# Define output prefix. 
prefix=$(basename ${bam/.rmdup.bam})_pilon

# Run pilon pipeline with default parameters. 
java -Xmx16G -jar ${PILON} --genome ${ref} --frags ${bam} --outdir ${VARS_DIR} --vcf --output ${prefix}

# Zip and index
bgzip -c ${prefix}.vcf > ${prefix}.vcf.gz
tabix -p vcf ${prefix}.vcf.gz

#Error handling
if [ "$?" != "0" ]; then
	echo "[Error]" $LINENO "pilon failed!" 1>&2
	exit 1
fi

#### PRINT OUTPUT ####
echo "Done====" >&2