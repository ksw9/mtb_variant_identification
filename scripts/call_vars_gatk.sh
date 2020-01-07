#!/bin/bash
# Variant calling with GATK. Requires: (1) reference genome, (2) bam file, (3) ploidy, (4) optionally, output directory.

# Read input arguments. 
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

# Set names of intermediate files.
prefix=${bam%.rmdup.bam}
prefix=$(basename $prefix)

# if no GATK dictionary, index reference genome
if [ ! -f ${ref%.*}".dict" ] ; then
echo "GATK dictionary for $ref" >&2
picard CreateSequenceDictionary \
	REFERENCE=${ref} \
	OUTPUT=${ref%.*}".dict" 
fi

# If BAM index does not exist, index BAM.
if [ ! -s ${bam}.bai ]; then 
  echo 'indexing bam' 
  samtools index ${bam}
fi
	 
# name gVCF file
gvcf=${prefix}"_gatk.g.vcf"

# name VCF file
vcf=${prefix}"_gatk.vcf.gz"
	  
# call variants with GATK 4.0 to GVCF. 
${GATK_41} --java-options "-Xmx50g" HaplotypeCaller \
-R ${ref} \
-ploidy ${ploidy} \
-I ${bam} \
-ERC GVCF \
-O ${gvcf}

# GVCF to VCF. 
# Use latest GATK version: 4.1.0.0 (allows outputing non-variant sites). 
${GATK_41} --java-options '-Xmx50g' GenotypeGVCFs \
-R ${ref} \
--variant ${gvcf} \
-ploidy ${ploidy} \
--include-non-variant-sites true \
--output ${vcf}
# min base quality score is 10 by default.

# Remove tmp files.
#rm ${gvcf} 
#rm ${gvcf}.idx

#Error handling
if [ "$?" != "0" ]; then
	echo "[Error]" $LINENO "GATK failed!" 1>&2
	exit 1
fi

#### PRINT OUTPUT ####
echo "Done====" >&2