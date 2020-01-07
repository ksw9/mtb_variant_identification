#!/bin/bash
# Variant calling with breseq. Requires: (1) reference genome, (2) paired-end read1, (3) paired-end read2, (4) variant output directory.

# Read input arguments. 
ref=$1 
p1=$2
p2=$3
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

# Work in special output directory so that calls against multiple ref genomes can be made. 
WORKING_DIR=tmp_$(basename ${ref/.fa})
rm -rf $WORKING_DIR
mkdir $WORKING_DIR
cd $WORKING_DIR

# Define output prefix. 
mapper=bowtie2
samp=$(basename ${p1/_1.fq*})
prefix=${samp}_${mapper}_$(basename ${ref%.*})_breseq

# Run breseq pipeline with default parameters. 
${BRESEQ} -r ${ref} ${p1} ${p2} -o ${VARS_DIR}/$WORKING_DIR -n ${prefix} -j 4

# Correctly format by adding genotype information
cat ${VARS_DIR}/$WORKING_DIR/data/output.vcf | grep '#' > ${prefix}_header.txt

# Update VCF header with sample name (allows later merging of multiple samples into one VCF)
sed -i.bak "s/FILTER\tINFO/FILTER\tINFO\tFORMAT\t$samp/" ${prefix}_header.txt 

# Add genotype information to reset of VCF
cat ${VARS_DIR}/$WORKING_DIR/data/output.vcf | grep -v '#' > ${prefix}_body.txt

awk '$0=$0"\tGT\t1"' ${prefix}_body.txt > ${prefix}_mod.txt

# Put header with body
cat ${prefix}_header.txt ${prefix}_mod.txt | bgzip > ${prefix}.vcf.gz
tabix -p vcf ${prefix}.vcf.gz

# Move vcf to vcf folder
mv ${prefix}.vcf.gz $VARS_DIR

# Remove intermediate files
rm ${prefix}_header.txt ${prefix}_mod.txt ${prefix}_body.txt 

# Move to VARS DIR
cd $VARS_DIR

# Remove tmp folder
rm -rf tmp_$(basename $ref)

#Error handling
if [ "$?" != "0" ]; then
	echo "[Error]" $LINENO "breseq failed!" 1>&2
	exit 1
fi

#### PRINT OUTPUT ####
echo "Done====" >&2