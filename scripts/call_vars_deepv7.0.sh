#!/bin/bash 
# Script to call variants with DeepVariant v.7.0. Requires (1) reference, (2) bam, (3) directory, (4) number of shards for parallelization. 
# Cals diploid variants. 

ref=$1
bam=$2
VARS_DIR=$3
NSHARDS=${4:-6} # Default NSHARDS is 6.

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

echo $ref

# Get prefix. 
prefix=$(basename ${bam%%.*})_deep

echo 'running script!'

## Running deepvariant on singularity container non-interactively.

## Stage1: Make examples

echo 'make examples'

time seq 0 $((NSHARDS-1)) | \
parallel --eta --halt 2 --joblog "log"  \
${DEEPVAR}make_examples \
--mode calling \
--ref $ref \
--reads $bam \
--examples ${prefix}@${NSHARDS}.tfrecord.gz \
--gvcf ${prefix}@${NSHARDS}.gvcf.tfrecord.gz \
--task {}

##Stage 2: Call variants

echo 'calling variants'

${DEEPVAR}call_variants \
--outfile ${prefix}.tfrecord.gz \
--examples ${prefix}@${NSHARDS}.tfrecord.gz \
--checkpoint ${DEEPVAR_MODEL}
 
##Stage 3: Post process variants

echo 'post processing variants'

${DEEPVAR}postprocess_variants \
--ref $ref \
--infile ${prefix}.tfrecord.gz  \
--nonvariant_site_tfrecord_path ${prefix}@${NSHARDS}.gvcf.tfrecord.gz \
--gvcf_outfile ${prefix}.g.vcf.gz \
--outfile ${prefix}.vcf.gz