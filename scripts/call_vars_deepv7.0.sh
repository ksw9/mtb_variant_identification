#!/bin/bash 
# Script to call variants with DeepVariant v.7.0. Requires (1) reference, (2) bam, (3) directory, (4) number of shards for parallelization. 
# Cals diploid variants. 

ref=$1
bam=$2
VARS_DIR=$3
NSHARDS=${4:-6}

# Path to trained model.
MODELS=/ifs/labs/andrews/walter/varcal/rui/data/DeepVariant-inception_v3-0.7.0+data-wgs_standard/model.ckpt

# Move to vars_dir if variable is set.
echo $VARS_DIR
if [ -z "$VARS_DIR" ]; then
  echo 'no output directory specified'
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
/opt/deepvariant/bin/make_examples \
--mode calling \
--ref $ref \
--reads $bam \
--examples ${prefix}@${NSHARDS}.tfrecord.gz \
--gvcf ${prefix}@${NSHARDS}.gvcf.tfrecord.gz \
--task {}

##Stage 2: Call variants

echo 'calling variants'

/opt/deepvariant/bin/call_variants \
--outfile ${prefix}.tfrecord.gz \
--examples ${prefix}@${NSHARDS}.tfrecord.gz \
--checkpoint ${MODELS}
 
##Stage 3: Post process variants

echo 'post processing variants'

/opt/deepvariant/bin/postprocess_variants \
--ref $ref \
--infile ${prefix}.tfrecord.gz  \
--nonvariant_site_tfrecord_path ${prefix}@${NSHARDS}.gvcf.tfrecord.gz \
--gvcf_outfile ${prefix}.g.vcf.gz \
--outfile ${prefix}.vcf.gz