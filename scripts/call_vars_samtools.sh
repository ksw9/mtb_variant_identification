#!/bin/bash
# Variant calling with Samtools. Requires: (1) reference genome, (2) bam file, (3) ploidy, (4) optionally, output directory.

# Source environment with all necessary software.
module load anaconda; source activate gatk_4.0.0.0_kwalter
module load samtools/1.9
module load htslib/1.9

# Read from command line: ref genome, fastq 1, fastq 2.
ref=$1 
bam=$2 
ploidy=$3
VARS_DIR=$4 

# Move to vars_dir if variable is set.
echo $VARS_DIR
if [ -z "$VARS_DIR" ]; then
  echo 'no output directory specified'
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
bcftools mpileup -Ou -f ${ref} --threads 8 -a INFO/ADF,INFO/ADR,FORMAT/AD,DP,SP ${bam}  | bcftools call --threads 8 --ploidy ${ploidy} -m -O z --gvcf 0 -o ${prefix}_samtools.vcf.gz 
	
# bgzip file and index
tabix -p vcf ${prefix}_samtools.vcf.gz
	
fi

if [ "$?" != "0" ]; then
	echo "[Error]" $LINENO "samtools variant calling failed!" 1>&2
	exit 1
fi