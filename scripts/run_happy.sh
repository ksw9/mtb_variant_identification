#!/bin/bash

# Run Illumina hap.py to calculate performance metrics. Requires (1) reference fasta, (2) query VCF, (3) truth VCF, 
# (4) genomic region and (5) output file to store summary of *unfiltered* SNP-based results.
# Do filtering before running hap.py.

# Read from command line
ref=$1
query=$2 # first input is query vcf (not g.vcf)
truth=$3
region=${4:-genome} # default is entire genome
out_file=${5} 

# Set up environment. Needs to use Python 2.
module load anaconda
source activate hap.py_0.3.10

# set file names
prefix=$(basename $query)
prefix=${prefix/.vcf*}_${region}
echo $prefix

# run hap.py
if [ "$region" == "genome" ]
  then
    echo "No filter applied, whole genome"
	# run hap.py
	hap.py -r ${ref} $truth $query -o ${prefix} --set-gt hom
	
# PPE genes
elif [ "$region" == "ppe" ]
  then
    echo "No filter applied, ppe"
	hap.py -r ${ref} $truth $query -o ${prefix} --set-gt hom -T 'ppe_genes_rename.bed.gz' 

# Non-PPE genes
elif [ "$region" == "noppe" ]
  then
    echo "No filter applied, ppe"
	hap.py -r ${ref} $truth $query -o ${prefix} --set-gt hom -T ppe_complement.bed.gz 
fi

# collect summary information, print file name along with prefix
echo 'printing out filtered variants'
echo $(basename $query) $(basename $ref) ${region} $(grep SNP,ALL ${prefix}'.summary.csv') >> ${out_file}