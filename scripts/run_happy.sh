#!/bin/bash
# Run Illumina hap.py to calculate performance metrics. Requires (1) reference fasta, (2) query VCF, (3) truth VCF, and optionally: (4) filter, (5) genomic region, (6) ROC field.

# Source environment with all necessary software.
module load anaconda; source activate gatk_4.0.0.0_kwalter
module purge
module load anaconda
source activate hap.py_0.3.10
REF_DIR=/ifs/labs/andrews/walter/varcal/data/refs/

# Read from command line
ref=$1
query=$2 # first input is query vcf (not g.vcf)
truth=$3
region=${4:-genome} # default is entire genome
out_file=${5:-QUAL} # output prefix

echo $query      
echo $region
echo $ref

# set file names
output=$(basename $query)
output=${output/.vcf*}_${region}
echo $output

# run hap.py
if [ "$region" == "genome" ]
  then
    echo "No filter applied, whole genome"
	# run hap.py
	hap.py -r ${ref} $truth $query -o ${output} --set-gt hom
	
# PPE genes
elif [ "$region" == "ppe" ]
  then
    echo "No filter applied, ppe"
		hap.py -r ${ref} $truth $query -o ${output} --set-gt hom -T ${REF_DIR}'ppe_genes_rename.bed.gz' 

# Non-PPE genes
elif [ "$region" == "noppe" ]
  then
    echo "No filter applied, ppe"
		hap.py -r ${ref} $truth $query -o ${output} --set-gt hom -T ${REF_DIR}ppe_complement.bed.gz 
fi

# collect summary information, print file name along with output
echo 'printing out filtered variants'
echo $(basename $query) $(basename $ref) ${region} $(grep SNP,ALL ${output}'.summary.csv') >> ${out_file}