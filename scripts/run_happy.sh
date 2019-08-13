#!/bin/bash
# Run Illumina hap.py to calculate performance metrics. Requires (1) reference fasta, (2) query VCF, (3) truth VCF, and optionally: (4) filter, (5) genomic region, (6) ROC field.

# Source environment with all necessary software.
module load anaconda; source activate gatk_4.0.0.0_kwalter
module purge
module load anaconda
source activate hap.py_0.3.10

# Read from command line
ref=$1
query=$2 # first input is query vcf (not g.vcf)
truth=$3
filter=$4
region=${5:-genome} # default is entire genome
roc_field=${6:-QUAL} # which annotation should ROC curve be calculated on 'INFO.VQSLOD'

# make perform directory
mkdir -p perform

# make file to hold stats
if [ ! -e "perform/genome_perform.csv" ] ; then
    touch "perform/genome_perform.csv"
fi

echo $query      
echo $filter
echo $region
echo $ref

# set file names
output=$(basename $query)
output=perform/${output/.vcf*}_${region}_${filter}
echo $output

# run hap.py
if [ "$filter" == "raw" ] && [ "$region" == "genome" ]
  then
    echo "No filter applied, whole genome"
	# run hap.py
	hap.py -r ${ref} $truth $query -o ${output} --set-gt hom #--no-leftshift --no-decompose --engine=vcfeval # --write-counts --roc ${roc_field} --set-gt hom --preserve-info  --write-vcf --restrict-regions ${confidence_bed}
fi

if [ "$filter" != "raw" ] && [ "$region" == "genome" ]
  then
    echo "Applying filter, whole genome"
		hap.py -r ${ref} $truth $query -o ${output} --set-gt hom --filters-only ${filter} #  --write-counts --roc ${roc_field} --set-gt hom --filters-only ${filter}  --preserve-info  --write-vcf --restrict-regions ${confidence_bed}
fi

if [ "$filter" == "raw" ] && [ "$region" == "ppe" ]
  then
    echo "No filter applied, ppe"
		hap.py -r ${ref} $truth $query -o ${output} --set-gt hom -T ${REF_DIR}'ppe_genes_rename.bed.gz' # --write-counts --roc ${roc_field}  --set-gt hom  --preserve-info  --write-vcf --restrict-regions ${confidence_bed}
fi

if [ "$filter" != "raw" ] && [ "$region" == "ppe" ]
  then
    echo "Applying filter, ppe"
		hap.py -r ${ref} $truth $query -o ${output} --set-gt hom --filters-only ${filter} -T ${REF_DIR}'ppe_genes_rename.bed.gz'  #  --write-counts --roc ${roc_field} --set-gt hom --filters-only ${filter} -T  'ppe_genes_rename.bed.gz'  --preserve-info  --write-vcf --restrict-regions ${confidence_bed}
fi

if [ "$filter" == "raw" ] && [ "$region" == "noppe" ]
  then
    echo "No filter applied, ppe"
		hap.py -r ${ref} $truth $query -o ${output} --set-gt hom -T ${REF_DIR}ppe_complement.bed.gz # --write-counts --roc ${roc_field}  --set-gt hom  --preserve-info  --write-vcf --restrict-regions ${confidence_bed}
fi

if [ "$filter" != "raw" ] && [ "$region" == "noppe" ]
  then
    echo "Applying filter, ppe"
		hap.py -r ${ref} $truth $query -o ${output} --set-gt hom --filters-only ${filter} -T ${REF_DIR}ppe_complement.bed.gz  #  --write-counts --roc ${roc_field} --set-gt hom --filters-only ${filter} -T  'ppe_genes_rename.bed.gz'  --preserve-info  --write-vcf --restrict-regions ${confidence_bed}
fi

# collect summary information, print file name along with output
echo 'printing out filtered variants'
echo $(basename $query) $(basename $ref) ${filter} ${region} $(grep SNP,ALL ${output}'.summary.csv') >> 'perform/genome_perform.csv'