#!/bin/bash
# Variant quality score recalibration (VQSR) with GATK. Takes (1) reference genome and (2) VCF file.

# read from command line
ref=$1
vcf=$2

# Local variables
GATK_PATH=/ifs/labs/andrews/walter/bin/gatk-4.1.0.0/gatk
BCFTOOLS=/ifs/labs/andrews/walter/repos/bcftools/bcftools

# set up environment
module load anaconda
source activate gatk_4.0.0.0_kwalter
module load gatk 

# get basename
base=${vcf%.vcf*}
base=$(basename $base)

# get mean QUAL, excluding QUAL of invariant sites, to select variants for training
qual=$(${BCFTOOLS} filter -i 'TYPE == "SNP"'  ${vcf}  | ${BCFTOOLS} query  -f '[%QUAL\n]' | grep -v inf | \
  awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' )

# echo qual
echo 'mean qual of variant sites': $qual

# Create truth set by selecting only high qual variants and only variant sites.
${BCFTOOLS} view --types snps ${vcf} | ${BCFTOOLS} filter -e 'QUAL == inf' | ${BCFTOOLS} filter -i " QUAL > $qual " > ${base}_training_sites.vcf

# Index vcf.
bgzip -f ${base}_training_sites.vcf > ${base}_training_sites.vcf.gz
tabix -f -p vcf ${base}_training_sites.vcf.gz

# Recalibrate variants.
if 
${GATK_PATH}  VariantRecalibrator \
-R ${ref} \
--variant ${vcf} \
--resource:truePos,known=false,training=true,truth=true,prior=15 ${base}_training_sites.vcf.gz \
-an DP \
-an QD \
-an MQRankSum \
-an ReadPosRankSum \
-an FS \
-an SOR \
-an MQ \
-mode SNP \
--max-gaussians 2 \
--output ${base}.recal \
--tranches-file ${base}".tranches"  \
-tranche 100.0 -tranche 99.0 -tranche 96.0 -tranche 93.0 -tranche 90.0 ; 
  then
 echo 'VQSR succeeded with MQRankSum'

elif
# If failed, remove MQRankSum which may not have sufficient variation.
${GATK_PATH}  VariantRecalibrator \
-R ${ref} \
--variant ${vcf} \
--resource:truePos,known=false,training=true,truth=true,prior=15 ${base}_training_sites.vcf.gz \
-an DP \
-an QD \
-an ReadPosRankSum \
-an FS \
-an SOR \
-an MQ \
-mode SNP \
--max-gaussians 2 \
--output ${base}.recal \
--tranches-file ${base}".tranches"  \
-tranche 100.0 -tranche 99.0 -tranche 96.0 -tranche 93.0 -tranche 90.0  
  then
    echo 'VQSR w/o MQRankSum'

elif
# If failed, remove ReadPosRankSum which may not have sufficient variation.
${GATK_PATH}  VariantRecalibrator \
-R ${ref} \
--variant ${vcf} \
--resource:truePos,known=false,training=true,truth=true,prior=15 ${base}_training_sites.vcf.gz \
-an DP \
-an QD \
-an MQRankSum \
-an FS \
-an SOR \
-an MQ \
-mode SNP \
--max-gaussians 2 \
--output ${base}.recal \
--tranches-file ${base}".tranches"  \
-tranche 100.0 -tranche 99.0 -tranche 96.0 -tranche 93.0 -tranche 90.0  
then 
  echo 'VQSR w/o ReadPosRankSum'
  
else 
# If failed, remove both MQRankSum and ReadPosRankSum which may not have sufficient variation.
${GATK_PATH}  VariantRecalibrator \
-R ${ref} \
--variant ${vcf} \
--resource:truePos,known=false,training=true,truth=true,prior=15 ${base}_training_sites.vcf.gz \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQ \
-mode SNP \
--max-gaussians 2 \
--output ${base}.recal \
--tranches-file ${base}".tranches"  \
-tranche 100.0 -tranche 99.0 -tranche 96.0 -tranche 93.0 -tranche 90.0  

echo 'VQSR w/o ReadPosRankSum or MQRankSum'

fi

# Set tranche filter level: this defines the sensitivity for the "truth variants."
ts_filter=99.0

# Apply VQSR.
${GATK_PATH}  ApplyVQSR \
-R ${ref} \
-mode SNP \
--variant ${vcf}  \
--recal-file ${base}.recal  \
--tranches-file ${base}".tranches"  \
--truth-sensitivity-filter-level ${ts_filter} \
--output ${base}_vqsr.vcf

# tabix index and zip all VCF files for vcf-merge to work
bgzip -f -c ${base}_vqsr.vcf >  ${base}_vqsr.vcf.gz
tabix -f -p vcf ${base}_vqsr.vcf.gz

# Test if this works in calculating VQSLOD. If the VQSLOD score is Nan/inf, rerun, setting lowering number of Gaussians.
echo 'VQSLOD contains inf/Nan at: ' $(${BCFTOOLS} query  -f '%INFO/VQSLOD\n' ${base}_vqsr.vcf | grep '\inf$\|\Nan$|\.$' | wc -l)

# Define test for rerunning 
failedSites=$(bcftools query  -f '%INFO/VQSLOD\n' ${base}_vqsr.vcf | grep '\inf$\|\Nan$|\.$' | wc -l)
if [ "$failedSites" -ne 0 ]; then
    echo 'rerunning with maxGaussians set to 1'
  
	${GATK_PATH} --java-options  "-Xmx20g"  VariantRecalibrator \
	-R ${ref} \
	--variant ${vcf} \
	--resource:truePos,known=false,training=true,truth=true,prior=15 ${base}_training_sites.vcf.gz \
	-an DP \
	-an QD \
	-an MQRankSum \
	-an ReadPosRankSum \
	-an FS \
	-an SOR \
	-an MQ \
	-mode SNP \
	--max-gaussians 1 \
	--output ${base}.recal \
	--tranches-file ${base}".tranches"  \
	-tranche 100.0 -tranche 99.0 -tranche 96.0 -tranche 93.0 -tranche 90.0 

	# define trance filter level
	ts_filter=99.0

	# recalibration to SNPs 
	${GATK_PATH} --java-options "-Xmx20g" ApplyVQSR \
	-R ${ref} \
	-mode SNP \
	--variant ${vcf}  \
	--recal-file ${base}.recal  \
	--tranches-file ${base}".tranches"  \
	--truth-sensitivity-filter-level ${ts_filter} \
	--output ${base}_vqsr.vcf

	# tabix index and zip all VCF files for vcf-merge to work
	bgzip -f -c ${base}_vqsr.vcf >  ${base}_vqsr.vcf.gz
	tabix -f -p vcf ${base}_vqsr.vcf.gz
    
    # test if this works in calculating VQSLOD -- issue may be # gaussians is to high - use this to rewrite VQSR dir - change # gaussians if the VQSLOD is Nan/inf
    echo 'VQSLOD contains inf/Nan at: ' $(${BCFTOOLS} query  -f '%INFO/VQSLOD\n' ${base}_vqsr.vcf | grep '\inf$\|\Nan$|\.$' | wc -l)

fi

# Remove intermediate files.
rm ${base}.recal
rm ${base}.recal.idx
rm ${base}_training_sites.vcf.gz
rm ${base}_training_sites.vcf.gz.tbi
rm ${base}".tranches"
rm ${base}_vqsr.vcf