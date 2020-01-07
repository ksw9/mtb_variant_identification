#!/bin/bash
# Soft-filter VCF files. Requires (1) reference fasta, (2) input VCF, 
# (3) if merging variants 'true', this will properly format single-sample variants for merge 
# so that QUAL score can be included for all samples.

# read from command line
ref=$1
vcf=$2
merge=$3

# Set up environment.
source config.txt

# Define prefix for naming subsequent VCF files. 
prefix=$(basename ${vcf/.vcf.gz})
if [[ "$vcf" == *deep.g.vcf.gz ]]; then 
  prefix=$(basename ${vcf/.g.vcf.gz})
fi

echo $prefix

# Define file names for intermediate/ filtered VCFs.  
allsites_vcf=${prefix}_allsites.vcf.gz
ann_vcf=${prefix}_ann.vcf.gz
filt_refcall=${prefix}_refcall.vcf.gz
filt_lowq=${prefix}_lowq.vcf.gz
filt_vqsr=${prefix}_vqsrfilt.vcf.gz
sort_vcf=${prefix}_sort.vcf.gz
vqsr_vcf=${prefix}_allsites_vqsr.vcf.gz
norm_vcf=${prefix}_norm.vcf.gz
hap_vcf=${prefix}_hap.vcf
tmp_vcf=${prefix}_tmp.vcf.gz
tmp_vqsr_vcf=${prefix}_tmp_vqsr.vcf.gz

# Tabix index
tabix -p vcf ${vcf}

#### Validate variants first. Ensures VCF is not truncated. Important for DeepVariant files, some of which are truncated. 
if  [[ "$vcf" != *"pilon.vcf.gz" ]]; 
then
if
  gatk ValidateVariants -V ${vcf}
then 
  echo 'continue, vcf validated'
else 
 echo ${vcf} truncated, exiting
 echo ${vcf} >> vcf_redo.txt
 #rm ${vcf}
 exit
fi
fi
# 
##### DeepVariant gVCF files: convert to allSites VCF, then haploidify, then annotate, then apply filter. ####
if [[ "$vcf" == *deep.g.vcf.gz ]]; then 
  
  #First haploidify DeepV calls within the gVCF.
  ~/.conda/envs/my_gatk_4.0.0.0_kwalter/bin/python ${SCRIPTS_DIR}haploidify.py ${vcf} ${hap_vcf}
  
  #Zip and index haploid VCF. 
  bgzip -f ${hap_vcf} > ${hap_vcf}.gz
  tabix -f -p vcf ${hap_vcf}.gz
  
  #Make an all sites samtools file before applying filters (filtering first on the GVCF version works strangely). 
  ${BCFTOOLS} convert --gvcf2vcf ${hap_vcf}.gz --fasta-ref ${ref} -O z -o ${allsites_vcf}

  #Index allsites haploid VCF.
  tabix -f -p vcf ${allsites_vcf}
  
  #Remove intermediate VCF files. 
  rm ${hap_vcf}.gz
    
  # Filter RefCalls from DeepVariant VCF.  
  ${BCFTOOLS} view  -e 'FILTER  = "RefCall" | TYPE == "INDEL" | QUAL == 0' ${allsites_vcf} -O z -o ${filt_refcall}

  # No RefCall/ exclude lowQ
  ${BCFTOOLS} view  -e 'FILTER  = "RefCall" | TYPE == "INDEL" | QUAL < 40' ${allsites_vcf} -O z -o ${filt_lowq}

  # Index filtered VCFs
  tabix -f -p vcf ${filt_refcall}
  tabix -f -p vcf ${filt_lowq}
  
  # For outbreak files,add filter and quality score fields to sample INFO fields before merging. 
  if [[ "$merge" == true ]]; then 
    ${SCRIPTS_DIR}format_single_vcf.sh ${allsites_vcf}
  fi
fi

#### Samtools files ####
if [[ "$vcf" == *samtools.vcf.gz ]]; then 
  echo 'filtering samtools vcf'
  
  # Index vcf
  tabix -f -p vcf ${vcf}
  
  # Define bam. Add .1 to ref name if GCA ref.
  if [[ "$vcf" == *"GCA"* ]]; then 
    bam=$(basename ${vcf/_samtools.vcf.gz}).1.rmdup.bam
  else
    # Set BAM file name. 
    bam=$(basename ${vcf/_samtools.vcf.gz}).rmdup.bam
  fi
  
  bam_dir=$(dirname ${vcf/vars/})/
  bam=${bam_dir}${bam}
  
  # Index bam
  ${BCFTOOLS} index ${bam}
  
  # Annotate variants.
  ${GATK_40} VariantAnnotator \
  -V ${vcf} \
  -O ${ann_vcf} \
  -I ${bam} \
  -A QualByDepth \
  -A MappingQualityRankSumTest \
  -A ReadPosRankSumTest \
  -A FisherStrand \
  -A StrandOddsRatio
   
  # Normalize alleles. 
  ${BCFTOOLS} norm --fasta-ref ${ref} ${ann_vcf} -O z -o ${norm_vcf}

  # Make an all sites samtools file before applying filters (filtering the GVCF version works strangely). 
  ${BCFTOOLS} convert --gvcf2vcf ${norm_vcf} --fasta-ref ${ref} -O z -o ${tmp_vcf}
 
  # Index allsites vcf - may not be able to index because of unsorted positions. 
  tabix -f -p vcf ${tmp_vcf}
  
  # Sort.
  ${BCFTOOLS} sort ${tmp_vcf} -O z -o ${allsites_vcf}

  # Index sorted vcf
  tabix -f -p vcf ${allsites_vcf}

  # VQSR on annotated, all sites variants.
  ${SCRIPTS_DIR}vqsr.sh ${ref} ${allsites_vcf}
  
  # Index filtered VCF 
  tabix -f -p vcf ${vqsr_vcf}

  # Filter VQSR
  ${BCFTOOLS} view  -e 'TYPE == "INDEL" | FILTER == "VQSRTrancheSNP99.00to100.00" ' ${vqsr_vcf} -O z -o ${filt_vqsr}

  # Filter lowQ
  ${BCFTOOLS} view  -e 'TYPE == "INDEL" | QUAL < 40' ${allsites_vcf} -O z -o ${filt_lowq}

  # Index filtered VCFs
  tabix -f -p vcf ${filt_lowq}
  tabix -f -p vcf ${filt_vqsr}
  
  # For outbreak files,add filter and quality score fields to sample INFO fields before merging. 
  if [[ "$merge" == true ]]; then 
    ${SCRIPTS_DIR}format_single_vcf.sh ${allsites_vcf}
    ${SCRIPTS_DIR}format_single_vcf.sh ${vqsr_vcf}
  
  fi

  # Remove intermediate samtools files. 
  rm ${tmp_vcf}* ${norm_vcf} ${ann_vcf}* 
    
fi

#### GATK files ####
 if [[ "$vcf" == *gatk.vcf.gz ]]; then 
  echo 'filtering gatk vcf'
  
  # Move raw VCF to filtered directory. 
  cp $vcf ${allsites_vcf}
  tabix -p vcf ${allsites_vcf}
  
  # VQSR on annotated variants.
  ${SCRIPTS_DIR}vqsr.sh ${ref} ${allsites_vcf}

  # Filter VQSR
  ${BCFTOOLS} view  -e 'TYPE == "INDEL" | FILTER == "VQSRTrancheSNP99.00to100.00" ' ${vqsr_vcf} -O z -o ${filt_vqsr}

  # Filter lowQ
  ${BCFTOOLS} view  -e 'TYPE == "INDEL" | QUAL < 40'  ${allsites_vcf} -O z -o ${filt_lowq}

  # Index filtered VCFs
  tabix -f -p vcf ${filt_lowq}
  tabix -f -p vcf ${filt_vqsr}
     
  # For outbreak files,add filter and quality score fields to sample INFO fields before merging. 
  if [[ "$merge" == true ]]; then 
    ${SCRIPTS_DIR}format_single_vcf.sh ${allsites_vcf}
    ${SCRIPTS_DIR}format_single_vcf.sh ${vqsr_vcf}
  
  fi

fi

#### Pilon files ####
 if [[ "$vcf" == *pilon.vcf.gz ]]; then 
  echo 'filtering pilon vcf'
  
  # Define bam. Add .1 to ref name if GCA ref
  if [[ "$vcf" == *"GCA"* ]]; then 
    bam=$(basename ${vcf/_samtools.vcf.gz}).1.rmdup.bam
  else
    # Set BAM file name. 
    bam=$(basename ${vcf/_pilon.vcf.gz}).rmdup.bam
  fi
  
  bam_dir=$(dirname ${vcf/vars/})/
  bam=${bam_dir}${bam}
  
  # Filter Pilon filtered LowCov/Amb/Del variants from Pilon VCF. 
  ${BCFTOOLS} view  -e 'FILTER = "Amb" | FILTER = "Del" | TYPE == "INDEL"' ${vcf}  -O z -o ${allsites_vcf}
  tabix -f -p vcf ${allsites_vcf}
  
  # Annotate variants - not easily applied to Pilon.
#   gatk --java-options "-Xmx60G" VariantAnnotator \
#   -V ${vcf} \
#   -O ${ann_vcf} \
#   -I ${bam} \
#   -A QualByDepth \
#   -A MappingQualityRankSumTest \
#   -A ReadPosRankSumTest \
#   -A FisherStrand \
#   -A StrandOddsRatio
  
  # VQSR on annotated variants.
  #${SCRIPTS_DIR}vqsr.sh ${ref} ${ann_vcf}

  # Filter Pilon filtered LowCov/Amb/Del variants from Pilon VCF. 
  ${BCFTOOLS} view  -e 'FILTER  = "LowCov" ' ${allsites_vcf}  -O z -o ${filt_refcall}
  tabix -f -p vcf ${filt_refcall}
    
  # Exclude lowQ and deletions (not allowed by hap.py). 
  ${BCFTOOLS} view  -e 'QUAL < 40' ${allsites_vcf}  -O z -o ${filt_lowq}
  tabix -f -p vcf ${filt_lowq}

  # For outbreak files,add filter and quality score fields to sample INFO fields before merging. 
  if [[ "$merge" == true ]]; then 
    echo 'formatting sample fields'
    ${SCRIPTS_DIR}format_single_vcf.sh ${allsites_vcf}
  
  fi
fi

#### Breseq files ####
  if [[ "$vcf" == *breseq*.vcf.gz ]]; then 
    echo 'filtering breseq vcf'
   # Define bam file
    bam=$(dirname $vcf)/data/reference.bam
 
  # Annotate variants not easily applied to Breseq.
#   gatk  --java-options "-Xmx60G" VariantAnnotator \
#   -V ${prefix}.vcf.gz \
#   -O ${ann_vcf} \
#   -I data/reference.bam \
#   -A MappingQuality \
#   -A QualByDepth \
#   -A MappingQualityRankSumTest \
#   -A ReadPosRankSumTest \
#   -A FisherStrand \
#   -A StrandOddsRatio
  
    # After reformatting, only adds MMQ, which has no variance. Only annotation available for VQSR is DP, skip this. 
 
    # Move raw VCF to filtered directory. 
    cp $vcf ${allsites_vcf}
    tabix -p vcf ${allsites_vcf}
  
    # Filter lowQ
    ${BCFTOOLS} view  -e 'TYPE == "INDEL" | QUAL < 40' ${allsites_vcf} -O z -o ${filt_lowq}
    
    # For outbreak files,add filter and quality score fields to sample INFO fields before merging. 
  if [[ "$merge" == true ]]; then 
    ${SCRIPTS_DIR}format_single_vcf.sh ${allsites_vcf}
  fi

  
fi
  
## If success
if [ $? == 0 ] 
then 
  echo $vcf 'finished filtering'
else
  echo $vcf 'failed'
fi