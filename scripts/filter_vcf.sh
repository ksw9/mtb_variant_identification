#!/bin/bash
# Soft-filter VCF files. Requires (1) reference fasta, (2) input VCF, (3) scripts directory.

# read from command line
ref=$1
vcf=$2
SCRIPTS_DIR=$3

# Define prefix for naming subsequent VCF files. 
prefix=$(basename ${vcf/.*})
echo $prefix

# Define file names. 
allsites_vcf=${prefix}_allsites.vcf.gz
ann_vcf=${prefix}_ann.vcf.gz
filt_vcf=${prefix}_filt.vcf.gz
sort_vcf=${prefix}_sort.vcf.gz
vqsr_vcf=${prefix}_allsites_vqsr.vcf.gz
norm_vcf=${prefix}_norm.vcf.gz
hap_vcf=${prefix}_hap.vcf
tmp_vcf=${prefix}_tmp.vcf.gz
tmp_vqsr_vcf=${prefix}_tmp_vqsr.vcf.gz

# Tabix index
tabix -p vcf ${vcf}

#### Validate variants first. Ensures VCF is not truncated. 
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

##### DeepVariant gVCF files: convert to allSites VCF, then haploidify, then annotate, then apply filter. ####
if [[ "$vcf" == *"deep.g.vcf.gz" ]]; then 
  
  # First haploidify DeepV calls within the gVCF.
  ${SCRIPTS_DIR}haploidify.py ${vcf} ${hap_vcf}
  
  # Zip and index haploid VCF. 
  bgzip -f ${hap_vcf} > ${hap_vcf}.gz
  tabix -f -p vcf ${hap_vcf}.gz
  
  # Make an all sites samtools file before applying filters (filtering first on the GVCF version works strangely). 
  bcftools convert --gvcf2vcf ${hap_vcf}.gz --fasta-ref ${ref} -O z -o ${allsites_vcf}

  # Index allsites haploid VCF.
  tabix -f -p vcf ${allsites_vcf}
  
  # Remove intermediate VCF files. 
  rm ${hap_vcf}.gz
    
fi
     
  # Next apply soft filter for QUAL < 40. Then soft filter for low D and Min D.
  $BCFTOOLS filter --mode + --soft-filter 'lowQ' -e 'QUAL != 0 & QUAL < 40' ${allsites_vcf} -O z | $BCFTOOLS filter --mode + --soft-filter 'lowD' -e 'MIN_DP < 5 | DP < 5' -O z -o ${filt_vcf}

  # Index filtered VCF
  tabix -f -p vcf ${filt_vcf}
  
  # Add filter and quality score fields to sample INFO fields before merging. 
  ${SCRIPTS_DIR}format_single_vcf.sh ${allsites_vcf}

#### Samtools files ####
if [[ "$vcf" == *"samtools.vcf.gz" ]]; then 
  echo 'filtering samtools vcf'
  
  # Index vcf
  tabix -f -p vcf ${vcf}
  
  # Define bam. Add .1 to ref name if GCA ref
  if [[ "$vcf" == *"GCA"* ]]; then 
    bam=$(basename ${vcf/_samtools.vcf.gz}).1.rmdup.bam
  else
    # Set BAM file name. 
    bam=$(basename ${vcf/_samtools.vcf.gz}).rmdup.bam
  fi
  
  bam_dir=$(dirname ${vcf/vars/})/
  bam=${bam_dir}${bam}
  
  # Index bam
  samtools index ${bam}
  
  # Annotate variants.
  gatk VariantAnnotator \
  -V ${vcf} \
  -O ${ann_vcf} \
  -I ${bam} \
  -A QualByDepth \
  -A MappingQualityRankSumTest \
  -A ReadPosRankSumTest \
  -A FisherStrand \
  -A StrandOddsRatio
   
  # Normalize alleles. 
  bcftools norm --fasta-ref ${ref} ${ann_vcf} -O z -o ${norm_vcf}

  # Make an all sites samtools file before applying filters (filtering the GVCF version works strangely). 
  bcftools convert --gvcf2vcf ${norm_vcf} --fasta-ref ${ref} -O z -o ${tmp_vcf}
 
  # Index allsites vcf - may not be able to index because of unsorted positions. 
  tabix -f -p vcf ${tmp_vcf}
  
  # Sort.
  bcftools sort ${tmp_vcf} -O z -o ${allsites_vcf}

  # Index sorted vcf
  tabix -f -p vcf ${allsites_vcf}

  # VQSR on annotated, all sites variants.
  ${SCRIPTS_DIR}vqsr.sh ${ref} ${allsites_vcf}
  
  # Index filtered VCF 
  tabix -f -p vcf ${vqsr_vcf}
  
  # Apply filters not to sorted file. There is no QUAL score for the expanded Ref blocks, ignore these when applying the filter. 
  $BCFTOOLS filter --mode + --soft-filter 'lowQ' -e 'QUAL != "." & QUAL < 40' ${sort_vcf} -O z | $BCFTOOLS filter  --mode + --soft-filter 'lowD' -e ' INFO/DP < 5 | MinDP < 5' -O z > ${filt_vcf}

  # Index filtered VCF
  tabix -f -p vcf ${filt_vcf}
  
  # Add filter and quality score fields to sample INFO fields before merging. 
  ${SCRIPTS_DIR}format_single_vcf.sh ${allsites_vcf}
  ${SCRIPTS_DIR}format_single_vcf.sh ${vqsr_vcf}
  
  # Remove intermediate samtools files. 
  rm ${tmp_vcf}* ${norm_vcf} ${ann_vcf}* ${allsites_vcf} 
    
fi

#### GATK files ####
 if [[ "$vcf" == *"gatk.vcf.gz" ]]; then 
  echo 'filtering gatk vcf'
  
  # Index filtered VCF
  tabix -f -p vcf ${norm_vcf}
  
  # Edit AD info for GATK only before merging. 
  zcat ${vcf} | sed 's/AD,Number=R/AD,Number=./g' |  bgzip -c > ${tmp_vcf}
  tabix -f -p vcf ${tmp_vcf}
  
  # VQSR on annotated variants.
  ${SCRIPTS_DIR}vqsr.sh ${ref} ${tmp_vcf}

  # Apply filters not to VQSR file. Don't use: RGQ < 20 filter (isn't applied to samtools files).
  bcftools filter --mode + --soft-filter 'lowQ' -e 'QUAL < 40' ${vcf} -O z | $BCFTOOLS filter --mode + --soft-filter 'lowD' -e 'INFO/DP < 5' -O z > ${filt_vcf}

  # Index filtered VCF
  tabix -f -p vcf ${filt_vcf}
  
  # Edit AD info for GATK only before merging. 
  zcat ${filt_vcf} | sed 's/AD,Number=R/AD,Number=./g' |  bgzip -c >  ${allsites_vcf}
  zcat ${tmp_vqsr_vcf} | sed 's/AD,Number=R/AD,Number=./g' |  bgzip -c > ${vqsr_vcf} 

  # Add filter and quality score fields to sample INFO fields before merging. 
  ${SCRIPTS_DIR}format_single_vcf.sh ${allsites_vcf}
  ${SCRIPTS_DIR}format_single_vcf.sh ${vqsr_vcf} 
  
  # Remove tmp files. 
  rm ${allsites_vcf} ${vqsr_vcf} ${tmp_vcf}
fi

## If success
if [ $? == 0 ] 
then 
  echo $vcf 'finished filtering'
else
  echo $vcf 'failed'
fi