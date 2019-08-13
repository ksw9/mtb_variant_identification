#!/bin/sh
# Requires a single sample VCF file. 

# Add single sample QUAL and FILTER INFO fields to the FORMAT field (before merging)in order to detect pairwise differences based on qual.
# Writes output to current directory. 

# read from command line.
vcf=$1

BCFTOOLS=/ifs/labs/andrews/walter/repos/bcftools/bcftools

prefix=$(basename ${vcf/.*})

# Print to output.
echo 'Adding QUAL and FILTER fields to $vcf'

# Index vcf
$BCFTOOLS index -f ${vcf}

# Extract qual scores and FILTER field (single sample)
$BCFTOOLS query -f'%CHROM\t%POS\t%FILTER\t%QUAL\n' ${vcf} | bgzip -c > ${prefix}_filter_qual.txt.gz

# Index qual_filts.txt.gz
tabix -s1 -b2 -e2 ${prefix}_filter_qual.txt.gz

# Create appropriate VCF header.
echo '##FORMAT=<ID=filt,Number=1,Type=String,Description="Per-sample FILTER">' > ${prefix}_hdr.txt
echo '##FORMAT=<ID=qual,Number=1,Type=Float,Description="Per-sample QUAL">' >> ${prefix}_hdr.txt

# Annotate single-sample VCF with header. 
$BCFTOOLS annotate -a ${prefix}_filter_qual.txt.gz -c CHROM,POS,FORMAT/filt,FORMAT/qual -h ${prefix}_hdr.txt ${vcf} | bgzip -c > ${prefix}_format.vcf.gz

# Index 
$BCFTOOLS index -f ${prefix}_format.vcf.gz

# Remove intermediate files. 
rm ${prefix}_filter_qual.txt.gz ${prefix}_hdr.txt