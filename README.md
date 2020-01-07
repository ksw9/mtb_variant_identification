# Genomic variant identification methods alter *Mycobacterium tuberculosis* transmission inference

Here we include scripts used in our investigation of variant calling approaches: https://www.biorxiv.org/content/10.1101/733642v1?rss=1. 

Pathogen genomic data are increasingly used to characterize global and local transmission patterns of important human pathogens and to inform public health interventions. Yet there is no current consensus on how to measure genomic variation. We investigated the effects of variant identification approaches on transmission inferences for M. tuberculosis by comparing variants identified by five different groups in the same sequence data from a clonal outbreak. We then measured the performance of commonly used variant calling approaches in recovering variation in a simulated tuberculosis outbreak and tested the effect of applying increasingly stringent filters on transmission inferences and phylogenies. 


## Scripts included. 

- map_reads.sh 
  - Maps short-read read data to a reference genome with BWA, Bowtie2, or SMALT.
- call_variants_gatk.sh 
  - Calls variants with GATK. 
- call_variants_samtools.sh 
  - Calls variants with Samtools. 
- call_variants_deepvariant.sh 
  - Calls variants with DeepVariant. 
- call_variants_pilon.sh 
  - Calls variants with Pilon.
- call_variants_breseq.sh 
  - Calls variants with Breseq.
- filter_vf.sh 
  - Filters VCF files. 
- haplidify.py
  - Converts diploid single-sample VCF file to haplid VCF file using allele depth information.  
- format_single_vcf.sh 
  - Adds single sample VCF INFO fields including QUAL and FILTER to the FORMAT field before merging in order to preserve single-sample information when merging several samples into a multi-sample VCF (i.e. for examining pairwise differences).
- run_happy.sh 
  - Tests performance of a pipeline in recovering true variants.
- vqsr.sh
  - Apply GATK's VQSR, using an internal training set of high-quality variants (defined here as variants with QUAL > the mean).
- config.txt. 
  - Set local paths to programs and directories to be sourced in config.txt.

## Data 
- Simulated data is available in a repository here: https://sdr.stanford.edu/items/druid:mr554nj9219. 

## Required tools are listed in the paper supplement.

- Bowtie2	http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

- BWA	http://bio-bwa.sourceforge.net/

- SMALT	https://www.sanger.ac.uk/science/tools/smalt-0

- Sambamba	https://lomereiter.github.io/sambamba/

- Bcftools/Samtools	https://samtools.github.io/bcftools/

- GATK	https://software.broadinstitute.org/gatk/

- DeepVariant	https://github.com/google/deepvariant

- Pilon https://github.com/broadinstitute/pilon

- Breseq https://github.com/barricklab/breseq/releases/tag/v0.34.1

- hap.py	https://github.com/Illumina/hap.py

