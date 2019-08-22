# Genomic variant identification methods alter *Mycobacterium tuberculosis* transmission inference

Here we include scripts to create the plots in our paper, available [here]. 

- Pathogen genomic data are increasingly used to characterize global and local transmission patterns of important human pathogens and to inform public health interventions. Yet there is no current consensus on how to measure genomic variation. We investigated the effects of variant identification approaches on transmission inferences for M. tuberculosis by comparing variants identified by five different groups in the same sequence data from a clonal outbreak. We then measured the performance of commonly used variant calling approaches in recovering variation in a simulated tuberculosis outbreak and tested the effect of applying increasingly stringent filters on transmission inferences and phylogenies. 


## Scripts included. 

- map_reads.sh 
  - Maps short-read read data to a reference genome with BWA, Bowtie2, or SMALT.
- call_variants_gatk.sh 
  - Calls variants with GATK. 
- call_variants_samtools.sh 
  - Calls variants with Samtools. 
- call_variants_deepvariant.sh 
  - Calls variants with DeepVariant. 
- filter_variants.sh 
  - Filters VCF files. 
- haplidify.py
  - Converts diploid single-sample VCF file to haplid VCF file using allele depth information.  
- format_single_vcf.sh 
  - Adds single sample VCF INFO fields including QUAL and FILTER to the FORMAT field before merging in order to preserve single-sample information when merging several samples into a multi-sample VCF (i.e. for examining pairwise differences).
- run_happy.sh 
  - Tests performance of a pipeline in recovering true variants.
- vqsr.sh
  - Apply GATK's VQSR, using an internal training set of high-quality variants (defined here as variants with QUAL > the mean).
- plot_roetzer_variants.R 
  - This script compares variants identified by 5 different groups within the same sequence data from a clonal M. tuberculosis outbreak. 
  - Sequence data from Roetzer et al. 2013 Plos Med, (ENA Study Accession: PRJEB6945).

## Data 
- A zipped archive including the submitted and formatted vcfs and fasta files from 5 different groups for the same sequence data from a clonal *M. tuberculosis* outbreak, that can be analyzed with plot_roetzer_variants.R.
- Includes vcfs, fastas, bootstrap replicate trees fit by raxml-ng, and accompanying sample metadata from the original study, Roetzer et al. 2013 Plos Med. 

## Required tools. 

- Bowtie2	http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

- BWA	http://bio-bwa.sourceforge.net/

- SMALT	https://www.sanger.ac.uk/science/tools/smalt-0

- Sambamba	https://lomereiter.github.io/sambamba/

- Bcftools/Samtools	https://samtools.github.io/bcftools/

- GATK	https://software.broadinstitute.org/gatk/

- DeepVariant	https://github.com/google/deepvariant

- hap.py	https://github.com/Illumina/hap.py

