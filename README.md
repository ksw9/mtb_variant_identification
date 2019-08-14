# Variant calling for *M. tuberculosis* transmission inference.

Here we include scripts to create the plots in our paper, Variant calling approaches alter transmission inferences in *Mycobacterium tuberculosis* genomic epidemiology studies, available [here] (https://www.nytimes.com/). 

- Add short abstract here. 

## Scripts included. 

- map_reads.sh 
  - Maps reads.
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
  - Apply GATK's VQSR, using an internal training set (defined here as variants with QUAL > the mean).

- plot_roetzer_variants.R 
  - includes the code for figures 1 and 2. The VCF, FASTA, and tree files associated with each pipeline are available in the Supplementary Files.  
- plot_outbreak_performance.R 
  - includes the code for figures 3-5. The data for those figures are available in the Supplementary Files.


## Required tools. 

- Bowtie2	http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

- BWA	http://bio-bwa.sourceforge.net/

- SMALT	https://www.sanger.ac.uk/science/tools/smalt-0

- Sambamba	https://lomereiter.github.io/sambamba/

- Bcftools/Samtools	https://samtools.github.io/bcftools/

- GATK	https://software.broadinstitute.org/gatk/

- DeepVariant	https://github.com/google/deepvariant

- hap.py	https://github.com/Illumina/hap.py

