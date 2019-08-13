# M. tuberculosis variant calling project.

Here we include scripts to create the plots in our paper, Variant calling approaches alter transmission inferences in Mycobacterium tuberculosis genomic epidemiology studies, available [here] (https://www.nytimes.com/). 

## Scripts included. 

- plot_roetzer_variants.R includes the code for figures 1 and 2. The VCF, FASTA, and tree files associated with each pipeline are available [here] (add URL). 
- plot_outbreak_performance.R includes the code for figures 3-6. The data for those figures are available here [here] (https://www.nytimes.com/).**update URL.
- map_reads.sh Maps reads.
- call_variants_gatk.sh Calls variants with GATK. 
- call_variants_samtools.sh Calls variants with Samtools. 
- call_variants_deepvariant.sh Calls variants with DeepVariant. 
- filter_variants.sh Filters VCF files. Haploidifies DeepVariant diploid VCF files. 
- run_happy.sh Tests performance of a pipeline in recovering true variants (requires a query VCF, truth VCF, and reference fasta).


