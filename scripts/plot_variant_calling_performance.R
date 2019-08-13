##########################
#### Plot performance ####
##########################
# Plot performance of variant calling tool combinations on both genome-wide and single genome variant calls. 

rm(list=ls())
library(dplyr)
library(ggplot2)
library(happyR)
library(magrittr)
library(gtools)
library(tidyr)
library(cowplot)
library(MASS)
library(plyr)
library(cowplot)
library(gridExtra)
library(grid)

##################################
## Plot genome-wide performance ##
##################################
# Read csv on genome-wide performance. All are mapped to the H37Rv reference. 
df <- read.csv('performance/genome_performance_H37Rv.csv', stringsAsFactors = F)

# Define error bars from genome-wide errors (so that errors correspond to entire genome).
error_bars <- df %>%
  filter(region == 'Genome') %>% 
  group_by(pipeline,filt) %>% 
  mutate(genome_fp_max = fp_mean[which(region == 'Genome')] + fp_std[which(region == 'Genome')] , 
         genome_fn_max = fn_mean[which(region == 'Genome')]  + fn_std[which(region == 'Genome')] ) %>%
  dplyr::select(pipeline,filt,genome_fp_max, genome_fn_max)

# Add error bars to summary. 
df_full <- merge(df, error_bars)

# Plot FPs, colored by genomic region.
g1 <- df_full %>%
  filter(region != 'Genome') %>% 
  ggplot(aes(x = filt, y = fp_mean, fill = region)) + #fill = pipeline,
  geom_errorbar(aes(ymin = fp_mean, ymax = genome_fp_max)) +
  geom_bar(stat = 'identity', color="black") + 
  theme_minimal() +
  facet_grid(cols = vars(pipeline)) +
  xlab('Filter') + 
  ylab('False positive errors (SNPs)') +
  guides(colour = guide_legend(order = 2)) + 
  theme(text = element_text(size=13), axis.text.x = element_text(angle = 90, size = 10), strip.text = element_text(size = 8)) + 
  scale_fill_discrete(name = 'Region')
g1

# Plot FNs
g2 <- df_full %>%
  filter(region != 'Genome') %>% 
  ggplot(aes(x = filt, y = fn_mean, fill = region)) + #fill = pipeline,
  geom_errorbar(aes(ymin = fn_mean, ymax = genome_fn_max)) +
  geom_bar(stat = 'identity', color="black") + 
  theme_minimal() +
  facet_grid(cols = vars(pipeline)) +
  xlab('Filter') + 
  ylab('False negative errors (SNPs)') +
  guides(colour = guide_legend(order = 2)) + 
  theme(text = element_text(size=13), axis.text.x = element_text(angle = 90, size = 10), 
        strip.text = element_text(size = 8)) + 
  scale_fill_discrete(name = 'Region')
g2

g3 <- plot_grid(g1,g2, ncol = 1, labels = c('a','b'))
g3

####################################
## Performance across ref genomes ##
####################################

# Read performance across ref genoems. 
df_dists <- read.csv('genome_performance_allrefs.csv', stringsAsFactors = F)

# Plot FPs vs. distance to ref.
g1 <- df_summary %>%
  filter(mapper == 'bwa' & caller == 'gatk') %>%
  ggplot(aes(x = dist, y = fp_mean)) +
  geom_errorbar(aes(ymin=fp_mean - fp_std, ymax=fp_mean + fp_std), alpha = .6) +
  geom_point() + 
  xlab('Distance to query genome (SNPs)') +
  ylab('False positive errors (SNPs)') +
  scale_x_continuous(trans='log10') +
  theme(legend.position="none") + 
  theme(axis.title.x=element_blank())
g1

# Plot FNs vs. distance to ref.
g2 <- df_dists %>%
  #filter(mapper == 'bwa') %>%
  filter(mapper == 'bwa' & caller == 'gatk') %>%
  ggplot(aes(x = dist, y = fn_mean)) +
  geom_errorbar(aes(ymin=fn_mean - fn_std, ymax=fn_mean + fn_std), alpha = .6) +
  geom_point() +
  xlab('Distance to query genome (SNPs)') +
  ylab('False negative errors (SNPs)') +
  scale_x_continuous(trans='log10') +
  theme(legend.position="none")+ 
  theme(axis.title.x=element_blank())
g2

# Plot F1 score vs. dist to ref. 
g3 <- df_dists %>%
  filter(mapper == 'bwa') %>%
  filter(mapper == 'bwa' & caller == 'gatk') %>%
  ggplot(aes(x = dist, y = f1_mean)) +
  geom_errorbar(aes(ymin=f1_mean - f1_std, ymax=f1_mean + f1_std), alpha = .6) +
  geom_point() +
  xlab('Distance to query genome (SNPs)') +
  ylab('F1 Score') +
  scale_x_continuous(trans='log10') +
  theme(legend.position="none")+ 
  theme(axis.title.x=element_blank())
g3

# Plot three together. Use single x-axis. 
g4 <- plot_grid(g1, g2, g3, labels = c('a','b', 'c'), ncol = 3)

# Create grouped x-axis. 
x.grob <- textGrob('Distance to query genome (SNPs)', 
                   gp=gpar( fontsize=15))

# Add x-axis to plot
g5 <- grid.arrange(arrangeGrob(g4, bottom = x.grob))

# Add margins.
g5 <- grid.arrange(grobs = g5, vp=viewport(width=0.9, height=0.9))
g5

#################################
## Plot pairwise performance  ##
################################

# Read csv. 
df_pairs <- read.csv("pairwise_performance_H37Rv.csv", stringsAsFactors = F)

# Rename regions
df_pairs$region <- mapvalues(df_pairs$region, from = c('genome','noppe','ppe'), to = c('Genome','Genome - PE/PPE','PE/PPE'))
df_pairs$region <- factor(df_pairs$region, levels =  c('Genome','PE/PPE','Genome - PE/PPE'))
df_pairs$filter <- mapvalues(df_pairs$filter, from = c('raw', 'depth','qual','vqsr'), to = c('Raw','Depth','Qual','VQSR'))
df_pairs$filter <- factor(df_pairs$filter, levels = c('Raw','Depth','Qual','VQSR'))
df_pairs$caller <- mapvalues(df_pairs$caller, from = c("deep", "gatk","samtools"), to = c("DeepVariant", "GATK", "Samtools"))
df_pairs$mapper <- mapvalues(df_pairs$mapper, from = c('bowtie2',"bwa",'smalt'), to = c("Bowtie 2","BWA",'SMALT'))

# Check for missing data. All data present
ddply(df_pairs, .(mapper,caller, filter), summarize, N=length(region))

# reshape for plotting.
df_pairs <- reshape(df_pairs,  idvar = c('mapper','caller','filter','region'), timevar = 'var',  direction = 'wide')

# Remove Depth filter
df_pairs <- df_pairs[-which(df_pairs$filter == 'Depth'),]

## Plot FPs, colored by genomic region.
df_pairs$pipeline <- paste(df_pairs$mapper,  "\n ",df_pairs$caller, sep = '')

# Add error bars from genome-wide errors.
error_bars <- df_pairs %>%
  filter(region == 'Genome') %>%
  group_by(pipeline,filter) %>%
  mutate(genome_fp_max = mean.fp[which(region == 'Genome')] + std.fp[which(region == 'Genome')] ,
         genome_fn_max = mean.fn[which(region == 'Genome')]  + std.fn[which(region == 'Genome')] ) %>%
  dplyr::select(pipeline,filter,genome_fp_max, genome_fn_max)

# Add error bars
df_pairs_errors <- merge(df_pairs, error_bars)

# False positive errors.
g1 <- df_pairs_errors %>%
  filter(region != 'Genome') %>% 
  ggplot(aes(x = filter, y = mean.fp, fill = region)) + #fill = pipeline,
  geom_errorbar(aes(ymin = mean.fp, ymax = genome_fp_max)) +
  geom_bar(stat = 'identity', color="black") + 
  theme_minimal() +
  facet_grid(cols = vars(pipeline)) +
  xlab('Filter') + 
  ylab('False positive pairwise errors (SNPs)') +
  guides(colour = guide_legend(order = 2)) + 
  theme(text = element_text(size=13), axis.text.x = element_text(angle = 90, size = 10), strip.text = element_text(size = 8)) + 
  scale_fill_discrete(name = 'Region') + 
  geom_hline(yintercept=5, color = 'red', lty = 'dotted')+
  geom_hline(yintercept=12, color = 'red', lty = 'dotted') 
g1

# False negatives. 
g2 <- df_pairs_errors %>%
  filter(region != 'Genome') %>% 
  ggplot(aes(x = filter, y = mean.fn, fill = region)) + 
  geom_errorbar(aes(ymin = mean.fn, ymax = genome_fn_max)) +
  geom_bar(stat = 'identity', color="black") + 
  theme_minimal() +
  facet_grid(cols = vars(pipeline)) +
  xlab('Filter') + 
  ylab('False negative pairwise errors (SNPs)') +
  guides(colour = guide_legend(order = 2)) + 
  theme(text = element_text(size=13), axis.text.x = element_text(angle = 90, size = 10), strip.text = element_text(size = 8)) + 
  scale_fill_discrete(name = 'Region') 
g2

g3 <- plot_grid(g1,g2, ncol = 1, labels = c('a','b'))
g3