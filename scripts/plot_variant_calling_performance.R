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

###############
## Read data ##
###############


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
df_dists <- read.csv('performance/genome_performance_allrefs.csv', stringsAsFactors = F)

# Plot FPs vs. distance to ref.
g1 <- df_dists %>%
  filter(mapper == 'bwa' & caller == 'gatk') %>%
  ggplot(aes(x = dist, y = fp_mean, color = mapper)) +
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