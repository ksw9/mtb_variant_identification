##############################
## Part A: Roetzer analysis ##
##############################

## This script compares variants identified by 5 different groups within the same sequence data from a clonal M. tuberculosis outbreak. 
## Sequence data from Roetzer et al. 2013 Plos Med, (ENA Study Accession: PRJEB6945).

rm(list = ls())

## Set up workspace.
library(ape)
library(ggtree)
library(treespace)
library(treeio)
library(vcfR)
library(pegas)
library(VennDiagram)
library(dplyr)
library(ggplot2)
library(plyr)
library(reshape2)
library(gridExtra)
library("adegenet")
library("adegraphics")
library("rgl")
library(treespace)
library(ips)
library(adephylo)
library(seqinr)
library(readxl)
library(r2d3)
library(gg3D)
library(cowplot)
library(viridis)
library(gridExtra)

# Set wd.
#setwd()

##############
## Organize ##
##############

## Read in multi-sample VCF files. And get SNP positions. 
Av <- as.numeric(read.vcfR('vcfs/A.vcf.gz', verbose = FALSE )@fix[,2])
Bv  <- as.numeric(read.vcfR('vcfs/B.vcf.gz', verbose = FALSE )@fix[,2])
Cv  <- as.numeric(read.vcfR('vcfs/C.vcf.gz', verbose = FALSE )@fix[,2])
Dv  <- as.numeric(read.vcfR('vcfs/D.vcf.gz', verbose = FALSE )@fix[,2])

# Roetzer sites (variants identified by original study and Sanger-confirmed).
Ev <- as.numeric(as.character(read.csv('metadata/roetz-tabS1.csv', header=F, stringsAsFactors=F)[1,2:86]))

## Read in FASTA files for each group with ape.
Afa <- read.dna('fastas/A.fa', format='fasta')
Bfa <- read.dna('fastas/B.fa', format = 'fasta')
Cfa <- read.dna('fastas/C.fa', format = 'fasta') 
Dfa <- read.dna('fastas/D.fa', format = 'fasta') 
Efa <- read.dna('fastas/E.fa', format = 'fasta') 

## Read in tree files (*raxml.support) for each group.
At <- unroot(read.tree(file = 'trees/A.fa_GTR.raxml.support'))
Bt <- unroot(read.tree(file = 'trees/B.fa_GTR.raxml.support'))
Ct <- unroot(read.tree(file = 'trees/C.fa_GTR.raxml.support'))
Dt <- unroot(read.tree(file = 'trees/D.fa_GTR.raxml.support'))
Et <- unroot(read.tree(file = 'trees/E.fa_GTR.raxml.support'))

## Read in bootstrapped trees. 
Abs <- unroot(read.tree(file = 'trees/A.fa_GTR.raxml.bootstraps'))
Bbs <- unroot(read.tree(file = 'trees/B.fa_GTR.raxml.bootstraps'))
Cbs <- unroot(read.tree(file = 'trees/C.fa_GTR.raxml.bootstraps'))
Dbs <- unroot(read.tree(file = 'trees/D.fa_GTR.raxml.bootstraps'))
Ebs <- unroot(read.tree(file = 'trees/E.fa_GTR.raxml.bootstraps'))

# Read in metadata from Roetzer et al. 2013, Plos Med. 
meta <- read.csv('metadata/roetz-tabS2.csv', stringsAsFactors = F)

#########################################
## Compare VCFs and pairwise distances ##
#########################################

#### Plot Venn diagram & SNP totals. ####
# make quintuple Venn diagram
area1 <- length(Av)
area2 <- length(Bv) + 2 # Add 2 snps called by pipeline B on  (could not be lifted over)
area3 <- length(Cv)
area4 <- length(Dv)
area5 <- length(Ev)

# Get total internal variants.
pipelines <-  c("A","B","C","D","E")
totals <- data.frame(pipelines, c( area1,area2,area3,area4,area5))
names(totals) <- c('pipeline', 'SNPs')
totals

# Define color palette.
pal <- viridis(5)
names(pal) <- pipelines

# Replot SNP totals with pipeline color.
g1 <- ggplot(totals, aes(x = pipeline, y = SNPs, fill = pipeline)) + 
  geom_col()  + 
  theme_minimal()+
  ylab('Total internal SNPs') + 
  xlab('Pipeline') + 
  scale_fill_manual(values = pal,  name = 'Pipeline') +
  theme_classic(base_size = 14)  
g1

## Plot Venn of shared SNP calls.
n12 <- length(intersect(Av, Bv))
n13 <- length(intersect(Av, Cv))
n14 <- length(intersect(Av,Dv))
n15 <- length(intersect(Av, Ev))
n23 <- length(intersect(Bv, Cv))
n24 <- length(intersect(Bv,Dv))
n25 <- length(intersect(Bv, Ev))
n34 <- length(intersect(Cv,Dv))
n35 <- length(intersect(Cv, Ev))
n45 <- length(intersect(Dv, Ev))
n123  <- length(Reduce(intersect,list(Av, Bv, Cv)))
n124 <- length(Reduce(intersect,list(Av, Bv,Dv)))
n125 <- length(Reduce(intersect,list(Av, Bv, Ev)))
n134 <- length(Reduce(intersect,list(Av, Cv,Dv)))
n135 <- length(Reduce(intersect,list(Av, Cv, Ev)))
n145 <- length(Reduce(intersect,list(Av,Dv, Ev)))
n234 <- length(Reduce(intersect,list(Cv, Bv,Dv)))
n235 <- length(Reduce(intersect,list(Cv, Bv, Ev)))
n245 <- length(Reduce(intersect,list(Bv,Dv, Ev)))
n345 <- length(Reduce(intersect,list(Cv,Dv, Ev)))
n1234 <- length(Reduce(intersect,list(Av, Bv, Cv,Dv)))
n1235 <- length(Reduce(intersect,list(Av, Bv, Cv, Ev)))
n1245 <- length(Reduce(intersect,list(Av, Bv,Dv, Ev)))
n1345 <- length(Reduce(intersect,list(Av, Cv,Dv, Ev)))
n2345 <- length(Reduce(intersect,list(Bv, Cv,Dv, Ev)))
n12345 <- length(Reduce(intersect,list(Av, Bv, Cv,Dv, Ev)))

# plot venn diagram. 
venn <- draw.quintuple.venn(area1 = area1, area2 = area2, area3 = area3, 
                            area4 = area4, 
                            area5=area5, 
                            n12 = n12, n13 = n13, n14 = n14, n15 = n15,
                            n23 = n23, n24 = n24, 
                            n25 = n25,
                            n34 = n34, 
                            n35 = n35, n45 = n45,
                            n123 = n123, n124 = n124, n125 = n125, n134 = n134, 
                            n135 = n135, n145 = n145, 
                            n234 = n234, 
                            n235 = n235, n245 = n245, n345 = n345, 
                            n1234 = n1234, 
                            n1235 = n1235, n1245 = n1245,n1345 = n1345, n2345 = n2345,
                            n12345 = n12345,
                            category = pipelines,
                            fill = pal,
                            lty = "dashed",
                            cex = 1,
                            cat.cex = 2,
                            cat.col = pal,
                            #cat.col = c("orange", "red", "green", "blue", "grey"), scaled = TRUE,
                            margin = 0.05, 
                            alpha = rep(.7,5))
# Plot.
grid.draw(venn)

#### Sensitivity in detecting Roetzer sites ####
# Length of Roetzer Sanger-confirmed SNPs.
n <- length(Ev)

# Sensitivity.
sens_A <- (n - length(setdiff(Ev, as.numeric(Av))))/ n
sens_B <- (n - length(setdiff(Ev, as.numeric(Bv))))/ n
sens_C <- (n - length(setdiff(Ev, as.numeric(Cv)))) /n
sens_D <- (n - length(setdiff(Ev, as.numeric(Dv)))) /n

# Get total internal variants.
totals$sensitivity <-  c(sens_A,sens_B,sens_C,sens_D, 1)*100

# plot sensitivity 
g2 <-  ggplot(totals[which(totals$pipeline != 'E'),], aes(x = pipeline, y = sensitivity, fill = pipeline)) + 
  geom_col()  + 
  theme_minimal() +
  ylab('Sensitivity (%) ') + 
  xlab('Pipeline') +
  scale_fill_manual(values = pal, name = 'Pipeline') +
  ylim(0,100)+
  theme_classic(base_size = 14) 
g2

# Get legend from g1
legend <- get_legend(
  # create some space to the left of the legend
  g1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# Arrange the three plots without legends.
p1 <- plot_grid(
  g1 + theme(legend.position="none"),
  venn, 
  g2 + theme(legend.position="none"),
  legend,
  labels = c('a','b','c'),
  ncol = 2, 
  hjust = -1
)
p1

# Save
#ggsave(p1, file = 'plots/fig1.pdf', width = 10, height = 10)

#### Pairwise distance plots. ####

# For each msa, get pairwise distances. 
Ad <- ape::dist.dna(Afa, model = 'N' )
Bd <- ape::dist.dna(Bfa, model = 'N' )
Cd <- ape::dist.dna(Cfa, model = 'N' )
Dd <- ape::dist.dna(Dfa, model = 'N' )
Ed <-  ape::dist.dna(Efa, model = 'N' )

# Create dataframe for plotting: pairwise distances and pipeline
d1 <- data.frame(c(Ad), 'A')
names(d1) <- c('distance', 'pipeline')
d2 <- data.frame(c(Bd), 'B')
names(d2) <- c('distance', 'pipeline')
d3 <- data.frame(c(Cd), 'C')
names(d3) <- c('distance', 'pipeline')
d4 <- data.frame(c(Dd), 'D')
names(d4) <- c('distance', 'pipeline')
d5 <- data.frame(c(Ed), 'E')
names(d5) <- c('distance', 'pipeline')

d <- bind_rows(d1,d2,d3,d4,d5)

# Violin plot of pairwise distances.
g3 <- ggplot(d, aes(y = distance, x = pipeline, fill = pipeline)) + 
  geom_violin() + 
  ylab('Pairwise SNP distances') +
  xlab('Pipeline') +
  geom_hline(yintercept=12, lty = 4) +
  geom_hline(yintercept=5, lty = 4) +
  scale_y_continuous(trans='log10') +
  theme_classic() +
  scale_fill_manual(values = pal )+
  theme(text = element_text(size=16)) + 
  theme(legend.position = "none") 
g3

# Count the number of pairs within transmission thresholds for each caller.
thresholds <- ddply(d, .(pipeline), summarise,
                    n = length(distance),
                    mean = round(mean(distance), digits = 2),
                    median = round(median(distance),digits = 2), 
                    identical = sum(distance == 0)/n * 100,
                    threshold5 = sum(distance <= 5)/n * 100, 
                    threshold12 = sum(distance <= 12)/n * 100)

# Add # of samples called to table. 
thresholds$samples <- c(86,86,68,86,86)

# Melt table for plotting. 
thresholds2 <- melt(thresholds[,c('pipeline','threshold5', 'threshold12')])
thresholds2$variable <- mapvalues(thresholds2$variable, from = c('threshold5', 'threshold12'), to = c('5 SNPs', '12 SNPs'))
thresholds2$variable <- factor(thresholds2$variable,levels(thresholds2$variable)[c(2,1)])

# order factor 
thresholds2$variable <- factor(thresholds2$variable, levels = c("5 SNPs", "12 SNPs"))

# Define 2-color palette. 
pal2 <- viridis(10)[c(2,8)]
names(pal2) <- unique(thresholds2$variable)

# Barplot of those distances falling under a threshold. 
g4 <- ggplot(thresholds2, aes(y = value, fill = variable, x = pipeline)) +
  geom_bar(stat = 'identity', position=position_dodge()) + 
  ylab('Potential transmission pairs (%)' ) +
  xlab('Pipeline' ) +
  theme_classic() +
  theme(legend.title=element_blank())  +
  theme(text = element_text(size=16)) + 
  scale_fill_manual(values = pal2)

g4

## Pairwise comparison of differences between A and B pipelines (all samples were called for both).
B_labels <- sapply(labels(Bfa), function(x) strsplit(x,'_')[[1]][1])

# Get ordering of B labels.
A_order <- sapply(B_labels, function(x) which(labels(Afa) == x)[1])

# Check order is correct: 
labels(Afa)[A_order] 

# Reorder A.
A_mat <- as.matrix(Ad)[A_order,A_order]

# Convert back to dist. 
A_mat <- as.dist(A_mat)

# Sort distances by sample. 
thresh12 <- 12
thresh5 <- 5
A_B_dists <- data.frame(cbind(as.vector(A_mat), as.vector(Bd)))

g5 <- ggplot(A_B_dists, aes(x = X1, y = X2)) + 
  geom_jitter(cex = .5) +
  xlab('Pipeline A pairwise distance (SNPs)') +
  ylab('Pipeline B pairwise distance (SNPs)') +
  geom_hline(yintercept=thresh5, color=pal2[1], lty = 4) + 
  geom_vline(xintercept=thresh5, color=pal2[1], lty = 4)+
  geom_hline(yintercept=thresh12, color=pal2[2], lty = 4) + 
  geom_vline(xintercept=thresh12, color=pal2[2], lty = 4) +
  # shade area in which inference is changed
  annotate("rect", ymin=0, ymax=thresh5, xmin = thresh5, xmax=15, alpha="0.5", fill= pal2[1])  +
  annotate("rect", ymin=thresh5, ymax=20, xmin = 0, xmax=thresh5, alpha="0.5", fill= pal2[1])  +
  annotate("rect", ymin=0, ymax=thresh12, xmin = thresh12, xmax=15, alpha="0.5", fill= pal2[2])  + 
  annotate("rect", ymin=thresh12, ymax=20, xmin = 0, xmax = thresh12, alpha="0.5", fill= pal2[2])  + 
  theme_classic() +
  theme(text = element_text(size=14))
g5

# What is correlation between pairwise distances for both pipelines?
A_B_cor <- cor.test(A_B_dists$X1,A_B_dists$X2) # 89.1 %
A_B_cor 

# How many pairs are linked by both callers?
length(which(A_mat <= 12 & Bd <= 12))
#3526

# How many discordant calls at 5 SNP threshold?
length(which(A_mat <= 5 & Bd > 5))
# 413

# How many discordant calls in other direction? 
length(which(A_mat >5 & Bd <=5))
# 14

# Cumulative discordant calls at 5 SNP threshold.
(length(which(A_mat <= 5 & Bd > 5)) +  length(which(A_mat >5 & Bd <=5))) / length(A_mat) # 11.68 % discordant calls. 

# Create table of important statistics.
t <- merge(totals, thresholds)

# Select rows for table
tt <- t[,c(1,10,2:3,5:9)]
colnames(tt) <- c('Pipeline','Samples', 'SNPs', 'Sensitivity','Mean pairwise', 'Median pairwise', 'Identical (%)', 
                  "<= 5 SNPs (%)", "<= 12 SNPs (%)") 

# Round to 1 decimal. 
tt[,c(4,5,7:9)] <- format(round(tt[,c(4,5,7:9)], digits=1), nsmall = 1) 
row.names(tt) <- NULL
tt

# Save table.
pdf(file ='plots/roetzer_table.pdf', height = 2, width = 10)
grid.table(tt, rows = NULL)
dev.off()

####################
## Plot tree sets ##
####################

# Name bootstrap treesets by pipeline. 
names(Abs) <- paste(rep('A', length(Abs)), 1:length(Abs), sep = '_')
names(Bbs) <- paste(rep('B', length(Bbs)), 1:length(Bbs), sep = '_')
names(Cbs) <- paste(rep('C', length(Cbs)), 1:length(Cbs), sep = '_')
names(Dbs) <- paste(rep('D', length(Dbs)), 1:length(Dbs), sep = '_')
names(Ebs) <- paste(rep('E', length(Ebs)), 1:length(Ebs), sep = '_')

# Rename tip labels for pipeline B. 
orig_labels <- Bbs[[1]]$tip.label
new_labels <- sub(pattern = '_.*',replacement = '', Bbs[[1]]$tip.label)

# Rename function. 
rename_tips <- function(x) { mapvalues(x, from = orig_labels, to = new_labels) }

# Apply function. 
for (nn in 1:length(Bbs)){
  print(nn)
  Bbs[[nn]]$tip.label <- rename_tips(Bbs[[nn]]$tip.label )
}

# Rename tip labels for pipeline D. 
orig_labels <- Dbs[[1]]$tip.label
new_labels <- sub(pattern = '_.*',replacement = '', Dbs[[1]]$tip.label)

# Apply function. 
for (nn in 1:length(Dbs)){
  print(nn)
  Dbs[[nn]]$tip.label <- rename_tips(Dbs[[nn]]$tip.label )
}

# Rename tip labels for pipeline E.
tmp <- sub(pattern = '-', replacement = '/', Ebs[[1]]$tip.label)
orig_labels <- Ebs[[1]]$tip.label
new_labels <- meta$run_accession[sapply(tmp, function(x) which(meta$KEY == x))]

# Apply function. 
for (nn in 1:length(Ebs)){
  print(nn)
  Ebs[[nn]]$tip.label <- rename_tips(Ebs[[nn]]$tip.label )
}

# Create set of trees with bootstrap (for memory: sample 50 bs replicates from each pipeline). 
tset_bs <- c(sample(Abs, 100), sample(Bbs,100), sample(Dbs,100), sample(Ebs,100))

# Tree names
names(tset_bs) <- c(names(Abs)[1:100], names(Bbs)[1:100], names(Cbs)[1:100], names(Dbs)[1:100])

# Root trees on earliest isolate from 1997. 
tset_rooted <- lapply(tset_bs, function(x) root(x, 'ERR552368', resolve.root = TRUE))
class(tset_rooted) <- 'multiPhylo'

# Get tree distances with RF distance for unrooted trees.
res_rf <- treespace(tset_bs, nf=5, 'RF') 

# Identify clusters in bs trees. 
rf_groves <- findGroves(res_rf, nclust=4)

# Table group assignments by pipeline.
rf_groves$pipeline <- c(rep('A',100), rep('B',100), rep('D',100), rep('E', 100))
table(rf_groves$groups, rf_groves$pipeline)

# Extract x,y,z coordinates.
scatter_rf <- cbind(rf_groves$treespace$pco$li, Cluster = rf_groves$groups)

# Define pipeline. 
scatter_rf$Pipeline <- c(rep('A',100), rep('B',100), rep('D',100), rep('E', 100))

# 3-D scatterplot: RF distances.
g9 <- ggplot(scatter_rf, aes(x=A1, y=A2, z=A3, shape=Cluster, color = Pipeline)) + 
  theme_void() +
  axes_3D() +
  scale_color_manual(values = c(pal)) +
  stat_3D() 
g9

# Save figure 2.
fig2 <- plot_grid(g3,g4,g5,g9, labels = c('a','b','c','d'), scale = 0.9)
fig2

# Save.
#ggsave(fig2, file = 'plots/fig2.pdf', width = 10, height = 10)