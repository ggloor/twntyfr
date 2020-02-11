# Load standard packages for analyses
# Jean M Macklaim

library(MicrobeR)	# Functions for microbiome analyses: filtering, plotting barplots
library(zCompositions) # Only needed for CLR zero estimates (CZM method) if NOT using ALDEx2 estimates
library(ALDEx2)	# For CLR transformations and differential expression analyses
library(ape)	# Phylogenetic plots (dendrograms)
library(tidyverse)
library(ggplot2)	# Figure generation
library(dplyr)	# Data manipulation in tidyverse
library(gplots)
library(RColorBrewer)	# Plotting colorschemes
library(phyloseq)	# Microbiome analyses
library(CoDaSeq)	# Functions for compositonal data anlyses
library(cowplot)	# For combining multiple ggplot figures


# [ ] Write quick descriptions for what each package is used for

