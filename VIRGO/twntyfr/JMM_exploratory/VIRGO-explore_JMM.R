# Exploratory for twntyfr metatranscriptomics
# Building ggplot2-style colored biplots. SEE FINAL PLOT AT END

# 10-Feb-2020
# Jean M Macklaim
#-------------------------------------------------------------------------------	

# See also
# /Users/jean.macklaim/Downloads/VIRGO-explore.R

source("common_functions/load_packages.r")

metadata<-read.table("metadata_full_cleaned.txt", sep="\t", quote="", header=T, check.names=F, row.names=1, comment.char="",stringsAsFactors=F)
gene <- read.table("summary.NR.abundance.txt", sep="\t", quote="", header=T, check.names=F, row.names=1, comment.char="", stringsAsFactors=F)

dim(gene)
# [1] 249233     24

d.1<-Confidence.Filter(gene, 2, 10, VERBOSE=T)
# Filtering features such that they are present in at least 2 samples with a total of at least 10 reads.
# ...There are 259489615 reads and 249233  features
# ...After filtering there are 257624889 reads and 154054 features

dim(d.1)
# [1] 154054     24

gene.filt <- codaSeq.filter(gene, min.occurrence=0.1, min.prop=0.00005, samples.by.row=F)
dim(gene.filt)
# [1] 23378    24

d.1<-gene.filt

#------------------------------------------------------------------------------
# CoDa
#------------------------------------------------------------------------------

d.czm <- cmultRepl(as.matrix(t(d.1)), label=0, method="CZM")

d.clr <- t(apply(d.czm, 1, function(x){log(x) - mean(log(x))}))

d.pca<-prcomp(d.clr)

# Keep only metadata rows that match the samples
# Some samples were dropped along the way
new.metadata <- metadata[rownames(metadata) %in% colnames(d.1),]
# Careful to keep the same order
new.metadata<-new.metadata[colnames(d.1),]

# https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
help(autoplot.prcomp)

# Plot by sample
autoplot(d.pca)
# plot featrues
autoplot(d.pca,  scale = 0)

# Select which coomponent to plot using x= and y=
# Careful, the metadata and the d.pca have to be in the same order
p1<-autoplot(d.pca, data = new.metadata, colour = 'nugent_status', label=TRUE, label.size = 3)
p2<-autoplot(d.pca, data = new.metadata, colour = 'nugent_status', label=TRUE, label.size = 3, x=1, y=3)

pdf("PCA_PC1PC2_donorcolor.pdf")
print(p1)
dev.off()

pdf("PCA_PC1PC3_donorcolor.pdf")
print(p2)
dev.off()

pdf("PCA_PC2PC3_donorcolor.pdf")
print(p3)
dev.off()

#-------------------------------------------------------------------------------
# Base R version
plot(d.pca$rotation[,1], d.pca$rotation[,2], pch=19, col=rgb(0,0,0,0.01))
text(d.pca$x[,1]/30000, d.pca$x[,2]/30000, labels=rownames(d))

# ggplot version/wrapper
# Make the plot, add points with transparancy
p<-ggplot(d.pca$rotation, aes(x=PC1, y=PC2)) +
 geom_point(alpha = 1/20)

# white background
p + theme_bw()

# Add sample points
p<-ggplot(d.pca$rotation, aes(x=PC1, y=PC2)) +
 geom_point(alpha = 1/20) +
 geom_point(data = d.pca$x, aes(x=PC1/30000, y=PC2/30000, colour = "red")
)
p + theme_bw()

#-------------------------------------------------
# Final version, using ggplot2
#-------------------------------------------------

# Get the % PCA for labels
d.mvar <- sum(d.pca$sdev^2)

PC1<-round(sum(d.pca$sdev[1]^2)/d.mvar, 4)*100
PC2<-round(sum(d.pca$sdev[2]^2)/d.mvar, 4)*100
PC3<-round(sum(d.pca$sdev[3]^2)/d.mvar, 4)*100


# Add sample labels, and color by metadata
# Need to put metadata on the same table
d2<-cbind(d.pca$x, new.metadata)

xlab<-paste("PC1: ", PC1, "%", sep="")
ylab<-paste("PC2: ", PC2, "%", sep="")

p<-ggplot() +
	geom_point(data = d.pca$rotation, aes(x=PC1, y=PC2), alpha = 1/20) +
	geom_point(data = d2, aes(x=PC1/30000, y=PC2/30000, colour = nugent_status)) +
	geom_text(data = d2, aes(x=PC1/30000, y=PC2/30000, label = rownames(d2), colour = nugent_status), hjust = 0, nudge_x = 0.0005
)
p + ggtitle("PCA") +
  xlab(xlab) + ylab(ylab) +
  theme_bw()

# NOTE: there is an arbitrary scaling factor (30000) that must be adjusted for the dataset

# Other parameters:
# check.overlap = TRUE to prevent label overlap. But labels will be removed if overlapping