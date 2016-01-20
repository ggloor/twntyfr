library(ALDEx2)
# you can get this from:
# https://github.com/DavidRLovell/propr
source("proprBayes/R/propr-functions.R")
library(igraph)

# this is the refseq set
d <- read.table("../twntyfr/analysis/merged_readcounts_subsys4.txt", header=T,
   row.names=2, check.names=F, sep="\t", comment.char="", quote="")

sub4 <- d$subsys4
d$subsys4 <- NULL

refseq <- d$refseqIDfaa
d$refseqIDfaa <- NULL

len <- d$length
d$length <- NULL

d.min <- d[which(apply(d,1,mean)> 5),]

# generate random DIR clr instances with ALDEx2
e.x <- aldex.clr(d.min)

# calculate phi divided by number of random instances
d.min.sma.df <- aldex.phi(e.x)
write.table(d.min.sma.df, file="aldex.phi.df.txt", sep="\t", quote=F, col.names=NA)
# find the set of connections with phi less than some value
# we choose an arbitrary cutoff, but it is higher after Bayesian estimation
# obviously
 phi.cutoff <- 0.03

d.min.sma.lo.phi <- subset(d.min.sma.df, phi < phi.cutoff)

# generate a graphical object
g <- graph.data.frame(d.min.sma.lo.phi, directed=FALSE)

# overview of all the proportional relationships
# this can take a long time!!!
plot(g, layout=layout.fruchterman.reingold.grid(g, weight=0.05/E(g)$phi), vertex.size=1,
    vertex.color="black", vertex.label=NA)

# # get the clusters from the graph object
 g.clust <- clusters(g)

# write them to a file

# # data frame containing the names and group memberships of each cluster
g.df <- data.frame(Systematic.name=V(g)$name, cluster=g.clust$membership,
    cluster.size=g.clust$csize[g.clust$membership])
write.table(g.df, file="g.df.txt", sep="\t", quote=F, col.names=NA)

max.clust <- induced.subgraph(g, which(g.clust$membership %in% which(g.clust$csize ==
    max(g.clust$csize))))
max.clust.names <- V(max.clust)$name


