#------------------
# HEATMAPS
#------------------

load("vsd.RData")
load("pvals.RData")
source("uniHeatmap.R")
# loading annotations - two-column tables
gg=read.table("Mcavernosa_Cladocopium_iso2geneName.tab",sep="\t", quote="", comment.char="")  # gene names
str(gg)
#gg=gg[-grep("Uncharacterized",gg)$V2,]
kogs=read.table("Mcavernosa_Cladocopium_iso2kogClass.tab",sep="\t", quote="", comment.char="")  # KOG classes
str(kogs)

design$time.site <- paste(design$time,design$site,sep=".")
design$time.site <- factor(design$time.site, levels = c("0.south","0.ledge","0.central","1.south","1.ledge","1.central","2.south","2.ledge","2.central","3.south","3.ledge","3.central","4.south","4.ledge","4.central","5.south","5.ledge","5.central","6.south","6.ledge","6.central"))
design[order(design$time.site),]

# reordering columns according to time.site
head(vsd)
vsd=vsd[,order(design$time.site)]
head(vsd)

#--------------
# all genes passing a cutoff, for treatment
# use for identifying large patterns, like genotypic differences
# install.packages("pheatmap")
library(pheatmap)

#--------------
# time pvals

pdf(file="heatmap_time_0.05.pdf", height=8, width=15)
uniHeatmap(vsd=vsd,gene.names=gg,
           metric=time$pvalue, # metric of gene significance
           cutoff=0.05, 
           sort=c(1:ncol(vsd)),   # overrides sorting of columns according to hierarchical clustering
           cex=0.8,
           pdf=F
)
dev.off()

pdf(file="heatmap_time_0.01.pdf", height=4.5, width=15)
uniHeatmap(vsd=vsd,gene.names=gg,
           metric=time$pvalue, # metric of gene significance
           cutoff=0.01, 
           sort=c(1:ncol(vsd)),   # overrides sorting of columns according to hierarchical clustering
           cex=0.8,
           pdf=F
)
dev.off()

#------------------------------
# site pvals

pdf(file="heatmap_site_0.05.pdf", height=3 , width=13)
uniHeatmap(vsd=vsd,gene.names=gg,
           metric=site$pvalue, # metric of gene significance
           cutoff=0.05, 
           sort=c(1:ncol(vsd)),   # overrides sorting of columns according to hierarchical clustering
           cex=0.8,
           pdf=F
)
dev.off()

#------------------------------
# int pvals

# no p < 0.1

#------------------------------