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

design$temp.water <- factor(design$temp.water, levels = c("control.offshore","control.discharge","elevated.offshore","elevated.discharge"))
design[order(design$temp.water),]

# reordering columns according to temp.water
head(vsd)
vsd.tw=vsd[,order(design$temp.water)]
head(vsd.tw)

#--------------
# all genes passing a cutoff, for treatment
# use for identifying large patterns, like genotypic differences
# install.packages("pheatmap")
library(pheatmap)

#--------------
# colony pvals

pdf(file="heatmap_colony_0.05.pdf", height=3, width=30)
uniHeatmap(vsd=vsd.tw,gene.names=gg,
           metric=colony$pvalue, # metric of gene significance
           cutoff=0.05, 
           sort=c(1:ncol(vsd)),   # overrides sorting of columns according to hierarchical clustering
           cex=0.8,
           pdf=F
)
dev.off()

#------------------------------
# temp pvals

pdf(file="heatmap_temp_0.05.pdf", height=5, width=30)
uniHeatmap(vsd=vsd.tw,gene.names=gg,
           metric=temp$pvalue, # metric of gene significance
           cutoff=0.05, 
           sort=c(1:ncol(vsd)),   # overrides sorting of columns according to hierarchical clustering
           cex=0.8,
           pdf=F
)
dev.off()

#------------------------------
# water pvals

pdf(file="heatmap_water_0.05.pdf", height=2, width=20)
uniHeatmap(vsd=vsd.tw,gene.names=gg,
	metric=water$pvalue, # metric of gene significance
	cutoff=0.05, 
	sort=c(1:ncol(vsd)),   # overrides sorting of columns according to hierarchical clustering
	cex=0.8,
	pdf=F
	)
dev.off()

#------------------------------
# temp:water pvals

pdf(file="heatmap_tw_0.05.pdf", height=1.5, width=20)
uniHeatmap(vsd=vsd.tw,gene.names=gg,
           metric=tw$pvalue, # metric of gene significance
           cutoff=0.05, 
           sort=c(1:ncol(vsd)),   # overrides sorting of columns according to hierarchical clustering
           cex=0.8,
           pdf=F
)
dev.off()

#------------------------------