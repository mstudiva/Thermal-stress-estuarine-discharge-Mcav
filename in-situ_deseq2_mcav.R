# run these once, then comment out
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.10")
# BiocManager::install("DESeq2",dependencies=T)
# BiocManager::install("arrayQualityMetrics",dependencies=T)  # requires Xquartz, xquartz.org
# BiocManager::install("BiocParallel")

# install.packages("pheatmap")
# install.packages("VennDiagram")
# install.packages("gplots")
# install.packages("vegan")
# install.packages("plotrix")
# install.packages("ape")
# install.packages("ggplot2")
# install.packages("rgl")
# install.packages("adegenet")

#---------------------
# assembling data, running outlier detection, and fitting models
# (skip this section if you don't need to remake models)

library(DESeq2)
library(arrayQualityMetrics)

#read in counts
counts = read.table("insitu_allcounts_mcav.txt")

# how many genes we have total?
nrow(counts) 
ncol(counts)

# how does the data look? 
head(counts)

#---------------------
keep <- rowSums(counts) >= 10
countData <- counts[keep,]
nrow(countData)
ncol(countData)
write.csv(countData, file="countData.csv")

# importing a design .csv file
design = read.csv("insitu_design.csv", head=TRUE)
design
design$time = as.factor(design$time)
str(design)

# make big dataframe including all factors and interaction, getting normalized data for outlier detection
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ time*site)

dds.season = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ site*wetdry)

# reorders treatment factor according to "control" vs "treatment" levels
dds$time <- factor(dds$time, levels = c("0","1","2","3","4","5","6"))
dds$site <- factor(dds$site, levels = c("south","ledge","central"))

dds.season$site <- factor(dds.season$site, levels = c("south","ledge","central"))
dds.season$wetdry <- factor(dds.season$wetdry, levels = c("dry","wet"))

# for large datasets, rlog may take too much time, especially for an unfiltered dataframe
# vsd is much faster and still works for outlier detection
Vsd=varianceStabilizingTransformation(dds)
Vsd.season=varianceStabilizingTransformation(dds.season)

library(Biobase)
e=ExpressionSet(assay(Vsd), AnnotatedDataFrame(as.data.frame(colData(Vsd))))

# running outlier detection
arrayQualityMetrics(e,intgroup=c("site"),force=T)
# open the directory "arrayQualityMetrics report for e" in your working directory and open index.html
# Array metadata and outlier detection overview gives a report of all samples, and which are likely outliers according to the 3 methods tested. I typically remove the samples that violate *1 (distance between arrays).
# Figure 2 shows a bar plot of array-to-array distances and an outlier detection threshold based on your samples. Samples above the threshold are considered outliers
# under Figure 3: Principal Components Analyses, look for any points far away from the rest of the sample cluster
# use the array number for removal in the following section

# if there were outliers:
outs=c(1,8,17)
countData=countData[,-outs]
Vsd=Vsd[,-outs]
Vsd.season=Vsd.season[,-outs]
design=design[-outs,]

# remaking model with outliers removed from dataset
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ time*site)
dds$time <- factor(dds$time, levels = c("0","1","2","3","4","5","6"))
dds$site <- factor(dds$site, levels = c("south","ledge","central"))

dds.season = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ site*wetdry)
dds.season$site <- factor(dds.season$site, levels = c("south","ledge","central"))
dds.season$wetdry <- factor(dds.season$wetdry, levels = c("dry","wet"))

# save all these dataframes as an Rdata package so you don't need to rerun each time
save(dds,dds.season,design,countData,Vsd,Vsd.season,file="initial.RData")

#---------------------
# generating normalized variance-stabilized data for PCoA, heatmaps, etc

load("initial.RData")
library(DESeq2)
library(BiocParallel)

# creating normalized dataframe
vsd=assay(Vsd)
# takes the sample IDs and factor levels from the design to create new column names for the dataframe
snames=paste(colnames(countData),design[,4],design[,7],sep=".")
# renames the column names
colnames(vsd)=snames

vsd.season=assay(Vsd.season)
snames.season=paste(colnames(countData),design[,7],design[,9],sep=".")
# renames the column names
colnames(vsd.season)=snames.season

save(vsd,vsd.season,design,file="vsd.RData")

#-------------------
# EXPLORING SIMILARITIES AMONG SAMPLES

# heatmap and hierarchical clustering:
load("vsd.RData")
library(pheatmap)
# similarity among samples
pdf(file="heatmap_insitu_mcav.pdf", width=10, height=10)
pheatmap(cor(vsd))
dev.off()

# Principal coordinates analysis
library(vegan)
library(rgl)
library(ape)

conditions=design
conditions$time <- factor(conditions$time, levels = c("0","1","2","3","4","5","6"))
conditions$site <- factor(conditions$site, levels = c("south","ledge","central"))

# creating a PCoA eigenvalue matrix
dds.pcoa=pcoa(dist(t(vsd),method="manhattan")/1000)
scores=dds.pcoa$vectors

# how many good PC's do we have? Compared to random ("broken stick") model
# plotting PCoA eigenvalues 
pdf(file="PCoA_Manhattan.pdf", width=6, height=6)
plot(dds.pcoa$values$Relative_eig)
points(dds.pcoa$values$Broken_stick,col="red",pch=3)
dev.off()
# the number of black points above the line of red crosses (random model) corresponds to the number of good PC's

# plotting PCoA by site
pdf(file="PCoA_insitu_site_mcav.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=c("#91bfdb","#fee090","#fc8d59","#e0f3f8","#ffffbf","#d73027","#4575b4")[as.numeric(as.factor(conditions$time))],pch=c(15,16,17)[as.numeric(as.factor(conditions$site))], xlab="Coordinate 1", ylab="Coordinate 2", main="Time")
# cluster overlay of time
ordiellipse(scores, conditions$time, label=F, draw= "polygon", col=c("#91bfdb","#fee090","#fc8d59","#e0f3f8","#ffffbf","#d73027","#4575b4"))
legend("bottomleft", legend=c("Oct 13","Sep 14","Jun 15","Oct 15","Mar 16","Jul 16","Nov 16"), fill = c("#91bfdb","#fee090","#fc8d59","#e0f3f8","#ffffbf","#d73027","#4575b4"), bty="n")
legend("topleft", legend=c("Central","Ledge","South"), pch=c(15,16,17), bty="n")
plot(scores[,1], scores[,2],col=c("#d8b365","#f6e8c3", "#5ab4ac")[as.numeric(as.factor(conditions$site))],pch=c(0,1,2,5,6,7,8)[as.numeric((as.factor(conditions$time)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Site")
# cluster overlay of site
ordiellipse(scores, conditions$site, label=F, draw= "polygon", col=c("#5ab4ac","#f6e8c3","#d8b365"))
legend("topleft", legend=c("Central","Ledge","South"), fill = c("#d8b365","#f6e8c3", "#5ab4ac"), bty="n")
legend("bottomleft", legend=c("Oct 13","Sep 14","Jun 15","Oct 15","Mar 16","Jul 16","Nov 16"), pch=c(0,1,2,5,6,7,8), bty="n")
dev.off()

conditions$site <- factor(conditions$site, levels = c("south","ledge","central"))
conditions$wetdry <- factor(conditions$wetdry, levels = c("dry","wet"))

# plotting PCoA by season
pdf(file="PCoA_insitu_season_mcav.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=c("#d8b365","#f6e8c3", "#5ab4ac")[as.numeric(as.factor(conditions$site))],pch=c(2,17)[as.numeric(as.factor(conditions$wetdry))], xlab="Coordinate 1", ylab="Coordinate 2", main="Site")
# cluster overlay of site
ordiellipse(scores, conditions$site, label=F, draw= "polygon", col=c("#5ab4ac","#f6e8c3","#d8b365"))
legend("topleft", legend=c("Central","Ledge","South"), fill = c("#d8b365","#f6e8c3", "#5ab4ac"), bty="n")
legend("bottomleft", legend=c("dry","wet"), pch=c(2,17), bty="n")
plot(scores[,1], scores[,2],col=c("#d73027","#4575b4")[as.numeric(as.factor(conditions$wetdry))],pch=c(15,16,17)[as.numeric(as.factor(conditions$site))], xlab="Coordinate 1", ylab="Coordinate 2", main="Season")
# cluster overlay of season
ordiellipse(scores, conditions$wetdry, label=F, draw= "polygon", col=c("#d73027","#4575b4"))
legend("bottomleft", legend=c("dry","wet"), fill = c("#d73027","#4575b4"), bty="n")
legend("topleft", legend=c("Central","Ledge","South"), pch=c(15,16,17), bty="n")
dev.off()

# neighbor-joining tree of samples (based on significant PCo's):
pdf(file="PCoA_tree.pdf", width=10, height=10)
tre=nj(dist(scores[,1:4]))
plot(tre,cex=0.8)
dev.off()

# formal analysis of variance in distance matricies: 
ad=adonis(t(vsd)~time*site,data=conditions,method="manhattan",permutations=9999)
ad

ad.season=adonis(t(vsd.season)~site*wetdry,data=conditions,method="manhattan",permutations=9999)
ad.season

# creating pie chart to represent ANOVA results
cols=c("blue","orange","lightblue","grey80")
pdf(file="ANOVA_pie.pdf", width=6, height=6)
pie(ad$aov.tab$R2[1:4],labels=row.names(ad$aov.tab)[1:4],col=cols,main="time vs site")
dev.off()

#----------------------
# FITTING GENE BY GENE MODELS

# with multi-factor, multi-level design - using LRT
load("initial.RData")
library(DESeq2)
library(BiocParallel)

# Running full model for contrast statements
dds=DESeq(dds, parallel=TRUE)

# Running DESeq with LRT for the effect of interaction (using full model dds = same as design argument on line 63, ~ time*site)
# then colony+time are removed, to look just at interaction term
dds.i=DESeq(dds,test="LRT",reduced=~time+site, parallel=TRUE)

# Creating dataset with design =~time+site (no interaction), to investigate importance of time and site
dds1=DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ time+site)

# model for the effect of site (>2 factor levels => LRT)
dds1$site <- factor(dds1$site, levels = c("south","ledge","central"))
dds.site=DESeq(dds1,parallel=TRUE)

# model for the effect of time: (>2 factor levels => LRT)
dds1$time <- factor(dds1$time, levels = c("0","1","2","3","4","5","6"))
dds.time=DESeq(dds1,test="LRT",reduced=~site, parallel=TRUE)

# full model for site and season
dds.season=DESeq(dds.season, parallel=TRUE)

dds.season.i=DESeq(dds.season,test="LRT",reduced=~site+wetdry, parallel=TRUE)

dds2=DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ site+wetdry)
dds2$wetdry <- factor(dds2$wetdry, levels = c("dry","wet"))
dds.season.season=DESeq(dds2, parallel=TRUE)
dds2$site <- factor(dds2$site, levels = c("south","ledge","central"))
dds.season.site=DESeq(dds2,test="LRT",reduced=~wetdry, parallel=TRUE)

# saving all models
save(dds,dds.time,dds.site,dds.i,dds.season,dds.season.site,dds.season.season,dds.season.i,file="realModels.RData")

#--------------
# GLM analysis

load("realModels.RData")
library(DESeq2)

# time treatment
time=results(dds.time) 
summary(time) 
degs.time=row.names(time)[time$padj<0.1 & !(is.na(time$padj))]

# site treatment
site=results(dds.site) 
summary(site) 
degs.site=row.names(site)[site$padj<0.1 & !(is.na(site$padj))]

# time:site interaction
int=results(dds.i)
summary(int)
degs.int=row.names(int)[int$padj<0.1 & !(is.na(int$padj))]

# site contrasts
ledge.south=results(dds.site,contrast=c("site","ledge","south"))
summary(ledge.south)
degs.ledge.south=row.names(ledge.south)[ledge.south$padj<0.1 & !(is.na(ledge.south$padj))]

central.south=results(dds.site,contrast=c("site","central","south"))
summary(central.south)
degs.central.south=row.names(central.south)[central.south$padj<0.1 & !(is.na(central.south$padj))]

central.ledge=results(dds.site,contrast=c("site","central","ledge"))
summary(central.ledge)
degs.central.ledge=row.names(central.ledge)[central.ledge$padj<0.1 & !(is.na(central.ledge$padj))]

# time contrasts
t1.0=results(dds.time,contrast=c("time","1","0"))
summary(t1.0)
degs.t1.0=row.names(t1.0)[t1.0$padj<0.1 & !(is.na(t1.0$padj))]

t2.0=results(dds.time,contrast=c("time","2","0"))
summary(t2.0)
degs.t2.0=row.names(t2.0)[t2.0$padj<0.1 & !(is.na(t2.0$padj))]

t3.0=results(dds.time,contrast=c("time","3","0"))
summary(t3.0)
degs.t3.0=row.names(t3.0)[t3.0$padj<0.1 & !(is.na(t3.0$padj))]

t4.0=results(dds.time,contrast=c("time","4","0"))
summary(t4.0)
degs.t4.0=row.names(t4.0)[t4.0$padj<0.1 & !(is.na(t4.0$padj))]

t5.0=results(dds.time,contrast=c("time","5","0"))
summary(t5.0)
degs.t5.0=row.names(t5.0)[t5.0$padj<0.1 & !(is.na(t5.0$padj))]

t6.0=results(dds.time,contrast=c("time","6","0"))
summary(t6.0)
degs.t6.0=row.names(t6.0)[t6.0$padj<0.1 & !(is.na(t6.0$padj))]

t2.1=results(dds.time,contrast=c("time","2","1"))
summary(t2.1)
degs.t2.1=row.names(t2.1)[t2.1$padj<0.1 & !(is.na(t2.1$padj))]

t3.1=results(dds.time,contrast=c("time","3","1"))
summary(t3.1)
degs.t3.1=row.names(t3.1)[t3.1$padj<0.1 & !(is.na(t3.1$padj))]

t4.1=results(dds.time,contrast=c("time","4","1"))
summary(t4.1)
degs.t4.1=row.names(t4.1)[t4.1$padj<0.1 & !(is.na(t4.1$padj))]

t5.1=results(dds.time,contrast=c("time","5","1"))
summary(t5.1)
degs.t5.1=row.names(t5.1)[t5.1$padj<0.1 & !(is.na(t5.1$padj))]

t6.1=results(dds.time,contrast=c("time","6","1"))
summary(t6.1)
degs.t6.1=row.names(t6.1)[t6.1$padj<0.1 & !(is.na(t6.1$padj))]

t3.2=results(dds.time,contrast=c("time","3","2"))
summary(t3.2)
degs.t3.2=row.names(t3.2)[t3.2$padj<0.1 & !(is.na(t3.2$padj))]

t4.2=results(dds.time,contrast=c("time","4","2"))
summary(t4.2)
degs.t4.2=row.names(t4.2)[t4.2$padj<0.1 & !(is.na(t4.2$padj))]

t5.2=results(dds.time,contrast=c("time","5","2"))
summary(t5.2)
degs.t5.2=row.names(t5.2)[t5.2$padj<0.1 & !(is.na(t5.2$padj))]

t6.2=results(dds.time,contrast=c("time","6","2"))
summary(t6.2)
degs.t6.2=row.names(t6.2)[t6.2$padj<0.1 & !(is.na(t6.2$padj))]

t4.3=results(dds.time,contrast=c("time","4","3"))
summary(t4.3)
degs.t4.3=row.names(t4.3)[t4.3$padj<0.1 & !(is.na(t4.3$padj))]

t5.3=results(dds.time,contrast=c("time","5","3"))
summary(t5.3)
degs.t5.3=row.names(t5.3)[t5.3$padj<0.1 & !(is.na(t5.3$padj))]

t6.3=results(dds.time,contrast=c("time","6","3"))
summary(t6.3)
degs.t6.3=row.names(t6.3)[t6.3$padj<0.1 & !(is.na(t6.3$padj))]

t5.4=results(dds.time,contrast=c("time","5","4"))
summary(t5.4)
degs.t5.4=row.names(t5.4)[t5.4$padj<0.1 & !(is.na(t5.4$padj))]

t6.4=results(dds.time,contrast=c("time","6","4"))
summary(t6.4)
degs.t6.4=row.names(t6.4)[t6.4$padj<0.1 & !(is.na(t6.4$padj))]

t6.5=results(dds.time,contrast=c("time","6","5"))
summary(t6.5)
degs.t6.5=row.names(t6.5)[t6.5$padj<0.1 & !(is.na(t6.5$padj))]

# site+season
# site treatment
season.site=results(dds.season.site) 
summary(season.site) 
degs.season.site=row.names(season.site)[season.site$padj<0.1 & !(is.na(season.site$padj))]

# season treatment
season.season=results(dds.season.season,contrast=c("wetdry","dry","wet")) 
summary(season.season) 
degs.season.season=row.names(season.season)[season.season$padj<0.1 & !(is.na(season.season$padj))]

# site:season interaction
season.int=results(dds.season.i)
summary(season.int)
degs.season.int=row.names(season.int)[season.int$padj<0.1 & !(is.na(season.int$padj))]

# site contrasts
season.ledge.south=results(dds.season.site,contrast=c("site","ledge","south"))
summary(season.ledge.south)
degs.season.ledge.south=row.names(season.ledge.south)[season.ledge.south$padj<0.1 & !(is.na(season.ledge.south$padj))]

season.central.south=results(dds.season.site,contrast=c("site","central","south"))
summary(season.central.south)
degs.season.central.south=row.names(season.central.south)[season.central.south$padj<0.1 & !(is.na(season.central.south$padj))]

season.central.ledge=results(dds.season.site,contrast=c("site","central","ledge"))
summary(season.central.ledge)
degs.season.central.ledge=row.names(season.central.ledge)[season.central.ledge$padj<0.1 & !(is.na(season.central.ledge$padj))]

save(time,site,int,ledge.south,central.south,central.ledge,t1.0,t2.0,t3.0,t4.0,t5.0,t6.0,t2.1,t3.1,t4.1,t5.1,t6.1,t3.2,t4.2,t5.2,t6.2,t4.3,t5.3,t6.3,t5.4,t6.4,t6.5,season.site,season.season,season.int,season.ledge.south,season.central.south,season.central.ledge,degs.time,degs.site,degs.int,degs.ledge.south,degs.central.south,degs.central.ledge,degs.t1.0,degs.t2.0,degs.t3.0,degs.t4.0,degs.t5.0,degs.t6.0,degs.t2.1,degs.t3.1,degs.t4.1,degs.t5.1,degs.t6.1,degs.t3.2,degs.t4.2,degs.t5.2,degs.t6.2,degs.t4.3,degs.t5.3,degs.t6.3,degs.t5.4,degs.t6.4,degs.t6.5,degs.season.site,degs.season.season,degs.season.int,degs.season.ledge.south,degs.season.central.south,degs.season.central.ledge,file="pvals.RData")

#-------------------
# density plots: are my DEGs high-abundant or low-abundant?

load("vsd.RData")
load("pvals.RData")

means=apply(vsd,1,mean)

pdf(file="DEG_density_time.pdf", height=5, width=5)
plot(density(means))
lines(density(means[degs.time]),col="blue")
lines(density(means[degs.site]),col="orange")
lines(density(means[degs.int]),col="lightblue")
legend("topright", title = "Factor", legend=c("time","site","interaction"), fill = c("blue","orange","lightblue"))
dev.off()

means.season=apply(vsd.season,1,mean)

pdf(file="DEG_density_season.pdf", height=5, width=5)
plot(density(means.season))
lines(density(means.season[degs.season.site]),col="blue")
lines(density(means.season[degs.season.season]),col="orange")
lines(density(means.season[degs.season.int]),col="lightblue")
legend("topright", title = "Factor", legend=c("site","season","interaction"), fill = c("blue","orange","lightblue"))
dev.off()

#-------------------
# venn diagrams

load("pvals.RData")
library(DESeq2)

candidates=list("time"=degs.time, "site"=degs.site, "interaction"=degs.int)

# install.packages("VennDiagram")
library(VennDiagram)

# overall factors, full model
fullmodel.venn=venn.diagram(
	x = candidates,
	filename=NULL,
	col = "transparent",
	fill = c("blue", "orange", "lightblue"),
	alpha = 0.5,
	label.col = c("darkblue", "white", "darkred", "white", "white", "white", "cornflowerblue"),
	cex = 5,
	fontfamily = "sans",
	fontface = "bold",
	cat.default.pos = "text",
	cat.col =c("darkblue", "darkred", "cornflowerblue"),
	cat.cex = 5,
	cat.fontfamily = "sans",
	cat.dist = c(0.06, 0.06, -0.06),
	cat.pos = 3
	)
pdf(file="Venn_insitu_mcav.pdf", height=12, width=12)
grid.draw(fullmodel.venn)
dev.off()

pairwise=list("ledge.south"=degs.ledge.south,"central.south"=degs.central.south, "central.ledge"=degs.central.ledge)

# overall factors, full model
pairwise.venn=venn.diagram(
  x = pairwise,
  filename=NULL,
  col = "transparent",
  fill = c("blue", "orange", "lightblue"),
  alpha = 0.5,
  label.col = c("darkblue", "white", "darkred", "white", "white", "white", "cornflowerblue"),
  cex = 5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("darkblue", "darkred", "cornflowerblue"),
  cat.cex = 5,
  cat.fontfamily = "sans",
  cat.dist = c(0.06, 0.06, -0.06),
  cat.pos = 3
)
pdf(file="Venn_insitu_mcav_pairwise.pdf", height=12, width=12)
grid.draw(pairwise.venn)
dev.off()

#-------------------
# saving data for GO and KOG analysis 
load("realModels.RData")
load("pvals.RData")

# time response
# signed log p-values: -log(pvalue)* direction:
source=time[!is.na(time$pvalue),]
time.p=data.frame("gene"=row.names(source))
time.p$lpv=-log(source[,"pvalue"],10)
time.p$lpv[source$stat<0]=time.p$lpv[source$stat<0]*-1
head(time.p)
write.csv(time.p,file="time_lpv.csv",row.names=F,quote=F)
save(time.p,file="time_lpv.RData")

# site response
# signed log p-values: -log(pvalue)* direction:
source=site[!is.na(site$pvalue),]
site.p=data.frame("gene"=row.names(source))
site.p$lpv=-log(source[,"pvalue"],10)
site.p$lpv[source$stat<0]=site.p$lpv[source$stat<0]*-1
head(site.p)
write.csv(site.p,file="site_lpv.csv",row.names=F,quote=F)
save(site.p,file="site_lpv.RData")

# time:site interaction 
# signed log p-values: -log(pvalue)* direction:
source=int[!is.na(int$pvalue),]
int.p=data.frame("gene"=row.names(source))
int.p$lpv=-log(source[,"pvalue"],10)
int.p$lpv[source$stat<0]=int.p$lpv[source$stat<0]*-1
head(int.p)
write.csv(int.p,file="int_lpv.csv",row.names=F,quote=F)
save(int.p,file="int_lpv.RData")

#--------------------------
# ledge.south
# log2 fold changes:
source=ledge.south[!is.na(ledge.south$pvalue),]
ledge.south.fc=data.frame("gene"=row.names(source))
ledge.south.fc$lfc=source[,"log2FoldChange"]
head(ledge.south.fc)
write.csv(ledge.south.fc,file="ledge.south_fc.csv",row.names=F,quote=F)
save(ledge.south.fc,file="ledge.south_fc.RData")

# signed log p-values: -log(pvalue)* direction:
ledge.south.p=data.frame("gene"=row.names(source))
ledge.south.p$lpv=-log(source[,"pvalue"],10)
ledge.south.p$lpv[source$stat<0]=ledge.south.p$lpv[source$stat<0]*-1
head(ledge.south.p)
write.csv(ledge.south.p,file="ledge.south_lpv.csv",row.names=F,quote=F)
save(ledge.south.p,file="ledge.south_lpv.RData")

# central.south
# log2 fold changes:
source=central.south[!is.na(central.south$pvalue),]
central.south.fc=data.frame("gene"=row.names(source))
central.south.fc$lfc=source[,"log2FoldChange"]
head(central.south.fc)
write.csv(central.south.fc,file="central.south_fc.csv",row.names=F,quote=F)
save(central.south.fc,file="central.south_fc.RData")

# signed log p-values: -log(pvalue)* direction:
central.south.p=data.frame("gene"=row.names(source))
central.south.p$lpv=-log(source[,"pvalue"],10)
central.south.p$lpv[source$stat<0]=central.south.p$lpv[source$stat<0]*-1
head(central.south.p)
write.csv(central.south.p,file="central.south_lpv.csv",row.names=F,quote=F)
save(central.south.p,file="central.south_lpv.RData")

# central.ledge
# log2 fold changes:
source=central.ledge[!is.na(central.ledge$pvalue),]
central.ledge.fc=data.frame("gene"=row.names(source))
central.ledge.fc$lfc=source[,"log2FoldChange"]
head(central.ledge.fc)
write.csv(central.ledge.fc,file="central.ledge_fc.csv",row.names=F,quote=F)
save(central.ledge.fc,file="central.ledge_fc.RData")

# signed log p-values: -log(pvalue)* direction:
central.ledge.p=data.frame("gene"=row.names(source))
central.ledge.p$lpv=-log(source[,"pvalue"],10)
central.ledge.p$lpv[source$stat<0]=central.ledge.p$lpv[source$stat<0]*-1
head(central.ledge.p)
write.csv(central.ledge.p,file="central.ledge_lpv.csv",row.names=F,quote=F)
save(central.ledge.p,file="central.ledge_lpv.RData")

#--------------------------
# t1.0
# log2 fold changes:
source=t1.0[!is.na(t1.0$pvalue),]
t1.0.fc=data.frame("gene"=row.names(source))
t1.0.fc$lfc=source[,"log2FoldChange"]
head(t1.0.fc)
write.csv(t1.0.fc,file="t1.0_fc.csv",row.names=F,quote=F)
save(t1.0.fc,file="t1.0_fc.RData")

# t2.0
# log2 fold changes:
source=t2.0[!is.na(t2.0$pvalue),]
t2.0.fc=data.frame("gene"=row.names(source))
t2.0.fc$lfc=source[,"log2FoldChange"]
head(t2.0.fc)
write.csv(t2.0.fc,file="t2.0_fc.csv",row.names=F,quote=F)
save(t2.0.fc,file="t2.0_fc.RData")

# t3.0
# log2 fold changes:
source=t3.0[!is.na(t3.0$pvalue),]
t3.0.fc=data.frame("gene"=row.names(source))
t3.0.fc$lfc=source[,"log2FoldChange"]
head(t3.0.fc)
write.csv(t3.0.fc,file="t3.0_fc.csv",row.names=F,quote=F)
save(t3.0.fc,file="t3.0_fc.RData")

# t4.0
# log2 fold changes:
source=t4.0[!is.na(t4.0$pvalue),]
t4.0.fc=data.frame("gene"=row.names(source))
t4.0.fc$lfc=source[,"log2FoldChange"]
head(t4.0.fc)
write.csv(t4.0.fc,file="t4.0_fc.csv",row.names=F,quote=F)
save(t4.0.fc,file="t4.0_fc.RData")

# t5.0
# log2 fold changes:
source=t5.0[!is.na(t5.0$pvalue),]
t5.0.fc=data.frame("gene"=row.names(source))
t5.0.fc$lfc=source[,"log2FoldChange"]
head(t5.0.fc)
write.csv(t5.0.fc,file="t5.0_fc.csv",row.names=F,quote=F)
save(t5.0.fc,file="t5.0_fc.RData")

# t6.0
# log2 fold changes:
source=t6.0[!is.na(t6.0$pvalue),]
t6.0.fc=data.frame("gene"=row.names(source))
t6.0.fc$lfc=source[,"log2FoldChange"]
head(t6.0.fc)
write.csv(t6.0.fc,file="t6.0_fc.csv",row.names=F,quote=F)
save(t6.0.fc,file="t6.0_fc.RData")

# t2.1
# log2 fold changes:
source=t2.1[!is.na(t2.1$pvalue),]
t2.1.fc=data.frame("gene"=row.names(source))
t2.1.fc$lfc=source[,"log2FoldChange"]
head(t2.1.fc)
write.csv(t2.1.fc,file="t2.1_fc.csv",row.names=F,quote=F)
save(t2.1.fc,file="t2.1_fc.RData")

# t3.1
# log2 fold changes:
source=t3.1[!is.na(t3.1$pvalue),]
t3.1.fc=data.frame("gene"=row.names(source))
t3.1.fc$lfc=source[,"log2FoldChange"]
head(t3.1.fc)
write.csv(t3.1.fc,file="t3.1_fc.csv",row.names=F,quote=F)
save(t3.1.fc,file="t3.1_fc.RData")

# t4.1
# log2 fold changes:
source=t4.1[!is.na(t4.1$pvalue),]
t4.1.fc=data.frame("gene"=row.names(source))
t4.1.fc$lfc=source[,"log2FoldChange"]
head(t4.1.fc)
write.csv(t4.1.fc,file="t4.1_fc.csv",row.names=F,quote=F)
save(t4.1.fc,file="t4.1_fc.RData")

# t5.1
# log2 fold changes:
source=t5.1[!is.na(t5.1$pvalue),]
t5.1.fc=data.frame("gene"=row.names(source))
t5.1.fc$lfc=source[,"log2FoldChange"]
head(t5.1.fc)
write.csv(t5.1.fc,file="t5.1_fc.csv",row.names=F,quote=F)
save(t5.1.fc,file="t5.1_fc.RData")

# t6.1
# log2 fold changes:
source=t6.1[!is.na(t6.1$pvalue),]
t6.1.fc=data.frame("gene"=row.names(source))
t6.1.fc$lfc=source[,"log2FoldChange"]
head(t6.1.fc)
write.csv(t6.1.fc,file="t6.1_fc.csv",row.names=F,quote=F)
save(t6.1.fc,file="t6.1_fc.RData")

# t3.2
# log2 fold changes:
source=t3.2[!is.na(t3.2$pvalue),]
t3.2.fc=data.frame("gene"=row.names(source))
t3.2.fc$lfc=source[,"log2FoldChange"]
head(t3.2.fc)
write.csv(t3.2.fc,file="t3.2_fc.csv",row.names=F,quote=F)
save(t3.2.fc,file="t3.2_fc.RData")

# t4.2
# log2 fold changes:
source=t4.2[!is.na(t4.2$pvalue),]
t4.2.fc=data.frame("gene"=row.names(source))
t4.2.fc$lfc=source[,"log2FoldChange"]
head(t4.2.fc)
write.csv(t4.2.fc,file="t4.2_fc.csv",row.names=F,quote=F)
save(t4.2.fc,file="t4.2_fc.RData")

# t5.2
# log2 fold changes:
source=t5.2[!is.na(t5.2$pvalue),]
t5.2.fc=data.frame("gene"=row.names(source))
t5.2.fc$lfc=source[,"log2FoldChange"]
head(t5.2.fc)
write.csv(t5.2.fc,file="t5.2_fc.csv",row.names=F,quote=F)
save(t5.2.fc,file="t5.2_fc.RData")

# t6.2
# log2 fold changes:
source=t6.2[!is.na(t6.2$pvalue),]
t6.2.fc=data.frame("gene"=row.names(source))
t6.2.fc$lfc=source[,"log2FoldChange"]
head(t6.2.fc)
write.csv(t6.2.fc,file="t6.2_fc.csv",row.names=F,quote=F)
save(t6.2.fc,file="t6.2_fc.RData")

# t4.3
# log2 fold changes:
source=t4.3[!is.na(t4.3$pvalue),]
t4.3.fc=data.frame("gene"=row.names(source))
t4.3.fc$lfc=source[,"log2FoldChange"]
head(t4.3.fc)
write.csv(t4.3.fc,file="t4.3_fc.csv",row.names=F,quote=F)
save(t4.3.fc,file="t4.3_fc.RData")

# t5.3
# log2 fold changes:
source=t5.3[!is.na(t5.3$pvalue),]
t5.3.fc=data.frame("gene"=row.names(source))
t5.3.fc$lfc=source[,"log2FoldChange"]
head(t5.3.fc)
write.csv(t5.3.fc,file="t5.3_fc.csv",row.names=F,quote=F)
save(t5.3.fc,file="t5.3_fc.RData")

# t6.3
# log2 fold changes:
source=t6.3[!is.na(t6.3$pvalue),]
t6.3.fc=data.frame("gene"=row.names(source))
t6.3.fc$lfc=source[,"log2FoldChange"]
head(t6.3.fc)
write.csv(t6.3.fc,file="t6.3_fc.csv",row.names=F,quote=F)
save(t6.3.fc,file="t6.3_fc.RData")

# t5.4
# log2 fold changes:
source=t5.4[!is.na(t5.4$pvalue),]
t5.4.fc=data.frame("gene"=row.names(source))
t5.4.fc$lfc=source[,"log2FoldChange"]
head(t5.4.fc)
write.csv(t5.4.fc,file="t5.4_fc.csv",row.names=F,quote=F)
save(t5.4.fc,file="t5.4_fc.RData")

# t6.4
# log2 fold changes:
source=t6.4[!is.na(t6.4$pvalue),]
t6.4.fc=data.frame("gene"=row.names(source))
t6.4.fc$lfc=source[,"log2FoldChange"]
head(t6.4.fc)
write.csv(t6.4.fc,file="t6.4_fc.csv",row.names=F,quote=F)
save(t6.4.fc,file="t6.4_fc.RData")

# t6.5
# log2 fold changes:
source=t6.5[!is.na(t6.5$pvalue),]
t6.5.fc=data.frame("gene"=row.names(source))
t6.5.fc$lfc=source[,"log2FoldChange"]
head(t6.5.fc)
write.csv(t6.5.fc,file="t6.5_fc.csv",row.names=F,quote=F)
save(t6.5.fc,file="t6.5_fc.RData")
