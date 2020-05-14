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
counts = read.table("exsitu_allcounts_sym.txt")

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

# for WCGNA: removing all genes with counts of <10 in more than 90 % of samples
counts4wgcna = counts[apply(counts,1,function(x) sum(x<10))<ncol(counts)*0.9,]
nrow(counts4wgcna)
ncol(counts4wgcna)
write.csv(counts4wgcna, file="counts4wgcna.csv")

# importing a design .csv file
design = read.csv("exsitu_design.csv", head=TRUE)
design
design$tank = as.factor(design$tank)
design$block = as.factor(design$block)
design$tag = as.factor(design$tag)
str(design)

# making dataframes
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ genotype*temp*water+block:genotype+block:temp+block:water)

dds.comb = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ temp.water)

# reorders treatment factor according to "control" vs "treatment" levels
dds$temp <- factor(dds$temp, levels = c("control", "elevated"))
dds$water <- factor(dds$water, levels = c("offshore","discharge"))

dds.comb$temp.water <- factor(dds.comb$temp.water, levels = c("control.offshore","control.discharge","elevated.offshore","elevated.discharge"))

# for large datasets, rlog may take too much time, especially for an unfiltered dataframe
# vsd is much faster and still works for outlier detection
Vsd=varianceStabilizingTransformation(dds)
Vsd.comb=varianceStabilizingTransformation(dds.comb)

library(Biobase)
e=ExpressionSet(assay(Vsd), AnnotatedDataFrame(as.data.frame(colData(Vsd))))

# running outlier detection
arrayQualityMetrics(e,intgroup=c("temp.water"),force=T)
# dev.off()
# open the directory "arrayQualityMetrics report for e" in your working directory and open index.html
# Array metadata and outlier detection overview gives a report of all samples, and which are likely outliers according to the 3 methods tested. I typically remove the samples that violate *1 (distance between arrays).
# Figure 2 shows a bar plot of array-to-array distances and an outlier detection threshold based on your samples. Samples above the threshold are considered outliers
# under Figure 3: Principal Components Analyses, look for any points far away from the rest of the sample cluster
# use the array number for removal in the following section

# if there were outliers, say, array 150:
outs=c(71,81,110,123,130)
countData=countData[,-outs]
Vsd=Vsd[,-outs]
Vsd.comb=Vsd.comb[,-outs]
counts4wgcna=counts4wgcna[,-outs]
design=design[-outs,]

# remaking model with outliers removed from dataset
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ genotype*temp*water+block:genotype+block:temp+block:water)

dds.comb = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ temp.water)

dds$temp <- factor(dds$temp, levels = c("control", "elevated"))
dds$water <- factor(dds$water, levels = c("offshore","discharge"))

dds.comb$temp.water <- factor(dds.comb$temp.water, levels = c("control.offshore","control.discharge","elevated.offshore","elevated.discharge"))

# save all these dataframes as an Rdata package so you don't need to rerun each time
save(dds,dds.comb,design,countData,Vsd,Vsd.comb,counts4wgcna,file="initial.RData")

#---------------------
# generating normalized variance-stabilized data for PCoA, heatmaps, etc

load("initial.RData")
library(DESeq2)
library(BiocParallel)

# creating normalized dataframes
vsd=assay(Vsd)
# takes the sample IDs and factor levels from the design to create new column names for the dataframe
snames=paste(colnames(countData),design[,5],design[,6],design[,7],sep=".")
# renames the column names
colnames(vsd)=snames

vsd.comb=assay(Vsd.comb)
colnames(vsd.comb)=snames

save(vsd,vsd.comb,design,file="vsd.RData")

# more reduced stabilized dataset for WGCNA
wg = DESeqDataSetFromMatrix(countData=counts4wgcna, colData=design, design=~ genotype*temp*water+block:genotype+block:temp+block:water)
vsd.wg=assay(varianceStabilizingTransformation(wg), blind=TRUE)
head(vsd.wg)
colnames(vsd.wg)=snames
save(vsd.wg,design,file="data4wgcna.RData")

#-------------------
# EXPLORING SIMILARITIES AMONG SAMPLES

# heatmap and hierarchical clustering:
load("vsd.RData")
library(pheatmap)
# similarity among samples
pdf(file="heatmap_exsitu_sym.pdf", width=25, height=25)
pheatmap(cor(vsd))
dev.off()

# Principal coordinates analysis
library(vegan)
library(rgl)
library(ape)

conditions=design
conditions$temp.water <- factor(conditions$temp.water, levels = c("control.offshore", "control.discharge","elevated.offshore","elevated.discharge"))
conditions$temp <- factor(conditions$temp, levels = c("control", "elevated"))
conditions$water <- factor(conditions$water, levels = c("offshore","discharge"))

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

# plotting PCoA by genotype and combined treatments
pdf(file="PCoA_exsitu_genotype_sym.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=c("#40004b","#762a83", "#9970ab","#c2a5cf","#e7d4e8","#f7f7f7", "#d9f0d3","#a6dba0","#a6dba0","#1b7837", "#00441b")[as.numeric((as.factor(conditions$colony)))],pch=c(0,15,2,17)[as.numeric(as.factor(conditions$temp.water))], xlab="Coordinate 1", ylab="Coordinate 2", main="Colony")
# cluster overlay of colony
ordiellipse(scores, conditions$colony, label=F, draw= "polygon", col=c("#40004b","#762a83", "#9970ab","#c2a5cf","#e7d4e8","#f7f7f7", "#d9f0d3","#a6dba0","#a6dba0","#1b7837", "#00441b"))
legend("bottomleft", legend=c("a","b","c","d","e","f","g","h","i","j","k"), fill = c("#40004b","#762a83", "#9970ab","#c2a5cf","#e7d4e8","#f7f7f7", "#d9f0d3","#a6dba0","#a6dba0","#1b7837", "#00441b"), bty="n")
legend("topleft", legend=c("control.offshore", "control.discharge","elevated.offshore","elevated.discharge"), pch=c(0,15,2,17), bty="n")
# cluster overlay of temp.water
plot(scores[,1], scores[,2],col=c("#92c5de","#0571b0","#f4a582","#ca0020")[as.numeric(as.factor(conditions$temp.water))],pch=c(1,3,4,5,6,7,8,9,10,16,18)[as.numeric(as.factor(conditions$colony))], xlab="Coordinate 1", ylab="Coordinate 2", main="Temp + Water")
ordiellipse(scores, conditions$temp.water, label=F, draw= "polygon", col=c("#92c5de","#0571b0","#f4a582","#ca0020"))
legend("topleft", legend=c("control.offshore", "control.discharge","elevated.offshore","elevated.discharge"), fill = c("#92c5de","#0571b0","#f4a582","#ca0020"), bty="n")
legend("bottomleft", legend=c("a","b","c","d","e","f","g","h","i","j","k"), pch=c(1,3,4,5,6,7,8,9,10,16,18), bty="n")
dev.off()

# plotting PCoA by treatments separately
pdf(file="PCoA_exsitu_sym.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=c("#0571b0","#ca0020")[as.numeric(as.factor(conditions$temp))],pch=c(0,15)[as.numeric(as.factor(conditions$water))], xlab="Coordinate 1", ylab="Coordinate 2", main="Temperature")
# cluster overlay of temp
ordiellipse(scores, conditions$temp, label=F, draw= "polygon", col=c("#0571b0","#ca0020"))
legend("bottomleft", legend=c("control","elevated"), fill = c("#0571b0","#ca0020"), bty="n")
legend("topleft", legend=c("offshore","discharge"), pch=c(0,15), bty="n")
# cluster overlay of water
plot(scores[,1], scores[,2],col=c("#018571","#a6611a")[as.numeric(as.factor(conditions$water))],pch=c(25,17)[as.numeric(as.factor(conditions$temp))], xlab="Coordinate 1", ylab="Coordinate 2", main="Water")
ordiellipse(scores, conditions$water, label=F, draw= "polygon", col=c("#018571","#a6611a"))
legend("topleft", legend=c("offshore","discharge"), fill = c("#018571","#a6611a"), bty="n")
legend("bottomleft", legend=c("control","elevated"), pch=c(25,17), bty="n")
dev.off()

# neighbor-joining tree of samples (based on significant PCo's):
pdf(file="PCoA_tree.pdf", width=15, height=40)
tre=nj(dist(scores[,1:4]))
plot(tre,cex=0.8)
dev.off()

# formal analysis of variance in distance matricies
ad=adonis(t(vsd)~genotype*temp*water+block:genotype+block:temp+block:water,data=conditions,method="manhattan",permutations=9999)
ad

# creating pie chart to represent ANOVA results
cols=c("#084081","#0868ac","#2b8cbe","#4eb3d3","#7bccc4","#a8ddb5","#ccebc5","#e0f3db","#f7fcf0","grey80")
pdf(file="ANOVA_pie.pdf", width=6, height=6)
pie(ad$aov.tab$R2[1:10],labels=row.names(ad$aov.tab)[1:10],col=cols,main="genotype*temp*water+block")
dev.off()

# DAPC analysis: creating discriminant function to tell treatments apart, based on controls
library(adegenet)

# runs simulations on randomly-chosen datasets of 90% of the total dataset to test the number of PCs to retain
set.seed(999)
# by temp.water, excluding elevated.discharge samples
xvalDapc(t(vsd.comb[,conditions$temp.water!="elevated.discharge"]),conditions$temp.water[conditions$temp.water!="elevated.discharge"], n.rep=100, parallel="multicore", ncpus= 8)
# 30 PCs
xvalDapc(t(vsd.comb[,conditions$temp.water!="elevated.discharge"]),conditions$temp.water[conditions$temp.water!="elevated.discharge"], n.rep=1000, n.pca=20:40, parallel="multicore", ncpus= 8)
# 32 PCs

# now running the dapc without elevated.discharge samples
dp.tw=dapc(t(vsd.comb[,conditions$temp.water!="elevated.discharge"]),conditions$temp.water[conditions$temp.water!="elevated.discharge"],n.pca=32, n.da=2)
dp.tw$ind.coord

# can we predict water for the elevated samples based on patterns among the control samples?
pred.ed=predict.dapc(dp.tw,newdata=(t(vsd.comb[,conditions$temp.water=="elevated.discharge"])))
pred.ed
# look at the posterior section for assignments by probability

# creating a new dapc object to add in mesophotic corals for plotting
dp.ed=dp.tw
dp.ed$ind.coord=pred.ed$ind.scores
dp.ed$posterior=pred.ed$posterior
dp.ed$assign=pred.ed$assign
dp.ed$grp<-as.factor(c("elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge","elevated.discharge"))

# now exporting side by side figures of control vs single-factor treatments and double-factor treatment
pdf(file="DAPC_exsitu_sym.pdf", width=12, height=6)
par(mfrow=c(1,2))
scatter(dp.tw, bg="white",scree.da=FALSE,legend=TRUE,solid=0.6, col= c("#92c5de","#0571b0","#f4a582","#ca0020"))
scatter(dp.ed, bg="white",scree.da=FALSE,legend=FALSE,solid=0.6, col= c("#ca0020"))
dev.off()

#----------------------------
# exporting for significance testing below
dpc.tw=data.frame(rbind(dp.tw$ind.coord,pred.ed$ind.scores))

write.csv(dpc.tw, "DAPC_exsitu_sym_DFA.csv", quote=F)
# modify the output CSV to add in a column for temp.water

# then reimport
dpc.tw<- read.csv("DAPC_exsitu_sym_DFA.csv")

# a little bit of rearranging
dpc.tw$temp.water<-factor(dpc.tw$temp.water, levels=c("control.offshore","elevated.offshore","control.discharge","elevated.discharge"))
dpc.tw[order(dpc.tw$temp.water),]
dpc.tw$id<-as.factor(dpc.tw$id)
str(dpc.tw)

# testing significance of DFA differences with MCMCglmm
# install.packages("MCMCglmm")
library(MCMCglmm)

# sets prior distribution and creates a glm
prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))
glm.tw <-MCMCglmm(LD1~temp.water,random=~id, family="gaussian", data=dpc.tw,prior=prior,nitt=75000,thin=25,burnin=5000)
summary(glm.tw)
# check to make sure you don't have autocorrelation with the reps (shown as "walks" in model traces)
plot(glm.tw)

# calculating difference in magnitudes of elevated.offshore and control.discharge using sampled sets of parameters
eo.cd.Delta=abs(glm.tw$Sol[,"temp.waterelevated.offshore"])-abs(glm.tw$Sol[,"temp.watercontrol.discharge"])
# 95% credible interval
HPDinterval(eo.cd.Delta)

# MCMC p-value 
if (is.na(table(eo.cd.Delta<0)[2])) {
  cat("p <",signif(1/length(eo.cd.Delta),1))
} else { cat("p =",signif(table(eo.cd.Delta<0)[2]/length(eo.cd.Delta),2)) }

#----------------------
# FITTING GENE BY GENE MODELS

# with multi-factor, multi-level design - using LRT
load("initial.RData")
library(DESeq2)
library(BiocParallel)

# Running full model for contrast statements
dds=DESeq(dds, parallel=TRUE)

# Using LRT for the effect of interactions (using full model dds = same as design argument on line 63)
# then genotype+temp+water are removed, to look just at interaction term
dds.ct=DESeq(dds,test="LRT",reduced=~genotype+temp+water+genotype:water+temp:water+block:genotype+block:temp+block:water, parallel=TRUE)
dds.cw=DESeq(dds,test="LRT",reduced=~genotype+temp+water+genotype:temp+temp:water+block:genotype+block:temp+block:water, parallel=TRUE)
dds.tw=DESeq(dds,test="LRT",reduced=~genotype+temp+water+genotype:temp+genotype:water+block:temp+block:water, parallel=TRUE)
dds.bc=DESeq(dds,test="LRT",reduced=~genotype+temp+water+genotype:temp+genotype:water+temp:water+block:temp+block:water, parallel=TRUE)
dds.bt=DESeq(dds,test="LRT",reduced=~genotype+temp+water+genotype:temp+genotype:water+temp:water+block:genotype+block:water, parallel=TRUE)
dds.bw=DESeq(dds,test="LRT",reduced=~genotype+temp+water+genotype:temp+genotype:water+temp:water+block:genotype+block:temp, parallel=TRUE)

# Creating dataset with design =~genotype+temp+water (no interaction), to investigate importance of genotype, temp, and water
dds1=DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ genotype+temp+water)

# model for the effect of water (2 factor levels => Wald)
dds1$water <- factor(dds1$water, levels = c("offshore", "discharge"))
dds.water=DESeq(dds1, parallel=TRUE)

# model for the effect of temp:
dds1$temp <- factor(dds1$temp, levels = c("control", "elevated"))
dds.temp=DESeq(dds1,test="LRT", reduced=~water, parallel=TRUE)

# model for the effect of genotype:
dds1$genotype <- factor(dds1$genotype, levels = c("a","b","c","d","e","f","g","h","i","j","k"))
dds.genotype=DESeq(dds1,test="LRT", reduced=~temp+water, parallel=TRUE)

dds.comb=DESeq(dds.comb, parallel=TRUE)

# saving all models
save(dds,dds.genotype,dds.temp,dds.water,dds.ct,dds.cw,dds.tw,dds.bc,dds.bt,dds.bw,dds.comb,file="realModels.RData")

#--------------
# GLM analysis

load("realModels.RData")
library(DESeq2)

# genotype treatment 
genotype=results(dds.genotype) 
summary(genotype) 
degs.genotype=row.names(genotype)[genotype$padj<0.1 & !(is.na(genotype$padj))]

# temp treatment - Wald test
temp=results(dds.temp,contrast=c("temp","elevated","control")) 
summary(temp) 
degs.temp=row.names(temp)[temp$padj<0.1 & !(is.na(temp$padj))]

# water treatment - Wald test
water=results(dds.water,contrast=c("water","discharge","offshore")) 
summary(water) 
degs.water=row.names(water)[water$padj<0.1 & !(is.na(water$padj))]

# genotype:temp interaction
ct=results(dds.ct)
summary(ct)
degs.ct=row.names(ct)[ct$padj<0.1 & !(is.na(ct$padj))]

# genotype:water interaction
cw=results(dds.cw)
summary(cw)
degs.cw=row.names(cw)[cw$padj<0.1 & !(is.na(cw$padj))]

# temp:water interaction
tw=results(dds.tw)
summary(tw)
degs.tw=row.names(tw)[tw$padj<0.1 & !(is.na(tw$padj))]

bc=results(dds.bc)
summary(bc)
degs.bc=row.names(bc)[bc$padj<0.1 & !(is.na(bc$padj))]

bt=results(dds.bt)
summary(bt)
degs.bt=row.names(bt)[bt$padj<0.1 & !(is.na(bt$padj))]

bw=results(dds.bw)
summary(bw)
degs.bw=row.names(bw)[bw$padj<0.1 & !(is.na(bw$padj))]

# For figuring out your contrast statements, the following vignette is helpful
# ?results

# combined model
comb=results(dds.comb)
summary(comb)

resultsNames(dds.comb)

cd.co=results(dds.comb,contrast=c("temp.water","control.discharge","control.offshore")) 
summary(cd.co) 
degs.cd.co=row.names(cd.co)[cd.co$padj<0.1 & !(is.na(cd.co$padj))]

eo.co=results(dds.comb,contrast=c("temp.water","elevated.offshore","control.offshore")) 
summary(eo.co) 
degs.eo.co=row.names(eo.co)[eo.co$padj<0.1 & !(is.na(eo.co$padj))]

ed.co=results(dds.comb,contrast=c("temp.water","elevated.discharge","control.offshore")) 
summary(ed.co) 
degs.ed.co=row.names(ed.co)[ed.co$padj<0.1 & !(is.na(ed.co$padj))]

eo.cd=results(dds.comb,contrast=c("temp.water","elevated.offshore","control.discharge")) 
summary(eo.cd) 
degs.eo.cd=row.names(eo.cd)[eo.cd$padj<0.1 & !(is.na(eo.cd$padj))]

ed.cd=results(dds.comb,contrast=c("temp.water","elevated.discharge","control.discharge")) 
summary(ed.cd) 
degs.ed.cd=row.names(ed.cd)[ed.cd$padj<0.1 & !(is.na(ed.cd$padj))]

ed.eo=results(dds.comb,contrast=c("temp.water","elevated.discharge","elevated.offshore")) 
summary(ed.eo) 
degs.ed.eo=row.names(ed.eo)[ed.eo$padj<0.1 & !(is.na(ed.eo$padj))]

save(genotype,temp,water,ct,cw,tw,bc,bt,bw,cd.co,eo.co,ed.co,eo.cd,ed.cd,ed.eo,degs.genotype,degs.temp,degs.water,degs.ct,degs.cw,degs.tw,degs.bc,degs.bt,degs.bw,degs.cd.co,degs.eo.co,degs.ed.co,degs.eo.cd,degs.ed.cd,degs.ed.eo,file="pvals.RData")

#-------------------
# density plots: are my DEGs high-abundant or low-abundant?

load("vsd.RData")
load("pvals.RData")

means=apply(vsd,1,mean)

pdf(file="DEG_density.pdf", height=5, width=5)
plot(density(means))
lines(density(means[degs.genotype]),col="blue")
lines(density(means[degs.temp]),col="orange")
lines(density(means[degs.water]),col="lightblue")
legend("topright", title = "Factor", legend=c("genotype","temp","water"), fill = c("blue","orange","lightblue"))
dev.off()

#-------------------
# venn diagrams

load("pvals.RData")
library(DESeq2)

candidates=list("genotype"=degs.genotype,"tw"=degs.tw, "temp"=degs.temp, "water"=degs.water)

# install.packages("VennDiagram")
library(VennDiagram)

# overall factors, full model
fullmodel.venn=venn.diagram(
	x = candidates,
	filename=NULL,
	col = "transparent",
	fill = c("blue","grey80", "orange", "lightblue"),
	alpha = 0.5,
	label.col = c("darkred","white","cornflowerblue","white","white","black","white", "white","darkblue","white","white","white","white","grey50","white"),
	cex = 3.5,
	fontfamily = "sans",
	fontface = "bold",
	cat.default.pos = "text",
	cat.col =c("darkblue", "grey50", "darkred","cornflowerblue"),
	cat.cex = 3.5,
	cat.fontfamily = "sans",
	cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
	)
pdf(file="Venn_exsitu_sym.pdf", height=10, width=12)
grid.draw(fullmodel.venn)
dev.off()

pairwise=list("cd.co"=degs.cd.co,"eo.co"=degs.eo.co, "ed.co"=degs.ed.co)

# overall factors, full model
pairwise.venn=venn.diagram(
  x = pairwise,
  filename=NULL,
  col = "transparent",
  fill = c("#0571b0","#f4a582","#ca0020"),
  alpha = 0.5,
  label.col = c("#0571b0","white","#f4a582","white","black","white","#ca0020"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("#0571b0","#f4a582","#ca0020"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0.25,0),c(0.75,0),c(0.5,1))
)
pdf(file="Venn_exsitu_sym_pairwise.pdf", height=6, width=8)
grid.draw(pairwise.venn)
dev.off()

#-------------------
# saving data for GO and KOG analysis 
load("realModels.RData")
load("pvals.RData")

# since DESeq2 was designed for 2 factor level designs (control vs treatment), only the treatment fold change data are relevant

# genotype response
# signed log p-values: -log(pvalue)* direction:
source=genotype[!is.na(genotype$pvalue),]
genotype.p=data.frame("gene"=row.names(source))
genotype.p$lpv=-log(source[,"pvalue"],10)
genotype.p$lpv[source$stat<0]=genotype.p$lpv[source$stat<0]*-1
head(genotype.p)
write.csv(genotype.p,file="genotype_lpv.csv",row.names=F,quote=F)
save(genotype.p,file="genotype_lpv.RData")

# temp response
# log2 fold changes:
head(temp)
source=temp[!is.na(temp$pvalue),]
temp.fc=data.frame("gene"=row.names(source))
temp.fc$lfc=source[,"log2FoldChange"]
head(temp.fc)
write.csv(temp.fc,file="temp_fc.csv",row.names=F,quote=F)
save(temp.fc,file="temp_fc.RData")

# signed log p-values: -log(pvalue)* direction:
temp.p=data.frame("gene"=row.names(source))
temp.p$lpv=-log(source[,"pvalue"],10)
temp.p$lpv[source$stat<0]=temp.p$lpv[source$stat<0]*-1
head(temp.p)
write.csv(temp.p,file="temp_lpv.csv",row.names=F,quote=F)
save(temp.p,file="temp_lpv.RData")

# water response
# log2 fold changes:
head(water)
source=water[!is.na(water$pvalue),]
water.fc=data.frame("gene"=row.names(source))
water.fc$lfc=source[,"log2FoldChange"]
head(water.fc)
write.csv(water.fc,file="water_fc.csv",row.names=F,quote=F)
save(water.fc,file="water_fc.RData")

# signed log p-values: -log(pvalue)* direction:
water.p=data.frame("gene"=row.names(source))
water.p$lpv=-log(source[,"pvalue"],10)
water.p$lpv[source$stat<0]=water.p$lpv[source$stat<0]*-1
head(water.p)
write.csv(water.p,file="water_lpv.csv",row.names=F,quote=F)
save(water.p,file="water_lpv.RData")

# temp:water interaction 
# signed log p-values: -log(pvalue)* direction:
source=tw[!is.na(tw$pvalue),]
tw.p=data.frame("gene"=row.names(source))
tw.p$lpv=-log(source[,"pvalue"],10)
tw.p$lpv[source$stat<0]=tw.p$lpv[source$stat<0]*-1
head(tw.p)
write.csv(tw.p,file="tw_lpv.csv",row.names=F,quote=F)
save(tw.p,file="tw_lpv.RData")

#--------------------------
# cd.co
# log2 fold changes:
source=cd.co[!is.na(cd.co$pvalue),]
cd.co.fc=data.frame("gene"=row.names(source))
cd.co.fc$lfc=source[,"log2FoldChange"]
head(cd.co.fc)
write.csv(cd.co.fc,file="cd.co_fc.csv",row.names=F,quote=F)
save(cd.co.fc,file="cd.co_fc.RData")

# signed log p-values: -log(pvalue)* direction:
cd.co.p=data.frame("gene"=row.names(source))
cd.co.p$lpv=-log(source[,"pvalue"],10)
cd.co.p$lpv[source$stat<0]=cd.co.p$lpv[source$stat<0]*-1
head(cd.co.p)
write.csv(cd.co.p,file="cd.co_lpv.csv",row.names=F,quote=F)
save(cd.co.p,file="cd.co_lpv.RData")

# eo.co
# log2 fold changes:
source=eo.co[!is.na(eo.co$pvalue),]
eo.co.fc=data.frame("gene"=row.names(source))
eo.co.fc$lfc=source[,"log2FoldChange"]
head(eo.co.fc)
write.csv(eo.co.fc,file="eo.co_fc.csv",row.names=F,quote=F)
save(eo.co.fc,file="eo.co_fc.RData")

# signed log p-values: -log(pvalue)* direction:
eo.co.p=data.frame("gene"=row.names(source))
eo.co.p$lpv=-log(source[,"pvalue"],10)
eo.co.p$lpv[source$stat<0]=eo.co.p$lpv[source$stat<0]*-1
head(eo.co.p)
write.csv(eo.co.p,file="eo.co_lpv.csv",row.names=F,quote=F)
save(eo.co.p,file="eo.co_lpv.RData")

# ed.co
# log2 fold changes:
source=ed.co[!is.na(ed.co$pvalue),]
ed.co.fc=data.frame("gene"=row.names(source))
ed.co.fc$lfc=source[,"log2FoldChange"]
head(ed.co.fc)
write.csv(ed.co.fc,file="ed.co_fc.csv",row.names=F,quote=F)
save(ed.co.fc,file="ed.co_fc.RData")

# signed log p-values: -log(pvalue)* direction:
ed.co.p=data.frame("gene"=row.names(source))
ed.co.p$lpv=-log(source[,"pvalue"],10)
ed.co.p$lpv[source$stat<0]=ed.co.p$lpv[source$stat<0]*-1
head(ed.co.p)
write.csv(ed.co.p,file="ed.co_lpv.csv",row.names=F,quote=F)
save(ed.co.p,file="ed.co_lpv.RData")

# eo.cd
# log2 fold changes:
source=eo.cd[!is.na(eo.cd$pvalue),]
eo.cd.fc=data.frame("gene"=row.names(source))
eo.cd.fc$lfc=source[,"log2FoldChange"]
head(eo.cd.fc)
write.csv(eo.cd.fc,file="eo.cd_fc.csv",row.names=F,quote=F)
save(eo.cd.fc,file="eo.cd_fc.RData")

# signed log p-values: -log(pvalue)* direction:
eo.cd.p=data.frame("gene"=row.names(source))
eo.cd.p$lpv=-log(source[,"pvalue"],10)
eo.cd.p$lpv[source$stat<0]=eo.cd.p$lpv[source$stat<0]*-1
head(eo.cd.p)
write.csv(eo.cd.p,file="eo.cd_lpv.csv",row.names=F,quote=F)
save(eo.cd.p,file="eo.cd_lpv.RData")

# ed.cd
# log2 fold changes:
source=ed.cd[!is.na(ed.cd$pvalue),]
ed.cd.fc=data.frame("gene"=row.names(source))
ed.cd.fc$lfc=source[,"log2FoldChange"]
head(ed.cd.fc)
write.csv(ed.cd.fc,file="ed.cd_fc.csv",row.names=F,quote=F)
save(ed.cd.fc,file="ed.cd_fc.RData")

# signed log p-values: -log(pvalue)* direction:
ed.cd.p=data.frame("gene"=row.names(source))
ed.cd.p$lpv=-log(source[,"pvalue"],10)
ed.cd.p$lpv[source$stat<0]=ed.cd.p$lpv[source$stat<0]*-1
head(ed.cd.p)
write.csv(ed.cd.p,file="ed.cd_lpv.csv",row.names=F,quote=F)
save(ed.cd.p,file="ed.cd_lpv.RData")

# ed.eo
# log2 fold changes:
source=ed.eo[!is.na(ed.eo$pvalue),]
ed.eo.fc=data.frame("gene"=row.names(source))
ed.eo.fc$lfc=source[,"log2FoldChange"]
head(ed.eo.fc)
write.csv(ed.eo.fc,file="ed.eo_fc.csv",row.names=F,quote=F)
save(ed.eo.fc,file="ed.eo_fc.RData")

# signed log p-values: -log(pvalue)* direction:
ed.eo.p=data.frame("gene"=row.names(source))
ed.eo.p$lpv=-log(source[,"pvalue"],10)
ed.eo.p$lpv[source$stat<0]=ed.eo.p$lpv[source$stat<0]*-1
head(ed.eo.p)
write.csv(ed.eo.p,file="ed.eo_lpv.csv",row.names=F,quote=F)
save(ed.eo.p,file="ed.eo_lpv.RData")