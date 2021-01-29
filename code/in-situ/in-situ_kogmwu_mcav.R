# install.packages("KOGMWU")
library(KOGMWU)

#---------------------------
# full model, assessing changes in KOG expression across colony, temp, and water

# loading KOG annotations
gene2kog=read.table("Mcavernosa_Cladocopium_iso2kogClass.tab",sep="\t")
head(gene2kog)

adt=load('time_lpv.RData')
adt # names of datasets in the package
lpv.t=kog.mwu(time.p,gene2kog) 
lpv.t 

ads=load('site_lpv.RData')
ads # names of datasets in the package
lpv.s=kog.mwu(site.p,gene2kog) 
lpv.s 

adi=load('int_lpv.RData')
adi # names of datasets in the package
lpv.i=kog.mwu(int.p,gene2kog) 
lpv.i

barplot(lpv.t$delta.rank,horiz=T)

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("Time"=lpv.t,"Site"=lpv.s,"Time:Site"=lpv.i))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)
  
# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_insitu_mcav_lpv.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
# creating a pub-ready corr plot
pdf(file="KOG_insitu_corr_mcav_lpv.pdf", width=7, height=5.2)
par(mfrow=c(1,3))
corrPlot(x="Time",y="Site",ktable)
corrPlot(x="Time",y="Time:Site",ktable)
corrPlot(x="Site",y="Time:Site",ktable)
dev.off()

#---------------------------
# pairwise comparisons of sites (lpv)

ledge.south=load('ledge.south_lpv.RData')
ledge.south # names of datasets in the package
lpv.ledge.south=kog.mwu(ledge.south.p,gene2kog) 
lpv.ledge.south

central.south=load('central.south_lpv.RData')
central.south # names of datasets in the package
lpv.central.south=kog.mwu(central.south.p,gene2kog) 
lpv.central.south

central.ledge=load('central.ledge_lpv.RData')
central.ledge # names of datasets in the package
lpv.central.ledge=kog.mwu(central.ledge.p,gene2kog) 
lpv.central.ledge

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("ledge.south"=lpv.ledge.south,"central.south"=lpv.central.south,"central.ledge"=lpv.central.ledge))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_insitu_site_mcav_lpv.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
# creating a pub-ready corr plot
pdf(file="KOG_insitu_site_corr_mcav_lpv.pdf", width=7, height=2.5)
par(mfrow=c(1,3))
corrPlot(x="ledge.south",y="central.south",ktable)
corrPlot(x="ledge.south",y="central.ledge",ktable)
corrPlot(x="central.south",y="central.ledge",ktable)
dev.off()

#---------------------------
# pairwise comparisons of sites (fc)

ledge.south=load('ledge.south_fc.RData')
ledge.south # names of datasets in the package
fc.ledge.south=kog.mwu(ledge.south.fc,gene2kog) 
fc.ledge.south

central.south=load('central.south_fc.RData')
central.south # names of datasets in the package
fc.central.south=kog.mwu(central.south.fc,gene2kog) 
fc.central.south

central.ledge=load('central.ledge_fc.RData')
central.ledge # names of datasets in the package
fc.central.ledge=kog.mwu(central.ledge.fc,gene2kog) 
fc.central.ledge

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("ledge.south"=fc.ledge.south,"central.south"=fc.central.south,"central.ledge"=fc.central.ledge))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_insitu_site_mcav_fc.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
# creating a pub-ready corr plot
pdf(file="KOG_insitu_site_corr_mcav_fc.pdf", width=7, height=2.5)
par(mfrow=c(1,3))
corrPlot(x="ledge.south",y="central.south",ktable)
corrPlot(x="ledge.south",y="central.ledge",ktable)
corrPlot(x="central.south",y="central.ledge",ktable)
dev.off()

#---------------------------
# pairwise comparisons of times (fc)

t1.0=load('t1.0_fc.RData')
t1.0 # names of datasets in the package
fc.t1.0=kog.mwu(t1.0.fc,gene2kog) 
fc.t1.0

t2.0=load('t2.0_fc.RData')
t2.0 # names of datasets in the package
fc.t2.0=kog.mwu(t2.0.fc,gene2kog) 
fc.t2.0

t3.0=load('t3.0_fc.RData')
t3.0 # names of datasets in the package
fc.t3.0=kog.mwu(t3.0.fc,gene2kog) 
fc.t3.0

t4.0=load('t4.0_fc.RData')
t4.0 # names of datasets in the package
fc.t4.0=kog.mwu(t4.0.fc,gene2kog) 
fc.t4.0

t5.0=load('t5.0_fc.RData')
t5.0 # names of datasets in the package
fc.t5.0=kog.mwu(t5.0.fc,gene2kog) 
fc.t5.0

t6.0=load('t6.0_fc.RData')
t6.0 # names of datasets in the package
fc.t6.0=kog.mwu(t6.0.fc,gene2kog) 
fc.t6.0

t2.1=load('t2.1_fc.RData')
t2.1 # names of datasets in the package
fc.t2.1=kog.mwu(t2.1.fc,gene2kog) 
fc.t2.1

t3.1=load('t3.1_fc.RData')
t3.1 # names of datasets in the package
fc.t3.1=kog.mwu(t3.1.fc,gene2kog) 
fc.t3.1

t4.1=load('t4.1_fc.RData')
t4.1 # names of datasets in the package
fc.t4.1=kog.mwu(t4.1.fc,gene2kog) 
fc.t4.1

t5.1=load('t5.1_fc.RData')
t5.1 # names of datasets in the package
fc.t5.1=kog.mwu(t5.1.fc,gene2kog) 
fc.t5.1

t6.1=load('t6.1_fc.RData')
t6.1 # names of datasets in the package
fc.t6.1=kog.mwu(t6.1.fc,gene2kog) 
fc.t6.1

t3.2=load('t3.2_fc.RData')
t3.2 # names of datasets in the package
fc.t3.2=kog.mwu(t3.2.fc,gene2kog) 
fc.t3.2

t4.2=load('t4.2_fc.RData')
t4.2 # names of datasets in the package
fc.t4.2=kog.mwu(t4.2.fc,gene2kog) 
fc.t4.2

t5.2=load('t5.2_fc.RData')
t5.2 # names of datasets in the package
fc.t5.2=kog.mwu(t5.2.fc,gene2kog) 
fc.t5.2

t6.2=load('t6.2_fc.RData')
t6.2 # names of datasets in the package
fc.t6.2=kog.mwu(t6.2.fc,gene2kog) 
fc.t6.2

t4.3=load('t4.3_fc.RData')
t4.3 # names of datasets in the package
fc.t4.3=kog.mwu(t4.3.fc,gene2kog) 
fc.t4.3

t5.3=load('t5.3_fc.RData')
t5.3 # names of datasets in the package
fc.t5.3=kog.mwu(t5.3.fc,gene2kog) 
fc.t5.3

t6.3=load('t6.3_fc.RData')
t6.3 # names of datasets in the package
fc.t6.3=kog.mwu(t6.3.fc,gene2kog) 
fc.t6.3

t5.4=load('t5.4_fc.RData')
t5.4 # names of datasets in the package
fc.t5.4=kog.mwu(t5.4.fc,gene2kog) 
fc.t5.4

t6.4=load('t6.4_fc.RData')
t6.4 # names of datasets in the package
fc.t6.4=kog.mwu(t6.4.fc,gene2kog) 
fc.t6.4

t6.5=load('t6.5_fc.RData')
t6.5 # names of datasets in the package
fc.t6.5=kog.mwu(t6.5.fc,gene2kog) 
fc.t6.5

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("1.0"=fc.t1.0,"2.0"=fc.t2.0,"3.0"=fc.t3.0,"4.0"=fc.t4.0,"5.0"=fc.t5.0,"6.0"=fc.t6.0,"2.1"=fc.t2.1,"3.1"=fc.t3.1,"4.1"=fc.t4.1,"5.1"=fc.t5.1,"6.1"=fc.t6.1,"3.2"=fc.t3.2,"4.2"=fc.t4.2,"5.2"=fc.t5.2,"6.2"=fc.t6.2,"4.3"=fc.t4.3,"5.3"=fc.t5.3,"6.3"=fc.t6.3,"5.4"=fc.t5.4,"6.4"=fc.t6.4,"6.5"=fc.t6.5))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_insitu_time_mcav_fc.pdf", width=10, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()