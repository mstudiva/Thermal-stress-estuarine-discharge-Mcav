# install.packages("KOGMWU")
library(KOGMWU)

#---------------------------
# full model, assessing changes in KOG expression across colony, temp, and water

# loading KOG annotations
gene2kog=read.table("Mcavernosa_Cladocopium_iso2kogClass.tab",sep="\t")
head(gene2kog)

adc=load('colony_lpv.RData')
adc # names of datasets in the package
lpv.c=kog.mwu(colony.p,gene2kog) 
lpv.c 

adt=load('temp_lpv.RData')
adt # names of datasets in the package
lpv.t=kog.mwu(temp.p,gene2kog) 
lpv.t 

adw=load('water_lpv.RData')
adw # names of datasets in the package
lpv.w=kog.mwu(water.p,gene2kog) 
lpv.w 

adtw=load('tw_lpv.RData')
adtw # names of datasets in the package
lpv.tw=kog.mwu(tw.p,gene2kog) 
lpv.tw

barplot(lpv.c$delta.rank,horiz=T)

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("Colony"=lpv.c,"Temperature"=lpv.t,"Water"=lpv.w,"Temperapture:Water"=lpv.tw))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)
  
# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_exsitu_sym_lpv.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
# creating a pub-ready corr plot
pdf(file="KOG_exsitu_corr_sym_lpv.pdf", width=7, height=5)
par(mfrow=c(2,3))
corrPlot(x="Colony",y="Temperature",ktable)
corrPlot(x="Colony",y="Water",ktable)
corrPlot(x="Colony",y="Temperapture:Water",ktable)
corrPlot(x="Temperature",y="Water",ktable)
corrPlot(x="Temperature",y="Temperapture:Water",ktable)
corrPlot(x="Water",y="Temperapture:Water",ktable)
dev.off()

#---------------------------
# pairwise comparisons of temp.water treatments (lpv)
cd.co=load('cd.co_lpv.RData')
cd.co # names of datasets in the package
lpv.cd.co=kog.mwu(cd.co.p,gene2kog) 
lpv.cd.co 

eo.co=load('eo.co_lpv.RData')
eo.co # names of datasets in the package
lpv.eo.co=kog.mwu(eo.co.p,gene2kog) 
lpv.eo.co 

ed.co=load('ed.co_lpv.RData')
ed.co # names of datasets in the package
lpv.ed.co=kog.mwu(ed.co.p,gene2kog) 
lpv.ed.co

eo.cd=load('eo.cd_lpv.RData')
eo.cd # names of datasets in the package
lpv.eo.cd=kog.mwu(eo.cd.p,gene2kog) 
lpv.eo.cd 

ed.cd=load('ed.cd_lpv.RData')
ed.cd # names of datasets in the package
lpv.ed.cd=kog.mwu(ed.cd.p,gene2kog) 
lpv.ed.cd

ed.eo=load('ed.eo_lpv.RData')
ed.eo # names of datasets in the package
lpv.ed.eo=kog.mwu(ed.eo.p,gene2kog) 
lpv.ed.eo

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("eo.co"=lpv.eo.co,"cd.co"=lpv.cd.co,"ed.co"=lpv.ed.co,"eo.cd"=lpv.eo.cd,"ed.cd"=lpv.ed.cd,"ed.eo"=lpv.ed.eo))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_exsitu_tw_sym_lpv.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
# creating a pub-ready corr plot
pdf(file="KOG_exsitu_tw_corr_sym_lpv.pdf", width=10, height=10)
par(mfrow=c(4,4))
corrPlot(x="eo.co",y="cd.co",ktable)
corrPlot(x="eo.co",y="ed.co",ktable)
corrPlot(x="eo.co",y="eo.cd",ktable)
corrPlot(x="eo.co",y="ed.cd",ktable)
corrPlot(x="eo.co",y="ed.eo",ktable)
corrPlot(x="cd.co",y="ed.co",ktable)
corrPlot(x="cd.co",y="eo.cd",ktable)
corrPlot(x="cd.co",y="ed.cd",ktable)
corrPlot(x="cd.co",y="ed.eo",ktable)
corrPlot(x="ed.co",y="eo.cd",ktable)
corrPlot(x="ed.co",y="ed.cd",ktable)
corrPlot(x="ed.co",y="ed.eo",ktable)
corrPlot(x="eo.cd",y="ed.cd",ktable)
corrPlot(x="eo.cd",y="ed.eo",ktable)
dev.off()

#---------------------------
# pairwise comparisons of temp.water treatments (fc)

cd.co=load('cd.co_fc.RData')
cd.co # names of datasets in the package
fc.cd.co=kog.mwu(cd.co.fc,gene2kog) 
fc.cd.co 

eo.co=load('eo.co_fc.RData')
eo.co # names of datasets in the package
fc.eo.co=kog.mwu(eo.co.fc,gene2kog) 
fc.eo.co 

ed.co=load('ed.co_fc.RData')
ed.co # names of datasets in the package
fc.ed.co=kog.mwu(ed.co.fc,gene2kog) 
fc.ed.co

eo.cd=load('eo.cd_fc.RData')
eo.cd # names of datasets in the package
fc.eo.cd=kog.mwu(eo.cd.fc,gene2kog) 
fc.eo.cd 

ed.cd=load('ed.cd_fc.RData')
ed.cd # names of datasets in the package
fc.ed.cd=kog.mwu(ed.cd.fc,gene2kog) 
fc.ed.cd

ed.eo=load('ed.eo_fc.RData')
ed.eo # names of datasets in the package
fc.ed.eo=kog.mwu(ed.eo.fc,gene2kog) 
fc.ed.eo

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("eo.co"=fc.eo.co,"cd.co"=fc.cd.co,"ed.co"=fc.ed.co,"eo.cd"=fc.eo.cd,"ed.cd"=fc.ed.cd,"ed.eo"=fc.ed.eo))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_exsitu_tw_sym_fc.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
# creating a pub-ready corr plot
pdf(file="KOG_exsitu_tw_corr_sym_fc.pdf", width=10, height=10)
par(mfrow=c(4,4))
corrPlot(x="eo.co",y="cd.co",ktable)
corrPlot(x="eo.co",y="ed.co",ktable)
corrPlot(x="eo.co",y="eo.cd",ktable)
corrPlot(x="eo.co",y="ed.cd",ktable)
corrPlot(x="eo.co",y="ed.eo",ktable)
corrPlot(x="cd.co",y="ed.co",ktable)
corrPlot(x="cd.co",y="eo.cd",ktable)
corrPlot(x="cd.co",y="ed.cd",ktable)
corrPlot(x="cd.co",y="ed.eo",ktable)
corrPlot(x="ed.co",y="eo.cd",ktable)
corrPlot(x="ed.co",y="ed.cd",ktable)
corrPlot(x="ed.co",y="ed.eo",ktable)
corrPlot(x="eo.cd",y="ed.cd",ktable)
corrPlot(x="eo.cd",y="ed.eo",ktable)
dev.off()

#---------------------------
# pairwise comparisons of temp.water treatments (lpv) to contrpol.offshore only
cd.co=load('cd.co_lpv.RData')
cd.co # names of datasets in the package
lpv.cd.co=kog.mwu(cd.co.p,gene2kog) 
lpv.cd.co 

eo.co=load('eo.co_lpv.RData')
eo.co # names of datasets in the package
lpv.eo.co=kog.mwu(eo.co.p,gene2kog) 
lpv.eo.co 

ed.co=load('ed.co_lpv.RData')
ed.co # names of datasets in the package
lpv.ed.co=kog.mwu(ed.co.p,gene2kog) 
lpv.ed.co

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("eo.co"=lpv.eo.co,"cd.co"=lpv.cd.co,"ed.co"=lpv.ed.co))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_exsitu_pairwise_sym_lpv.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
# creating a pub-ready corr plot
pdf(file="KOG_exsitu_pairwise_corr_sym_lpv.pdf", width=7, height=2.5)
par(mfrow=c(1,3))
corrPlot(x="eo.co",y="cd.co",ktable)
corrPlot(x="eo.co",y="ed.co",ktable)
corrPlot(x="cd.co",y="ed.co",ktable)
dev.off()

#---------------------------
# pairwise comparisons of temp.water treatments (fc)

cd.co=load('cd.co_fc.RData')
cd.co # names of datasets in the package
fc.cd.co=kog.mwu(cd.co.fc,gene2kog) 
fc.cd.co 

eo.co=load('eo.co_fc.RData')
eo.co # names of datasets in the package
fc.eo.co=kog.mwu(eo.co.fc,gene2kog) 
fc.eo.co 

ed.co=load('ed.co_fc.RData')
ed.co # names of datasets in the package
fc.ed.co=kog.mwu(ed.co.fc,gene2kog) 
fc.ed.co

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("eo.co"=fc.eo.co,"cd.co"=fc.cd.co,"ed.co"=fc.ed.co))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_exsitu_pairwise_sym_fc.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
# creating a pub-ready corr plot
pdf(file="KOG_exsitu_pairwise_corr_sym_fc.pdf", width=7, height=2.5)
par(mfrow=c(1,3))
corrPlot(x="eo.co",y="cd.co",ktable)
corrPlot(x="eo.co",y="ed.co",ktable)
corrPlot(x="cd.co",y="ed.co",ktable)
dev.off()