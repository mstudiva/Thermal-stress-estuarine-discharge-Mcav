# run these once, then comment out
# install.packages("vegan")
# install.packages("rcompanion")
# install.packages("nparcomp")

library(vegan)
library(rcompanion)
library(nparcomp)

#--------------------------------------------
# Data import
sym <- read.csv(file="ex-situ_design.csv",header=TRUE,sep=",")

# displays the format of all variables in the dataframe
str(sym)
# converts variables from integer to factor for statistical analyses
sym$tank <- as.factor(sym$tank)
sym$block <- as.factor(sym$block)
sym$tag <- as.factor(sym$tag)
str(sym)

# Data transformations for large numbers
# performs linear transformation on variables with numbers in the millions to make it easier to graph
sym$zoox_cm <- sym$zoox_cm/1000000
sym$chla_cm <- sym$chla_cm/1000000
sym$chlc_cm <- sym$chlc_cm/1000000
str(sym)

#--------------------------------------------
# Normality tests

# creates a PDF for all the histograms
pdf("sym_hist.pdf", height=6, width=12)
# creates a plot window with 2 columns and 3 rows
par(mfrow=c(2,3))
# creates a histogram of data values for each variable
hist(sym$zoox_cm, breaks=15, xlab= expression(10^6 ~ cells ~ cm^-2), main=expression(Symbionts ~ cm^-2), col="grey50")
hist(sym$chla_cm, breaks=15, xlab= expression(ug ~ cm^-2), main=expression(Chl ~ a ~ cm^-2), col="grey50")
hist(sym$chlc_cm, breaks=15, xlab= expression(ug ~ cm^-2), main=expression(Chl ~ c[2] ~ cm^-2), col="grey50")
hist(sym$chl_ac, breaks=15, xlab= NA, main=expression(Chl ~ a  : c[2]), col="grey50")
hist(sym$chla_cell, breaks=15, xlab= expression(pg ~ cell^-1), main=expression(Chl ~ a ~ cell^-1), col="grey50")
hist(sym$chlc_cell, breaks=15, xlab= expression(pg ~ cell^-1), main=expression(Chl ~ c[2] ~ cell^-1), col="grey50")
# stops plotting, final save of the file
dev.off()

# runs a statistical test (Shapiro test) to determine if the data violates the assumptions of normality
# p < 0.05 means non-normal data
shapiro.test(sym$zoox_cm)
shapiro.test(sym$chla_cm)
shapiro.test(sym$chlc_cm)
shapiro.test(sym$chl_ac)
shapiro.test(sym$chla_cell)
shapiro.test(sym$chlc_cell)
# copy the output from each of these tests into a .txt file

#----------------------------------------------
# creating a distance matrix with Euclidean distance
symDist <- vegdist(sym[9:14], method="euclidean")

# testing for homogeneity of variance using betadisper
# setting seed allows you to reproduce exact results in iteration-based tests
set.seed(799)

colony.disp = betadisper(symDist, sym$colony)
anova(colony.disp)
# not significant

temp.disp = betadisper(symDist, sym$temp)
anova(temp.disp)
# not significant

water.disp = betadisper(symDist, sym$water)
anova(water.disp)
# not significant

# univariate PERMANOVAs
# run one for each response variable
zoox_cm.perm <- adonis(zoox_cm ~ colony * temp * water, data=sym, method= "euclidean", permutations = 9999, parallel=8)
chla_cm.perm <- adonis(chla_cm ~ colony * temp * water, data=sym, method= "euclidean", permutations = 9999, parallel=8)
chlc_cm.perm <- adonis(chlc_cm ~ colony * temp * water, data=sym, method= "euclidean", permutations = 9999, parallel=8)
chl_ac.perm <- adonis(chl_ac ~ colony * temp * water, data=sym, method= "euclidean", permutations = 9999, parallel=8)
chla_cell.perm <- adonis(chla_cell ~ colony * temp * water, data=sym, method= "euclidean", permutations = 9999, parallel=8)
chlc_cell.perm <- adonis(chlc_cell ~ colony * temp * water, data=sym, method= "euclidean", permutations = 9999, parallel=8)

# copy the outputs of each of these objects into the same .txt file as before
zoox_cm.perm
chla_cm.perm
chlc_cm.perm
chl_ac.perm
chla_cell.perm
chlc_cell.perm

#----------------------------------------------
# pairwise comparisons for compact letter display
zoox_cm.pair <- mctp(zoox_cm ~ temp.water, data=sym, type="Tukey", asy.method="fisher")
chla_cm.pair <- mctp(chla_cm ~ temp.water, data=sym, type="Tukey", asy.method="fisher")
chlc_cm.pair <- mctp(chlc_cm ~ temp.water, data=sym, type="Tukey", asy.method="fisher")
chl_ac.pair <- mctp(chl_ac ~ temp.water, data=sym, type="Tukey", asy.method="fisher")
chla_cell.pair <- mctp(chla_cell ~ temp.water, data=sym, type="Tukey", asy.method="fisher")
chlc_cell.pair <- mctp(chlc_cell ~ temp.water, data=sym, type="Tukey", asy.method="fisher")

summary(zoox_cm.pair)
summary(chla_cm.pair)
summary(chlc_cm.pair)
summary(chl_ac.pair)
summary(chla_cell.pair)
summary(chlc_cell.pair)

# export each of the Analysis outputs to a .csv file
write.csv(zoox_cm.pair$Analysis, file="zoox_cm_pair.csv")
write.csv(chla_cm.pair$Analysis, file="chla_cm_pair.csv")
write.csv(chlc_cm.pair$Analysis, file="chlc_cm_pair.csv")
write.csv(chl_ac.pair$Analysis, file="chl_ac_pair.csv")
write.csv(chla_cell.pair$Analysis, file="chla_cell_pair.csv")
write.csv(chlc_cell.pair$Analysis, file="chlc_cell_pair.csv")

# combine all .csv files into one sheet and reformat in a way that has comparisons as a single column and each response variable as columns, then reimport
# make sure to rename factor levels in the comparisons to match what you're using
sym.pair <- read.csv(file="PERMANOVA_outputs.csv", head=T)
str(sym.pair)

zoox_cm.list <- cldList(zoox_cm ~ Comparison, data = sym.pair, threshold = 0.05)
chla_cm.list <- cldList(chla_cm ~ Comparison, data = sym.pair, threshold = 0.05)
chlc_cm.list <- cldList(chlc_cm ~ Comparison, data = sym.pair, threshold = 0.05)

# commenting out due to no significant pairwise comparisons
# chl_ac.list <- cldList(chl_ac ~ Comparison, data = sym.pair, threshold = 0.05)
# chla_cell.list <- cldList(chla_cell ~ Comparison, data = sym.pair, threshold = 0.05)
# chlc_cell.list <- cldList(chlc_cell ~ Comparison, data = sym.pair, threshold = 0.05)

#----------------------------------------------
# 3-way PERMANOVA
# non-parametric multivariate test to identify effects of colony, temp, and water on the 6 response variables
permanova <- adonis(symDist ~ colony * temp * water, data=sym, permutations=9999, parallel=8)
# copy the outputs of this object into the same .txt file as before
permanova

# pairwise comparisons for colony
colony.pair <- pairwise.perm.manova(symDist, sym$colony, nperm=9999)
# copy the outputs of each of these objects into the same .txt file as before
colony.pair

#----------------------------------------------
#----------------------------------------------
# figures

# install.packages("ggplot2")
# install.packages("gridExtra")
# install.packages("ggpubr")
# install.packages("FSA")
# install.packages("rcompanion")

library(ggplot2)
library(gridExtra)
library(ggpubr)
library(FSA)
library(rcompanion)

# generates citations for theses, manuscripts, etc.
# citation(package="FSA")
# citation(package="ggpubr")

#----------------------------------------------
# statistical tests for plot overlay

# maintains temp.water order as imported
sym$temp.water=factor(sym$temp.water, levels=unique(sym$temp.water)) 

#-------------------
# plotting means of each response variable across time then depth using boxplots with data point scatter overlay

# creates a manual color palette for graphs
colors= c(control.offshore="#92c5de", control.discharge="#0571b0", elevated.offshore="#f4a582", elevated.discharge="#ca0020")

# plots boxplot of zooxcm with sample points
# palette forces colors from above color palette
# theme() removes x axis label
zooxcm <-
  ggboxplot(
    sym,
    x = "temp.water",
    y = "zoox_cm",
    color = "grey30",
    fill = "temp.water",
    palette = colors,
    add = "jitter",
    width = 0.7,
    size = 0.75) + 
  labs(x = "",
       y = "Symbiont Density \n(106 cells cm-2)",
       title = "Symbiont Density",
       fill = 'temp.water') + 
  ylim(0.25,4.25) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") + 
  annotate("text",
           x = 0.5,
           y = 0.56,
           hjust = 0,
           label = "geno: p=0.0001\ntemp: p=0.0001\nwater: p=0.3771")  + 
  geom_text(data = zoox_cm.list, 
            aes(x = Group, 
                y = 4.25, 
                label =Letter)) 

zooxcm

# has y axis label only
chlacm <-
  ggboxplot(
    sym,
    x = "temp.water",
    y = "chla_cm",
    color = "grey30",
    fill = "temp.water",
    palette = colors,
    add = "jitter",
    width = 0.7,
    size = 0.75) + 
  labs(x = "",
       y = "Chlorophyll Density \n(ug cm-2)",
       title = "Areal Chlorophyll a",
       fill = 'temp.water') + 
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")  + 
  annotate("text",
           x = 0.5,
           y = 2.6,
           hjust = 0,
           label = "geno: p=0.0001\ntemp: p=0.0001\nwater: p=0.0612")  + 
  geom_text(data = chla_cm.list, 
            aes(x = Group, 
                y = 14, 
                label =Letter)) 

chlacm

# has y axis label only
chlccm <-
  ggboxplot(
    sym,
    x = "temp.water",
    y = "chlc_cm",
    color = "grey30",
    fill = "temp.water",
    palette = colors,
    add = "jitter",
    width = 0.7,
    size = 0.75) + 
  labs(x = "",
       y = "Chlorophyll Density \n(ug cm-2)",
       title = "Areal Chlorophyll c2",
       fill = 'temp.water') + 
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")  + 
  annotate("text",
           x = 0.5,
           y = 0.6,
           hjust = 0,
           label = "geno: p=0.0001\ntemp: p=0.0001\nwater: p=0.3788")  + 
  geom_text(data = chlc_cm.list, 
            aes(x = Group, 
                y = 4.5, 
                label =Letter)) 

chlccm

# has both axes labels
chlac <-
  ggboxplot(
    sym,
    x = "temp.water",
    y = "chl_ac",
    color = "grey30",
    fill = "temp.water",
    palette = colors,
    add = "jitter",
    width = 0.7,
    size = 0.75) + 
  labs(x = "",
       y = "Chlorophyll a:c2\n",
       title = "Chlorophyll a:c2",
       fill = 'temp.water') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")  + 
  annotate("text",
           x = 0.5,
           y = 7.1,
           hjust = 0,
           label = "geno: p=0.0068\ntemp: p=0.0197\nwater: p=0.2855") +
  scale_x_discrete(labels = c("Control\nOffshore","Control\nDischarge","Elevated\nOffshore","Elevated\nDischarge"))

chlac

# has both axis labels
chlacell <-
  ggboxplot(
    sym,
    x = "temp.water",
    y = "chla_cell",
    color = "grey30",
    fill = "temp.water",
    palette = colors,
    add = "jitter",
    width = 0.7,
    size = 0.75) + 
  labs(x = "",
       y = "Chlorophyll cell-1 \n(pg cell-1)",
       title = "Cellular Chlorophyll a",
       fill = 'temp.water') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust =0.5),
        legend.position = "none")  + 
  annotate("text",
           x = 0.5,
           y = 7.7,
           hjust = 0,
           label = "geno: p=0.0401\ntemp: p=0.8461\nwater: p=0.4129") +
  scale_x_discrete(labels = c("Control\nOffshore","Control\nDischarge","Elevated\nOffshore","Elevated\nDischarge"))

chlacell

# has x axis label only
chlccell <-
  ggboxplot(
    sym,
    x = "temp.water",
    y = "chlc_cell",
    color = "grey30",
    fill = "temp.water",
    palette = colors,
    add = "jitter",
    width = 0.7,
    size = 0.75) + 
  labs(x = "",
       y = "Chlorophyll cell-1 \n(pg cell-1)",
       title = "Cellular Chlorophyll c2",
       fill = 'temp.water') + 
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")  + 
  annotate("text",
           x = 0.5,
           y = 3,
           hjust = 0,
           label = "geno: p=0.0145\ntemp: p=0.5062\nwater: p=0.6886") +
  scale_x_discrete(labels = c("Control\nOffshore","Control\nDischarge","Elevated\nOffshore","Elevated\nDischarge"))

chlccell

# creates plot object with a 3x2 grid of plots
# forces layout to be zooxcm, chlacm, chlccm, chlac, chlacell, chlccell
plot=grid.arrange(zooxcm, chlac, chlacm, chlacell, chlccm, chlccell, ncol=3, nrow=2, layout_matrix=cbind(c(1,2),c(3,4),c(5,6)), widths=c(3,3.05,2.75), heights=c(3,3.4))

#saves plot as PDF and EPS
ggsave("sym_boxplot.pdf", plot= plot, width=12, height=7, units="in", dpi=300)
ggsave("sym_boxplot.eps", plot= plot, width=12, height=7, units="in", dpi=300)

#-------------------