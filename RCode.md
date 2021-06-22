# Relevant R code for the statistical analysis and visualization of sequencing data as performed for "Selective feeding in Southern Ocean key grazers – diet composition of krill and salps" by Pauli et al.

## Compositional analyses of sequencing data
````
### Filtering refined dataset and replacing zero counts
merged.df <- read.csv("FILENAME.CSV", header=T, row.names=1, sep=";",check.names=F)

df.filt <- merged.df[rowSums(merged.df)>100,colSums(merged.df)>300]
dim(df.filt)

require(zCompositions)
df.n0 <- t(cmultRepl(t(df.filt), method = "CZM", label = 0))
head(df.n0)
colSums(df.n0)# sanity check of data orientation

### Perform centered-log-ratio transformation
df.clr <- t(apply(df.n0, 2, function(x){log(x) - mean(log(x))}))
head(df.clr)

#alternatively use codaseq package
require(CoDaSeq)
coda.clr <- codaSeq.clr(df.n0, samples.by.row = FALSE)
head(coda.clr)

### Perform and plot principal component analyses
svd1 <- prcomp(df.clr)

pc1.1 <- round(svd1$sdev[1]^2/sum(svd1$sdev^2),2)
pc2.1 <- round(svd1$sdev[2]^2/sum(svd1$sdev^2),2)
xlab1 <- paste("PC1: ", pc1.1, sep="")
ylab1 <- paste("PC2: ", pc2.1, sep="")

par(mfrow=c(1,1))
biplot(svd1, cex=c(0.6,0.4), var.axes=F, scale=1, xlab=xlab1, ylab=ylab1)

summary(svd1)
screeplot(svd1)#inertia of single components

## Create groups for five sampling groups
kp <- df.filt[,1:14]#krill fecal pellets
sp <- df.filt[,15:25]#salp fecal pellets
kg <- df.filt[,26:86]#krill stomach content
sg <- df.filt[,87:146]#salp stomach content
p <- df.filt[,147:156]#water column

## Create nicer graph for PCA
rn<-rownames(df.clr)
svd2<- prcomp(df.clr[rn,])
svd2.mvar <- sum(svd2$sdev^2)

# Calculate the PC1 and PC2 variance
df0.PC1 <- paste("PC1: ", round(sum(svd2$sdev[1]^2)/svd2.mvar, 3))
df0.PC2 <- paste("PC2: ", round(sum(svd2$sdev[2]^2)/svd2.mvar, 3))

# Create colros for ellipses of groups
names = c("kp", "sp", "kg", "sg", "p")
# custom colours
rbcol <- c("coral1", "coral4", "deepskyblue3", "deepskyblue4","mediumseagreen")
cols <- c(rep(rbcol[2], ncol(kp)),rep(rbcol[4], ncol(sp)), 
          rep(rbcol[1], ncol(kg)),rep(rbcol[3], ncol(sg)),
          rep(rbcol[5], ncol(p)))
leg.col <-  c(rbcol[2],rbcol[4],rbcol[1],rbcol[3],rbcol[5])
rbcol.idx <- c(2,4,1,3,5)

pdf("FILENAME.pdf")

par(mfrow=c(1,2))
par(mfrow=c(1,1))

plot(svd2$x[,1], svd2$x[,2], pch=19, xlab=df0.PC1, ylab=df0.PC2,
     col=rgb(0,0,0,0.1),  main="Samples", cex.lab=1.2,
     xlim=c(min(svd2$x[,1]) -5, max(svd2$x[,1]) +5    ),
     ylim=c(min(svd2$x[,2]) -15, max(svd2$x[,2]) +15    )
)
abline(h=0, lty=2, lwd=2, col=rgb(0,0,0,0.3))
abline(v=0, lty=2, lwd=2, col=rgb(0,0,0,0.3))


for(i in 1:length(names)){
  jnk <- colnames(get(names[i]))[colnames(get(names[i])) %in% rn]

  dataEllipse(svd2$x[jnk,1], svd2$x[jnk,2],
              levels=c(0.75), center.cex=FALSE, plot.points=TRUE, add=TRUE, col=rbcol[rbcol.idx[i]],
              fill = TRUE, fill.alpha = 0.4)
}

legend.names = c("Krill Pellets", "Salp Pellets", "Krill Stomachs", "Salp Stomachs", "Plankton")
legend("topright", legend.names, col=leg.col, pch=19, cex=1.2, bg=rgb(0,0,0,0), box.lty = 0)

dev.off()

### Unsupervised clustering and plotting dendrogram
df.dist <- dist(df.clr, method="euclidian")
hc <- hclust(df.dist, method = "ward.D2")#method="average"|"ward.D2"

# Plot Dendrogram
par(mfrow=c(1,1))
plot(hc, cex=0.5)

# Plot using custom colurs for pre-defined groups
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    #labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    labCol <- if (grepl("KP", a$label) == TRUE ) "coral4" else if (grepl("SP", a$label) == TRUE ) "deepskyblue4"
    else if (grepl("KG", a$label) == TRUE ) "coral1" else if (grepl("SG", a$label) == TRUE ) "deepskyblue3"
    else if (grepl("F0", a$label) == TRUE ) "mediumseagreen" else "black" 
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}

hcd = as.dendrogram(hc)
clusDendro <- dendrapply(hcd, colLab)
plot(clusDendro)


### Permoanova using adonis
names = c("kp", "sp", "kg", "sg", "p")

for(i in 1:(length(names) -1)){
  names.kin <- colnames(eval(parse(text=names[i])))
  names.pup <- colnames(eval(parse(text=names[i+1])))
  
  kvsp.clr <- data.frame(rbind(coda.clr[names.kin,], coda.clr[names.pup,]))
  conds <- data.frame(c(rep("K",length(names.kin)), rep("P",length(names.pup))))
  colnames(conds) <- "grp"
  
  test <- adonis(kvsp.clr~grp, data=conds, method="euclidean", permutations=999)
  print(c(names[i], names[i+1], round(test$aov.tab$R2[1],3), round(test$aov.tab$F.Model[1],3),test$aov.tab[["Pr(>F)"]][1]))
}

# Anosim between groups using Aitchison distance
dist.clr <- dist(df.clr)
conds <- c(rep("kp", 14), rep("sp", 11),rep("kg", 61),rep("sg", 60),rep("p", 10))
ano <- anosim(dist.clr, conds, permutations=999)
summary(ano)

````
### Other multivariate analyses
For additional compositional data analyses  using the ALDEx2 package in R that was performed in this sudy, pelase see available scripts in the following publications:

- Bian, G. et al. The gut microbiota of healthy aged chinese is similar to that of the healthy young. mSphe- re 2, e00327-00317 (2017).

- Gloor, G. B. & Reid, G. Compositional analysis: a valid approach to analyze microbiome high-throughput sequencing data. Can. J. Microbiol. 62, 692-703 (2016).

## Selectivity analysis
```
require(dietr)
require(dplyr)
require(ggplot2)
require(reshape2)

PieChart <- read.csv("Realtive_Abundances.CSV", sep=";", header=TRUE)
str(PieChart)

PieChart_Available <- subset(PieChart, Group=="P")
PieChart_Available <- PieChart_Available[,-1]
str(PieChart_Available)

Piechart_Diet <- subset(PieChart, Group=="KG"|Group=="SG")
str(Piechart_Diet)

PieChart_Indices <- Electivity(Diet = Piechart_Diet, Available = PieChart_Available,
                         Indices = c("ForageRatio","Ivlev","Strauss","JacobsQ","JacobsD","Chesson",
                                     "VanderploegScavia"),LogQ = TRUE, Depleting = FALSE)


PlotElectivity(Electivity.Calcs = PieChart_Indices, Indices = "Ivlev")

write.csv(PieChart_Indices$Ivlev,"Selectivity_PieChartData_Ivlev.csv")

AllStations<-read.csv("Realtive_Abundances_PieChartData_Ivlev.CSV", sep=";", header=TRUE,check.names = "F")#slightly modified version of "Selectivity_PieChartData_Ivlev.csv"

str(AllStations)
colnames(AllStations)

AllStations_PieChart <- AllStations[,c(1:10,12:19)]#select most important taxa

AllSt_PieCharts_l<- melt(AllStations_PieChart, id.vars =c("Diet","Index"))
str(AllSt_PieCharts_l)

col <- c("firebrick3", "darkblue")

ggplot(AllSt_PieCharts_l, aes(y=value, x=variable, fill=Diet))+
  geom_bar(stat = "identity", position=position_dodge())+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, size=12))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.y = element_text(size=14))+
  theme(legend.text = element_text(size = 13))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="", y="Ivlev's selectivity Index")+
  theme(legend.title=element_blank())+
  scale_fill_manual(values=c("coral3","dodgerblue4"))+
  scale_y_continuous(breaks = seq(-1, 1, by=0.3))
````
