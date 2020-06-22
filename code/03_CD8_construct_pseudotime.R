phate <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_t/data/chris/cd8_phate2.csv')
meta <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/meta.rds')
active.ident <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/active.ident.rds')
meta = cbind(meta, active.ident = active.ident[rownames(meta)])
meta <- meta[meta$Tumor.Grade == 'IV' & meta$Treatment == 'Untreated', ]
selectcell <- rownames(meta[meta$active.ident %in% c('CD8 Naive', 'CD8 EFF', 'CD8 TEFF', 'CD8 RM', 'CD8 EM', 'CD8 EX'), ])
library(ggplot2)
pd = data.frame(x = phate[match(selectcell, as.character(phate[,1])),2], y = phate[match(selectcell, as.character(phate[,1])), 3], ct = active.ident[selectcell])
ggplot() + geom_point(data = pd, aes(x = x, y = y, col = ct), size = 0.1, alpha = 0.1) + theme_classic() +
  guides(color=guide_legend(override.aes = list(size=5, alpha = 1)))
set.seed(10)
clu <- kmeans(phate[match(selectcell, as.character(phate[,1])),2:3],6)$cluster
pd = data.frame(x = phate[match(selectcell, as.character(phate[,1])),2], y = phate[match(selectcell, as.character(phate[,1])), 3], ct = factor(clu))
ggplot() + geom_point(data = pd, aes(x = x, y = y, col = ct), size = 0.5) + theme_classic() +
  guides(color=guide_legend(override.aes = list(size=5, alpha = 1)))

library(TSCAN)
d = phate[match(selectcell, as.character(phate[,1])),2:3]
d = cbind(phate1 = d[,1], phate2 = d[,2])
rownames(d) = selectcell

mc <- exprmclust(t(d),cluster=clu,reduce=F)
ggplot() + geom_point(data = pd, aes(x = x, y = y, col = as.character(clu)), size = 0.5) + theme_classic() +
  guides(color=guide_legend(override.aes = list(size=5, alpha = 1)))
plotmclust(mc,show_full_tree=T)
ord <- TSCANorder(mc,MSTorder=c(4,1),orderonly=T)

####### slingshot
expr = readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/CD8_log2norm.rds')
dm <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_t/data/chris/cd8_diffusion.csv')
library(SingleCellExperiment)
sim <- SingleCellExperiment(assays = List(data = expr[, selectcell]))
rd = cbind(DC1 = dm$DC_1, DC2 = dm$DC_2)
rd = rd[match(selectcell, dm[,1]),]
rownames(rd) = selectcell
reducedDims(sim) <- SimpleList(PHATE = d, DiffMap = rd)

### cluster
library(mclust, quietly = TRUE)
cl1 <- Mclust(d)$classification
colData(sim)$GMM <- cl1
library(RColorBrewer)
plot(d, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

cl2 <- kmeans(d, centers = 20)$cluster
colData(sim)$kmeans <- cl2
plot(d, col = colorRampPalette(brewer.pal(9,"Set1"))(max(cl2))[cl2], pch=16, asp = 1, cex = 0.2)

### slingshot
library(princurve)
library(slingshot, quietly = TRUE)
sim <- slingshot(sim, clusterLabels = 'kmeans', reducedDim = 'PHATE')

### plot
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sim$slingPseudotime_1, breaks=100)]
plot(reducedDims(sim)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, col='black')

##### TSCAN
cm <- aggregate(d,list(cl2),mean)
colnames(cm)=c('cluster','phate1','phate2')
ggplot() + geom_point(data=data.frame(d,cl2=factor(cl2)),aes(x=phate1,y=phate2,col=cl2), size = 0.1) + geom_text(data=cm,aes(x=phate1,y=phate2,label=cluster)) + theme_bw() 

id <- which(cl2 %in% c(11,4,14,19,18,13))
names(cl2) <- row.names(d)
tmp <- cl2[id]
n <- names(tmp)
tmp <- as.numeric(as.factor(tmp))
names(tmp) <- n
mc <- exprmclust(t(d[id,]),cluster=tmp,reduce=F)
ord <- TSCANorder(mc,orderonly=T)
plotmclust(mc,  show_full_tree = T,
       cell_point_size = 0.2)
ggplot(data.frame(ct=active.ident[ord],pt=1:length(ord)),aes(x=pt,y=ct)) + geom_point(size = 0.1)

psn = cbind(Cell = ord, Pseudotime = 1:length(ord), Celltype = as.character(active.ident)[match(ord, names(active.ident))])
saveRDS(psn, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_t/result/cd8_EM_EX_pseudotime.rds')



########### check only CD8 EM and EX
selectcell2 <- rownames(meta[meta$active.ident %in% c('CD8 EM', 'CD8 EX'), ])
pd = data.frame(x = phate[match(selectcell2, as.character(phate[,1])),2], y = phate[match(selectcell2, as.character(phate[,1])), 3], ct = active.ident[selectcell2])
ggplot() + geom_point(data = pd, aes(x = x, y = y, col = ct), size = 0.1, alpha = 0.1) + theme_classic() +
  guides(color=guide_legend(override.aes = list(size=5, alpha = 1)))

set.seed(1)
clu <- kmeans(phate[match(selectcell2, as.character(phate[,1])),2:3],9)$cluster

pd = data.frame(x = phate[match(selectcell2, as.character(phate[,1])),2], y = phate[match(selectcell2, as.character(phate[,1])), 3], ct = factor(clu))
ggplot() + geom_point(data = pd, aes(x = x, y = y, col = ct), size = 0.5) + theme_classic() +
  guides(color=guide_legend(override.aes = list(size=5, alpha = 1)))


