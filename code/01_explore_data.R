suppressMessages(library(Seurat))
obj <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/seuratObject/Lymph_ser.RDS')
cnt <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/count.rds')
cmeta <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/meta.rds')

active.ident = obj@active.ident
pca <- obj@reductions$pca@cell.embeddings
umap <- obj@reductions$umap@cell.embeddings
phate <- obj@reductions$phate@cell.embeddings
saveRDS(active.ident, '/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/active.ident.rds')
saveRDS(pca, '/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/pca.rds')
saveRDS(umap, '/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/umap.rds')
saveRDS(phate, '/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/phate.rds')

DimPlot(obj,label=T, pt.size=0.01)


meta = cbind(meta, active.ident = active.ident[rownames(meta)])
meta <- meta[meta$Grade == 'High' & meta$Treatment == 'Untreated', ]
selectcell <- rownames(meta[meta$active.ident %in% c('CD8 Naive', 'CD8 EFF', 'CD8 TEFF', 'CD8 RM', 'CD8 EM', 'CD8 EX'), ])

library(umap)
set.seed(12345)
u <- umap(d = pca[selectcell, 1:10])$layout


library(ggplot2)
ggplot() + geom_point(data = data.frame(x = pca[selectcell, 1], y = pca[selectcell, 2], ct = as.factor(active.ident[selectcell])), aes(x = x, y = y, col=ct), size=0.01, alpha=0.5) + theme_classic()+xlab('PC1')+ylab('PC2') + guides(color=guide_legend(override.aes = list(size=5, alpha=1)))

ggplot() + geom_point(data = data.frame(x = pca[selectcell, 1], y = pca[selectcell, 3], ct = as.factor(active.ident[selectcell])), aes(x = x, y = y, col=ct), size=0.01, alpha=0.5) + theme_classic()+xlab('PC1')+ylab('PC3') + guides(color=guide_legend(override.aes = list(size=5, alpha=1)))


ggplot() + geom_point(data = data.frame(x = pca[selectcell, 3], y = pca[selectcell, 2], ct = as.factor(active.ident[selectcell])), aes(x = x, y = y, col=ct), size=0.01, alpha=0.5) + theme_classic()+xlab('PC3')+ylab('PC2') + guides(color=guide_legend(override.aes = list(size=5, alpha=1)))

library("scatterplot3d")
library(RColorBrewer)
# mycolor <- factor(active.ident[selectcell], levels = brewer.pal(length(unique(active.ident[selectcell])), 'Set3'))
scatterplot3d(pca[selectcell, 1:3], pch=20)


pca <- pca[, 1:10]
library(umap)
set.seed(12345)
u <- umap(d = pca[selectcell, 1:10])$layout
ggplot() + geom_point(data = data.frame(x = u[,1], y = u[,2], ct = as.factor(active.ident[rownames(u)])), aes(x = x, y = y, col = ct), size = 0.01, alpha=0.1) + xlab('UMAP1') + ylab('UMAP2') + guides(color=guide_legend(override.aes = list(size=5,alpha=1))) + theme_classic() + facet_wrap(~ct)
saveRDS(u, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_t/data/order/subset_redo_umap.rds')

library(TSCAN)
set.seed(12345)
clu <- kmeans(pca[selectcell,], 6)$cluster
ggplot() + geom_point(data = data.frame(x = u[,1], y = u[,2], clu = as.factor(clu[rownames(u)])), aes(x = x, y = y, col = clu), size = 0.01) + xlab('UMAP1') + ylab('UMAP2') + guides(color=guide_legend(override.aes = list(size=5))) + theme_classic() 

mc <- exprmclust(t(pca[,1:3]),cluster=clu,reduce=F)
tapply(active.ident[selectcell],list(clu),table)

order <- rev(TSCANorder(exprmclust(t(pca[selectcell,1:3]),cluster=clu,reduce=F),orderonly = T, startcluster =2, listbranch = TRUE))

plotmclust(exprmclust(t(pca[selectcell,1:3]),cluster=clu,reduce=F),cell_point_size = 0.5,show_full_tree =T)

ggplot() + geom_point(data=data.frame(x = u[order,1], y = u[order, 2], pseudotime = seq(1, length(order))), aes(x = x, y = y, col=pseudotime), size = 0.1, alpha = 0.5) + scale_color_gradient(low = 'grey', high = 'red') + xlab('UMAP1') + ylab('UMAP2') + theme_classic()

ggplot() + geom_point(data=data.frame(x = pca[order,1], y = pca[order, 2], pseudotime = seq(1, length(order))), aes(x = x, y = y, col=pseudotime), size = 0.1, alpha = 0.5) + scale_color_gradient(low = 'grey', high = 'red') + xlab('PC1') + ylab('PC2') + theme_classic()


