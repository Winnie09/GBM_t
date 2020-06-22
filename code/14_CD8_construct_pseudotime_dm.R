### GBM CD8 T cell
dm = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_t/data/chris/cd8_diffusion.csv')
rn = as.character(dm[,1])
dm = cbind(DC1 = dm[,2], DC2 = dm[,3])
rownames(dm) = rn
meta <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/meta.rds')
active.ident <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/active.ident.rds')
meta = cbind(meta, active.ident = active.ident[rownames(meta)])
meta <- meta[meta$Tumor.Grade == 'IV' & meta$Treatment == 'Untreated', ]

selectcell <- rownames(meta[meta$active.ident %in% c('CD8 Naive', 'CD8 EFF', 'CD8 TEFF', 'CD8 RM', 'CD8 EM', 'CD8 EX'), ])
dm = dm[selectcell, ]
library(ggplot2)
pd = data.frame(x = dm[,1], y = dm[,2], ct = active.ident[selectcell])
ggplot() + geom_point(data = pd, aes(x = x, y = y, col = ct), size = 0.1, alpha = 0.1) + theme_classic() +guides(color=guide_legend(override.aes = list(size=5, alpha = 1))) + xlab('DM1') + ylab('DM2') +
  # scale_color_brewer(palette = 'Dark2') + 
  facet_wrap(~ct)
set.seed(10)
clu <- kmeans(dm,6)$cluster
pd = data.frame(x = dm[,1], y = dm[,2], ct = factor(clu))
ggplot() + geom_point(data = pd, aes(x = x, y = y, col = ct), size = 0.1, alpha = 0.1) + theme_classic() +
  guides(color=guide_legend(override.aes = list(size=5, alpha = 1)))

library(TSCAN)
library(RColorBrewer)
clu = as.character(active.ident[rownames(dm)])
clu = as.numeric(as.factor(clu))
names(clu) = rownames(dm)
pd = data.frame(x = dm[,1], y = dm[,2], ct = factor(clu))
ggplot() + geom_point(data = pd, aes(x = x, y = y, col = ct), size = 0.1, alpha = 0.1) + theme_classic() +
  guides(color=guide_legend(override.aes = list(size=5, alpha = 1))) + xlab('DM1') + ylab('DM2')


id = which(clu %in% c(2,3))
newclu <- clu[id]
n <- names(newclu)
newclu <- as.numeric(as.factor(newclu))
names(newclu) <- n
mc <- exprmclust(t(dm[id, ]),cluster=newclu,reduce=F)
plotmclust(mc,show_full_tree=T, cell_point_size = 0.01)
ord <- TSCANorder(mc,MSTorder=c(1,2),orderonly=T)
ggplot(data.frame(ct=active.ident[ord],pt=1:length(ord)),aes(x=pt,y=ct)) + geom_point(size = 0.1)
ggplot(data.frame(x = dm[ord,1], y = dm[ord,2], t = seq(1, length(ord)))) + 
  geom_point(aes(x = x, y = y, col = t), size = 0.1) + 
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(11, 'RdYlBu'))(length(ord))) + 
  theme_classic() + xlab('DM1') + ylab('DM2')   ## useful
psn = data.frame(Cell = ord, Pseudotime = 1:length(ord), Celltype = as.character(active.ident)[match(ord, names(active.ident))], stringsAsFactors = FALSE)
saveRDS(psn, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_t/data/order/dm_CD8EM_CD8EX_pseudotime.rds')

#######
id = which(clu %in% c(5,3))
newclu <- clu[id]
n <- names(newclu)
newclu <- as.numeric(as.factor(newclu))
names(newclu) <- n
mc <- exprmclust(t(dm[id, ]),cluster=newclu,reduce=F)
plotmclust(mc,show_full_tree=T, cell_point_size = 0.01)
ord <- TSCANorder(mc,MSTorder=c(2,1),orderonly=T)
ggplot(data.frame(ct=active.ident[ord],pt=1:length(ord)),aes(x=pt,y=ct)) + geom_point(size = 0.1)
ggplot(data.frame(x = dm[ord,1], y = dm[ord,2], t = seq(1, length(ord)))) + 
  geom_point(aes(x = x, y = y, col = t), size = 0.1) + 
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(11, 'RdYlBu'))(length(ord))) + 
  theme_classic() + xlab('DM1') + ylab('DM2')   ## useful
psn = data.frame(Cell = ord, Pseudotime = 1:length(ord), Celltype = as.character(active.ident)[match(ord, names(active.ident))], stringsAsFactors = FALSE)
saveRDS(psn, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_t/data/order/dm_CD8RM_CD8EX_pseudotime.rds')

