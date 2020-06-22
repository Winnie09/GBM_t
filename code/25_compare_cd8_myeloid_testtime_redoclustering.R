library(Matrix)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
path1 = 'EMDSC_MMDSC_PMNMDSC'
path2 = 'dm_CD8EM_CD8EX'
rdir1 = '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/testtime/M_MDSC/'
rdir2 = '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_t/result/testtime/M_MDSC/'
ddir1 = '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/'
ddir2 = '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_t/data/order/'
sdir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_t/result/testtime/overlap/', path1, '_and_', path2,'/')
dir.create(sdir, recursive = T, showWarnings = F)
### myeloid
r1 <- read.csv(paste0(rdir1, path1, '/fdr_foldchange.csv'))
r1[,1] = as.character(r1[,1])
o1 <- readRDS(paste0(ddir1, path1, '_pseudotime.rds'))
res1 <- readRDS(paste0(rdir1, path1, '/final.rds'))

### t cell
r2 <- read.csv(paste0(rdir2, path2, '/fdr_foldchange.csv'))
r2[,1] = as.character(r2[,1])
o2 <- readRDS(paste0(ddir2, path2, '_pseudotime.rds'))
res2 <- readRDS(paste0(rdir2, path2, '/final.rds'))

### get overlap genes
g = intersect(as.character(r1[,1]), as.character(r2[,1]))
gmat = cbind(r1[match(g, r1[,1]), 2:3], r2[match(g, r2[,1]), 2:3])
colnames(gmat) = c(paste0(path1, '_', colnames(gmat)[1:2]), paste0(path2, '_', colnames(gmat)[3:4]))
saveRDS(g, paste0(sdir, 'DEG.rds'))
write.csv(gmat, paste0(sdir,'DEG_fdr_foldchange.csv'))

###### path1
cnt.bak1 <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/M.rds')
meta <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/meta.rds')
cnt <- cnt.bak1[, o1[,1]]
cellanno <- data.frame(cell = colnames(cnt), sample = sapply(colnames(cnt), function(i) sub('_.*','',sub('.*-','',i))), stringsAsFactors = FALSE)
mdsc <- read.csv('/home-4/whou10@jhu.edu/data2/whou10/GBM/meta/mdsc_proportions.csv', header = T)
design <- data.frame(MdscProp = mdsc[,8])  ## 3 is E-MDSC, 8 is M-MDSC
rownames(design) <- as.character(mdsc[,2])
cellanno <- cellanno[cellanno[,2] %in% rownames(design),]
cnt <- cnt[, cellanno[,1]]
design = cbind(1, design)
p <- plotGene(testptObj = res1, Gene=g, Mat=cnt[g, ], Pseudotime = o1, Cellanno=cellanno, Design=design,  Alpha=0.5, Size=0.1, PlotPoints = F, FreeScale = TRUE, BySample = F, type = 'Time')
ggsave(paste0(sdir, path1, '_DEG.png'), p, width=22, height=20, dpi = 100)
###### path2
cnt.bak <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/L.rds')
meta <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/meta.rds')
cnt <- cnt.bak[, o2[,1]]
cellanno <- data.frame(cell = colnames(cnt), sample = sapply(colnames(cnt), function(i) sub('_.*','',sub('.*-','',i))), stringsAsFactors = FALSE)
mdsc <- read.csv('/home-4/whou10@jhu.edu/data2/whou10/GBM/meta/mdsc_proportions.csv', header = T)
design <- data.frame(MdscProp = mdsc[,8])  ## 3 is E-MDSC, 8 is M-MDSC
rownames(design) <- as.character(mdsc[,2])
cellanno <- cellanno[cellanno[,2] %in% rownames(design),]
cnt <- cnt[, cellanno[,1]]
design = cbind(1, design)
p <- plotGene(testptObj = res2, Gene=g, Mat=cnt[g, ], Pseudotime = o2, Cellanno=cellanno, Design=design,  Alpha=0.5, Size=0.1, PlotPoints = F, FreeScale = TRUE, BySample = F, type = 'Time')
ggsave(paste0(sdir, path2, '_DEG.png'), p, width=22, height=20, dpi = 100)


#####  use intersect DEG to redoclustering
### get fitted values
fit = get_fitted_values(testObj = res1, Gene = g, Pseudotime = o1[,2])
fit_scale <- t(apply(fit, 1, scale))

### elbow's method, determine number of clsuters
bttss <- sapply(1:20, function(i){
  set.seed(12345)
  clures = kmeans(fit_scale,  centers = i)
  (clures$totss - sum(clures$withinss))/clures$totss
})
plot(bttss ~ seq(1,20), pch=20, xlab='Number of clusters', ylab='BTSS/TOTSS')

### clustering
set.seed(12345)
clu = kmeans(fit_scale,  centers = 10)$cluster

### plot
par(mfrow=c(2,5))
for (i in  seq(1, max(clu))){
  plot(colMeans(fit_scale[names(clu[clu==i]), ,drop=F]), pch = 20, xlab='Pseudotime', ylab='Scaled Expression', main = paste0('clu', i, ':', length(clu[clu==i])))
}
clu1 = clu

##### repeat the process for the second trajectory
res = res2
r = r2
o = o2
fit = get_fitted_values(testObj = res, Gene = g, Pseudotime = as.numeric(o[,2]))
fit_scale <- t(apply(fit, 1, scale))

bttss <- sapply(1:20, function(i){
  set.seed(12345)
  clures = kmeans(fit_scale,  centers = i)
  (clures$totss - sum(clures$withinss))/clures$totss
})
plot(bttss ~ seq(1,20), pch=20, xlab='Number of clusters', ylab='BTSS/TOTSS')

set.seed(12345)
clu = kmeans(fit_scale,  centers = 5)$cluster

par(mfrow=c(2,5))
for (i in  seq(1, max(clu))){
  plot(colMeans(fit_scale[names(clu[clu==i]), ]), pch = 20, xlab='Pseudotime', ylab='Scaled Expression', main = paste0('clu', i, ':', length(clu[clu==i])))
}
clu2 = clu

########## myeloid up, t cell down
v <- sapply(1:max(clu2), function(i){
  length(intersect(names(clu1[clu1==7]), names(clu2[clu2==i])))
})
v

df <- sapply(1:max(clu2), function(i){
  intersect(names(clu1[clu1==7]), names(clu2[clu2==i]))
})
names(df) = paste0(path1,'_clu7_', path2, '_clu', 1:max(clu2))
write.csv(as.matrix(df), paste0(sdir, path1,'_clu7.csv'))

v <- sapply(1:10, function(i){
  length(intersect(names(clu1[clu1==10]), names(clu2[clu2==i])))
})
v




