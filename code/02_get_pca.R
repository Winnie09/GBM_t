cntnorm <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/CD8_log2norm.rds')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
pr <- PCA(cntnorm, maxVariableGenes = 2000)
saveRDS(pr, '/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/CD8_pc.rds')

library(umap)
set.seed(13)
u <- umap(d = pr[, 1:10])$layout
saveRDS(u, 'umap.rds')

