expr <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/L.rds')
meta <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/meta.rds')
active.ident <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/active.ident.rds')
meta = cbind(meta, active.ident = active.ident[rownames(meta)])
meta <- meta[meta$Grade == 'High' & meta$Treatment == 'Untreated', ]
selectcell <- rownames(meta[meta$active.ident %in% c('CD8 Naive', 'CD8 EFF', 'CD8 TEFF', 'CD8 RM', 'CD8 EM', 'CD8 EX'), ])
cnt <- expr[, selectcell]
saveRDS(cnt, '/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/Lymph/CD8_saver.rds') ## finished

cnt <- cnt[rowMeans(cnt>0.1)>0.01, ]
cv = apply(cnt, 1, sd)/rowMeans(cnt)
library(destiny)
dm2 <- DiffusionMap(cnt[names(sort(cv, decreasing = T)[1:3e3]),])
saveRDS(dm2, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_t/data/diffusionmap/CD8_3e3cvgene_diffusionmap.rds')
