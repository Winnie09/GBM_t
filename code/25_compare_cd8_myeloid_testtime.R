library(Matrix)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
r1 <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/testtime/M_MDSC/EMDSC_MMDSC_PMNMDSC/fdr_foldchange.csv')
r2 <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_t/result/testtime/E_MDSC/phate_CD8EM_CD8EX/phate_fdr_foldchange.csv')
# M_MDSC/dm_CD8RM_CD8EX

o1 <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/EMDSC_MMDSC_PMNMDSC_pseudotime.rds')
o2 <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_t/data/order/cd8_EM_EX_pseudotime.rds')

res1 <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/testtime/M_MDSC/EMDSC_MMDSC_PMNMDSC/final.rds')
res2 <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_t/result/testtime/E_MDSC/phate_CD8EM_CD8EX/final.rds')

### get fitted values
get_fitted_values <- function(testObj, Gene, Pseudotime=NULL){
  ### testObj: test result from testpt_Time
  mat = matrix(0, nrow = length(Gene),ncol=length(Pseudotime),dimnames = list(Gene,NULL))
  num.knot = testObj$knotnum[Gene]
  if (is.null(Pseudotime)){
    Pseudotime = testObj$pseudotime
  }
  for (i in unique(num.knot)){
    x <- kronecker(diag(i + 4), 1)
    genes = names(num.knot[num.knot == i])
    if (i==0) {
      phi <- cbind(1,bs(Pseudotime))
    } else {
      knots = seq(min(Pseudotime),max(Pseudotime),length.out=i+2)[2:(i+1)]
      phi <- cbind(1,bs(Pseudotime,knots = knots))  
    }
    beta = sapply(genes,function(i) testObj$parameter[[i]]$beta)
    mat[genes,] <- t(phi %*% x %*% beta)
  } 
  return(mat)
}

fit = get_fitted_values(testObj = res1, Gene = as.character(r1$X), Pseudotime = o1[,2])
fit_scale <- t(apply(fit, 1, scale))

### elbow's method, determine number of clsuters
bttss <- sapply(1:20, function(i){
  set.seed(12345)
  clures = kmeans(fit_scale,  centers = i)
  (clures$totss - sum(clures$withinss))/clures$totss
})

### clustering
set.seed(12345)
clu = kmeans(fit_scale,  centers = 10)$cluster

### plot
par(mfrow=c(2,5))
for (i in  seq(1, max(clu))){
  plot(colMeans(fit_scale[names(clu[clu==i]), ]), pch = 20, xlab='Pseudotime', ylab='Scaled Expression', main = paste0('clu', i, ':', length(clu[clu==i])))
}
clu1 = clu

### repeat the process for the second trajectory
res = res2
r = r2
o = o2
fit = get_fitted_values(testObj = res, Gene = as.character(r$X), Pseudotime = as.numeric(o[,2]))
fit_scale <- t(apply(fit, 1, scale))

bttss <- sapply(1:20, function(i){
  set.seed(12345)
  clures = kmeans(fit_scale,  centers = i)
  (clures$totss - sum(clures$withinss))/clures$totss
})
plot(bttss ~ seq(1,20), pch=20, xlab='Number of clusters', ylab='BTSS/TOTSS')

set.seed(12345)
clu = kmeans(fit_scale,  centers = 10)$cluster

par(mfrow=c(2,5))
for (i in  seq(1, max(clu))){
  plot(colMeans(fit_scale[names(clu[clu==i]), ]), pch = 20, xlab='Pseudotime', ylab='Scaled Expression', main = paste0('clu', i, ':', length(clu[clu==i])))
}
clu2 = clu

########## myeloid up, t cell down
v <- sapply(1:10, function(i){
  length(intersect(names(clu1[clu1==2]), names(clu2[clu2==i])))
})



