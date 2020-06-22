Mclu = list()
Mclu[['EMDSC_MMDSC_PMNMDSC']] = c(1,3)
Mclu[['EMDSC_MMDSC']] = c(1,2)
Mclu[['EMDSC_MMDSC_MAC1']] = c(1,2)
Mclu[['EMDSC_MAC1']] = c(4)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
for (path1 in names(Mclu)){
  print(path1)
  for (path2 in c('dm_CD8EM_CD8EX', 'dm_CD8RM_CD8EX')){
    print(path2)
    rdir1 <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/testtime/M_MDSC/', path1, '/')
    rdir2 = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_t/result/testtime/M_MDSC/', path2, '/')
    rdir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_t/result/testtime/M_MDSC/', path1,'_and_', path2, '/')
    dir.create(rdir, recursive = T)
    
    res1 <- read.csv(paste0(rdir1, 'fdr_foldchange_cluster.csv'), row.names = 1, stringsAsFactors = F)
    res2 <- read.csv(paste0(rdir2, 'fdr_foldchange_cluster.csv'), row.names = 1, stringsAsFactors = F)
    
    for (i in Mclu[[path1]]){
      df <- sapply(seq(1, max(res2[,3])), function(j){
        g = intersect(rownames(res1[res1[,'Cluster']==i, ]), rownames(res2[res2[,'Cluster']==j, ]))
      })
        
      mlen = max(sapply(df, length))
      df <- sapply(seq(1, max(res2[,3])), function(j){
        g = intersect(rownames(res1[res1[,'Cluster']==i, ]), rownames(res2[res2[,'Cluster']==j, ]))
        c(g, rep("", mlen-length(g)))
      })
      colnames(df) <- paste0(path1, '_clu', i, '_and_', path2, '_clu', 1:ncol(df))
      write.csv(df, paste0(rdir, 'intersect_DEG.csv'))
    }
  }
}
    
