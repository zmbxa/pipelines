add_doublets_info <- function(se) {
  # find doublet, this step should be done on each sample, if on merged data, could cause some issue
  # https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html
  # https://www.jianshu.com/p/6770c6a05287
  # https://mp.weixin.qq.com/s/V4IjVBpXA3jw9_HKrgq0ew
  
  se[["percent_mt"]] <- PercentageFeatureSet(se, pattern = "^mt-")
  se <-  SCTransform(se, method = "glmGamPoi", vars.to.regress = c("nFeature_RNA", "percent_mt"))
  se <- RunPCA(se, verbose = F)
  pc.num <- 1:15
  se <- RunUMAP(se, dims = pc.num)
  se <- FindNeighbors(se, dims = pc.num) %>% FindClusters(resolution = 0.5)
  
  # find best pk, it takes long time 
  # pK_bcmvn <- 0.09 
  sweep.res.list <- paramSweep_v3(se, PCs = pc.num, sct = T)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  paste("pK value:", pK_bcmvn)
  
  ## https://uofuhealth.utah.edu/huntsman/shared-resources/gba/htg/single-cell/genomics-10x.php
  ## 10000 cells with 7.6% doublet rate with homotypic adjustment
  ## TODO: try different rate and see if it matters
  doubletrate = ncol(se) * 7.6 * 1e-6
  homotypic.prop <- modelHomotypic(se$seurat_clusters)
  nExp <- round(doubletrate * ncol(se))
  nExp.adj <- round(nExp * (1-homotypic.prop))
  
  se <- doubletFinder_v3(se, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, nExp = nExp.adj, sct = T)
  
  # name of the DF prediction can change, so extract the correct column name.
  tissue <- se@project.name
  DF.name = colnames(se@meta.data)[grepl("DF.classification", colnames(se@meta.data))]
  p1 = DimPlot(se, group.by = DF.name) + NoAxes()
  ggsave(filename=paste0("plots/doublet_dimplot_", tissue,".png"), plot=p1)
  p1 = VlnPlot(se, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
  ggsave(paste0("plots/doublet_vlnplot_", tissue,".png"), plot=p1)
  
  colnames(se@meta.data)[grep("DF.classification", colnames(se@meta.data))] <- "DF.classification"
  colnames(se@meta.data)[grep("pANN", colnames(se@meta.data))] <- "pANN"
  
  return(se)
}
save(add_doublets_info,file = "~/pipelines/Paired-Tag-Pipeline/analysis/doubletFinder_add.RData")
