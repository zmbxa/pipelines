## have reps
get_EdgeR_DEout=function(counts,group,fdr=0.05,abs_FC=0,TMM_norm=T){
  design = model.matrix(~1+group)
  library(edgeR)
  
  y = DGEList(counts = counts,group = group)
  y = y[which(rowSums(edgeR::cpm(y)>=1)>=1),]
  
  ### TMM normalize
  if(TMM_norm){y<-calcNormFactors(y)}
  ### estimate dispersion
  y <- estimateGLMCommonDisp(y,design)
  y <- estimateGLMTrendedDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  # plotBCV(y)
  
  ### DE gene 
  fit_tag<-glmFit(y,design)  # fit GLM model
  lrt.tagwise<-glmLRT(fit_tag,coef=2)  # statistical analysis
  pvals_tag <- lrt.tagwise$table$PValue
  FDR_tag<- p.adjust(pvals_tag, method="BH")
  
  out = cbind(edgeR::cpm(y),lrt.tagwise$table,FDR_tag)
  out = out[order(out$PValue),]
  out$gene_name = rownames(out)
  out$change = ifelse(out$FDR_tag <fdr & abs(out$logFC)>abs_FC,
                      ifelse(out$logFC>abs_FC,"Up","Down"),"Stable")
  return(out)
}
## no replicates
get_EdgeR_DEout_bcv = function(counts,group,bcv=0.1,fdr=0.05,abs_FC=0,TMM_norm=T){
  library(edgeR)
  
  y = DGEList(counts = counts,group = group)
  keep <- rowSums(cpm(y)>1) >= 1
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  ### TMM normalize
  if(TMM_norm){y<-calcNormFactors(y)}
  
  ### DE gene 
  res <- exactTest(y, dispersion = bcv ^ 2)
  out <- cbind(cpm(y),res$table,FDR_tag = p.adjust(res$table$PValue,method = "BH"))
  out$gene_name = rownames(out)
  out$change = ifelse(out$FDR_tag <fdr & abs(out$logFC)>abs_FC,
                      ifelse(out$logFC>abs_FC,"Up","Down"),"Stable")
  return(arrange(out,FDR_tag))
}

save(get_EdgeR_DEout,get_EdgeR_DEout_bcv,file = "~/pipelines/RNAseq/analysis/get_edgeR_DE.RData")

## for DAVID chart GO enrichment
getDAVID = function(geneENS,readable=FALSE,DEout=NULL){
  ibrary(reticulate)
  source_python("/storage/zhangyanxiaoLab/niuyuxiao/pipelines/miniscripts/DAVIDChart4Rreticulate.py")
  inputG=paste(geneENS,collapse = ",")
  out_DAVID=david_gene_enrichment(gene_list = inputG)
  if(readable==TRUE){
    for (i in 1:nrow(out_DAVID)) {
      out_DAVID$geneSYMBOLs[i] = paste(DEout[match(as.character(limma::strsplit2(out_DAVID[i,"geneIds"],", ")),DEout$gene_name),"SYMBOL"],collapse = ", ")
    }
  }
}
save(getDAVID,file = "~/pipelines/RNAseq/analysis/get_DAVID_GOchart.RData")









