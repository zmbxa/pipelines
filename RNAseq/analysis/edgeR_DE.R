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
  library(reticulate)
  source_python("/storage/zhangyanxiaoLab/niuyuxiao/pipelines/miniscripts/DAVIDChart4Rreticulate.py")
  inputG=paste(geneENS,collapse = ",")
  out_DAVID=david_gene_enrichment(gene_list = inputG)
  if(readable==TRUE){
    for (i in 1:nrow(out_DAVID)) {
      out_DAVID$geneSYMBOLs[i] = paste(DEout[match(as.character(limma::strsplit2(out_DAVID[i,"geneIds"],", ")),DEout$gene_name),"SYMBOL"],collapse = ", ")
    }
  }
  return(out_DAVID)
}
save(getDAVID,file = "~/pipelines/RNAseq/analysis/get_DAVID_GOchart.RData")

### PCA Plot and batch effect
edgeR_PCA = function(counts = counts, group = group, batch=NULL,batch_rm=FALSE,hollow=FALSE,sex=sex){
  ### check parameters
  if(is.null(batch) & (batch_rm)){
    stop("Can not remove batch if 'batch' is empty!\n")
  }
  if(hollow){
    if(sum(grep("WT",group))==0 ) stop("You should give a 'WT' group for solid points!\n")
    if(is.null(sex)) stop("You should specify sex for coloring when hollow is True!\n")
  }
  design = model.matrix(~1+group)
  library(edgeR)
  y = DGEList(counts = counts,group = group)
  y = y[which(rowSums(edgeR::cpm(y)>1)>=2),]
  y=calcNormFactors(y)
  ## not given batch or sex
  
  if(!is.null(batch)){
    y$samples$group = group
    y$samples$batch = batch
  }
  if(!is.null(sex)){
    y$samples$group = group
    y$samples$sex = sex
    y$samples$batch = batch
  }
  if (batch_rm) {
    batch <- factor(y$samples$batch)
    logCPMs <- cpm(y, log = TRUE)
    rv <- apply(logCPMs, 1, var)  # Calculate rowwise variance
    logCPMs_corrected <- limma::removeBatchEffect(logCPMs, batch = batch)
    #  PCA using top 1k genes
    logCPM_corrected_top1000 <- logCPMs_corrected[head(order(rv, decreasing=TRUE), 1000),]
    pca <- prcomp(t(logCPM_corrected_top1000))
  }else{
    logCPMs <- cpm(y, log = TRUE)
    rv <- apply(logCPMs, 1, var)  # Calculate rowwise variance
    # From the logCPMs subset for the top-1000
    logCPM_top1000 <- logCPMs[head(order(rv, decreasing=TRUE), 1000),]
    # Run PCA
    pca <- prcomp(t(logCPM_top1000))
  }
  # Combine PCA coordinates with the metadata from the DGEList
  to_plot <- data.frame(pca$x, y$samples)
  # Calculate how many % of total variance is explained by each principal component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  # We focus here on PC1 and PC2
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  if(hollow){
    return(ggplot(to_plot, aes(x=PC1, y=PC2, color=tolower(sex), shape=batch))+theme_bw()+  xlab(labs[1]) + ylab(labs[2])+
             geom_point(data=filter(to_plot,group == "WT"),size=3)+scale_color_manual(values = c("brown2","royalblue"))+
             geom_point(data=filter(to_plot,group != "WT"),size=2,shape=ifelse(filter(to_plot,group != "WT")$batch==unique(batch)[1],2,1),stroke=1.3)+
             ggrepel::geom_text_repel(aes(label=rownames(to_plot)),size=3)+labs(color="gender"))
  }else{
    return(ggplot(to_plot, aes(x=PC1, y=PC2, color=group, shape=batch)) + 
             geom_point(size=3) +  xlab(labs[1]) + ylab(labs[2])+labs(color = NULL,shape=NULL))
  }
  
}
save(edgeR_PCA,file = "~/pipelines/RNAseq/analysis/PCAplot_EdgeR.RData")





