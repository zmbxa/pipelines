## have reps
### add batch adjust
get_EdgeR_DEout=function(counts,group,fdr=0.05,abs_FC=0,TMM_norm=T,filter=T,IDtrans=T,batch_rm = F,batch=NULL){
  if(batch_rm & is.null(batch)){
    stop("Can not remove batch if 'batch' is empty!\n")
  }
  library(edgeR)

  design = model.matrix(~1+group)
  y = DGEList(counts = counts,group = group)
  if (filter) {y = y[which(rowSums(edgeR::cpm(y)>=1)>=1),]}
  
  if(batch_rm){
    y$samples$batch = batch
    design <- model.matrix(~batch+group, y$samples)
    
  }
  
  ### TMM normalize
  if(TMM_norm){y<-calcNormFactors(y)}
  ### estimate dispersion
  # y <- estimateCommonDisp(y,design)
  y<-estimateCommonDisp(y)
  # y <- estimateGLMTrendedDisp(y, design)
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
  if(IDtrans){
    if(sum(startsWith(out$gene_name,"ENSMUSG"))) {
      print("mouse, do transfering")
      gene_SYMBOL = read.csv("/storage/zhangyanxiaoLab/niuyuxiao/annotations/genes_SYMBOL_mm10.csv")
      out = merge(out,gene_SYMBOL,by = "gene_name",all.x = T) %>% arrange(PValue)
      out$SYMBOL = ifelse(is.na(out$SYMBOL),out$gene_name,out$SYMBOL)
      out$TE = ifelse(grepl(":",out$gene_name),"TE","gene")
    }
    if(sum(startsWith(out$gene_name,"ENSG"))) {
      print("human, do transfering")
      gene_SYMBOL = read.csv("/storage/zhangyanxiaoLab/niuyuxiao/annotations/genes_SYMBOL_hg38.csv")
      out = merge(out,gene_SYMBOL,by = "gene_name",all.x = T) %>% arrange(PValue)
      out$SYMBOL = ifelse(is.na(out$SYMBOL),out$gene_name,out$SYMBOL)
      out$TE = ifelse(grepl(":",out$gene_name),"TE","gene")
    }
  }
  return(out)
}


## no replicates
get_EdgeR_DEout_bcv = function(counts,group,bcv=0.1,fdr=0.05,abs_FC=0,TMM_norm=T,filter=T,IDtrans=T){
  library(edgeR)
  
  y = DGEList(counts = counts,group = group)
  if(filter){keep <- rowSums(cpm(y)>1) >= 1
  y <- y[keep, , keep.lib.sizes=FALSE]}
  
  ### TMM normalize
  if(TMM_norm){y<-calcNormFactors(y)}
  
  ### DE gene 
  res <- exactTest(y, dispersion = bcv ^ 2)
  out <- cbind(cpm(y),res$table,FDR_tag = p.adjust(res$table$PValue,method = "BH"))
  out$gene_name = rownames(out)
  out$change = ifelse(out$FDR_tag <fdr & abs(out$logFC)>abs_FC,
                      ifelse(out$logFC>abs_FC,"Up","Down"),"Stable")
  
  if(IDtrans){
    if(sum(startsWith(out$gene_name,"ENSMUSG"))) {
      print("mouse, do transfering")
      gene_SYMBOL = read.csv("/storage/zhangyanxiaoLab/niuyuxiao/annotations/genes_SYMBOL_mm10.csv")
      out = merge(out,gene_SYMBOL,by = "gene_name",all.x = T) %>% arrange(PValue)
      out$SYMBOL = ifelse(is.na(out$SYMBOL),out$gene_name,out$SYMBOL)
      out$TE = ifelse(grepl(":",out$gene_name),"TE","gene")
    }
    if(sum(startsWith(out$gene_name,"ENSG"))) {
      print("human, do transfering")
      gene_SYMBOL = read.csv("/storage/zhangyanxiaoLab/niuyuxiao/annotations/genes_SYMBOL_hg38.csv")
      out = merge(out,gene_SYMBOL,by = "gene_name",all.x = T) %>% arrange(PValue)
      out$SYMBOL = ifelse(is.na(out$SYMBOL),out$gene_name,out$SYMBOL)
      out$TE = ifelse(grepl(":",out$gene_name),"TE","gene")
    }
  }
  return(arrange(out,FDR_tag))
}

### modified volcano plot
plot_volcano = function(out,FDR_c="FDR_tag",logFC_c="logFC",pointLabel=TRUE,label_c="SYMBOL",label_num = NULL,pointSize_c=NULL,color_c = "change",colors=NULL,numLabel=TRUE,title="volcano Plot",subtitle=NULL){
  library(ggplot2);library(ggsci);library(ggpubr)
  if(is.null(colors)){
    if(nlevels(factor(out[,color_c]))==3)
      colors = c(Down = "royalblue",Stable="gray",Up="firebrick")
    else{
      if(nlevels(factor(out[,color_c]))==5)
        colors=c(`Down gene`="royalblue1",`Down TE`="navy",`Stable`="grey", `Up gene`="brown1",`Up TE`="firebrick4")
    }
  }
  
  # make basic plot
  p=ggplot(out,aes(x=get(logFC_c),y=-log(get(FDR_c)),color=get(color_c)))+theme_bw()+labs(color=color_c)+
    theme(plot.title = element_text(hjust = 0.5,face="bold"))+xlab(logFC_c)+ylab(paste0("-log(",FDR_c,")"))
  if(is.null(pointSize_c))
    p=p+geom_point(alpha=0.65)
  else
    p=p+geom_point(aes(size=get(pointSize_c)),alpha=0.65)+labs(size=pointSize_c)
  
  # mod color
  if(!is.null(colors)){
    p=p+scale_color_manual(values = colors)
  }
  # subtitle
  if(!is.null(subtitle))
    p=p+ggtitle(title,subtitle)
  else
    p=p+ggtitle(title)
  # add number label
  if(nlevels(factor(out[,color_c]))==5){
    p=p+annotate("label", x = min(out[,logFC_c]), y = max(-log(out[,FDR_c])), label = nrow(out[which(out[,color_c]=="Down gene"),]), 
                 vjust = 1, hjust = 0,colour="royalblue1",size=4,fontface = "bold")+
      annotate("label", x = max(out[,logFC_c]), y = max(-log(out[,FDR_c])), label = nrow(out[which(out[,color_c]=="Up gene"),]), 
               vjust = 1, hjust = 1,colour="brown1",size=4,fontface = "bold")+
      annotate("label", x = 0, y = 0, label = nrow(out[which(out[,color_c]=="Stable"),]), 
               vjust = 2, hjust = 0,colour="grey66",size=4,fontface = "bold")+
      annotate("label", x = min(out[,logFC_c]), y = max(-log(out[,FDR_c])), label = nrow(out[which(out[,color_c]=="Down TE"),]), 
               vjust = 2, hjust = 0,colour="navy",size=4,fontface = "bold")+
      annotate("label", x = max(out[,logFC_c]), y = max(-log(out[,FDR_c])), label = nrow(out[which(out[,color_c]=="Up TE"),]), 
               vjust = 2, hjust = 1,colour="firebrick4",size=4,fontface = "bold")
  }
  if (nlevels(factor(out[,color_c]))==3) {
    p=p+annotate("label", x = min(out[,logFC_c]), y = max(-log(out[,FDR_c])), label = nrow(out[which(out[,color_c]=="Down"),]), 
                 vjust = 1, hjust = 0,colour="royalblue1",size=4,fontface = "bold")+
      annotate("label", x = 0, y = 0, label = nrow(out[which(out[,color_c]=="Stable"),]), 
               vjust = 0, hjust = 0.5,colour="grey33",size=4,fontface = "bold")+
      annotate("label", x = max(out[,logFC_c]), y = max(-log(out[,FDR_c])), label = nrow(out[which(out[,color_c]=="Up"),]), 
               vjust = 1, hjust = 1,colour="brown1",size=4,fontface = "bold")
  }
  # add point label
  if(pointLabel){
    if(is.null(label_num))
      p=p+ggrepel::geom_text_repel(data=filter(out,get(color_c) != "Stable") %>% head(label_num),size=3,mapping = aes(label=get(label_c)),show.legend = F)
    else
      p=p+ggrepel::geom_text_repel(data=head(out,label_num),size=3,mapping = aes(label=get(label_c)),show.legend = F)
  }
  return(p)
}
### modified MA plot
plot_MA = function(out,logCPM_c="logCPM",logFC_c="logFC",pointLabel=TRUE,label_c="SYMBOL",label_num = NULL,pointSize_c=NULL,color_c = "change",colors=NULL,numLabel=TRUE,title="volcano Plot",subtitle=NULL){
  library(ggplot2);library(ggsci);library(ggpubr)
  if(is.null(colors)){
    if(nlevels(factor(out[,color_c]))==3)
      colors = c(Down = "royalblue",Stable="gray",Up="firebrick")
    else{
      if(nlevels(factor(out[,color_c]))==5)
        colors=c(`Down gene`="royalblue1",`Down TE`="navy",`Stable`="grey", `Up gene`="brown1",`Up TE`="firebrick4")
    }
  }
  
  # make basic plot
  p=ggplot(out,aes(y=get(logFC_c),x=get(logCPM_c),color=get(color_c)))+theme_bw()+labs(color=color_c)+
    theme(plot.title = element_text(hjust = 0.5,face="bold"))+ylab(logFC_c)+xlab(logCPM_c)
  if(is.null(pointSize_c))
    p=p+geom_point(alpha=0.65)
  else
    p=p+geom_point(aes(size=get(pointSize_c)),alpha=0.65)+labs(size=pointSize_c)
  
  # mod color
  if(!is.null(colors)){
    p=p+scale_color_manual(values = colors)
  }
  # subtitle
  if(!is.null(subtitle))
    p=p+ggtitle(title,subtitle)
  else
    p=p+ggtitle(title)
  # add number label
  if(nlevels(factor(out[,color_c]))==5){
    p=p+annotate("label", y = min(out[,logFC_c]), x = max(out[,logCPM_c]), label = nrow(out[which(out[,color_c]=="Down gene"),]), 
                 hjust = 1, vjust = 0,colour="royalblue1",size=4,fontface = "bold")+
      annotate("label", y = max(out[,logFC_c]), x = max(out[,logCPM_c]), label = nrow(out[which(out[,color_c]=="Up gene"),]), 
               hjust = 1, vjust = 1,colour="brown1",size=4,fontface = "bold")+
      annotate("label", x = 0, y = 0, label = nrow(out[which(out[,color_c]=="Stable"),]), 
               hjust = 2, vjust = 0,colour="grey66",size=4,fontface = "bold")+
      annotate("label", y = min(out[,logFC_c]), x = max(out[,logCPM_c]), label = nrow(out[which(out[,color_c]=="Down TE"),]), 
               hjust = 2, vjust = 0,colour="navy",size=4,fontface = "bold")+
      annotate("label", y = max(out[,logFC_c]), x = max(out[,logCPM_c]), label = nrow(out[which(out[,color_c]=="Up TE"),]), 
               hjust = 2, vjust = 1,colour="firebrick4",size=4,fontface = "bold")
  }
  if (nlevels(factor(out[,color_c]))==3) {
    p=p+annotate("label", y = min(out[,logFC_c]), x = max(out[,logCPM_c]), label = nrow(out[which(out[,color_c]=="Down"),]), 
                 hjust = 1, vjust = 0,colour="royalblue1",size=4,fontface = "bold")+
      annotate("label", x = 0, y = 0, label = nrow(out[which(out[,color_c]=="Stable"),]), 
               hjust = 0, vjust = 0.5,colour="grey33",size=4,fontface = "bold")+
      annotate("label", y = max(out[,logFC_c]), x = max(out[,logCPM_c]), label = nrow(out[which(out[,color_c]=="Up"),]), 
               hjust = 1, vjust = 1,colour="brown1",size=4,fontface = "bold")
  }
  # add point label
  if(pointLabel){
    if(is.null(label_num))
      p=p+ggrepel::geom_text_repel(data=filter(out,get(color_c) != "Stable") %>% head(label_num),size=3,mapping = aes(label=get(label_c)),show.legend = F)
    else
      p=p+ggrepel::geom_text_repel(data=head(out,label_num),size=3,mapping = aes(label=get(label_c)),show.legend = F)
  }
  return(p)
}

save(get_EdgeR_DEout,get_EdgeR_DEout_bcv,plot_volcano,plot_MA,file = "~/pipelines/RNAseq/analysis/get_edgeR_DE.RData")



### PCA Plot and batch effect
edgeR_PCA = function(counts = counts, group = group, batch=NULL,batch_rm=FALSE,label=T,hollow=FALSE,sex=NULL,sourceTab=FALSE){
  load("~/pipelines/RNAseq/analysis/colorset.RData")
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
  }
  if(!is.null(sex) & !is.null(batch)){
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
    if(!is.null(batch)){
      graph = (ggplot(to_plot, aes(x=PC1, y=PC2, color=tolower(sex), shape=batch))+theme_bw()+  xlab(labs[1]) + ylab(labs[2])+
               geom_point(data=to_plot,size=3)+ggsci::scale_color_startrek()+
               geom_point(data=filter(to_plot,group != "WT"),size=1.2,color="white",show.legend = F)+
               labs(color="gender"))
    }else{
      graph = (ggplot(to_plot, aes(x=PC1, y=PC2, color=tolower(sex)))+theme_bw()+  xlab(labs[1]) + ylab(labs[2])+
                 geom_point(data=filter(to_plot,group == "WT"),size=3)+geom_point(data=filter(to_plot,group != "WT"),size=2,stroke=1.3,shape=1)+
                 ggsci::scale_color_startrek()+
                 labs(color="gender"))
    }
  }else{
    graph = (ggplot(to_plot, aes(x=PC1, y=PC2, color=group, shape=batch)) + 
               geom_point(size=3) +  xlab(labs[1]) + ylab(labs[2])+labs(color = NULL,shape=NULL))
  }
  if(label){
    graph = graph+ggrepel::geom_text_repel(aes(label=rownames(to_plot)),size=3)
  }
  if(sourceTab){
    return(list(plot=graph,data=to_plot))
  }else 
    return(graph)
  
}
save(edgeR_PCA,file = "~/pipelines/RNAseq/analysis/PCAplot_EdgeR.RData")





