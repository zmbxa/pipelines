## GSEA by fgsea
# library(fgsea)
# library(org.Mm.eg.db)
# library(msigdbr)
# library(tidyverse)

msigdb.mm = msigdbr(species = "Mus musculus",category = "C5",subcategory = "BP")
msigdb.hs = msigdbr(category = "C5",subcategory = "BP")

pathways = sapply(unique(msigdb.mm$gs_name), function(pathway_name){
  as.character(filter(msigdb.mm,gs_name==pathway_name)$"ensembl_gene")
})

rank_df = setNames(as.numeric(tmp$logFC),nm = gsub("\\..*","",tmp$gene_name)) %>% sort(decreasing = T)

out = tmp

get_fgsea_mm10 = function(out=NULL,rank_df=NULL,species="Mus musculus",msigdb = NULL,pathways = NULL,by="ensembl_gene",readable=T,colN_logFC = "logFC",colN_by = "gene_name"){
  library(fgsea)
  library(org.Mm.eg.db)
  library(msigdbr)
  library(tidyverse)
  
  # check mtx
  if(is.null(out) & is.null(rank_df))
    stop("Please give a outTable with logFC and ENSEMBL_id or a prepared rank_df!")
  if(is.null(rank_df)){
    if(length(grep(colN_logFC,colnames(out)))!=1)
      stop("Please check DE out matrix or colName! It should have only 1 logFC column T^T")
    if(length(grep(colN_by,colnames(out)))!=1)
      stop("Please check DE out matrix or colName! It should have only 1 gene_name column T^T")
    }
  
  # check database
  if(is.null(msigdb)) {
    print(paste0("Didn't provide msigdb, extracting from package for ",species," ---",date()))
    msigdb = msigdbr(species = species,category = "C5",subcategory = "BP")
  }
  if(is.null(pathways)){
    pathways = sapply(unique(msigdb$gs_name), function(pathway_name){
      as.character(filter(msigdb,gs_name==pathway_name)$"ensembl_gene")
    })
  }

  # rank_df not specified, use out mtx
  if(is.null(rank_df)){
    print(paste0("Begin GSEA analysis using given DEout mtx --- ",date()))
    colN_logFC = colnames(out)[grep(colN_logFC,colnames(out))]
    colN_by = colnames(out)[grep(colN_by,colnames(out))]
    rank_df=setNames(as.numeric(out[,colN_logFC]),nm = gsub("\\..*","",out[,colN_by])) %>% sort(decreasing = T)
  }
  # given rank_df
  else{
    print(paste0("Begin GSEA analysis using given rank_df --- ",date()))
  }
  names(rank_df) = gsub("\\..*","",names(rank_df))
  rank_df = sort(rank_df[na.omit(names(rank_df))],decreasing = T)
  GSEA_out <- fgsea(pathways = pathways, stats = rank_df)
  GSEA_out = dplyr::arrange(GSEA_out,padj)
  GSEA_out$Count = unlist(lapply(GSEA_out$leadingEdge, length))
  if(readable){
    to_map = data.frame(unique(msigdb[,c("gene_symbol", by)]))
    GSEA_out$geneIds = unlist(lapply(GSEA_out$leadingEdge, function(x){
      paste(to_map$gene_symbol[match(x,to_map[,by])],collapse = ", ")
    }))
  }
  return(GSEA_out)
}

save(msigdb.mm,pathways,get_fgsea_mm10,msigdb.hs,file = "~/pipelines/RNAseq/analysis/get_GSEA_mm10.RData")


#####
## for DAVID chart GO enrichment
getDAVID = function(geneENS,readable=FALSE,DEout=NULL){
  library(reticulate)
  source_python("/storage/zhangyanxiaoLab/niuyuxiao/pipelines/miniscripts/DAVIDChart4Rreticulate.py")
  inputG=paste(geneENS,collapse = ",")
  out_DAVID=david_gene_enrichment(gene_list = inputG)
  if(readable==TRUE){
    DEout$gene_name = gsub("\\..*","",DEout$gene_name)
    for (i in 1:nrow(out_DAVID)) {
      out_DAVID$geneSYMBOLs[i] = paste(DEout[match(as.character(limma::strsplit2(out_DAVID[i,"geneIds"],", ")),DEout$gene_name),"SYMBOL"],collapse = ", ")
    }
  }
  return(out_DAVID)
}
save(getDAVID,file = "~/pipelines/RNAseq/analysis/get_DAVID_GOchart.RData")


