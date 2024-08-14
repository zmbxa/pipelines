#!/usr/bin/Rscript
##----------------------------
# @Author: Yuxiao Niu
##----------------------------
# QC plots
library(ggplot2)

argSet <- commandArgs(trailingOnly = TRUE)
setwd(argSet[1])
### quality check barplot
mapSum = read.csv("reportSummary.csv")
mapSum$subLib = gsub(".Report.*","",basename(mapSum$File))
mapSum$Fully.ligated=gsub(".*\\((.*)%\\).*", "\\1", mapSum$Fully.ligated)
mapSum$Useful.reads=gsub(".*\\((.*)%\\).*", "\\1", mapSum$Useful.reads)
mapSum$Duplication=gsub(".*\\((.*)%\\).*", "\\1", mapSum$Duplication)

bar_data = tidyr::gather(mapSum[,c(8,3:6)],key = "index",value = "stats",-'subLib',-'library')
bar_data$stats = as.numeric(bar_data$stats)

ggplot(bar_data,mapping = aes(x=index,y=stats,fill = library))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+xlab("subLibrary")+
  geom_text(aes(x=index,y=stats,label = paste0(stats,"%")), vjust = -0.3,position = position_dodge2(0.9))+
  geom_label(aes(label = subLib), vjust = 2.7,position = position_dodge2(0.9),
            color = "mediumpurple4",size = 3,family="sans",fontface = "bold",fill="ivory",label.size = NA,alpha=0.7)+
  ggtitle(paste0(mapSum$subLib[1],"-",mapSum$subLib[nrow(mapSum)],", Quality overview"))+ggsci::scale_fill_d3()
ggsave('MapQuality_Bar.pdf',width = 18,height = 8)

### barcodes reads and cells
library(Seurat)
library(dplyr)

mouse = CreateSeuratObject(counts = Read10X("merge/RNA_matrix/"),project = "scPairedTag",min.cells = 3,min.features = 20)
bc1 <- sub('.*:', '', colnames(mouse))
mouse$bc1 = bc1

### reads-bc
ReadsPerc_box_data = data.frame(bc1 = sub('.*:', '', colnames(mouse)),subLib = sub("^([^:]+):.*", "\\1", colnames(mouse)),nReads = mouse$nCount_RNA) %>% group_by(bc1,subLib) %>% summarise(totalReads = sum(nReads))
ReadsPerc_box_data$perc = ReadsPerc_box_data$totalReads / sum(ReadsPerc_box_data$totalReads)
### cell-bc
CellPerc_box_data = prop.table(table(sub("^([^:]+):.*", "\\1", colnames(mouse)),sub('.*:', '', colnames(mouse)))) %>% as.data.frame()
colnames(CellPerc_box_data) = c("subLib","barcodes","perc")
ggplot(ReadsPerc_box_data,aes(x=bc1,y=perc,fill = bc1))+geom_boxplot(outlier.size=0)+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4,show.legend = F)+
  geom_point(shape=16, show.legend = F)+scale_fill_brewer(palette="Set3")+
  ggtitle("RNA reads-barcode percentage")+
ggplot(CellPerc_box_data,aes(x=barcodes,y=perc,fill = barcodes))+geom_boxplot()+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4,show.legend = F)+
  geom_jitter(shape=16, position = position_jitter(0.2),show.legend = F)+scale_fill_brewer(palette="Set3")+
  ggtitle("RNA cells-barcode percentage")
ggsave("RNA_barcodesDistribution.pdf",width = 12,height = 4)


# DNA
mouse = CreateSeuratObject(counts = Read10X("merge/DNA_matrix/"),project = "scPairedTag",min.cells = 3,min.features = 20,assay="DNA")
bc1 <- sub('.*:', '', colnames(mouse))
mouse$bc1 = bc1

### reads-bc
ReadsPerc_box_data = data.frame(bc1 = sub('.*:', '', colnames(mouse)),subLib = sub("^([^:]+):.*", "\\1", colnames(mouse)),nReads = mouse$nCount_DNA) %>% group_by(bc1,subLib) %>% summarise(totalReads = sum(nReads))
ReadsPerc_box_data$perc = ReadsPerc_box_data$totalReads / sum(ReadsPerc_box_data$totalReads)
### cell-bc
CellPerc_box_data = prop.table(table(sub("^([^:]+):.*", "\\1", colnames(mouse)),sub('.*:', '', colnames(mouse)))) %>% as.data.frame()
colnames(CellPerc_box_data) = c("subLib","barcodes","perc")
ggplot(ReadsPerc_box_data,aes(x=bc1,y=perc,fill = bc1))+geom_boxplot(outlier.size=0)+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4,show.legend = F)+
  geom_point(shape=16, show.legend = F)+scale_fill_brewer(palette="Set3")+
  ggtitle("DNA reads-barcode percentage")+
  ggplot(CellPerc_box_data,aes(x=barcodes,y=perc,fill = barcodes))+geom_boxplot()+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4,show.legend = F)+
  geom_jitter(shape=16, position = position_jitter(0.2),show.legend = F)+scale_fill_brewer(palette="Set3")+
  ggtitle("DNA cells-barcode percentage")
ggsave("DNA_barcodesDistribution.pdf",width = 12,height = 4)

