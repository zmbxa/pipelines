##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/1/17
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)
library(stringr)

argSet <- commandArgs(trailingOnly = TRUE)
# file1-Feature-Cell matrix
counts <- Read10X(data.dir=argSet[1])

# Perform standard Seurat analysis for RNA
rna <- CreateSeuratObject(counts=counts, min.cells=3)
rna <- subset(rna, subset=nFeature_RNA > 1000)
rna <- NormalizeData(rna, normalization.method='LogNormalize')
rna <- FindVariableFeatures(rna, nfeatures=1000)

# remove cell cycle influence, use hg38 genes as reference
s.genes <- str_to_title(Seurat::cc.genes$s.genes)
g2m.genes <- str_to_title(Seurat::cc.genes$g2m.genes)
rna <- CellCycleScoring(rna, s.features=s.genes, g2m.features=g2m.genes, set.ident=T)

#-----------------------------
rna <- ScaleData(rna, vars.to.regress=c('S.Score', 'G2M.Score'))
rna <- RunPCA(rna, features=VariableFeatures(rna), verbose=F)

# ElbowPlot(rna, ndims=50)
# ggsave('PCA_dimension.pdf', width=4, height=4)

dimUse <- 2:30

rna <- FindNeighbors(rna, dims=dimUse)
rna <- FindClusters(rna, resolution=.5)
rna <- RunUMAP(rna, dims=dimUse)

nNuclei <- paste(gsub(' ', '', format(ncol(rna), big.mark=',')), 'nuclei')
posX <- (max(rna@reductions$umap@cell.embeddings[, 'UMAP_1']))*0.88+(min(rna@reductions$umap@cell.embeddings[, 'UMAP_1']))*0.12
posY <- max(rna@reductions$umap@cell.embeddings[, 'UMAP_2'])

DimPlot(rna, reduction='umap', label=TRUE, label.size=5, pt.size=.1)+
	labs(title=NULL, x='UMAP dimension 1', y='UMAP dimension 2')+
	theme_classic()+scale_color_igv()+
	annotate('text', x=posX, y=posY, label=nNuclei, color='#555555')+
	scale_x_continuous(expand=expansion(mult=c(.05, .05), add=0))+
#	scale_y_continuous(limits=c(-13, 8.3), expand=c(0, 0))+
	theme(plot.margin=unit(c(.01, .8, 0, .01), 'inches'),
		plot.title=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.75, color='black'),
		legend.title=element_blank(),
		legend.text=element_text(size=11),
		legend.key=element_blank(),
		legend.key.size=unit(.2, 'in'),
		legend.spacing.x=unit(0, 'inches'),
		legend.spacing.y=unit(0, 'inches'),
		legend.position=c(1.12, .5),
		legend.background=element_blank(),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.title.x=element_text(colour='black', size=12),
		axis.title.y=element_text(colour='black', size=12),
		axis.text.x=element_blank(),
		axis.text.y=element_blank())
ggsave('01_UMAP.pdf', width=4.6, height=4)

featureGenes <- c('Pdgfra', 'Gad2', 'Slc17a7', 'Apoe', 'Mog', 'Sst')
FeaturePlot(rna, features=featureGenes,
	cols=c('#D3D3D3', '#FF8C00'), ncol=3)
ggsave('02_Features.pdf', width=8, height=4)


clusters <- rna@meta.data
clusters$Barcode <- rownames(clusters)
write.table(clusters[, c('Barcode', 'seurat_clusters')], 'Cell_Clusters.txt', row.names=F, col.names=F, sep='\t', quote=F)
save(rna, file='RNA_Clusters_Seurat.RData')


# markers <- FindAllMarkers(rna, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
# top10Markers <- markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
# write.table(top10Markers, 'top10Markers.txt', row.names=F, sep='\t', quote=F)

# DoHeatmap(rna, features=top10Markers$gene, angle=0, group.colors=pal_igv('default')(9))+
# 	theme_classic()+
# 	scale_fill_gradient2(low='#0073c2', mid='#FFFFF0', high='#efc000')+
# 	theme(plot.margin=unit(c(.05, .01, 0, .01), 'inches'),
# 		plot.title=element_blank(),
# 		legend.title=element_text(size=15),
# 		legend.text=element_text(size=15),
# 		legend.key=element_blank(),
# 		legend.key.size=unit(.3, 'in'),
# 		legend.spacing.x=unit(0, 'inches'),
# 		legend.spacing.y=unit(0, 'inches'),
# 		legend.background=element_blank(),
# 		axis.line.x=element_blank(),
# 		axis.line.y=element_blank(),
# 		axis.ticks.x=element_blank(),
# 		axis.ticks.y=element_blank(),
# 		axis.title.x=element_blank(),
# 		axis.title.y=element_blank(),
# 		axis.text.x=element_blank(),
# 		axis.text.y=element_text(colour='black', size=12))
# ggsave('UMAP_heatmap.pdf', width=15, height=20)
