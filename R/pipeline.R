library(Seurat)
library(stringi)
library(dplyr)
library(SingleR)
library(celldex)
library(readxl)
library(ggplot2)
library(SeuratDisk)
library(dittoSeq)
library(gridExtra)
library(scales)

# Reading aggregated H5 file ####
adata = Read10X_h5("count/filtered_feature_bc_matrix.h5")
so = CreateSeuratObject(adata)

# Annotating batches and treatment
batches = c()
for (y in so@assays[["RNA"]]@counts@Dimnames[[2]]){
  if (stri_sub(y,-1) == 1){batches <- append(batches, 1)}
  if (stri_sub(y,-1) == 2){batches <- append(batches, 2)}
  if (stri_sub(y,-1) == 3){batches <- append(batches, 3)}
  if (stri_sub(y,-1) == 4){batches <- append(batches, 4)}
}
so[["Batches"]] = batches

Batches = so[["Batches"]]
Batches = replace(Batches,Batches == "1", "N-Ex")
Batches = replace(Batches,Batches == "2", "N-Sed")
Batches = replace(Batches,Batches == "3", "O-Ex")
Batches = replace(Batches,Batches == "4", "O-Sed")
so[["Batches"]] = Batches

Exercised = c()
for (i in so[["Batches"]]){
  Exercised = append(Exercised,i)
}
Exercised = replace(Exercised,Exercised == 1|Exercised == 3, "Ex")
Exercised = replace(Exercised,Exercised == 2|Exercised == 4, "Sed")
so[["Exercised"]] = Exercised

saveRDS(so, file = "output/pig_Batches.rds")

# Checkpoint ####
so = readRDS("output/pig_Batches.rds")

#Pre-processing graphs
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Normalization
so <- NormalizeData(so)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(so), 10)

# Scaling and Dimensionality Reduction
all.genes <- rownames(so)
so <- ScaleData(so, features = all.genes)
so <- RunPCA(so, features = VariableFeatures(object = so))
DimPlot(so, reduction = "pca")

#Clustering
so <- FindNeighbors(so, dims = 1:10)
so <- FindClusters(so, resolution = 0.5)

#UMAP
so <- RunUMAP(so, dims = 1:10)
DimPlot(so, reduction = "umap")

## By batches
jpeg('output/Pig_Cluster_Batches.jpeg',width = 1200, height = 1200, quality = 100,pointsize = 100)
DimPlot(so, reduction = "umap", label = F,
        pt.size = 0.5, group.by = 'Batches')+
  ggtitle('UMAP by Batches')+
  geom_text(aes(x = 3, y = -4.5, label = "N-Ex"),size = 8)+
  geom_text(aes(x = 7, y = 1, label = "N-Sed"),size = 8)+
  geom_text(aes(x = -3.5, y = -6, label = "O-Ex"),size = 8)+
  geom_text(aes(x = 0.7, y = 2, label = "O-Sed"),size = 8)+
  theme(axis.title = element_text(size = 30),
        plot.title = element_text(size=30),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.5, 'cm'))
dev.off()
saveRDS(so, file = "output/pig_UMAP_new.rds")

# Checkpoint ####
so = readRDS("output/pig_UMAP_new.rds")

# Manual annotation of orthologs of some important genes
geneold = c('ENSSSCG00045037762',
            'ENSSSCG00045008562',
            'ENSSSCG00045011271',
            'ENSSSCG00045009046',
            'ENSSSCG00045038246',
            'ENSSSCG00045019298',
            'ENSSSCG00045010998',
            'ENSSSCG00045031665',
            'ENSSSCG00045028985',
            'ENSSSCG00045026854',
            'ENSSSCG00045008986',
            'ENSSSCG00045003554',
            'ENSSSCG00045012558',
            'ENSSSCG00045008615',
            'ENSSSCG00045019911',
            'ENSSSCG00045034516',
            'ENSSSCG00045022800',
            'ENSSSCG00045016758',
            'ENSSSCG00045018110',
            'ENSSSCG00045011067',
            'ENSSSCG00045022911',
            'ENSSSCG00045000930',
            'ENSSSCG00045031502',
            'ENSSSCG00045020561',
            'ENSSSCG00045002945',
            'ENSSSCG00045033061',
            'ENSSSCG00045017315',
            'ENSSSCG00045002110',
            'ENSSSCG00045038038',
            'ENSSSCG00045009746',
            'ENSSSCG00045002714')

genenew = c('LRRC43',
            'FN1',
            'EGF7',
            'RYR2',
            'RBFOX1',
            'LRP1B',
            'KCNB2',
            'TMEM132D',
            'ADGRB3',
            'PLSCR',
            'MAP3K1',
            'DENND5A',
            'ADIPOQ',
            'PTP',
            'FRMD',
            'TTN',
            'FAT3',
            'TNFRSF21',
            'PDE3B',
            'TERF',
            'ARMC2',
            'PKHD1L1',
            'LOC110258364',
            'CST9',
            'ADAMTSL3',
            'CSF1',
            'ZNF395',
            'MACF1',
            'CNKSR2',
            'CORIN',
            'CALD1')

rnames = row.names(so@assays$RNA@counts)
for (i in 1:length(geneold)){
  row.names(so@assays$RNA@counts) =
    replace(row.names(so@assays$RNA@counts),
            row.names(so@assays$RNA@counts)==geneold[i],genenew[i])
  print(genenew[i] %in% row.names(so@assays$RNA@counts))
}
for (i in 1:length(geneold)){
  row.names(so@assays$RNA@data) =
    replace(row.names(so@assays$RNA@data),
            row.names(so@assays$RNA@data)==geneold[i],genenew[i])
  print(genenew[i] %in% row.names(so@assays$RNA@data))
}
for (i in 1:length(geneold)){
  row.names(so@assays$RNA@meta.features) =
    replace(row.names(so@assays$RNA@meta.features),
            row.names(so@assays$RNA@meta.features)==geneold[i],genenew[i])
  print(genenew[i])
  print(genenew[i] %in% row.names(so@assays$RNA@data))
}

features = so@assays$RNA@var.features
for (i in length(geneold)){
  features[features == geneold[i]] = genenew[i]
}

row.names(so@assays$RNA@counts) = rnames
row.names(so@assays$RNA@data) = rnames


# Manual annotation of orthologs of some important genes - SECOND TIME
manann = read.csv("manann3.csv",header = T)
geneold = manann$Gene
genenew = manann$Gene_name

rnames = row.names(so@assays$RNA@counts)
for (i in 1:length(geneold)){
  row.names(so@assays$RNA@counts) =
    replace(row.names(so@assays$RNA@counts),
            row.names(so@assays$RNA@counts)==geneold[i],genenew[i])
}
for (i in 1:length(geneold)){
  row.names(so@assays$RNA@data) =
    replace(row.names(so@assays$RNA@data),
            row.names(so@assays$RNA@data)==geneold[i],genenew[i])
}
for (i in 1:length(geneold)){
  row.names(so@assays$RNA@meta.features) =
    replace(row.names(so@assays$RNA@meta.features),
            row.names(so@assays$RNA@meta.features)==geneold[i],genenew[i])
}


features = so@assays$RNA@var.features
for (i in length(geneold)){
  features[features == geneold[i]] = genenew[i]
}
saveRDS(so, file = "output/pig_UMAP_clustered_new.rds")

# Checkpoint ####
so = readRDS("output/pig_UMAP_clustered_new.rds")

#Find marker genes
cluster0.markers <- FindMarkers(so, ident.1 = 0, min.pct = 0.25)
Marker.genes = head(cluster0.markers, n = 10)
Marker.genes[['Cluster']] = 0

for (i in 1:16){
  cluster.markers <- FindMarkers(so, ident.1 = i, min.pct = 0.25)
  Mgenes = head(cluster.markers, n = 10)
  Mgenes[['Cluster']] = i
  Marker.genes = bind_rows(Marker.genes,Mgenes)
}

write.csv(Marker.genes,"PigCellTypeMarkers.csv")
marker_genes = read.csv("PigMarkers.csv")

## Celltype annotation
labels = read_xls("PigCellTypeMarkers_11232022.xls")
labels = labels$Annotations
labels = labels[!is.na(labels)]
cell_types = c()
for (i in so[["seurat_clusters"]]){
  cell_types = append(cell_types,labels[i])
}
so[["cell_types"]] = cell_types

#Removing 1 and 11 Clusters
Idents(so) = "cell_types"
so = subset(so, idents = "Unknown", invert=TRUE)

saveRDS(so, file = "output/pig_UMAP_clustered_new.rds")
saveRDS(so, file = "output/pig_v3.rds")


# Checkpoint ####
so = readRDS('output/pig_v3.rds')

#Batch Visualization
jpeg('output/Pig_batches_final.jpeg',width = 720, height = 540)
DimPlot(so, reduction = "umap", label = TRUE, repel = TRUE, label.size = 5, pt.size = 0.5, group.by = 'Batches')
dev.off()

#Cluster Visualization
jpeg('Results/9-18-23_Publication_Figures//Pig_cluster_Celltype.jpeg',width = 2000, height = 2000, quality = 100,pointsize = 10)
DimPlot(so, reduction = "umap", pt.size = 3, group.by = 'cell_types',)+
  geom_text(aes(x = 9, y = -3.7, label = "B cells"),size = 20)+
  geom_text(aes(x = 1, y = 7.3, label = "Beige adipocytes"),size = 20)+
  geom_text(aes(x = -6, y = -1, label = "Endothelial cells"),size = 20)+
  geom_text(aes(x = 8.2, y = 2, label = "Erythroid cells"),size = 20)+
  geom_text(aes(x = -5, y = -7, label = "Fibroblasts"),size = 20)+
  geom_text(aes(x = 6, y = -8, label = "Macrophages"),size = 20)+
  geom_text(aes(x = -6.3, y = 6, label = "Mesenchymal"),size = 20)+
  geom_text(aes(x = -7.6, y = 5, label = "cells"),size = 20,)+
  geom_text(aes(x = 8, y = -6, label = "Smooth muscle cells"),size = 20)+
  geom_text(aes(x = 1, y = -2.6, label = "T cells"),size = 20)+
  theme(axis.title = element_text(size = 60),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 40),
        plot.title = element_blank(),
        legend.text = element_text(size = 60),
        legend.key.size = unit(2, 'cm'),
        legend.position=c(.65,.83))+
  guides(color = guide_legend(override.aes=list(shape = 15,size = 10)))
dev.off()


jpeg('Results/9-18-23_Publication_Figures/Pig_batches_final.jpeg',width = 2000, height = 2000, quality = 100,pointsize = 10)
DimPlot(so, reduction = "umap", pt.size = 3, group.by = 'Batches')+
  theme(axis.title = element_text(size = 60),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 40),
        plot.title = element_blank(),
        legend.text = element_text(size = 80),
        legend.key.size = unit(3, 'cm'),
        legend.position=c(.77,.85))+
  guides(color = guide_legend(override.aes=list(shape = 15,size = 10)))
dev.off()


So_endo = subset(so, idents = c('Endothelial cells'))
jpeg('Results/Other_Figures/Pig_Endothelial.jpeg',width = 2000, height = 2000, quality = 100,pointsize = 10)
DimPlot(So_endo, reduction = "umap", label.size = 5, pt.size = 5, group.by = 'Batches')+
  theme(axis.title = element_text(size = 40),
        plot.title = element_blank(),
        legend.text = element_text(size = 60),
        legend.key.size = unit(3, 'cm'),
        legend.position=c(.8,.9),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())+
  guides(color = guide_legend(override.aes=list(shape = 15,size = 15)))
dev.off()

#Violin plots showing nFeature_RNA of cells
jpeg('output/Pig_vln_celltype.jpeg',width = 1000, height = 800, quality = 100,pointsize = 100)
VlnPlot(so, features = c("nFeature_RNA"), split.by = "Exercised",pt.size = 0)
dev.off()

#Box plots showing number of Feature genes detected in each cell type
jpeg('output/Pig_Boxplot_nFeatures.jpeg',width = 1200, height = 800, quality = 100,pointsize = 100)
Percelltype = ggplot(so@meta.data, aes(x = nFeature_RNA, y = cell_types, fill=Exercised)) +
  coord_flip() +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=0, notch=FALSE)+
  scale_color_brewer(palette="Dark2")+
  ggtitle('nFeature Genes vs Cell Type')+
  xlab('nFeature Genes')+
  ylab('Cell Types')+
  theme(axis.title = element_text(size = 30),
        plot.title = element_text(color="red", size=30),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, 'cm'))
Percelltype
dev.off()

#Histogram Graph cell type Percentage
celltypes = rep(unique(so@meta.data$cell_types),2)
celltypedf = as.data.frame(celltypes)
colnames(celltypedf) = 'Cell_Types'
celltypedf[['Cells']] = 0
celltypedf[['Type']] = 0
num=1
for (i in unique(celltypes)){
  celltypedf[['Type']][num] = 'Ex'
  celltypedf[['Cells']][num] = 100*(length(so@meta.data[so[["Exercised"]]=="Ex" & so[["cell_types"]]==i])/length(so@meta.data[so[["Exercised"]]=="Ex"]))
  num = num+1
}
for (i in unique(celltypes)){
  celltypedf[['Type']][num] = 'Sed'
  celltypedf[['Cells']][num] = 100*(length(so@meta.data[so[["Exercised"]]=="Sed" & so[["cell_types"]]==i])/length(so@meta.data[so[["Exercised"]]=="Sed"]))
  num = num+1
}

# Bar plots
jpeg('Results/2-6-23_Poster_Figures/Pig_Barplot_celltype.jpeg',width = 2000, height = 3000, quality = 100,pointsize = 10)
gg1 = dittoBarPlot(so,var = "Batches", group.by = "cell_types",scale = "count")+
  scale_fill_manual("Batches",values=c("#F68282", "#7CAE00", "#00BFC4","#C77CFF"))+
  xlab("Cell Type")+
  theme(text = element_text(size = 50),
        legend.text = element_text(size = 50),
        legend.key.size = unit(3,"cm"),
        legend.position = c(.9,.85),
        plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 50, face = "bold"))
gg2 = dittoBarPlot(so,var = "Batches", group.by = "cell_types")+
  scale_fill_manual("Batches",values=c("#F68282", "#7CAE00", "#00BFC4","#C77CFF"))+
  xlab("Cell Type")+
  theme(text = element_text(size = 50),
        legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 50, angle = 90, vjust = 0.5, hjust=1),
        axis.title = element_text(size = 50, face = "bold"))
grid.arrange(gg1,gg2, ncol=1)
dev.off()

#Stacked Violin plot
features <- c("REEP1", "FLT1", "SPON1", "ADGRD1",
              "SYNPO2", "TPT1", "BLNK", "SYK","PFKFB3")

jpeg('Results/2-6-23_Poster_Figures/Pig_Stack_VlnMarker.jpeg',width = 2000,
     height = 3000, quality = 100, pointsize = 10)
VlnPlot(so, features, stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(text = element_text(size = 50),
        legend.position = "none",
        axis.text.y.right = element_text(size = 50),
        axis.text.x = element_text(size = 50, angle = 90, vjust = 0.5, hjust=1),
        axis.title = element_text(size = 70))
dev.off()



# Checkpoint ####
#Sub-clustering endothelial cells
so = readRDS('output/pig_v3.rds')
So_endo = subset(so, idents = c('Endothelial cells'))
So_endo <- ScaleData(object = So_endo, features = rownames(So_endo))

emb = Embeddings(So_endo,reduction = "umap")
kmeans.so = kmeans(emb, 5)

So_endo@meta.data['kmeans'] = kmeans.so$cluster

jpeg('Results/Other_Figures/Pig_Endothelial_kmeans.jpeg',width = 2000, height = 2000, quality = 100,pointsize = 10)
DimPlot(So_endo, reduction = "umap", label.size = 5, pt.size = 5, group.by = 'kmeans')+
  theme(axis.title = element_text(size = 40),
        plot.title = element_blank(),
        legend.text = element_text(size = 60),
        legend.key.size = unit(3, 'cm'),
        legend.position=c(.85,.85),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())+
  guides(color = guide_legend(override.aes=list(shape = 15,size = 15)))
dev.off()

Idents(object = So_endo) <- "kmeans"
so.markers <- FindAllMarkers(So_endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so.markers,'Endothelial_Subcluster_Markers.csv')

so.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

jpeg('Results/Other_Figures/Pig_Endothelial_Heatmap.jpeg',
     width = 2000, height = 800, quality = 100,pointsize = 1)
DoHeatmap(So_endo, features = top5$gene,size = 10)+
  theme(text = element_text(size = 30),
        legend.key.size = unit(1.5, 'cm'),
        group.bar = NULL)+
  guides(color = FALSE)
dev.off()

#Sub-clustering Mesenchymal cells
so = readRDS("output/pig_v3.rds")
so_foranno = subset(so, subset = cell_types == 'Mesenchymal cells')
Batches = so_foranno$Batches
source("../Codes/scType_anno.R")
so_foranno = Annotate_cells(so_foranno@assays$RNA@counts,
                            res = 1, tissue = 'Heart', Scale = T)
so_foranno$Batches = factor(Batches)
Idents(object = so_foranno) <- "scType_anno"
DimPlot(so_foranno, reduction = "umap", pt.size = 1, label = F, repel = TRUE, group.by = 'scType_anno')

levels(so_foranno)
new_names = c("Erythroblasts","Stromal cells","Epicardial fat cells","Epicardial fat cells",
              "Endocardial cells")
names(new_names) <- levels(so_foranno)
so_renamed <- RenameIdents(so_foranno, new_names)
so_renamed$new_names = so_renamed@active.ident
#so_renamed = subset(so_renamed, subset = new_names == c("Erythroblasts","Stromal cells","Epicardial fat cells","Endocardial cells"))
DimPlot(so_renamed, reduction = "umap", pt.size = 1, label = F, repel = TRUE)

png('Results/10-10-23_Publication_Figures/Mesenchymal Bar PLot.png', height = 4, width = 8,
    units = 'in', res = 300)
g1 = dittoBarPlot(so_renamed,var = "Batches", group.by = "new_names",scale = 'count')+
  scale_fill_manual("Batches",values=hue_pal()(4)) + ggtitle("Counts")+ xlab("Cell Types")
g2 = dittoBarPlot(so_renamed,var = "Batches", group.by = "new_names")+
  scale_fill_manual("Batches",values=hue_pal()(4),drop= FALSE) + ggtitle("Percentages") + xlab("Cell Types")
grid.arrange(g1, g2, ncol=2)
dev.off()

png('Results/10-10-23_Publication_Figures/Mesenchymal UMAP.png', height = 4, width = 7,
    units = 'in', res = 300)
DimPlot(so_renamed, group.by = 'new_names', reduction = "umap", pt.size = 1, label = F, repel = TRUE,)+
  ggtitle("")
dev.off()

