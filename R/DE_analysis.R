library(edgeR)
library(Seurat)
library(dplyr)
library(ggplot2)
library(readxl)

runedger = function(x,y){
  #x = Seurat Object
  #y = Cell Type
  soforDGE = subset(x,idents = c(y))
  dge = DGEList(counts=soforDGE@assays$RNA@counts,
                group = factor(soforDGE@meta.data$Exercised))

  ## Filtering the data
  dge.full <- dge
  apply(dge$counts, 2, sum)
  keep <- rowSums(cpm(dge)>100) >= 2
  dge <- dge[keep,]
  dim(dge)
  dge$samples$lib.size <- colSums(dge$counts)
  dge$samples
  dge <- calcNormFactors(dge)
  d1 <- estimateCommonDisp(dge, verbose=T)
  d1 <- estimateTagwiseDisp(d1)
  et <- exactTest(d1, pair=c("Ex","Sed"))
  et$table = et$table %>% mutate(PValue_Adj = p.adjust (PValue, method='BH')) %>%
    arrange(PValue_Adj)
  name = paste("EdgeR/EdgeR_",y,"_new.csv",sep="")
  write.csv(et,name)
  topTags(et, n=10)
}

so = readRDS("output/pig_v3.rds")

# Between Exercised and Sedentary ####
y='Fibroblasts'
runedger(so,'T cells')
runedger(so,'Beige adipocytes')
runedger(so,'B cells')
runedger(so,'Endothelial cells')
runedger(so,'Smooth muscle cells')
runedger(so,'Macrophages')
runedger(so,'Fibroblasts')
runedger(so,'Erythroid cells')
runedger(so,'Mesenchymal cells')


for (i in unique(so$cell_types)){
  print(i)
  name = paste0("EdgeR/DE_Ex_Sed/EdgeR_",i,"_new.csv",sep="")
  et = read.csv(name)
  et = et %>% arrange(PValue_Adj) %>% filter( PValue_Adj < 0.05)
  et = et[1:10,]
  soforDE_dot = so
  n1 = paste0("EdgeR/DE_Ex_Sed/EdgeR_Dotplot_",i,".png")

  png(n1,width = 6, height = 7, res = 300, units = 'in')
  g = DotPlot(object = soforDE_dot, features = unique(et$X),group.by = "Exercised",
              scale = T, cols = c("lightgrey", "#FF68A1"), dot.scale = 10)+
    ggtitle(paste(i))+
    theme(axis.title = element_text(face="bold"),
          plot.title = element_text(face="bold", hjust = 0.5))+
    coord_flip()
  plot(g)
  dev.off()
}


#Marker genes
soformarkers = subset(so,idents = c("Erythroid cells","Mesenchymal cells","Fibroblasts"), invert=TRUE)
so.markers <- FindAllMarkers(soformarkers, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

marker_genes = read_xls("PigCellTypeMarkers_11232022.xls")
labels = read_xls("PigCellTypeMarkers_11232022.xls")
labels = labels$Annotations
labels = labels[!is.na(labels)]
marker_genes$Cluster = marker_genes$Cluster+1
cell_types = c()
for (i in marker_genes$Cluster){
  cell_types = append(cell_types,labels[i])
}
marker_genes$Cluster = cell_types
marker_genes = subset(marker_genes,marker_genes$Cluster != "Unknown")
marker_genes = marker_genes %>%
  group_by(Cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

jpeg('output/Pig_Dotplot_DEGs.jpeg',width = 800, height = 1200, quality = 100,pointsize = 100)
DotPlot(object = soformarkers, features = unique(marker_genes$M.Genesnew))+
  ggtitle('nFeature Genes vs Cell Type')+
  theme(axis.title = element_text(size = 25),
        plot.title = element_text(color="red", size=30),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'))+
  coord_flip()
dev.off()

write.csv(marker_genes,"PigMarkers.csv")

Idents(soformarkers)
levels(soformarkers)
as.character(levels(soformarkers))




# Between Occluded and Non-occluded ####

## Exercised
runedger = function(x,y){
  #x = Seurat Object
  #y = Cell Type
  soforDGE = subset(x,idents = c(y), subset = Exercised == "Ex")
  dge = DGEList(counts=soforDGE@assays$RNA@counts,
                group = factor(soforDGE@meta.data$Batches))

  ## Filtering the data
  dge.full <- dge
  apply(dge$counts, 2, sum)
  keep <- rowSums(cpm(dge)>100) >= 2
  dge <- dge[keep,]
  dim(dge)
  dge$samples$lib.size <- colSums(dge$counts)
  dge$samples
  dge <- calcNormFactors(dge)
  d1 <- estimateCommonDisp(dge, verbose=T)
  d1 <- estimateTagwiseDisp(d1)
  et <- exactTest(d1, pair=c("O-Ex","N-Ex"))
  et$table = et$table %>% mutate(PValue_Adj = p.adjust (PValue, method='BH')) %>%
    arrange(PValue_Adj)
  name = paste("EdgeR/Ex/Revised/EdgeR_",y,"_Ex_O-N.csv",sep="")
  write.csv(et,name)
  topTags(et, n=10)
}

runedger(so,'T cells')
runedger(so,'Beige adipocytes')
runedger(so,'B cells')
runedger(so,'Endothelial cells')
runedger(so,'Smooth muscle cells')
runedger(so,'Macrophages')
runedger(so,'Fibroblasts')
runedger(so,'Erythroid cells')
runedger(so,'Mesenchymal cells')


i = 'Mesenchymal cells'
for (i in unique(so$cell_types)){
  if ( i != 'Mesenchymal cells'){
    print(i)
    name = paste0("EdgeR/Ex/Revised/EdgeR_",i,"_Ex_O-N.csv",sep="")
    et = read.csv(name)
    et = et %>% arrange(PValue_Adj) %>% filter( PValue_Adj < 0.05)
    et = et[1:10,]
    soforDE_dot = subset(so,idents = c(i), subset = Batches == c("O-Ex","N-Ex"))
    n1 = paste0("EdgeR/Ex/Revised/Dotplot_",i,".png")

    png(n1,width = 6, height = 7, res = 300, units = 'in')
    g = DotPlot(object = soforDE_dot, features = unique(et$X),group.by = "Batches",
                scale = T, cols = c("lightgrey", "#FF68A1"), dot.scale = 10)+
      ggtitle(paste(i))+
      theme(axis.title = element_text(face="bold"),
            plot.title = element_text(face="bold", hjust = 0.5))+
      coord_flip()
    plot(g)
    dev.off()
  }
}

## Sedentary

so = readRDS("output/pig_v3.rds")

runedger = function(x,y){
  #x = Seurat Object
  #y = Cell Type
  soforDGE = subset(x,idents = c(y), subset = Exercised == "Sed")
  dge = DGEList(counts=soforDGE@assays$RNA@counts,
                group = factor(soforDGE@meta.data$Batches))

  ## Filtering the data
  dge.full <- dge
  apply(dge$counts, 2, sum)
  keep <- rowSums(cpm(dge)>100) >= 2
  dge <- dge[keep,]
  dim(dge)
  dge$samples$lib.size <- colSums(dge$counts)
  dge$samples
  dge <- calcNormFactors(dge)
  d1 <- estimateCommonDisp(dge, verbose=T)
  d1 <- estimateTagwiseDisp(d1)
  et <- exactTest(d1, pair=c("O-Sed","N-Sed"))
  et$table = et$table %>% mutate(PValue_Adj = p.adjust (PValue, method='BH')) %>%
    arrange(PValue_Adj)
  name = paste("EdgeR/Sed/Revised/EdgeR_",y,"_Sed_O-N.csv",sep="")
  write.csv(et,name)
  topTags(et, n=10)
}


for (i in unique(so$cell_types)){
  if ( i != 'Fibroblasts'){
    print(i)
    name = paste0("EdgeR/Sed/Revised/EdgeR_",i,"_Sed_O-N.csv",sep="")
    et = read.csv(name)
    et = et %>% arrange(PValue_Adj) %>% filter( PValue_Adj < 0.05)
    et = et[1:10,]
    soforDE_dot = subset(so,idents = c(i), subset = Batches == c("O-Sed","N-Sed"))
    n1 = paste0("EdgeR/Sed/Revised/Dotplot_",i,".png")

    png(n1,width = 6, height = 7, res = 300, units = 'in')
    g = DotPlot(object = soforDE_dot, features = unique(et$X),group.by = "Batches",
                scale = T, cols = c("lightgrey", "#FF68A1"), dot.scale = 10)+
      ggtitle(paste(i))+
      theme(axis.title = element_text(face="bold"),
            plot.title = element_text(face="bold", hjust = 0.5))+
      coord_flip()
    plot(g)
    dev.off()
  }
}

# Between Exercised and Sedentary ####

## Non-Occ

runedger = function(x,y){
  #x = Seurat Object
  #y = Cell Type
  soforDGE = subset(x,idents = c(y), subset = Batches == c("N-Ex","N-Sed"))
  dge = DGEList(counts=soforDGE@assays$RNA@counts,
                group = factor(soforDGE@meta.data$Batches))

  ## Filtering the data
  dge.full <- dge
  apply(dge$counts, 2, sum)
  keep <- rowSums(cpm(dge)>100) >= 2
  dge <- dge[keep,]
  dim(dge)
  dge$samples$lib.size <- colSums(dge$counts)
  dge$samples
  dge <- calcNormFactors(dge)
  d1 <- estimateCommonDisp(dge, verbose=T)
  d1 <- estimateTagwiseDisp(d1)
  et <- exactTest(d1, pair=c("N-Ex","N-Sed"))
  et$table = et$table %>% mutate(PValue_Adj = p.adjust (PValue, method='BH')) %>%
    arrange(PValue_Adj)
  name = paste("EdgeR/NOcc/Revised/EdgeR_",y,"_NOcc_E-S.csv",sep="")
  write.csv(et,name)
  topTags(et, n=10)
}

runedger(so,'T cells')
runedger(so,'Beige adipocytes')
runedger(so,'B cells')
runedger(so,'Endothelial cells')
runedger(so,'Smooth muscle cells')
runedger(so,'Macrophages')
runedger(so,'Fibroblasts')
runedger(so,'Erythroid cells')
runedger(so,'Mesenchymal cells')


for (i in unique(so$cell_types)){
  name = paste0("EdgeR/NOcc/Revised/EdgeR_",i,"_NOcc_E-S.csv",sep="")
  et = read.csv(name)
  et = et %>% arrange(PValue_Adj) %>% filter( PValue_Adj < 0.05)
  et = et[1:10,]
  soforDE_dot = subset(so,idents = c(i), subset = Batches == c("N-Ex","N-Sed"))
  n1 = paste0("EdgeR/NOcc/Revised/Dotplot_",i,".jpeg",sep="")

  png(n1,width = 6, height = 7, res = 300, units = 'in')
  g= DotPlot(object = soforDE_dot, features = unique(et$X),group.by = "Batches",
             scale = F, cols = c("lightgrey", "#FF68A1"), dot.scale = 10)+
    ggtitle(paste(i))+
    theme(axis.title = element_text(face="bold"),
          plot.title = element_text(face="bold", hjust = 0.5))+
    coord_flip()
  print(g)
  dev.off()
}


## Occ

so = readRDS("output/pig_v3.rds")

runedger = function(x,y){
  #x = Seurat Object
  #y = Cell Type
  soforDGE = subset(x,idents = c(y), subset = Batches == c("O-Ex","O-Sed"))
  dge = DGEList(counts=soforDGE@assays$RNA@counts,
                group = factor(soforDGE@meta.data$Batches))

  ## Filtering the data
  dge.full <- dge
  apply(dge$counts, 2, sum)
  keep <- rowSums(cpm(dge)>100) >= 2
  dge <- dge[keep,]
  dim(dge)
  dge$samples$lib.size <- colSums(dge$counts)
  dge$samples
  dge <- calcNormFactors(dge)
  d1 <- estimateCommonDisp(dge, verbose=T)
  d1 <- estimateTagwiseDisp(d1)
  et <- exactTest(d1, pair=c("O-Ex","O-Sed"))
  et$table = et$table %>% mutate(PValue_Adj = p.adjust (PValue, method='BH')) %>%
    arrange(PValue_Adj)
  name = paste("EdgeR/Occ/Revised/EdgeR_",y,"_Occ_E-S.csv",sep="")
  write.csv(et,name)
  topTags(et, n=10)
}

runedger(so,'T cells')
runedger(so,'Beige adipocytes')
runedger(so,'B cells')
runedger(so,'Endothelial cells')
runedger(so,'Smooth muscle cells')
runedger(so,'Macrophages')
runedger(so,'Fibroblasts')
runedger(so,'Erythroid cells')
runedger(so,'Mesenchymal cells')

for (i in unique(so$cell_types)){
  name = paste0("EdgeR/Occ/Revised/EdgeR_",i,"_Occ_E-S.csv",sep="")
  et = read.csv(name)
  et = et %>% arrange(PValue_Adj) %>% filter( PValue_Adj < 0.05)
  et = et[1:10,]
  soforDE_dot = subset(so,idents = c(i), subset = Batches == c("O-Ex","O-Sed"))
  n1 = paste0("EdgeR/Occ/Revised/Dotplot_",i,".jpeg",sep="")

  png(n1,width = 6, height = 7, res = 300, units = 'in')
  g= DotPlot(object = soforDE_dot, features = unique(et$X),group.by = "Batches",
             scale = F, cols = c("lightgrey", "#FF68A1"), dot.scale = 10)+
    ggtitle(paste(i))+
    theme(axis.title = element_text(face="bold"),
          plot.title = element_text(face="bold", hjust = 0.5))+
    coord_flip()
  print(g)
  dev.off()
}

# Number of DEGs plot ####

# Occ and Non-Occ
y = c('Beige adipocytes','T cells','B cells','Endothelial cells',
      'Smooth muscle cells','Macrophages','Fibroblasts','Erythroid cells',
      'Mesenchymal cells')
nrows = length(y)
ncols = 4
DEG_Heat = data.frame(matrix(ncol = ncols, nrow = nrows),row.names = y)
colnames(DEG_Heat) = c('Upregulated','Downregulated','Upregulated','Downregulated')

for(i in 1:length(y)){
  cnt = 1
  name = paste("EdgeR/NOcc/Revised/EdgeR_",y[i],"_NOcc_E-S.csv",sep="")
  if (file.exists(name)){
    et = read.csv(name)
    et = et %>% arrange(PValue_Adj) %>%
      filter(PValue_Adj<0.05) %>%
      filter(abs(logFC)>1)
    upreg = nrow(subset(et,logFC>0))
    DEG_Heat[i,1] = upreg
    downreg = nrow(subset(et,logFC<0))
    DEG_Heat[i,2] = downreg
  }

  name = paste("EdgeR/Occ/Revised/EdgeR_",y[i],"_Occ_E-S.csv",sep="")
  if (file.exists(name)){
    et = read.csv(name)
    et = et %>% arrange(PValue_Adj) %>%
      filter(PValue_Adj<0.05) %>%
      filter(abs(logFC)>1)
    upreg = nrow(subset(et,logFC>0))
    DEG_Heat[i,3] = upreg
    downreg = nrow(subset(et,logFC<0))
    DEG_Heat[i,4] = downreg
  }
}

DEG_Heat = na.omit(DEG_Heat)

library(ComplexHeatmap)
library(colorRamp2)

col_fun = colorRamp2(c(0, max(DEG_Heat)), c("white", "#F68282"))

Heatmap(DEG_Heat, name = "DEGs",col = col_fun, cluster_columns = F, cluster_rows = F,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%d", DEG_Heat[i, j]), x, y, gp = gpar(fontsize = 15))
        },
        border = TRUE, column_dend_reorder = FALSE,
        column_split = c("Non-Occluded","Non-Occluded","Occluded", "Occluded"), cluster_column_slices = F,
        right_annotation = rowAnnotation(nDEGs = anno_barplot(rowSums(DEG_Heat))))

# Ex and Sed
y = c('Beige adipocytes','T cells','B cells','Endothelial cells',
      'Smooth muscle cells','Macrophages','Fibroblasts','Erythroid cells',
      'Mesenchymal cells')
nrows = length(y)
ncols = 4
DEG_Heat = data.frame(matrix(ncol = ncols, nrow = nrows),row.names = y)
colnames(DEG_Heat) = c('Upregulated','Downregulated','Upregulated','Downregulated')

for(i in 1:length(y)){
  cnt = 1
  name = paste("EdgeR/Ex/Revised/EdgeR_",y[i],"_Ex_O-N.csv",sep="")
  if (file.exists(name)){
    et = read.csv(name)
    et = et %>% arrange(PValue_Adj) %>%
      filter(PValue_Adj<0.05) %>%
      filter(abs(logFC)>1)
    upreg = nrow(subset(et,logFC>0))
    DEG_Heat[i,1] = upreg
    downreg = nrow(subset(et,logFC<0))
    DEG_Heat[i,2] = downreg
  }

  name = paste("EdgeR/Sed/Revised/EdgeR_",y[i],"_Sed_O-N.csv",sep="")
  if (file.exists(name)){
    et = read.csv(name)
    et = et %>% arrange(PValue_Adj) %>%
      filter(PValue_Adj<0.05) %>%
      filter(abs(logFC)>1)
    upreg = nrow(subset(et,logFC>0))
    DEG_Heat[i,3] = upreg
    downreg = nrow(subset(et,logFC<0))
    DEG_Heat[i,4] = downreg
  }
}

DEG_Heat = na.omit(DEG_Heat)

col_fun = colorRamp2(c(0, max(DEG_Heat)), c("white", "#F68282"))

Heatmap(DEG_Heat, name = "DEGs",col = col_fun, cluster_columns = F, cluster_rows = F,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%d", DEG_Heat[i, j]), x, y, gp = gpar(fontsize = 15))
        },
        border = TRUE, column_dend_reorder = FALSE,
        column_split = c("Exercised","Exercised","Sedentary", "Sedentary"), cluster_column_slices = F,
        right_annotation = rowAnnotation(nDEGs = anno_barplot(rowSums(DEG_Heat))))
