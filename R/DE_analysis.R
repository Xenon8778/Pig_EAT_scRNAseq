library(edgeR)
library(Seurat)
library(dplyr)
library(ggplot2)
library(readxl)
library(ComplexHeatmap)
library(colorRamp2)

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
  name = paste("EdgeR/EdgeR_",y,"_new.csv",sep="")
  write.csv(et,name)
  topTags(et, n=10)
}

so = readRDS("output/pig_v3.rds")

# DE Between Exercised and Sedentary ####
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

y = 'Beige adipocytes'
y = 'T cells'
y = 'B cells'
y = 'Endothelial cells'
y = 'Smooth muscle cells'
y = 'Macrophages'
y = 'Fibroblasts'
y = 'Erythroid cells'
y = 'Mesenchymal cells'

# Dotplots 
name = paste("EdgeR/EdgeR_",y,"_new.csv",sep="")
et = read.csv(name)
et = et[order(et$PValue),]
et = et[1:10,]
soforDE_dot = subset(so,idents = c(y))
n1 = paste("Results/2-6-23_Poster_Figures/EdgeR/EdgeR_Dotplot_",y,".jpeg",sep="")
jpeg(n1,width = 1500, height = 2000, quality = 100,pointsize = 10)
DotPlot(object = soforDE_dot, features = unique(et$X),group.by = "Exercised",
        dot.scale = 50,cols = c("lightgrey", "#FF68A1"))+
  ggtitle(paste(y))+
  theme(axis.title = element_text(size = 50, face="bold"),
        plot.title = element_text(size = 70, face="bold", hjust = 0.5),
        axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 40),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.key.size = unit(2, 'cm'))+
  coord_flip()
dev.off()

y = 'Beige adipocytes'
y = 'B cells'
y = 'Endothelial cells'
y = 'Smooth muscle cells'
y = 'Macrophages'


# DE Between Occluded and Non-occluded ####

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
  et <- exactTest(d1, pair=c("N-Ex","O-Ex"))
  name = paste("EdgeR/Ex/EdgeR_",y,"_new_n-o.csv",sep="")
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

x=so
y='Fibroblasts'

y = 'Beige adipocytes'
y = 'T cells'
y = 'B cells'
y = 'Endothelial cells'
y = 'Smooth muscle cells'
y = 'Macrophages'
y = 'Fibroblasts'
y = 'Erythroid cells'
y = 'Mesenchymal cells'

# Dotplots 
name = paste("EdgeR/Ex/EdgeR_",y,"_new_n-o.csv",sep="")
et = read.csv(name)
et = et[order(et$PValue),]
et = et[1:10,]
soforDE_dot = subset(so,idents = c(y), subset = Exercised == "Ex")
n1 = paste("Results/Finalised_EdgeR/Ex_Dotplot_",y,".jpeg",sep="")
jpeg(n1,width = 1500, height = 2000, quality = 100,pointsize = 10)
DotPlot(object = soforDE_dot, features = unique(et$X),group.by = "Batches",
        dot.scale = 50,cols = c("lightgrey", "#FF68A1"))+
  ggtitle(paste(y))+
  theme(axis.title = element_text(size = 50, face="bold"),
        plot.title = element_text(size = 70, face="bold", hjust = 0.5),
        axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 40),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.key.size = unit(2, 'cm'))+
  coord_flip()
dev.off()

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
  et <- exactTest(d1, pair=c("N-Sed","O-Sed"))
  name = paste("EdgeR/Sed/EdgeR_",y,"_n-o.csv",sep="")
  write.csv(et,name)
  topTags(et, n=10)
}

y = 'Beige adipocytes'
y = 'T cells'
y = 'B cells'
y = 'Endothelial cells'
y = 'Smooth muscle cells'
y = 'Macrophages'
y = 'Fibroblasts'
y = 'Erythroid cells'
y = 'Mesenchymal cells'

# Dotplots 
name = paste("EdgeR/sed/EdgeR_",y,"_n-o.csv",sep="")
et = read.csv(name)
et = et[order(et$PValue),]
et = et[1:10,]
soforDE_dot = subset(so,idents = c(y), subset = Exercised == "Sed")
n1 = paste("Results/Finalised_EdgeR/Sed/Sed_Dotplot_",y,".jpeg",sep="")
jpeg(n1,width = 1500, height = 2000, quality = 100,pointsize = 10)
DotPlot(object = soforDE_dot, features = unique(et$X),group.by = "Batches",
        dot.scale = 50,cols = c("lightgrey", "#FF68A1"))+
  ggtitle(paste(y))+
  theme(axis.title = element_text(size = 50, face="bold"),
        plot.title = element_text(size = 70, face="bold", hjust = 0.5),
        axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 40),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.key.size = unit(2, 'cm'))+
  coord_flip()
dev.off()

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
  name = paste("EdgeR/NOcc/EdgeR_",y,"_NOcc.csv",sep="")
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

# Dotplots 
for (i in unique(so$cell_types)){
  name = paste0("EdgeR/NOcc/EdgeR_",i,"_NOcc.csv",sep="")
  et = read.csv(name)
  et = et[order(et$PValue),]
  et = et[1:10,]
  soforDE_dot = subset(so,idents = c(i), subset = Batches == c("N-Ex","N-Sed"))
  n1 = paste0("EdgeR/NOcc/Dotplot_",i,".jpeg",sep="")

  png(n1,width = 6, height = 7, res = 300, units = 'in')
  g= DotPlot(object = soforDE_dot, features = unique(et$X),group.by = "Batches",
             scale = F, cols = c("yellow", "red"), dot.scale = 10)+
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
  name = paste("EdgeR/Occ/EdgeR_",y,"_Occ.csv",sep="")
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

# Dotplots 
for (i in unique(so$cell_types)){
  if ( i != 'Mesenchymal cells'){
    print(i)
    name = paste0("EdgeR/Occ/EdgeR_",i,"_Occ.csv",sep="")
    et = read.csv(name)
    et = et[order(et$PValue),]
    et = et[1:10,]
    soforDE_dot = subset(so,idents = c(i), subset = Batches == c("O-Ex","O-Sed"))
    n1 = paste0("EdgeR/Occ/Dotplot_",i,".jpeg")

    png(n1,width = 6, height = 7, res = 300, units = 'in')
    g= DotPlot(object = soforDE_dot, features = unique(et$X),group.by = "Batches",
               scale = F, cols = c("yellow", "red"), dot.scale = 10)+
      ggtitle(paste(i))+
      theme(axis.title = element_text(face="bold"),
            plot.title = element_text(face="bold", hjust = 0.5))+
      coord_flip()
    print(g)
    dev.off()
  }
}

# Number of DEGs Heatmap
y = c('Beige adipocytes','T cells','B cells','Endothelial cells',
      'Smooth muscle cells','Macrophages','Fibroblasts','Erythroid cells',
      'Mesenchymal cells')
nrows = length(y)
ncols = 4
DEG_Heat = data.frame(matrix(ncol = ncols, nrow = nrows),row.names = y)
colnames(DEG_Heat) = c('Upregulated','Downregulated','Upregulated','Downregulated')

for(i in 1:length(y)){
  cnt = 1
  name = paste("EdgeR/NOcc/EdgeR_",y[i],"_NOcc.csv",sep="")
  if (file.exists(name)){
    et = read.csv(name)
    et = et[order(et$PValue),]
    et  = subset(et,PValue<0.01)
    upreg = nrow(subset(et,logFC>0))
    DEG_Heat[i,1] = upreg
    downreg = nrow(subset(et,logFC<0))
    DEG_Heat[i,2] = downreg
  }

  name = paste("EdgeR/Occ/EdgeR_",y[i],"_Occ.csv",sep="")
  if (file.exists(name)){
    et = read.csv(name)
    et = et[order(et$PValue),]
    et  = subset(et,PValue<0.01)
    upreg = nrow(subset(et,logFC>0))
    DEG_Heat[i,3] = upreg
    downreg = nrow(subset(et,logFC<0))
    DEG_Heat[i,4] = downreg
  }
}

# Heatmap
DEG_Heat = na.omit(DEG_Heat)
col_fun = colorRamp2(c(0, 5000), c("white", "#F68282"))
Heatmap(DEG_Heat, name = "DEGs",col = col_fun, cluster_columns = F, cluster_rows = F,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%d", DEG_Heat[i, j]), x, y, gp = gpar(fontsize = 15))
        },
        border = TRUE, column_dend_reorder = FALSE,
        column_split = c("Non-Occluded","Non-Occluded","Occluded", "Occluded"), cluster_column_slices = F,
        right_annotation = rowAnnotation(nDEGs = anno_barplot(rowSums(DEG_Heat))))
