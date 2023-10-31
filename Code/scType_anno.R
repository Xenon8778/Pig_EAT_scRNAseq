## Load Libs
library(Seurat)
library(openxlsx)
library(Matrix)

Annotate_cells <- function(data = x, res = 0.8, tissue = NULL, Scale = F, annot_only = F,
                           filter_conf = T){
  #> data - is a count matrix
  #> res - resolution for graph based clustering
  #> tissue - Tissue type for which we are looking markers for
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
  # db_ = "../codes/data/ScTypeDB_full.xlsx"
  dataDB = read.xlsx(db_)

  if (is.null(tissue)){
    print(unique(dataDB$tissueType))
    tissue <- readline("Select Tissue Type for annotation?\n")
  }
  data = CreateSeuratObject(data)
  data <- NormalizeData(data)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  if (isTRUE(Scale)){data <- ScaleData(data, features = rownames(data))}
  data <- RunPCA(data, features = NULL,verbose = F)
  data <- FindNeighbors(data, dims = 1:10)
  data <- FindClusters(data, resolution = res)
  data <- RunUMAP(data, dims = 1:10)

  # load libraries and functions
  lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  # load gene set preparation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  # load cell type annotation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

  # prepare gene sets
  gs_list = gene_sets_prepare(db_,tissue)

  # get cell-type by cell matrix
  es.max = sctype_score(scRNAseqData = data[["RNA"]]@scale.data, scaled = TRUE,
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  # merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(data@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(data@meta.data[data@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(data@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

  # set low-confident (low ScType score) clusters to "unknown"
  if (filter_conf == T){
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
    }
  print(sctype_scores[,1:3])

  data@meta.data$scType_anno = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,];
    data@meta.data$scType_anno[data@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  if (annot_only == T){
    return(unlist(data$scType_anno))
  }
  else{return(data)}
}

#> library(ggplot2)
#> library(dittoSeq)
#> library(gridExtra)
#> DimPlot(so, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'scType_anno')
#> dittoBarPlot(so,var = "Batches", group.by = "scType_anno")+
#>  scale_fill_manual("Batches",values=c("#F68282", "#7CAE00", "#00BFC4","#C77CFF"))
