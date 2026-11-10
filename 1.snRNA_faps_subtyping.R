library(Seurat)
library(harmony)
library(SeuratDisk)
library(dplyr)

################################################################
#Prepare separate objects for control and duchenne samples 
################################################################

#1. CONTROL SAMPLES

ctrl_obj <- list(Control1, Control2, Control3, Control4, Control5)

ctrl.features <- SelectIntegrationFeatures(object.list = ctrl_obj, nfeatures = 3000)

control_obj <- PrepSCTIntegration(object.list = ctrl_obj, anchor.features = ctrl.features, 
                                  verbose = FALSE)

#identify anchors and integrate the datasets
ctrl_anchors <- FindIntegrationAnchors(object.list = control_obj, normalization.method = "SCT", 
                                       anchor.features = ctrl.features, verbose = FALSE)

ctrl.integrated <- IntegrateData(anchorset = ctrl_anchors, normalization.method = "SCT", 
                                 verbose = FALSE)

set.seed(1234)

ctrl.integrated <- RunPCA(ctrl.integrated, verbose = FALSE)
ctrl.integrated <- RunUMAP(ctrl.integrated, dims = 1:20, reduction = "pca")
ctrl.integrated <- FindNeighbors(ctrl.integrated, reduction = "pca", dims = 1:20)
ctrl.integrated <- FindClusters(ctrl.integrated, resolution = c(0.3, 0.5, 0.7,1.0))

Idents(ctrl.integrated) <- "integrated_snn_res.0.5"
DimPlot(ctrl.integrated, reduction = "umap", label = T)
DimPlot(ctrl.integrated, reduction = "umap",split.by= "orig.ident", label = T)

#first normalize the RNA assay
DefaultAssay(ctrl.integrated) <- "RNA"
ctrl.integrated <- NormalizeData(ctrl.integrated, assay = "RNA") %>% ScaleData()

DefaultAssay(ctrl.integrated)
ctrl_markers <- FindAllMarkers(ctrl.integrated, only.pos = T)
write.csv(ctrl_markers, file = "ctrl_markers.csv")

##annotating the clusters based on the marker gene expression
ctrl.integrated <- RenameIdents(ctrl.integrated,
                                "0" = "Slow Fibres",
                                "1" = "Fast Fibres",
                                "2" = "Fast Fibres",
                                "3" = "Slow Fibres",
                                "4" = "FAPs",
                                "5" = "Slow Fibres",
                                "6" = "Fast Fibres",
                                "7" = "Slow Fibres",
                                "8" = "Slow Fibres",
                                "9" = "Endothelial Cells",
                                "10" = "Smooth Muscle Cells",
                                "11"= "Satellite Cells",
                                "12"= "Inflammatory Cells",
                                "13"= "Smooth Muscle Cells",
                                "14"= "Slow Fibres"
)

#2. DMD SAMPLES

dmd_obj <- list(DMD1, DMD2, DMD3, DMD4, DMD5, DMD6, DMD7)

dmd.features <- SelectIntegrationFeatures(object.list = dmd_obj, nfeatures = 3000)

dmd_obj <- PrepSCTIntegration(object.list = dmd_obj, anchor.features = dmd.features, 
                              verbose = FALSE)

#identify anchors and integrate the datasets
dmd_anchors <- FindIntegrationAnchors(object.list = dmd_obj, normalization.method = "SCT", 
                                      anchor.features = dmd.features, verbose = FALSE)

dmd.integrated <- IntegrateData(anchorset = dmd_anchors, normalization.method = "SCT", 
                                verbose = FALSE)

set.seed(123456)

dmd.integrated <- RunPCA(dmd.integrated, verbose = FALSE)
dmd.integrated <- RunUMAP(dmd.integrated, dims = 1:20, reduction = "pca")
dmd.integrated <- FindNeighbors(dmd.integrated, reduction = "pca", dims = 1:20)
dmd.integrated <- FindClusters(dmd.integrated, resolution = c(0.3, 0.5, 0.7,1.0))

Idents(dmd.integrated) <- "integrated_snn_res.0.5"
DimPlot(dmd.integrated, reduction = "umap", label = T)
DimPlot(dmd.integrated, reduction = "umap",group.by= "orig.ident", label = T)

#first normalize the RNA assay
dmd.integrated <- NormalizeData(dmd.integrated, assay = "RNA") %>% ScaleData()

DefaultAssay(dmd.integrated)
dmd_markers <- FindAllMarkers(dmd.integrated, only.pos = T)
write.csv(dmd_markers, file = "dmd_markers.csv")

dmd.integrated <- RenameIdents(dmd.integrated,
                               "0" = "FAPs",
                               "1" = "Slow Fibres",
                               "4" = "Slow Fibres",
                               "14"= "Slow Fibres",
                               "2" = "Fast Fibres",
                               "3" = "Fast Fibres",
                               "5" = "Regenerative Fibres",
                               "6" = "Fast Fibres",
                               "7" = "Inflammatory Cells",
                               "8" = "Mixed Fibres",
                               "9" = "Mixed Fibres",
                               "10"= "Satellite Cells",
                               "11"= "Smooth Muscle Cells",
                               "12"= "Endothelial Cells",
                               "13"= "Adipocytes",
                               "15"= "Regenerative Fibres",
                               "16" = "FAPs"
)

################################################################
##Subset FAPs from control and DMD objects separately
################################################################

#1. Control
DefaultAssay(ctrl.integrated) <- "RNA"

#subset FAPS
ctrl_faps <- subset(ctrl.integrated, idents = "FAPs")
DefaultAssay(ctrl_faps) <- "RNA"

#recluster
ctrl_faps_list <- SplitObject(ctrl_faps, split.by = "orig.ident")

#Merge normalized samples
merged_seurat <- merge(x = ctrl_faps_list[[1]],
                       y = ctrl_faps_list[2:length(ctrl_faps_list)],
                       merge.data = TRUE)

ctrl_faps_final <- merged_seurat %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA( npcs = 20)

print(ctrl_faps_final@reductions$pca)


#run harmony 
ctrl_faps_final <- RunHarmony(ctrl_faps_final, group.by.vars= c("orig.ident"),
                              reduction = "pca", reduction.save = "harmony")

#Check the new harmony object
ctrl_faps_final@reductions$harmony

#run umap
ctrl_faps_final <- ctrl_faps_final %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5)

DimPlot(ctrl_faps_final, group.by = "orig.ident")
DimPlot(ctrl_faps_final, label = T)

#cluster specific marker genes
fapsMark_ctrl <- FindAllMarkers(ctrl_faps_final, only.pos = T)

#2. DMD 

DefaultAssay(dmd.integrated) <- "RNA"

#subset FAPS
dmd_faps <- subset(dmd.integrated, idents = "FAPs")
DefaultAssay(dmd_faps) <- "RNA"

###new recluster workflow
dmd_faps_list <- SplitObject(dmd_faps, split.by = "orig.ident")


#Merge normalized samples
merged_seurat <- merge(x = dmd_faps_list[[1]],
                       y = dmd_faps_list[2:length(dmd_faps_list)],
                       merge.data = TRUE)

dmd_faps_final <- merged_seurat %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA( npcs = 20)

print(ctrl_faps_final@reductions$pca)


#run harmony 
dmd_faps_final <- RunHarmony(dmd_faps_final, group.by.vars= c("orig.ident"),
                             reduction = "pca", reduction.save = "harmony")

#Check the new harmony object
dmd_faps_final@reductions$harmony

#run umap
dmd_faps_final <- dmd_faps_final %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5)

DimPlot(dmd_faps_final, group.by = "orig.ident")
DimPlot(dmd_faps_final, label = T)

#cluster specific marker genes
fapsMark_dmd <- FindAllMarkers(dmd_faps_final, only.pos = T)


##Add chosen faps subtypes to the main DMD seurat object

dmd.integrated$new_celltype <- NA

#for the cells that were in the subset, assign the new sub-cluster labels
dmd.integrated$new_celltype[Cells(dmd_faps_final)] <- as.character(dmd_faps_final$cell_type)

dmd.integrated$merge_type <- ifelse(
  is.na(dmd.integrated$new_celltype),
  as.character(dmd.integrated$cell_type),                   
  as.character(dmd.integrated$new_celltype)      
)

################################################################
##Convert to anndata
################################################################

DefaultAssay(dmd.integrated) <- "RNA"

#to make the data compatible for cell2location deconvolution
dmd.integrated@assays$RNA@scale.data <- matrix(0, 0, 0)

#slim down a Seurat object
seu = DietSeurat(
  dmd.integrated,
  counts = TRUE, 
  data = FALSE, 
  scale.data = FALSE,
  features = rownames(dmd.integrated), 
  assays = "RNA",
  dimreducs = c("pca","umap"),
  graphs = c("RNA_nn", "RNA_snn"), 
  misc = TRUE
)

#all the columns in metadata are integers which might cause issue when converting to anndata, so run the below steps
i <- sapply(seu@meta.data, is.factor)
seu@meta.data[i] <- lapply(seu@meta.data[i], as.character)

#convert for dmd
SaveH5Seurat(seu, filename = "dmd_withfaps.h5Seurat", overwrite = T)
Convert("dmd_withfaps.h5Seurat", dest = "h5ad", assay = "RNA", overwrite = T)

