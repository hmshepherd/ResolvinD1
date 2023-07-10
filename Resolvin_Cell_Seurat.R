library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(harmony)
library(glmGamPoi)
library(sctransform)
library(tidyr)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(ggpubr)


#Load data
etohr.data = Read10X(data.dir = "./filtered_feature_bc_matrix")
etohd.data = Read10X(data.dir = "./filtered_feature_bc_matrix")
rvd1r.data = Read10X(data.dir = "./filtered_feature_bc_matrix")
rvd1d.data = Read10X(data.dir = "./filtered_feature_bc_matrix")

#Create Suerat Objects
etohr = CreateSeuratObject(counts = etohr.data)
etohd = CreateSeuratObject(counts = etohd.data)
rvd1r = CreateSeuratObject(counts = rvd1r.data)
rvd1d = CreateSeuratObject(counts = rvd1d.data)

#QC
etohr[["percent.mt"]] <- PercentageFeatureSet(etohr, pattern = "^mt-")
etohd[["percent.mt"]] <- PercentageFeatureSet(etohd, pattern = "^mt-")
rvd1r[["percent.mt"]] <- PercentageFeatureSet(rvd1r, pattern = "^mt-")
rvd1d[["percent.mt"]] <- PercentageFeatureSet(rvd1d, pattern = "^mt-")

VlnPlot(etohr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(etohd, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(rvd1r, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(rvd1d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

etohr=subset(etohr, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
etohd=subset(etohd, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
rvd1r=subset(rvd1r, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
rvd1d=subset(rvd1d, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)

#Add relevant metadata
etohr$sample <- "etohr"
etohd$sample <- "etohd"
rvd1r$sample <- "rvd1r"
rvd1d$sample <- "rvd1d"

#Merge
all.merged <- merge(etohr, y = c(etohd, rvd1r, rvd1d))

#Pipeline
sample <- SCTransform(all.merged, method = "glmGamPoi", return.only.var.genes=FALSE)
s.genes <- c("Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung", "Gins2", "Mcm6", "Cdca7", "Dtl", "Prim1", "Uurf1", "Cenpu", "Hells", "Rfc2", "Rpa2", "Nasp", "Rad51ap1", "Gmnn", "Wdr76", "Slbp", "Ccne2", "Ubr7", "Pold3", "Msh2", "Atad2", "Rad51", "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm", "Casp8ap2", "Usp1", "Clspn", "Pola1", "Chaf1b", "Brip1", "E2f8")
g2m.genes <- c("Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", "Ndc80", "Cks2", "Nuf2", "Cks1b", "Mki67", "Tmpo", "Cenpf", "Tacc3", "Pimreg", "Smc4", "Ccnb2", "Ckap2l", "Ckap2", "Aurkb", "Bub1", "Kif11", "Anp32e", "Tubb4b", "Gtse1", "Kif20b", "Hjurp", "Cdca3", "Jpt1", "Cdc20", "Ttk", "Cdc25c", "Kif2c", "Rangap1", "Ncapd2", "Dlgap5", "Cdca2", "Cdca8", "Ect2", "Kif23", "Hmmr", "Aurka", "Psrc1", "Anln", "Lbr", "Ckap5", "Cenpe", "Ctcf", "Nek2", "G2e3", "Gas2l3", "Cbx5", "Cenpa")
sample <- CellCycleScoring(sample, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
sample$CC.Difference <- sample$S.Score-sample$G2M.Score
DefaultAssay(sample) <- 'RNA'
sample <- SCTransform(sample, vars.to.regress = c("percent.mt", "CC.Difference"), method = "glmGamPoi", return.only.var.genes=FALSE)
sample <- RunPCA(sample, features = VariableFeatures(object = sample), npcs=100, verbose=TRUE)
sample <- RunHarmony(object = sample, assay.use = "SCT", reduction = "pca", dims.use = 1:50, group.by.vars = "sample", plot_convergence = TRUE)
sample <- RunUMAP(object = sample, assay = "SCT", reduction = "harmony", dims = 1:50)
sample <- FindNeighbors(object = sample, assay = "SCT", reduction = "harmony", dims = 1:50)
sample <- FindClusters(sample, graph.name = "SCT_snn", algorithm = 3, resolution = c(0.1,0.2,0.3,0.4,0.5), verbose = FALSE)

#Add other metadata
metadata <- sample@meta.data
metadata$split <- NA
metadata$split[which(str_detect(metadata$sample, "etohr"))] <- "etoh"
metadata$split[which(str_detect(metadata$sample, "etohd"))] <- "etoh"
metadata$split[which(str_detect(metadata$sample, "rvd1r"))] <- "rvd1"
metadata$split[which(str_detect(metadata$sample, "rvd1d"))] <- "rvd1"
metadata$type <- NA
metadata$type[which(str_detect(metadata$sample, "etohr"))] <- "recipient"
metadata$type[which(str_detect(metadata$sample, "etohd"))] <- "donor"
metadata$type[which(str_detect(metadata$sample, "rvd1r"))] <- "recipient"
metadata$type[which(str_detect(metadata$sample, "rvd1d"))] <- "donor"
sample@meta.data <- metadata

#Set Default
DefaultAssay(sample) <- 'SCT'
Idents(sample) <- "SCT_snn_res.0.5"

#Plots
Idents(rvd1)<-"SCT_snn_res.0.3"
rvd1<-subset(rvd1, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "15", "16", "17"))
fun <- function(x) {
 if (x == "0") {"Neutrophils"} 
 else if (x == "1") {"Monocytes/Macrophages"}
 else if (x == "2") {"NKT cells/T cells"}
 else if (x == "3") {"NK cells"}
 else if (x == "4") {"B cells"}
 else if (x == "5") {"NKT cells/T cells"}
 else if (x == "6") {"Monocytes/Macrophages"}
 else if (x == "7") {"Neutrophils"}
 else if (x == "8") {"Mast cells"}
 else if (x == "9") {"ILC"}
 else if (x == "10") {"Dendritic cells"}}
rvd1$celltype <- mapply(fun, rvd1$SCT_snn_res.0.1)
Idents(rvd1)<-"celltype"
DimPlot(rvd1, label = TRUE)
VlnPlot(rvd1, features = "Retnlg", pt.size = 1)+ theme(legend.position = 'none')
VlnPlot(rvd1, features = "Lyz2", pt.size = 1)+ theme(legend.position = 'none')
VlnPlot(rvd1, features = "Cd3e", pt.size = 1)+ theme(legend.position = 'none')
VlnPlot(rvd1, features = "Cd19", pt.size = 1)+ theme(legend.position = 'none')
VlnPlot(rvd1, features = "Klrb1c", pt.size = 1)+ theme(legend.position = 'none')
VlnPlot(rvd1, features = "Siglech", pt.size = 1)+ theme(legend.position = 'none')
VlnPlot(rvd1, features = "Gata3", pt.size = 1)+ theme(legend.position = 'none')
VlnPlot(rvd1, features = "Ms4a2", pt.size = 1)+ theme(legend.position = 'none')
fun <- function(x) {
 if (x == "0") {"Neutrophil"} 
 else if (x == "1") {"T Cell"}
 else if (x == "2") {"NK Cell"}
 else if (x == "3") {"B Cell"}
 else if (x == "4") {"Neutrophil"}
 else if (x == "5") {"CM"}
 else if (x == "6") {"NCM"}
 else if (x == "7") {"Macrophage"}
 else if (x == "8") {"Neutrophil"}
 else if (x == "9") {"T Cell"}
 else if (x == "10") {"T Cell"}
 else if (x == "11") {"Neutrophil"}
 else if (x == "12") {"AM"}
 else if (x == "13") {"Mast Cell"}
 else if (x == "15") {"ILC"}
 else if (x == "16") {"Neutrophil"}
 else if (x == "17") {"Dendritic Cell"}}
rvd1$cell.type <- mapply(fun, rvd1$SCT_snn_res.0.3)
VlnPlot(rvd1, features = "Fpr2", pt.size = 1, group.by = "cell.type", sort = TRUE)+ theme(legend.position = 'none')
celltype.markers <- FindAllMarkers(rvd1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
celltype.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top10
rvd1.small <- subset(rvd1, downsample = 300)
DoHeatmap(rvd1.small, features = top10$gene) + NoLegend()
expdata <- GetAssayData(rvd1)
Pop1 <- c("Retnlg", "Ly6g", "Fcgr3b", "Ifitm2")
pops<-list(Pop1)
z_scores<-NULL
for (i in 1:length(pops)) {
  genes <- pops[[i]]
  zz <- which(tolower(rownames(expdata)) %in% tolower(genes))
  av <- numeric(ncol(expdata))
  geneExp <- as.matrix(expdata[zz, ])
  geneExp <- t(scale(t(geneExp)))
  geneExp[is.nan(geneExp)] <- 0
  z_scores <- rbind(z_scores,(av + colSums(geneExp) / length(zz)))}
rvd1@meta.data$pop_z<-z_scores[1,]
FeaturePlot(object=rvd1, features = "pop_z",pt.size=.3, reduction = "umap") + scale_color_gradientn(colors=c("blue","turquoise2","yellow","red","red4"), oob=scales::squish, limits=c(0,1.0))


#Monocyte/Macrophage Subset
Idents(sample)<-"SCT_snn_res.0.5"
mono.mac.rvd1<-subset(sample, idents = c("5", "6","8", "13", "18"))
DefaultAssay(mono.mac.rvd1) <- "RNA"
mono.mac.rvd1 <- NormalizeData(mono.mac.rvd1)
all.genes <- rownames(mono.mac.rvd1)
mono.mac.rvd1 <- ScaleData(mono.mac.rvd1, features = all.genes)
mono.mac.rvd1 <- FindVariableFeatures(mono.mac.rvd1, selection.method = "vst", nfeatures = 2000)
mono.mac.rvd1 <- RunPCA(mono.mac.rvd1, features = VariableFeatures(object = mono.mac.rvd1))
mono.mac.rvd1 <- FindNeighbors(mono.mac.rvd1, dims = 1:20)
mono.mac.rvd1 <- FindClusters(mono.mac.rvd1, resolution = c(0.2, 0.4, 0.5, 0.8, 1.0))
mono.mac.rvd1 <- RunUMAP(mono.mac.rvd1, dims = 1:20)
DefaultAssay(mono.mac.rvd1) <- "RNA"
mono.mac.rvd1 <- RunHarmony(object = mono.mac.rvd1, assay.use = "RNA", reduction = "pca", dims.use = 1:20, group.by.vars = "sample", plot_convergence = TRUE)
mono.mac.rvd1 <- RunUMAP(object = mono.mac.rvd1, assay = "RNA", reduction = "harmony", dims = 1:20)
mono.mac.rvd1 <- FindNeighbors(object = mono.mac.rvd1, assay = "RNA", reduction = "harmony", dims = 1:20)
mono <- FindClusters(mono.mac.rvd1, graph.name = "RNA_snn", algorithm = 3, resolution = c(0.1,0.2,0.3,0.4,0.5), verbose = FALSE)

#Myeloid subset plots
Idents(mono)<-"RNA_snn_res.0.5"
mono <- subset(mono, idents = c("0", "1", "2", "4", "5", "6", "8", "9", "10", "11"))
fun <- function(x) {
if (x == "0") {"Classical Monocytes"} 
else if (x == "1") {"Non-classical Monocytes"}
else if (x == "2") {"Macrophage 3"}
else if (x == "4") {"Macrophage 2"}
else if (x == "5") {"Macrophage 1"}
else if (x == "6") {"Alveolar Macrophage"}
else if (x == "8") {"cDC1"}
else if (x == "9") {"pDC"}
else if (x == "10") {"cDC2"}
else if (x == "11") {"Non-classical Monocytes"}}
mono$celltype <- mapply(fun, mono$RNA_snn_res.0.5)
DimPlot(mono, group.by = "celltype", label = TRUE)
my_comparisons <- list(c("rvd1", "etoh"))
c<- VlnPlot(mono, features = "Cxcl2", group.by = "split", y.max = 9, pt.size = 0)+stat_compare_means(comparisons = my_comparisons, label = "p.format",  label.y = c(7))
c+scale_fill_manual(values = c('#F8766D', '#00BFC4'))
c<- VlnPlot(mono, features = "Tnf", group.by = "split", y.max = 9, pt.size = 0)+stat_compare_means(comparisons = my_comparisons, label = "p.format",  label.y = c(7))
c+scale_fill_manual(values = c('#F8766D', '#00BFC4'))

