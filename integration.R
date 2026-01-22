getwd();library(Seurat);library(Matrix);suppressMessages(library(dplyr));library(ggplot2);suppressMessages(library('pagoda2'))#;suppressMessages(library("velocyto.R"));suppressMessages(library(monocle));
suppressMessages(library(cowplot));suppressMessages(library("patchwork"));suppressMessages(library(harmony))
options(repr.plot.width=7, repr.plot.height=6); library("RColorBrewer");library(magrittr);library(scales)
setwd("/sdc/xxjiang/rapl")
getwd()
packageVersion('Seurat')

pbmc1=readRDS("rapl114572_new.rds");pbmc1;names(pbmc1@meta.data)

pbmc2=readRDS("rapl148060n-pca.rds");pbmc2;names(pbmc2@meta.data)

pbmc1$cx=pbmc1$C25FX
pbmc2$cx=pbmc2$stage
pbmc1$type='SC'
pbmc2$type='SN'
pbmc1$Annotation=pbmc1$Name
pbmc2$Annotation=pbmc2$stage

pbmc <- merge(x =pbmc1, y = pbmc2);      pbmc

pbmc <- CreateSeuratObject(pbmc$RNA@counts, meta.data = pbmc@meta.data ,min.cells = 3, project = "trophoblast");pbmc

names(pbmc@meta.data);pbmc;table(pbmc$sample)

table(Idents(pbmc))

Idents(pbmc, cells = WhichCells(pbmc,idents = "GD26n-2")) <- "35"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD26n-1")) <- "34"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD26")) <- "33"

Idents(pbmc, cells = WhichCells(pbmc,idents = "GD23n-2")) <- "32"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD23n-1")) <- "31"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD23")) <- "30"

Idents(pbmc, cells = WhichCells(pbmc,idents = "GD20n-2")) <- "29"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD20n-1")) <- "28"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD20-2")) <- "27"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD20-1")) <- "26"

Idents(pbmc, cells = WhichCells(pbmc,idents = "GD17n-2")) <- "25"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD17n-1")) <- "24"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD17-2")) <- "23"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD17-1")) <- "22"

Idents(pbmc, cells = WhichCells(pbmc,idents = "GD15n-2")) <- "21"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD15n-1")) <- "20"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD15-2")) <- "19"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD15-1")) <- "18"

Idents(pbmc, cells = WhichCells(pbmc,idents = "GD13n-2")) <- "17"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD13n-1")) <- "16"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD13")) <- "15"

Idents(pbmc, cells = WhichCells(pbmc,idents = "GD11n")) <- "14"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD11-2")) <- "13"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD11-1")) <- "12"

Idents(pbmc, cells = WhichCells(pbmc,idents = "GD10n-2")) <- "11"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD10n-1")) <- "10"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD10-2")) <- "9"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD10-1")) <- "8"

Idents(pbmc, cells = WhichCells(pbmc,idents = "GD9n-2")) <- "7"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD9n-1")) <- "6"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD9-2")) <- "5"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD9-1")) <- "4"

Idents(pbmc, cells = WhichCells(pbmc,idents = "GD8n")) <- "3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD8-2")) <- "2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "GD8-1")) <- "1"

pbmc$orig2=Idents(pbmc)

Idents(pbmc)="orig2"
new.cluster.ids <- c("GD8-1","GD8-2", "GD8n",
                     "GD9-1", "GD9-2","GD9n-1","GD9n-2",
                     "GD10-1","GD10-2" ,"GD10n-1","GD10n-2",
                     "GD11-1","GD11-2","GD11n",
                     "GD13","GD13n-1","GD13n-2",
                     "GD15-1","GD15-2","GD15n-1", "GD15n-2",
                     "GD17-1","GD17-2","GD17n-1", "GD17n-2",
                     "GD20-1","GD20-2","GD20n-1", "GD20n-2" ,
    "GD23", "GD23n-1", "GD23n-2",
    "GD26","GD26n-1" , "GD26n-2"
                     )
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
pbmc$sample2=Idents(pbmc)
table(Idents(pbmc))

Idents(pbmc)="orig2"
new.cluster.ids <- c("GD8","GD8", "GD8",
                     "GD9", "GD9","GD9","GD9",
                     "GD10","GD10" ,"GD10","GD10",
                     "GD11","GD11","GD11",
                     "GD13","GD13","GD13",
                     "GD15","GD15","GD15", "GD15",
                     "GD17","GD17","GD17", "GD17",
                     "GD20","GD20","GD20", "GD20" ,
    "GD23","GD23","GD23",
    "GD26","GD26","GD26"  )
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
pbmc$stage2=Idents(pbmc)
table(Idents(pbmc))

options(repr.plot.width=12, repr.plot.height=8)
 VlnPlot(pbmc,pt.size=0, features = c("nFeature_RNA","percent.mt"),raster=FALSE,ncol=1) +  theme(axis.title=element_blank(),axis.text=element_text(size=16),plot.title=element_text(size=16))+NoLegend()

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(pbmc))

pancreas=pbmc
pancreas.list <- SplitObject(object = pancreas, split.by = "type")

pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, anchor.features = 8000, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(object = pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(object = pancreas.integrated)
pancreas.integrated <- RunPCA(object = pancreas.integrated, npcs = 30, verbose = FALSE)
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(pancreas.integrated);  ElbowPlot(pancreas.integrated)


options(repr.plot.width=10, repr.plot.height=6)
pbmc <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:13)
DimPlot(object = pbmc, label =T,label.size = 8)+ theme(legend.text=element_text(size=16,face="plain"),text=element_text(size=20,face="plain"))

Idents(pbmc)=pbmc$sample;table(Idents(pbmc))

options(repr.plot.width=18, repr.plot.height=7)
DimPlot(pbmc,label = T,label.size = 8,group.by=c('type','Annotation'),raster=FALSE,ncol=2)& theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=14, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,group.by=c('cx','origin2'),raster=FALSE,ncol=2)& theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))



DefaultAssay(pbmc) <- "integrated"

pbmc <- FindNeighbors(pbmc,dims = 1:13) 
pbmc <- FindClusters(object = pbmc, algorithm = 1, resolution =0.3) 
options(repr.plot.width=7.5, repr.plot.height=6)
DimPlot(object = pbmc,label = T,label.size = 8)+ theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))


options(repr.plot.width=7, repr.plot.height=6)
pbmc <- ReorderIdent(pbmc, c("PC_1"), reorder.numeric = T,reverse=F)
DimPlot(pbmc,label = T,label.size = 8,raster=FALSE)+ theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.29,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

pbmc$C22m=Idents(pbmc)
table(pbmc$C22m)






DefaultAssay(pbmc) <- "RNA"

options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("nFeature_RNA","nCount_RNA","score","percent.mt"), min.cutoff = "q01", max.cutoff = "q90",raster=FALSE,ncol=4) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


