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

names(pbmc@meta.data)

Idents(pbmc)=pbmc$sample;table(Idents(pbmc))

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

###Idents(pbmc, cells = WhichCells(pbmc,idents = "GD11n-2")) <- "14"
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

table(Idents(pbmc))

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

head(colnames(pbmc))

pbmc3=pbmc;pbmc3

names(pbmc@meta.data)

table(Idents(pbmc))

options(repr.plot.width=12, repr.plot.height=8)
 VlnPlot(pbmc,pt.size=0.01, features = c("nFeature_RNA","percent.mt"),raster=FALSE,ncol=1) +  theme(axis.title=element_blank(),axis.text=element_text(size=16),plot.title=element_text(size=16))+NoLegend()

options(repr.plot.width=12, repr.plot.height=8)
 VlnPlot(pbmc,pt.size=0, features = c("nFeature_RNA","percent.mt"),raster=FALSE,ncol=1) +  theme(axis.title=element_blank(),axis.text=element_text(size=16),plot.title=element_text(size=16))+NoLegend()

options(repr.plot.width=12, repr.plot.height=8)
 VlnPlot(pbmc,pt.size=0.01, group.by = "sample2", features = c("nFeature_RNA","percent.mt"),raster=FALSE,ncol=1) +  theme(axis.title=element_blank(),axis.text=element_text(size=16),plot.title=element_text(size=16))+NoLegend()

options(repr.plot.width=12, repr.plot.height=8)
 VlnPlot(pbmc,pt.size=0, group.by = "sample2", features = c("nFeature_RNA","percent.mt"),raster=FALSE,ncol=1) +  theme(axis.title=element_blank(),axis.text=element_text(size=16),plot.title=element_text(size=16))+NoLegend()

options(repr.plot.width=12, repr.plot.height=8)
 VlnPlot(pbmc,pt.size=0, group.by = "sample2",split.by="type", features = c("nFeature_RNA","percent.mt"),raster=FALSE,ncol=1) +  theme(axis.title=element_blank(),axis.text=element_text(size=16),plot.title=element_text(size=16))+NoLegend()

options(repr.plot.width=12, repr.plot.height=8)
 VlnPlot(pbmc,pt.size=0, group.by = "stage2",split.by="type", features = c("nFeature_RNA","percent.mt"),raster=FALSE,ncol=1) +  theme(axis.title=element_blank(),axis.text=element_text(size=16),plot.title=element_text(size=16))+NoLegend()



options(repr.plot.width=12, repr.plot.height=8)
 VlnPlot(pbmc,pt.size=0.01,features = c("nFeature_RNA","nCount_RNA"),ncol=1,raster=FALSE) & theme(axis.title=element_blank(),axis.text=element_text(size=16),plot.title=element_text(size=16))+NoLegend()

options(repr.plot.width=12, repr.plot.height=8)
 VlnPlot(pbmc,pt.size=0,features = c("nFeature_RNA","nCount_RNA"),ncol=1,raster=FALSE) & theme(axis.title=element_blank(),axis.text=element_text(size=16),plot.title=element_text(size=16))+NoLegend()

options(repr.plot.width=12, repr.plot.height=8)
 VlnPlot(pbmc,pt.size=0, group.by = "sample2",features = c("nFeature_RNA","nCount_RNA"),ncol=1,raster=FALSE) & theme(axis.title=element_blank(),axis.text=element_text(size=16),plot.title=element_text(size=16))+NoLegend()

colnames(pbmc@meta.data)

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

###Idents(pbmc, cells = WhichCells(pbmc,idents = "GD11n-2")) <- "14"
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

table(Idents(pbmc))

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

Idents(pbmc)="orig2"
new.cluster.ids <- c("GD8","GD8", "GD8n",
                     "GD9", "GD9","GD9n","GD9n",
                     "GD10","GD10" ,"GD10n","GD10n",
                     "GD11","GD11","GD11n",
                     "GD13","GD13n","GD13n",
                     "GD15","GD15","GD15n", "GD15n",
                     "GD17","GD17","GD17n", "GD17n",
                     "GD20","GD20","GD20n", "GD20n" ,
    "GD23","GD23n","GD23n",
    "GD26","GD26n","GD26n"  )
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
pbmc$stage3=Idents(pbmc)
table(Idents(pbmc))

Idents(pbmc)=pbmc$sample
table(Idents(pbmc))
Idents(pbmc)=pbmc$stage
table(Idents(pbmc))

Idents(pbmc)=pbmc$sample
table(Idents(pbmc))
pbmc$sample=Idents(pbmc)
table(Idents(pbmc))

Idents(pbmc)=pbmc$stage
table(Idents(pbmc))
pbmc$stage=Idents(pbmc)
table(Idents(pbmc))

options(repr.plot.width=15, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,group.by=c('sample'),ncol=2,raster=FALSE)& theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=15, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,group.by=c('sample2'),ncol=2,raster=FALSE)& theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=18, repr.plot.height=7)
DimPlot(pbmc,label = T,label.size = 8,group.by=c('type','Annotation'),raster=FALSE,ncol=2)& theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=14, repr.plot.height=6)
DimPlot(pbmc,label = F,label.size = 8,split.by=c('type'),raster=FALSE,ncol=2)& theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=14, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,split.by=c('type'),group.by='Annotation',raster=FALSE,ncol=2)& theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=14, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,split.by=c('type'),group.by='cx',raster=FALSE,ncol=2)& theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=16, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,group.by=c('cx','origin'),raster=FALSE,ncol=2)& theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=14, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,group.by=c('CX','origin2'),raster=FALSE,ncol=2)& theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

PBMC=pbmc;pbmc

options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("nFeature_RNA","nCount_RNA","score","PTPRC"), min.cutoff = "q01", max.cutoff = "q90",raster=FALSE,ncol=4) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("nFeature_RNA","nCount_RNA","score","percent.mt"), min.cutoff = "q01", max.cutoff = "q90",raster=FALSE,ncol=4) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=10, repr.plot.height=16)
FeaturePlot(pbmc, c("nFeature_RNA","nCount_RNA","score","percent.mt"), split.by='type',min.cutoff = "q01", max.cutoff = "q90",raster=FALSE,ncol=4) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("TFAP2C","VIM","EGFR","GATA3"), min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("GRIP1","IL20","RB1","LPL"), min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=10, repr.plot.height=16)
FeaturePlot(pbmc, c("TFAP2C","VIM","EGFR","GATA3"),split.by='type', min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=60, repr.plot.height=4)
FeaturePlot(pbmc, split.by='sample2',c("VIM","GATA3"), min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


pbmc

pbmc2=pbmc

saveRDS(pbmc, file = "rapl262632m-pca.rds")



PBMC=pbmc

median(pbmc$nFeature_RNA)
median(pbmc$percent.mt)
range(pbmc$nFeature_RNA)
range(pbmc$percent.mt)



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

options(repr.plot.width=8, repr.plot.height=6)
DimPlot(object = pbmc,label = T,label.size = 8,raster=FALSE)+ theme(legend.key.height=unit(0.3,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=18, repr.plot.height=20)
DimPlot(object = pbmc,label = T,label.size = 8,split.by="sample2",ncol=5,raster=FALSE)+ theme(axis.text=element_blank(),axis.ticks=element_blank(),
       legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.5,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20,face="plain"))#+guides(col = guide_legend(override.aes=list(size=3),ncol =4)) 

options(repr.plot.width=18, repr.plot.height=15)
DimPlot(object = pbmc,label = T,label.size = 8,split.by="stage",ncol=5,raster=FALSE)+ theme(axis.text=element_blank(),axis.ticks=element_blank(),
       legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.5,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20,face="plain"))#+guides(col = guide_legend(override.aes=list(size=3),ncol =4)) 

options(repr.plot.width=18, repr.plot.height=10)
DimPlot(object = pbmc,label = T,label.size = 8,split.by="stage2",group.by='C22m',ncol=5,raster=FALSE)+ theme(axis.text=element_blank(),axis.ticks=element_blank(),
       legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.5,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20,face="plain"))#+guides(col = guide_legend(override.aes=list(size=3),ncol =4)) 







DefaultAssay(pbmc) <- "RNA"

options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("nFeature_RNA","nCount_RNA","score","percent.mt"), min.cutoff = "q01", max.cutoff = "q90",raster=FALSE,ncol=4) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=10, repr.plot.height=16)
FeaturePlot(pbmc, c("nFeature_RNA","nCount_RNA","score","percent.mt"),split.by='type', min.cutoff = "q01", max.cutoff = "q90",raster=FALSE,ncol=4) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("CDH1","KDR","PECAM1","KRT7"), min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("GATA3","TFAP2C","GCM1","IL20"), min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("SOX11","HMGA2","DLK1","EPCAM"), min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("KLRD1","NKG7","CD79B","FCER1A"), min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("SATB1","PLCL1","TUBB1","MZB1"), min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("S100A12","CD14","TUBB1","MZB1"), min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("TTR","APOA1","LMO2","RUNX1"), min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("TTR","APOA1","AFP","RUNX1"), min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=10, repr.plot.height=8)
FeaturePlot(pbmc, c("APOA1","AFP"),split.by='type', min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("GJB3","RHOX6","HBB-BH1","PTPRC"), min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("DCN","CITED4","ACTA2","RGS5"), min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("CDX2","ACTA2","HBB-BH1","PTPRC"), min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))

options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("ALAS2","PLAC8","DCN","KDR"), min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("VIM","EGFR","GATA3","CDH1"), min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))


options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(pbmc, c("EGFR","TFAP2C","GATA3","POU5F1"), min.cutoff = "q01", max.cutoff = "q90",ncol=4,raster=FALSE) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))



Idents(pbmc)='origin2'
options(repr.plot.width=7, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,ncol=2,raster=FALSE)& theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

Idents(pbmc)='origin'
options(repr.plot.width=7, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,ncol=2,raster=FALSE)& theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

pbmc2=pbmc;pbmc2

saveRDS(pbmc, file = "rapl262632m-inte.rds")

pbmc=readRDS("rapl262632m-inte.rds");pbmc;names(pbmc@meta.data);pbmc2=pbmc



options(repr.plot.width=18, repr.plot.height=4)
FeaturePlot(pbmc, c("PTPRC", "CD14", "EGFR", "CDH1"),raster=FALSE, min.cutoff = "q01", max.cutoff = "q90",ncol=4) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"))

options(repr.plot.width=18, repr.plot.height=4)
FeaturePlot(pbmc, c("GATA3", "STAB1", "S100A12", "S100A8"),raster=FALSE, min.cutoff = "q01", max.cutoff = "q90",ncol=4) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"))

options(repr.plot.width=15.6, repr.plot.height=3.5)
FeaturePlot(pbmc, c("ENSOCUG00000005133","TOP2A","MCM5","PCNA"),raster=FALSE, min.cutoff = "q01", max.cutoff = "q90",ncol=5) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.size=unit(0.4,'cm'))


pbmc2=pbmc;pbmc


Idents(pbmc)='C22m'

pbmc4=pbmc
table(Idents(pbmc))

pbmc=pbmc4

Idents(pbmc, cells = WhichCells(pbmc,idents =22))<-"Unk."
Idents(pbmc, cells = WhichCells(pbmc,idents =16))<-"Mac and DC"
Idents(pbmc, cells = WhichCells(pbmc,idents =14))<-"B and Mast"
Idents(pbmc, cells = WhichCells(pbmc,idents =15))<-"NK and T"
Idents(pbmc, cells = WhichCells(pbmc,idents =20))<-"Neu"
Idents(pbmc, cells = WhichCells(pbmc,idents =13)) <- "Ery"
Idents(pbmc, cells = WhichCells(pbmc,idents =8)) <- "Endo"
Idents(pbmc, cells = WhichCells(pbmc,idents =c(7,9)))<- "Epi_m"
Idents(pbmc, cells = WhichCells(pbmc,idents =c(1))) <- "SMC"
Idents(pbmc, cells = WhichCells(pbmc,idents =c(4))) <- "PV_f"
Idents(pbmc, cells = WhichCells(pbmc,idents =c(6))) <- "SC_f"
Idents(pbmc, cells = WhichCells(pbmc,idents =c(12,10,5,3,2))) <- "SC_m"
Idents(pbmc, cells = WhichCells(pbmc,idents =c(11,17,18,19,21))) <- "Tr"
pbmc$Name2=Idents(pbmc)

options(repr.plot.width=9, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,raster=FALSE)+ theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.35,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=18, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,raster=FALSE,split.by='origin2')+ theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

pbmc

saveRDS(pbmc, file = "rapl262632m-inte.rds")

pbmc=readRDS("rapl262632m-inte.rds");pbmc;names(pbmc@meta.data)

DefaultAssay(pbmc) <- "RNA"





set.seed(123)  
Idents(pbmc) <- "Remained"
Idents(object = pbmc, cells = sample(colnames(pbmc), 10000)) <- "Random"
pbmcx<-subset(pbmc, idents ="Random", invert = F)
options(repr.plot.width=15, repr.plot.height=6)
DimPlot(pbmcx,label = T,label.size = 8)+ theme(legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

head(colnames(pbmcx))

Idents(pbmcx)='C22m'

table(Idents(pbmcx))

pbmc.markers <- FindAllMarkers(object = pbmcx, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dim(pbmc.markers)
pbmc

write.csv(pbmc.markers, "DEG-rapl262632m-inte-22.csv")

pbmcy <- ScaleData(pbmcx, features = rownames(pbmcx))

dim(pbmcx$RNA@scale.data)

pbmc.markers %>% group_by(cluster) %>% top_n(n = -10,  wt = p_val_adj) -> top10
top10 %>% group_by(cluster) %>% top_n(n = 10,  wt = avg_log2FC) -> top10
length(top10$gene)
options(repr.plot.width=30, repr.plot.height=45)
#pdf("MF248018_heat_14.pdf",width=10, height=15)
DoHeatmap(object = pbmcy,features = as.character(top10$gene),cells=sample(colnames(pbmcx), 10000)) + NoLegend()+theme(text=element_text(size=20,face="italic"))
#dev.off()),cells=sample(colnames(pbmc), 10000)

options(repr.plot.width=7, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,raster=FALSE)+ theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

table((pbmc$C22m))


Idents(pbmc)='Name2'
table(Idents(pbmc))

pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dim(pbmc.markers)
pbmc

write.csv(pbmc.markers, "DEG-rapl262632m-inte-22-name.csv")

pbmcy <- ScaleData(pbmc,  features = rownames(pbmc))

dim(pbmc$RNA@scale.data)

pbmc.markers %>% group_by(cluster) %>% top_n(n = -10,  wt = p_val_adj) -> top10
top10 %>% group_by(cluster) %>% top_n(n = 10,  wt = avg_log2FC) -> top10
length(top10$gene)
options(repr.plot.width=20, repr.plot.height=35)
#pdf("MF248018_heat_14.pdf",width=10, height=15)
DoHeatmap(object = pbmcy,features = as.character(top10$gene),cells=sample(colnames(pbmc), 10000)) + NoLegend()+theme(text=element_text(size=22,face="italic"))
#dev.off()

options(repr.plot.width=9, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,raster=FALSE)+ theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.35,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=9, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,raster=FALSE,group.by='type')+ theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.35,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))



pbmc=readRDS("rapl262632m-inte.rds");pbmc;names(pbmc@meta.data);pbmc2=pbmc

y=c("#635547","#90E3BF","#FF891C","#FF4A46","#f0e400" ,"#FF9C90" , "#0000A6","#63FFAC","#8EC792","#8B8A55",
"#6B2075","#E35F82","#12ED4C","#549E79","#127D4C",   "#BF9DA8","#BF9D88","#772600","#4A6798","#BD84B0",
"#BD3400","#C19F70","#0F4A9C","#FACB12", "#FCFF00",  "#FFD731","#dea93c",  "#2FDA00",   "#FFAAFF",  "#FF5C00"  ,
"#2F9A00","#2F4A60","#DABE99","#89C1F5",             "#C594BF","#5581CA","#005579","#00BFC4","#E3CB3A",
"#E3CB9A","#D5E839","#E85639","#5ADBE4","#F79083",   "#A64D7E","#B51D8D","#3F84AA","#456722","#F9DFE6",
"#683ED8","#685ED8","#683EA8","#FF7F9C","#AAFFAA",   "#532C8A","#532C5A","#DDAA22", "#CC5522","#CC7818",
"#9D0049","#FBBE92","#7C2A47","#C72228","#EF4E22",   "#95E1FF",'#FFB7FF','#FF00B2','#F397C0',"#7295e2",
"#46c5bf","#1A1A1A","#97BAD3","#7F6874","#9d506e",   "#F6BFCB","#5B4534","#989898","#7766db","#153B3D")
options(repr.plot.width=16, repr.plot.height=3)
barplot(rep(2,79),col=y)

z=c("#635547","#90E3BF","#FF891C","#BBCCBB","#89C1F5","#127D4C",
"#12ED4C",  "#BD84B0",  "#C19F70",   "#0F4A9C" , 
"#FACB12",  "#FCFF00",  "#dea93c",   "#2F9A00",    "#DABE99",
"#C594BF",  "#5581CA",  "#E85639",   'lightblue',   "#F79083",   
"#A64D7E",  "#B51D8D",  "#3F84AA",   "#456722",    "#F9DFE6",
"#685ED8",  "#AAFFAA",  "#683EA8",   "#FF7F9C",    "#532C8A",  
"#CC7818",  "#9D0049",  "#FBBE92",   "#C72228",    "#EF4E22", 
'#FFB7FF',  "#EF5A9D",  "#7F6874",   "#989898",    "#153B3D")
length(z)
options(repr.plot.width=16, repr.plot.height=3)
barplot(rep(2,40),col=z)

pbmc1=readRDS("rapl262632m-inte.rds");pbmc1;names(pbmc1@meta.data)

pbmc2=readRDS("sn148053-sc114572-2nd.rds");pbmc2;names(pbmc2@meta.data)

pbmc3=readRDS("/sdc/xxjiang/rapl/rapl114572_new.rds");pbmc3;names(pbmc3@meta.data)

pbmc3$CX3=pbmc3$CX2
pbmc3$Type3=pbmc3$Type2
pbmc3$cx=pbmc1$C25FX
pbmc2$cx=pbmc2$C16H
pbmc3$type='SC'
pbmc2$type='SN'
pbmc3$Annotation=pbmc1$Name
pbmc2$Annotation=pbmc2$stage
names(pbmc3@meta.data);names(pbmc2@meta.data)

options(repr.plot.width=8, repr.plot.height=6)
DimPlot(object = pbmc3,label = T,label.size = 8,raster=FALSE,group.by='origin2')+ theme(axis.line=element_blank(),axis.title=element_blank(),legend.key.height=unit(0.4,"inch"),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

pbmcx=pbmc3
Idents(pbmcx)="CX3"
options(repr.plot.width=18, repr.plot.height=6)
DimPlot(pbmcx,label = T,label.size = 6, raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.35,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=y[as.numeric(names(table(Idents(pbmcx))))],labels =rownames(table(pbmcx$Type3))) +guides(col = guide_legend(override.aes=list(size=3),ncol =4)) 
Idents(pbmcx)="CX1"
options(repr.plot.width=9.6, repr.plot.height=6)
DimPlot(pbmcx,label = T,label.size = 6, ncol=1,raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=z[as.numeric(names(table(Idents(pbmcx))))],labels =rownames(table(pbmcx$Type1))) +guides(col = guide_legend(override.aes=list(size=3),ncol =2)) 


pbmcx=pbmc2
Idents(pbmcx)="CX3"
options(repr.plot.width=15, repr.plot.height=6)
DimPlot(pbmcx,label = T,label.size = 6, raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.35,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=y[as.numeric(names(table(Idents(pbmcx))))],labels =rownames(table(pbmcx$Type3))) +guides(col = guide_legend(override.aes=list(size=3),ncol =4)) 
Idents(pbmcx)="CX1"
options(repr.plot.width=9.6, repr.plot.height=6)
DimPlot(pbmcx,label = T,label.size = 6, ncol=1,raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=z[as.numeric(names(table(Idents(pbmcx))))],labels =rownames(table(pbmcx$Type1))) +guides(col = guide_legend(override.aes=list(size=3),ncol =2)) 


pbmc4 <- merge(x =pbmc3, y = pbmc2);      pbmc4;names(pbmc4@meta.data)

pancreas.query <- AddMetaData(pbmc1, metadata = pbmc4@meta.data)

names(pancreas.query @meta.data)

pancreas.query 

table(Idents(pancreas.query))



other2=pancreas.query
Idents(other2)="Type3"
levels(other2)=names(table(Idents(other2)))[order( as.numeric(substr(names(table(Idents(other2))), 1, 2)))]
names(table(Idents(other2)))
other2$Type3=Idents(other2)
Idents(other2)="Type1"
levels(other2)=names(table(Idents(other2)))[order( as.numeric(substr(names(table(Idents(other2))), 1, 2)))]
names(table(Idents(other2)))
other2$Type1=Idents(other2)

Idents(other2)="CX3"
levels(other2)=names(table(Idents(other2)))[order( as.numeric(substr(names(table(Idents(other2))), 1, 2)))]
names(table(Idents(other2)))
other2$CX3=Idents(other2)
Idents(other2)="CX1"
levels(other2)=names(table(Idents(other2)))[order( as.numeric(substr(names(table(Idents(other2))), 1, 2)))]
names(table(Idents(other2)))
other2$CX1=Idents(other2)

options(repr.plot.width=12, repr.plot.height=6)
DimPlot(object = other2,split.by='type',group.by='CX3',label = T,label.size = 8,raster=FALSE)+ theme(legend.key.height=unit(0.3,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=12, repr.plot.height=6)
DimPlot(object = other2,split.by='type',group.by='CX1',label = T,label.size = 8,raster=FALSE)+ theme(legend.key.height=unit(0.3,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))



Idents(other2)='CX1'

other2 <-subset(other2, idents =c(1,2,3,4,5,6,7,8,9,10,11,12,13), invert = F)

other2

y=c("#635547","#90E3BF","#FF891C","#FF4A46","#f0e400" ,"#FF9C90" , "#0000A6","#63FFAC","#8EC792","#8B8A55",
"#6B2075","#E35F82","#12ED4C","#549E79","#127D4C",   "#BF9DA8","#BF9D88","#772600","#4A6798","#BD84B0",
"#BD3400","#C19F70","#0F4A9C","#FACB12", "#FCFF00",  "#FFD731","#dea93c",  "#2FDA00",   "#FFAAFF",  "#FF5C00"  ,
"#2F9A00","#2F4A60","#DABE99","#89C1F5",             "#C594BF","#5581CA","#005579","#00BFC4","#E3CB3A",
"#E3CB9A","#D5E839","#E85639","#5ADBE4","#F79083",   "#A64D7E","#B51D8D","#3F84AA","#456722","#F9DFE6",
"#683ED8","#685ED8","#683EA8","#FF7F9C","#AAFFAA",   "#532C8A","#532C5A","#DDAA22", "#CC5522","#CC7818",
"#9D0049","#FBBE92","#7C2A47","#C72228","#EF4E22",   "#95E1FF",'#FFB7FF','#FF00B2','#F397C0',"#7295e2",
"#46c5bf","#1A1A1A","#97BAD3","#7F6874","#9d506e",   "#F6BFCB","#5B4534","#989898","#7766db","#153B3D")
options(repr.plot.width=16, repr.plot.height=3)
barplot(rep(2,79),col=y)

z=c("#635547","#90E3BF","#FF891C","#BBCCBB","#89C1F5","#127D4C",
"#12ED4C",  "#BD84B0",  "#C19F70",   "#0F4A9C" , 
"#FACB12",  "#FCFF00",  "#dea93c",   "#2F9A00",    "#DABE99",
"#C594BF",  "#5581CA",  "#E85639",   'lightblue',   "#F79083",   
"#A64D7E",  "#B51D8D",  "#3F84AA",   "#456722",    "#F9DFE6",
"#685ED8",  "#AAFFAA",  "#683EA8",   "#FF7F9C",    "#532C8A",  
"#CC7818",  "#9D0049",  "#FBBE92",   "#C72228",    "#EF4E22", 
'#FFB7FF',  "#EF5A9D",  "#7F6874",   "#989898",    "#153B3D")
length(z)
options(repr.plot.width=16, repr.plot.height=3)
barplot(rep(2,40),col=z)

pbmcx=other2
Idents(pbmcx)="CX3"
options(repr.plot.width=18, repr.plot.height=6)
DimPlot(pbmcx,label = T,label.size = 6, raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.35,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=y[as.numeric(names(table(Idents(pbmcx))))],labels =rownames(table(pbmcx$Type3))) +guides(col = guide_legend(override.aes=list(size=3),ncol =4)) 
Idents(pbmcx)="CX1"
options(repr.plot.width=9.8, repr.plot.height=6)
DimPlot(pbmcx,label = T,label.size = 6, ncol=1,raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=z[as.numeric(names(table(Idents(pbmcx))))],labels =rownames(table(pbmcx$Type1))) +guides(col = guide_legend(override.aes=list(size=3),ncol =2)) 


pbmcx=other2
Idents(pbmcx)="CX3"
options(repr.plot.width=24, repr.plot.height=6)
DimPlot(pbmcx,label = T,label.size = 6, raster=FALSE,repel=F,split.by='type')+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.35,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=y[as.numeric(names(table(Idents(pbmcx))))],labels =rownames(table(pbmcx$Type3))) +guides(col = guide_legend(override.aes=list(size=3),ncol =4)) 
Idents(pbmcx)="CX1"
options(repr.plot.width=14, repr.plot.height=6)
DimPlot(pbmcx,label = T,label.size = 6, ncol=2,raster=FALSE,repel=F,split.by='type')+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=z[as.numeric(names(table(Idents(pbmcx))))],labels =rownames(table(pbmcx$Type1))) +guides(col = guide_legend(override.aes=list(size=3),ncol =2)) 


other2

saveRDS(other2, file = "rapl262625-inte-ref.rds")

names(pbmc@meta.data)

options(repr.plot.width=9, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,group.by='Name2',raster=FALSE)+ theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.35,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=13, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,group.by='Name2',split.by='type',raster=FALSE)+ theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.35,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=9, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,raster=FALSE,group.by='type')+ theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.35,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=8, repr.plot.height=6)
DimPlot(object = other2,label = T,label.size = 8,raster=FALSE,group.by='origin2')+ theme(axis.line=element_blank(),axis.title=element_blank(),legend.key.height=unit(0.4,"inch"),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

table(pbmc$type)



pbmcx=other2
Idents(pbmcx)="CX3"
options(repr.plot.width=32, repr.plot.height=24)
DimPlot(pbmcx,label = T,label.size = 6, raster=FALSE,repel=F,split.by='sample2',ncol=5)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.35,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=y[as.numeric(names(table(Idents(pbmcx))))],labels =rownames(table(pbmcx$Type3))) +guides(col = guide_legend(override.aes=list(size=3),ncol =4)) 
Idents(pbmcx)="CX1"
options(repr.plot.width=14, repr.plot.height=16)
DimPlot(pbmcx,label = T,label.size = 6, ncol=5,raster=FALSE,repel=F,split.by='sample2')+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=z[as.numeric(names(table(Idents(pbmcx))))],labels =rownames(table(pbmcx$Type1))) +guides(col = guide_legend(override.aes=list(size=3),ncol =2)) 


Idents(other2)="CX3"
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(other2,cells.highlight = CellsByIdentities(object = other2, idents =c(17)),raster=FALSE,reduction="umap", label=F,label.size = 8) & theme(legend.key.height=unit(0.3,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))

Idents(other2)="CX3"
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(other2,split.by='type',cells.highlight = CellsByIdentities(object = other2, idents =c(17)),raster=FALSE,reduction="umap", label=F,label.size = 8) & theme(legend.key.height=unit(0.3,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))

other2



DefaultAssay(other2)='RNA';pbmc=other2;pbmc

Idents(pbmc)='Type1'
table(Idents(pbmc))

set.seed(123) 
Idents(pbmc) <- "Remained"
Idents(object = pbmc, cells = sample(colnames(pbmc), 10000)) <- "Random"
pbmcx<-subset(pbmc, idents ="Random", invert = F)
options(repr.plot.width=15, repr.plot.height=6)
DimPlot(pbmcx,label = T,label.size = 8)+ theme(legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

head(colnames(pbmcx))

Idents(pbmcx)='Type1'

table(Idents(pbmcx))

pbmc.markers <- FindAllMarkers(object = pbmcx, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dim(pbmc.markers)
pbmc

write.csv(pbmc.markers, "DEG-rapl262632m-inte-Type1.csv")

pbmcy <- ScaleData(pbmcx, features = rownames(pbmcx))

dim(pbmcx$RNA@scale.data)

pbmc.markers %>% group_by(cluster) %>% top_n(n = -10,  wt = p_val_adj) -> top10
top10 %>% group_by(cluster) %>% top_n(n = 10,  wt = avg_log2FC) -> top10
length(top10$gene)
options(repr.plot.width=17, repr.plot.height=27)
#pdf("MF248018_heat_14.pdf",width=10, height=15)
DoHeatmap(object = pbmcy,features = as.character(top10$gene),cells=sample(colnames(pbmcx), 10000)) + NoLegend()+theme(text=element_text(size=20,face="italic"))
#dev.off()),cells=sample(colnames(pbmc), 10000)







options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(other2, c("PTX3","ARG1","VIM","ENSOCUG00000017253"), min.cutoff = "q01", max.cutoff = "q90",raster=FALSE,ncol=4) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))

options(repr.plot.width=17, repr.plot.height=4)
FeaturePlot(other2, c("VIM","LUM","DLK1","ACTA2"), min.cutoff = "q01", max.cutoff = "q90",raster=FALSE,ncol=4) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))

options(repr.plot.width=17, repr.plot.height=28)
FeaturePlot(other2, c('ENSOCUG00000023508',
'CRISP3',
'MMRN1',
'ENSOCUG00000027987',
'CXCR4',
'ENSOCUG00000038559',
'CD52',
'ACTG1',
'GNG11',
'LGALS1',
'ESM1',
'SNORA65',
'PTP4A3',
'MARCKSL1',
'ENSOCUG00000025758',
'ENSOCUG00000027378',
'UNC5B',
'ENSOCUG00000009767',
'DUSP5',
'FABP5',
'BTG2',
'ENSOCUG00000035524',
'LXN',
'FXYD6',
'COX4I1',
'ENSOCUG00000021938',
'ENSOCUG00000011176',
'ENSOCUG00000005867'

), min.cutoff = "q01", max.cutoff = "q90",raster=FALSE,ncol=4) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))

options(repr.plot.width=17, repr.plot.height=28)
FeaturePlot(other2, c('ASRGL1',
'MECOM',
'ENSOCUG00000017253',
'PAX8',
'IYD',
'RNF144B',
'PCYT1B',
'PTHLH',
'ENSOCUG00000027736',
'ARHGEF26',
'ARG1',
'GRM8',
'PKHD1',
'F3',
'TEX15',
'ACSS1',
'CDC42EP3',
'ATP7B',
'ADCY5',
'B3GNT2',
'CLDN8',
'RIMKLB',
'USP43',
'TRPM6',
'ENSOCUG00000004655',
'BCL9',
'PRR15',
'PAX2'
), min.cutoff = "q01", max.cutoff = "q90",raster=FALSE,ncol=4) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.key.height=unit(0.4,'inch'))





getwd();library(Seurat);library(Matrix);suppressMessages(library(dplyr));library(ggplot2);suppressMessages(library('pagoda2'))#;suppressMessages(library("velocyto.R"));suppressMessages(library(monocle));
suppressMessages(library(cowplot));suppressMessages(library("patchwork"));suppressMessages(library(harmony))
options(repr.plot.width=7, repr.plot.height=6); library("RColorBrewer");library(magrittr);library(scales)
setwd("/sdc/xxjiang/rapl")
getwd()
packageVersion('Seurat')

pbmc=readRDS("rapl262625-inte-ref.rds");pbmc;names(pbmc@meta.data);pbmcx=pbmc

DefaultAssay(pbmc)='RNA'
pbmc

df1 =as.data.frame(rownames(pbmc))
colnames(df1) <- 'key'
df1[1:10,]

PBMC <- Read10X(data.dir = "/sdc/xxjiang/rapl/ncbref/GD11/outs/filtered_feature_bc_matrix");dim(PBMC)
pbmcx <- CreateSeuratObject(counts = PBMC)  #建立字典
pbmcx

table(duplicated(rownames(PBMC)))

df2=read.table("/sdc/xxjiang/rapl/rabbit-gene-cellranger.tsv",sep="\t",check.names = FALSE, header=F,row.names=1)
df2[9:12,]
df2$key=rownames(pbmcx)#rownames(df2)
row.names(df2) <- NULL
head(df2);dim(df2)

length(setdiff(df2$key,df1$key))
length(setdiff(df1$key,df2$key))

library(dplyr)
# 使用 left_join() 根据 'Key' 列合并数据框
merged_df <- left_join(df1, df2, by = "key")
head(merged_df)

dim(merged_df)
pbmc
head(merged_df)
head(rownames(pbmc))
head(merged_df$V4)
length(merged_df$V4)

pbmc

rownames(pbmc$RNA@counts)=merged_df$V4
rownames(pbmc$RNA@data)=merged_df$V4
head(rownames(pbmc$RNA@counts))

options(repr.plot.width=20, repr.plot.height=6)
FeaturePlot(pbmc, c('CD3E','IFIT5',"TREM2","CD209"),raster=FALSE,  min.cutoff = "q01", max.cutoff = "q90",ncol=4) & theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(size=16),plot.title =element_text(face = "bold.italic"),legend.position=c(0.8,0.85),legend.key.size=unit(0.3,'cm'))

pbmc

DefaultAssay(pbmc)='integrated'

# 检查参考数据集 HS 的变量基因
ref_features <- VariableFeatures(pbmc)
print(paste("参考数据集变量基因数:", length(ref_features)))

ref_features





Idents(pbmc)='type'
pbmc1=subset(pbmc,idents='SC');pbmc1
pbmc2=subset(pbmc,idents='SN');pbmc2

pbmc1 <- CreateSeuratObject(pbmc1$RNA@counts, meta.data = pbmc1@meta.data ,min.cells = 3, project = "trophoblast");pbmc1
pbmc2 <- CreateSeuratObject(pbmc2$RNA@counts, meta.data = pbmc2@meta.data ,min.cells = 3, project = "trophoblast");pbmc2



pbmc <- merge(x =pbmc1, y = pbmc2);      pbmc

pbmc <- CreateSeuratObject(pbmc$RNA@counts, meta.data = pbmc@meta.data ,min.cells = 3, project = "trophoblast");pbmc

names(pbmc@meta.data);pbmc;table(pbmc$sample)

table(Idents(pbmc))

names(pbmc@meta.data)

Idents(pbmc)=pbmc$sample;table(Idents(pbmc))

names(pbmc@meta.data)

table(Idents(pbmc))

colnames(pbmc@meta.data)

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

###Idents(pbmc, cells = WhichCells(pbmc,idents = "GD11n-2")) <- "14"
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

table(Idents(pbmc))

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

Idents(pbmc)="orig2"
new.cluster.ids <- c("GD8","GD8", "GD8n",
                     "GD9", "GD9","GD9n","GD9n",
                     "GD10","GD10" ,"GD10n","GD10n",
                     "GD11","GD11","GD11n",
                     "GD13","GD13n","GD13n",
                     "GD15","GD15","GD15n", "GD15n",
                     "GD17","GD17","GD17n", "GD17n",
                     "GD20","GD20","GD20n", "GD20n" ,
    "GD23","GD23n","GD23n",
    "GD26","GD26n","GD26n"  )
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
pbmc$stage3=Idents(pbmc)
table(Idents(pbmc))

Idents(pbmc)=pbmc$sample
table(Idents(pbmc))
Idents(pbmc)=pbmc$stage
table(Idents(pbmc))

Idents(pbmc)=pbmc$sample
table(Idents(pbmc))
pbmc$sample=Idents(pbmc)
table(Idents(pbmc))

Idents(pbmc)=pbmc$stage
table(Idents(pbmc))
pbmc$stage=Idents(pbmc)
table(Idents(pbmc))

options(repr.plot.width=15, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8,group.by=c('sample'),ncol=2,raster=FALSE)& theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

pbmcx=readRDS("rapl262625-inte-ref.rds");pbmcx;names(pbmcx@meta.data)

DefaultAssay(pbmcx)='RNA'

pbmcx

sp=pbmc

sp$umap@cell.embeddings=pbmcx$umap@cell.embeddings
sp$umap@misc$model$embedding=pbmcx$umap@cell.embeddings

options(repr.plot.width=15, repr.plot.height=6)
DimPlot(sp,label = T,label.size = 8,group.by=c('sample'),ncol=2,raster=FALSE)& theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=18, repr.plot.height=6)
DimPlot(sp,label = T,label.size = 8,group.by=c('Typex'),ncol=2,raster=FALSE)& theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

pbmc=sp;pbmc

saveRDS(pbmc, file = "rapl262625-inte-ref-V4.rds",compress=F)














