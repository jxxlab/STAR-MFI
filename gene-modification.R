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
pbmcx <- CreateSeuratObject(counts = PBMC)  
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


ref_features <- VariableFeatures(pbmc)
print(paste("genes:", length(ref_features)))

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














