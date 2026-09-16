getwd();library(Seurat);library(Matrix);suppressMessages(library(dplyr));library(ggplot2);suppressMessages(library('pagoda2'));suppressMessages(library("velocyto.R"));suppressMessages(library(monocle));
suppressMessages(library(cowplot));suppressMessages(library("patchwork"));suppressMessages(library(harmony))
options(repr.plot.width=7, repr.plot.height=6); library("RColorBrewer");library(magrittr);library(scales)
setwd("/sdc/xxjiang/rapl")
getwd()
packageVersion('Seurat')

#mouse  rabbit

pbmc1=readRDS("tr79843m-inte.rds");pbmc1;names(pbmc1@meta.data)

#Idents(pbmc1)='C5t'
options(repr.plot.width=8.5, repr.plot.height=6)
DimPlot(pbmc1 , label=T,label.size = 8) & theme(legend.key.height=unit(0.35,"inch"),axis.text=element_blank(),axis.ticks=element_blank(),text=element_text(size=20),legend.text=element_text(size=20,face = "italic"))

names(table(pbmc1$stage2))

Idents(pbmc1)="stage2"
pbmc1=subset(pbmc1,idents=c('GD8',  'GD9' ,'GD10' ,'GD11' ,'GD13','GD15',"GD17"));pbmc1



pbmc2=readRDS("other/MTR32518INTE-R4.rds");pbmc2;names(pbmc2@meta.data)

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

Idents(pbmc2)='C21'
options(repr.plot.width=14, repr.plot.height=6)
DimPlot(object = pbmc2,label = T,raster=FALSE,split.by='origin',label.size =8,repel=F,reduction='umap')+  theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.35,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20))  +
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.2,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=z,labels =rownames(table(pbmc2$type3)))+guides(col = guide_legend(override.aes=list(size=3),ncol =1)) 



Idents(pbmc2)='type3'
names(table(pbmc2$type3))

pbmc=pbmc2


Idents(pbmc, cells = WhichCells(pbmc,idents = c('18 EPC Migratory Cell'))) <- "9   EPC Migratory Cell"
Idents(pbmc, cells = WhichCells(pbmc,idents = '17 SpA-TGC')) <- '8   SpA-TGC'
Idents(pbmc, cells = WhichCells(pbmc,idents = c('14-15   SpT','E1 SpT Precursor','P1-2-3 SpT&S-TGC-P Precursor'))) <- '7   SpT'
Idents(pbmc, cells = WhichCells(pbmc,idents = c('10 Primary P-TGC','11 Secondary P-TGC','12 Secondary P-TGC Precursor'))) <- '6   P-TGC'
Idents(pbmc, cells = WhichCells(pbmc,idents = c('6-7-8 S-TGC Precursor','9   S-TGC','E2 S-TGC-P Precursor'))) <- '5   S-TGC'
Idents(pbmc, cells = WhichCells(pbmc,idents = '16 Gly-T')) <- '4   Gly-T'
Idents(pbmc, cells = WhichCells(pbmc,idents = c('19 SynTII Precursor','20 SynTII'))) <- '3   SynTII'
Idents(pbmc, cells = WhichCells(pbmc,idents = c('4-5 SynTI Precursor','21 SynTI'))) <- '2   SynTI'
Idents(pbmc, cells = WhichCells(pbmc,idents = c('1   TSC','2   LaTP','3   LaTP 2'))) <- "1   LaTP"
pbmc$type3x=Idents(pbmc)
table(Idents(pbmc))


#Idents(pbmc1)='C5t'
options(repr.plot.width=8.5, repr.plot.height=6)
DimPlot(pbmc , label=T,label.size = 8) & theme(legend.key.height=unit(0.35,"inch"),axis.text=element_blank(),axis.ticks=element_blank(),text=element_text(size=20),legend.text=element_text(size=20,face = "italic"))

new.cluster.ids <- 1:9
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
pbmc$cx3=Idents(pbmc)  

pbmc2=pbmc;pbmc2

saveRDS(pbmc2,"other/MTR32518INTE-R4.rds");pbmc2;names(pbmc2@meta.data)





Idents(pbmc2)='cx'
options(repr.plot.width=14, repr.plot.height=6)
DimPlot(object = pbmc2,label = T,raster=FALSE,split.by='origin',label.size =8,repel=F,reduction='umap')+  theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20))  +
scale_colour_manual(values=z[as.numeric(names(table(Idents(pbmc2))))],labels =rownames(table(pbmc2$typex)))+guides(col = guide_legend(override.aes=list(size=3),ncol =1)) 

Idents(pbmc2)='cx3'
options(repr.plot.width=14, repr.plot.height=6)
DimPlot(object = pbmc2,label = T,raster=FALSE,split.by='origin',label.size =8,repel=F,reduction='umap')+  theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20))  +
scale_colour_manual(values=z[as.numeric(names(table(Idents(pbmc2))))],labels =rownames(table(pbmc2$type3x)))+guides(col = guide_legend(override.aes=list(size=3),ncol =1)) 

Idents(pbmc2)='cx3'
options(repr.plot.width= 18, repr.plot.height=6)
DimPlot(object = pbmc2,label = T,raster=FALSE,split.by='stage',label.size =8,repel=F,reduction='umap')+  theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.35,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20))  +
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.2,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=z,labels =rownames(table(pbmc2$type3x)))+guides(col = guide_legend(override.aes=list(size=3),ncol =1)) 

table(pbmc2$origin)

table(pbmc2$stage)



pbmc1;pbmc2



Idents(pbmc1)="sample";table(Idents(pbmc1))

Idents(pbmc1)="stage";table(Idents(pbmc1))

Idents(pbmc1, cells = WhichCells(pbmc1,idents = c( 'GD8-1', 'GD8-2' ,'GD9-1', 'GD9-2'))) <- "GD8"
Idents(pbmc1, cells = WhichCells(pbmc1,idents = c( 'GD10-1','GD10-2','GD11-1','GD13'))) <- "GD8" 

Idents(pbmc1, cells = WhichCells(pbmc1,idents = c(  'GD8n', 'GD9n-1','GD9n-2','GD10n-1','GD10n-2'))) <- "GD11n" 


pbmc1$samplex=Idents(pbmc1)
table(Idents(pbmc1))

Idents(pbmc2)="stage";table(Idents(pbmc2))

pbmc2$samplex=Idents(pbmc2)
table(Idents(pbmc2))

pbmc1$typex=pbmc1$Typex1
pbmc2$typex=pbmc2$type3x

pbmc1$cx=pbmc1$CX1
pbmc2$cx=pbmc2$cx3

pbmc1$species='Rabbit'
pbmc2$species='Mouse'

pbmc1$type=pbmc1$type
pbmc2$type=pbmc2$origin

table(Idents(pbmc1));table(Idents(pbmc2))

pbmc1;pbmc2

#兔子与human整合第一种 取交集

#rownames(pbmc_ht$RNA@counts)=toupper(rownames(pbmc_ht$RNA@counts))
#rownames(pbmc_r$RNA@counts)=toupper(rownames(pbmc_r$RNA@counts))

#pbmc_ht=pbmc_ht1;pbmc_r=pbmc_r1

library('homologene')
transformHomoloGene <- function(exp_sc_mat, inTaxID = 10090, outTaxID = 9606) {
    library(homologene)
    genes.in <- rownames(exp_sc_mat)
    res.home <- homologene(genes.in, inTax = inTaxID, outTax = outTaxID)
    res.home <- res.home[!duplicated(res.home[, 1]),]
    res.home <- res.home[!duplicated(res.home[, 2]),]
    genes.out <- res.home[, 1]
    genes.homo <- res.home[, 2]
    exp.out <- exp_sc_mat[genes.out,]
    rownames(exp.out) <- genes.homo
    return(exp.out)
}

packageVersion('homologene')

pbmc2
count=transformHomoloGene(pbmc2$RNA@counts)
dim(count)
pbmc2 <- CreateSeuratObject(count, meta.data = pbmc2@meta.data, min.cells = 3 );pbmc2



pbmc1 <- CreateSeuratObject(pbmc1$RNA@counts, meta.data = pbmc1@meta.data, min.cells = 3, project = "trophoblast");pbmc1
pbmc2 <- CreateSeuratObject(pbmc2$RNA@counts, meta.data = pbmc2@meta.data, min.cells = 3, project = "trophoblast");pbmc2

pbmc <- merge(x =pbmc1, y =pbmc2)
pbmc;table(pbmc$sample);table(Idents(pbmc))

Idents(pbmc)="samplex"
table(Idents(pbmc))



strh=rownames(pbmc1);length(strh)
strm=rownames(pbmc2$RNA@counts);length(strm)
str=intersect(strh,strm);length(str);head(str)
count=(pbmc$RNA@counts[str,]);dim(count)

pbmc <- CreateSeuratObject(count, meta.data = pbmc@meta.data, min.cells = 3, project = "trophoblast");pbmc



names(pbmc@meta.data);pbmc;table(pbmc$samplex,pbmc$species)



pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(pbmc))

pancreas=pbmc
pancreas.list <- SplitObject(object = pancreas, split.by = "samplex")

pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, anchor.features = 5000, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(object = pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(object = pancreas.integrated)
pancreas.integrated <- RunPCA(object = pancreas.integrated, npcs = 30, verbose = FALSE)
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(pancreas.integrated);  ElbowPlot(pancreas.integrated)

options(repr.plot.width=10, repr.plot.height=6)
pbmc <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:8)
DimPlot(object = pbmc, label =T,label.size = 8)+ theme(legend.text=element_text(size=16,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=18, repr.plot.height=6)
DimPlot(object = pbmc,split.by='species',group.by='typex',label = T,label.size = 8,raster=FALSE)+ theme(legend.key.height=unit(0.25,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=18, repr.plot.height=6)
DimPlot(object = pbmc,split.by='species',group.by='typex',label = T,label.size = 8,raster=FALSE)+ theme(legend.key.height=unit(0.3,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=7, repr.plot.height=6)
DimPlot(object = pbmc,group.by='species',label = T,label.size = 8,raster=FALSE)+ theme(legend.key.height=unit(0.3,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))


z=rbind(as.matrix(pbmc1$Typex1),as.matrix(pbmc2$type3x));length(z);head(z)
z=z[order(rownames(z)),];length(z);head(z)
pbmc$typex=z



names(table(pbmc$typex))


Idents(pbmc, cells = WhichCells(pbmc,idents = "8   STB2")) <- "8   JZT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "7   STB1")) <- "7   STB"

pbmc$typex=Idents(pbmc)

Idents(pbmc)="typex"
levels(pbmc)=names(table(Idents(pbmc)))[order( as.numeric(substr(names(table(Idents(pbmc))), 1, 2)))]
names(table(Idents(pbmc)))
pbmc$typex=Idents(pbmc)

levels(pbmc)=c('1   CTB1','2   CTB2','3   ColumnT','4   GlyT1','5   GlyT2','6   GlyT3',
               '7   STB','8   JZT','9   XTB',
               '1   LaTP','2   SynTI','3   SynTII','4   Gly-T','5   S-TGC','6   P-TGC',
               '7   SpT','8   SpA-TGC','9   EPC Migratory Cell' )

pbmc$typex=Idents(pbmc)

options(repr.plot.width=18, repr.plot.height=6)
DimPlot(object = pbmc,split.by='species',group.by='typex',label = T,label.size = 8,raster=FALSE)+ theme(legend.key.height=unit(0.25,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

new.cluster.ids <- c("1","2","3","4",'5','6','7','8','9',
                    '1 ',"2 ","3 ","4 ","5 ",'6 ','7 ','8 ','9 ')
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
pbmc$cx=Idents(pbmc)
names(table(Idents(pbmc)))

options(repr.plot.width=18, repr.plot.height=6)
DimPlot(object = pbmc,split.by='species',label = T,label.size = 8,raster=FALSE)+ theme(legend.key.height=unit(0.3,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

pbmc
saveRDS(pbmc, file = "inte/GD17-tr74522-scsn-mr-inte.rds",compress=F)
