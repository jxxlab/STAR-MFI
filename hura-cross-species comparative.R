getwd();library(Seurat);library(Matrix);suppressMessages(library(dplyr));library(ggplot2);suppressMessages(library('pagoda2'));suppressMessages(library("velocyto.R"));suppressMessages(library(monocle));
suppressMessages(library(cowplot));suppressMessages(library("patchwork"));suppressMessages(library(harmony))
options(repr.plot.width=7, repr.plot.height=6); library("RColorBrewer");library(magrittr);library(scales)
setwd("/sdc/xxjiang/rapl")
getwd()
packageVersion('Seurat')

#human  rabbit

pbmc1=readRDS("tr79843m-inte.rds");pbmc1;names(pbmc1@meta.data)

#Idents(pbmc1)='C5t'
options(repr.plot.width=8.5, repr.plot.height=6)
DimPlot(pbmc1 , label=T,label.size = 8) & theme(legend.key.height=unit(0.35,"inch"),axis.text=element_blank(),axis.ticks=element_blank(),text=element_text(size=20),legend.text=element_text(size=20,face = "italic"))

#Idents(pbmc1)='C5t'
options(repr.plot.width=8.5, repr.plot.height=6)
DimPlot(pbmc1, label=T,label.size = 8,group.by="Typex1") & theme(legend.key.height=unit(0.35,"inch"),axis.text=element_blank(),axis.ticks=element_blank(),text=element_text(size=20),legend.text=element_text(size=20,face = "italic"))

Idents(pbmc1)="stage2"
pbmc1=subset(pbmc1,idents=c('GD8','GD9','GD10' ,'GD11' ,'GD13','GD15',"GD17"));pbmc1

Idents(pbmc1)="Typex1"
pbmc1$Typex1=Idents(pbmc1)







pbmc2=readRDS("other/hu-early-tr-14083.rds");pbmc2;names(pbmc2@meta.data)

Idents(pbmc2)="stage2";table(Idents(pbmc2))

options(repr.plot.width=12, repr.plot.height=6)
DimPlot(pbmc2,  label=T,label.size = 8,group.by=c('cx','C20A')) & theme(legend.key.height=unit(0.3,"inch"),axis.text=element_blank(),axis.ticks=element_blank(),text=element_text(size=20),legend.text=element_text(size=20,face = "italic"))

options(repr.plot.width=12, repr.plot.height=6)
DimPlot(pbmc2,  label=T,label.size = 8,group.by=c('typex2','typex')) & theme(legend.key.height=unit(0.3,"inch"),axis.text=element_blank(),axis.ticks=element_blank(),text=element_text(size=20),legend.text=element_text(size=20,face = "italic"))

Idents(pbmc2)='cx'
table(Idents(pbmc2))





pbmc3=readRDS("other/snenEVT-Ashly_P13.rds");pbmc3;names(pbmc3@meta.data)

Idents(pbmc3)="age";table(Idents(pbmc3))

mycol=c(  "#E31A1C","#1F78B4","#A6CEE3","#B2DF8A","#33A02C",
          "#FB9A99","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A",
          "#FFFF99","#B15928","#66C2A5","#FC8D62","#8DA0CB",
          "#B3B3B3","#A6D854","#FFD92F","#E5C494","#E78AC3")
options(repr.plot.width=16, repr.plot.height=3)
barplot(rep(2,20),col=mycol)

Idents(pbmc3)='final_annot_all_troph_corrected'
options(repr.plot.width=7.9, repr.plot.height=6)
DimPlot(object = pbmc3,label = T,raster=FALSE,label.size =8,repel=F)+  theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.35,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20))  +
scale_colour_manual(values=mycol,labels =rownames(table(pbmc3$final_annot_all_troph_corrected)))+guides(col = guide_legend(override.aes=list(size=3),ncol =1)) 

set.seed(123)  # 这里的123是一个示例值，可以是任意整数
Idents(pbmc5) <- "Remained"
Idents(object = pbmc5, cells = sample(colnames(pbmc5), 15000)) <- "Random"
pbmcx<-subset(pbmc5, idents ="Random", invert = F)
options(repr.plot.width=15, repr.plot.height=6)
DimPlot(pbmcx,label = T,label.size = 8)+ theme(legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))



pbmc1;pbmc2;pbmc3





Idents(pbmc1)="stage2";table(Idents(pbmc1))

Idents(pbmc1)="sample";table(Idents(pbmc1))

Idents(pbmc1, cells = WhichCells(pbmc1,idents = c( 'GD8-1', 'GD8-2' ,'GD9-1', 'GD9-2'))) <- "GD8"
Idents(pbmc1, cells = WhichCells(pbmc1,idents = c( 'GD10-1','GD10-2','GD11-1','GD13'))) <- "GD8" 

Idents(pbmc1, cells = WhichCells(pbmc1,idents = c(  'GD8n', 'GD9n-1','GD9n-2','GD10n-1','GD10n-2'))) <- "GD11n" 


pbmc1$samplex=Idents(pbmc1)
table(Idents(pbmc1))





Idents(pbmc2)="sample";table(Idents(pbmc2))

Idents(pbmc2, cells = WhichCells(pbmc2,idents = c('D6' , 'D7'))) <- "D9" 
pbmc2$samplex=Idents(pbmc2)
table(Idents(pbmc2))





Idents(pbmc3)="sample";table(Idents(pbmc3))

#Idents(pbmc3, cells = WhichCells(pbmc3,idents = c('villi_24',"villi_17.6"    ))) <- "sc_24" 

pbmc3$samplex=Idents(pbmc3)
table(Idents(pbmc3))




pbmc1$typex=pbmc1$Typex1
pbmc2$typex=pbmc2$C20A
pbmc3$typex=pbmc3$final_annot_all_troph_corrected


pbmc1$species='Rabbit'
pbmc2$species='Human'
pbmc3$species='Human'

pbmc1$type=pbmc1$type
pbmc2$type='Human-sc'
pbmc3$type='Human-sn'

pbmc1$Species='Rabbit'
pbmc2$Species='Human'
pbmc3$Species='Human'

Idents(pbmc2)="typex";pbmc2$typex=Idents(pbmc2)
Idents(pbmc3)="final_annot_all_troph_corrected";pbmc3$typex=Idents(pbmc3)

table(pbmc1$typex);table(pbmc2$typex);table(pbmc3$typex)#;table(pbmc4$typex);table(pbmc5$final_annot_all_troph_corrected)

table(pbmc1$samplex);table(pbmc2$samplex);table(pbmc3$samplex)#;table(pbmc4$samplex);table(pbmc5$samplex)

pbmc1;pbmc2;pbmc3


pbmc <- merge(x =pbmc1, y =c(pbmc2,pbmc3))
pbmc;table(pbmc$sample);table(Idents(pbmc))

Idents(pbmc)="samplex"
table(Idents(pbmc))

strh=rownames(pbmc1);length(strh)
strm=rownames(pbmc2);length(strm)
strr=rownames(pbmc3);length(strr)
str=Reduce(intersect,list(strh,strm,strr));length(str);head(str)
count=(pbmc$RNA@counts[str,]);dim(count)


pbmc <- CreateSeuratObject(count, meta.data = pbmc@meta.data, min.cells = 3, project = "trophoblast");pbmc

names(pbmc@meta.data);pbmc;table(pbmc$sample,pbmc$species)

Idents(pbmc)="samplex"
table(Idents(pbmc))

names(table(Idents(pbmc)))

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(pbmc))

pancreas=pbmc
pancreas.list <- SplitObject(object = pancreas, split.by = "samplex")

pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, anchor.features = 4000, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(object = pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(object = pancreas.integrated)
pancreas.integrated <- RunPCA(object = pancreas.integrated, npcs = 30, verbose = FALSE)
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(pancreas.integrated);  ElbowPlot(pancreas.integrated)

options(repr.plot.width=10, repr.plot.height=6)
pbmc <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:11)
DimPlot(object = pbmc, label =T,label.size = 8)+ theme(legend.text=element_text(size=16,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=10, repr.plot.height=6)
pbmc <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:8)
DimPlot(object = pbmc, label =T,label.size = 8)+ theme(legend.text=element_text(size=16,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=10, repr.plot.height=6)
pbmc <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:17)
DimPlot(object = pbmc, label =T,label.size = 8)+ theme(legend.text=element_text(size=16,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=16, repr.plot.height=8)
DimPlot(object = pbmc,split.by='species',group.by='typex',label = T,label.size = 8,raster=FALSE)+ theme(legend.key.height=unit(0.3,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))





z=rbind(as.matrix(pbmc1$Typex1),as.matrix(pbmc2$annotation2),as.matrix(pbmc3$final_annot_all_troph_corrected));length(z);head(z)
z=z[order(rownames(z)),];length(z);head(z)
pbmc$typex=z



Idents(pbmc)="typex"
levels(pbmc)=names(table(Idents(pbmc)))[order( as.numeric(substr(names(table(Idents(pbmc))), 1, 2)))]
names(table(Idents(pbmc)))
pbmc$typex=Idents(pbmc)


Idents(pbmc, cells = WhichCells(pbmc,idents = 'SCT')) <- "7   STB "
Idents(pbmc, cells = WhichCells(pbmc,idents = 'GC')) <- '6   GC '
Idents(pbmc, cells = WhichCells(pbmc,idents = 'eEVT')) <- '5   eEVT '
Idents(pbmc, cells = WhichCells(pbmc,idents = c('iEVT','EVT_2','EVT_1'))) <- '4   iEVT '
Idents(pbmc, cells = WhichCells(pbmc,idents = c('VCT_CCC'))) <- "3   Column "
Idents(pbmc, cells = WhichCells(pbmc,idents = c('VCT_fusing'))) <- "2   FC "
Idents(pbmc, cells = WhichCells(pbmc,idents = c('VCT_p','VCT'))) <- "1   CTB "

Idents(pbmc, cells = WhichCells(pbmc,idents = "XTB")) <- "8   XTB "
Idents(pbmc, cells = WhichCells(pbmc,idents = "STB")) <- "7   STB "
Idents(pbmc, cells = WhichCells(pbmc,idents = "EVT")) <- "4   iEVT "
Idents(pbmc, cells = WhichCells(pbmc,idents = "FC")) <- "2   FC "
Idents(pbmc, cells = WhichCells(pbmc,idents = "CTB")) <- "1   CTB "
               

Idents(pbmc, cells = WhichCells(pbmc,idents = "9   XTB")) <- "9   XTB"
Idents(pbmc, cells = WhichCells(pbmc,idents = "8   STB2")) <- "8   JZT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "7   STB1")) <- "7   STB"
Idents(pbmc, cells = WhichCells(pbmc,idents = "6   GlyT3")) <- "6   GlyT3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "5   GlyT2")) <- "5   GlyT2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "4   GlyT1")) <- "4   GlyT1"
Idents(pbmc, cells = WhichCells(pbmc,idents = "3   ColumnT")) <- "3   ColumnT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "2   CTB2")) <- "2   CTB2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "1   CTB1")) <- "1   CTB1"


pbmc$typex2=Idents(pbmc)
table(Idents(pbmc))


options(repr.plot.width=13, repr.plot.height=6)
DimPlot(object = pbmc,split.by='species',group.by='typex2',label = T,label.size = 8,raster=FALSE)+ theme(legend.key.height=unit(0.2,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

levels(pbmc)=c('1   CTB1','2   CTB2','3   ColumnT','4   GlyT1','5   GlyT2','6   GlyT3','7   STB','8   JZT','9   XTB',
               '1   CTB ','2   FC ','3   Column ','4   iEVT ','5   eEVT ','6   GC ','7   STB ','8   XTB '
              )

pbmc$typex2=Idents(pbmc)

Idents(pbmc)=pbmc$typex2

new.cluster.ids <- c("1","2","3","4",'5','6','7','8','9',
                     "1 ","2 ","3 ","4 ",'5 ','6 ','7 ','8 '
                    )
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
pbmc$cx2=Idents(pbmc)
names(table(Idents(pbmc)))

options(repr.plot.width=12, repr.plot.height=6)
DimPlot(object = pbmc,split.by='Species',group.by='typex2',label = T,label.size = 8,raster=FALSE)+ theme(legend.key.height=unit(0.2,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))



Idents(pbmc)="typex"
levels(pbmc)=names(table(Idents(pbmc)))[order( as.numeric(substr(names(table(Idents(pbmc))), 1, 2)))]
names(table(Idents(pbmc)))
pbmc$typex=Idents(pbmc)


Idents(pbmc, cells = WhichCells(pbmc,idents = 'SCT')) <- "4   STB "
Idents(pbmc, cells = WhichCells(pbmc,idents = c('iEVT','EVT_2','EVT_1','eEVT','GC'))) <- '3   EVT '
Idents(pbmc, cells = WhichCells(pbmc,idents = c('VCT_CCC'))) <- "2   ColumnT "
Idents(pbmc, cells = WhichCells(pbmc,idents = c('VCT_p','VCT','VCT_fusing'))) <- "1   CTB "

               
Idents(pbmc, cells = WhichCells(pbmc,idents = "XTB")) <- "5   XTB "
Idents(pbmc, cells = WhichCells(pbmc,idents = "STB")) <- "4   STB "
Idents(pbmc, cells = WhichCells(pbmc,idents = "EVT")) <- "3   EVT "
Idents(pbmc, cells = WhichCells(pbmc,idents = "FC")) <- "1   CTB "
Idents(pbmc, cells = WhichCells(pbmc,idents = "CTB")) <- "1   CTB "


Idents(pbmc, cells = WhichCells(pbmc,idents = "9   XTB")) <- "9   XTB"
Idents(pbmc, cells = WhichCells(pbmc,idents = "8   STB2")) <- "8   JZT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "7   STB1")) <- "7   STB"
Idents(pbmc, cells = WhichCells(pbmc,idents = "6   GlyT3")) <- "6   GlyT3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "5   GlyT2")) <- "5   GlyT2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "4   GlyT1")) <- "4   GlyT1"
Idents(pbmc, cells = WhichCells(pbmc,idents = "3   ColumnT")) <- "3   ColumnT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "2   CTB2")) <- "2   CTB2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "1   CTB1")) <- "1   CTB1"


pbmc$typex3=Idents(pbmc)
table(Idents(pbmc))


options(repr.plot.width=12, repr.plot.height=6)
DimPlot(object = pbmc,split.by='species',group.by='typex3',label = T,label.size = 8,raster=FALSE)+ theme(legend.key.height=unit(0.2,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

levels(pbmc)=c('1   CTB1','2   CTB2','3   ColumnT',
               '4   GlyT1','5   GlyT2','6   GlyT3','7   STB','8   JZT','9   XTB',
               '1   CTB ',"2   ColumnT ",'3   EVT ', '4   STB ','5   XTB '               
              )

pbmc$typex3=Idents(pbmc)

Idents(pbmc)=pbmc$typex3

new.cluster.ids <- c("1","2","3","4",'5','6','7','8','9',
                     "1 ","2 ","3 ","4 ",'5'
                    )
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
pbmc$cx3=Idents(pbmc)
names(table(Idents(pbmc)))

options(repr.plot.width=12, repr.plot.height=6)
DimPlot(object = pbmc,split.by='Species',group.by='typex3',label = T,label.size = 8,raster=FALSE)+ theme(legend.key.height=unit(0.2,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))





Idents(pbmc)="typex"
levels(pbmc)=names(table(Idents(pbmc)))[order( as.numeric(substr(names(table(Idents(pbmc))), 1, 2)))]
names(table(Idents(pbmc)))
pbmc$typex=Idents(pbmc)

Idents(pbmc, cells = WhichCells(pbmc,idents = "XTB")) <- "4   XTB "
Idents(pbmc, cells = WhichCells(pbmc,idents = "STB")) <- "3   STB "
Idents(pbmc, cells = WhichCells(pbmc,idents = "EVT")) <- "2   EVT "
Idents(pbmc, cells = WhichCells(pbmc,idents = "FC")) <- "1   CTB "
Idents(pbmc, cells = WhichCells(pbmc,idents = "CTB")) <- "1   CTB "


Idents(pbmc, cells = WhichCells(pbmc,idents = 'SCT')) <- "3   STB "
Idents(pbmc, cells = WhichCells(pbmc,idents = c('iEVT','EVT_2','EVT_1','eEVT','GC'))) <- '2   EVT '
Idents(pbmc, cells = WhichCells(pbmc,idents = c('VCT_p','VCT','VCT_fusing','VCT_CCC'))) <- "1   CTB "


Idents(pbmc, cells = WhichCells(pbmc,idents = "9   XTB")) <- "9   XTB"
Idents(pbmc, cells = WhichCells(pbmc,idents = "8   STB2")) <- "8   JZT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "7   STB1")) <- "7   STB"
Idents(pbmc, cells = WhichCells(pbmc,idents = "6   GlyT3")) <- "6   GlyT3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "5   GlyT2")) <- "5   GlyT2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "4   GlyT1")) <- "4   GlyT1"
Idents(pbmc, cells = WhichCells(pbmc,idents = "3   ColumnT")) <- "3   ColumnT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "2   CTB2")) <- "2   CTB2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "1   CTB1")) <- "1   CTB1"


pbmc$typex=Idents(pbmc)
table(Idents(pbmc))

options(repr.plot.width=12, repr.plot.height=6)
DimPlot(object = pbmc,split.by='species',group.by='typex',label = T,label.size = 8,raster=FALSE)+ theme(legend.key.height=unit(0.2,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

levels(pbmc)=c('1   CTB1','2   CTB2','3   ColumnT',
               '4   GlyT1','5   GlyT2','6   GlyT3','7   STB','8   JZT','9   XTB',
               
               '1   CTB ','2   EVT ', '3   STB ','4   XTB '
              )

pbmc$typex=Idents(pbmc)

Idents(pbmc)='typex'

new.cluster.ids <- c("1","2","3","4",'5','6','7','8','9',
                     "1 ","2 ","3 ","4 "
                    )
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
pbmc$cx=Idents(pbmc)
names(table(Idents(pbmc)))

options(repr.plot.width=12, repr.plot.height=6)
DimPlot(object = pbmc,split.by='species',group.by='typex',label = T,label.size = 8,raster=FALSE)+ theme(legend.key.height=unit(0.2,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

pbmc

saveRDS(pbmc, file = "inte/GD17-tr93762-hr-snsc-inte.rds",compress=F)

