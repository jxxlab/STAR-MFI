getwd();library(Seurat);library(Matrix);suppressMessages(library(dplyr));library(ggplot2);suppressMessages(library(monocle));
suppressMessages(library(cowplot));suppressMessages(library('pagoda2'));suppressMessages(library("velocyto.R"));suppressMessages(library("patchwork"))
options(repr.plot.width=7, repr.plot.height=6); library("RColorBrewer");library(scales)
setwd("/sdc/xxjiang/rapl/")
getwd()



pbmc1=readRDS("/sdc/xxjiang/rapl-sp/rapl-sn-148060.rds");pbmc1;names(pbmc1@meta.data)

table(Idents(pbmc1),pbmc1$C16H)

Idents(pbmc1)='C16H'

pbmc=pbmc1
options(repr.plot.width=7, repr.plot.height=6)
DimPlot(object = pbmc,label = T,raster=FALSE,label.size =8,repel=F,reduction='umap')+  theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.35,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20))  +
scale_colour_manual(values=c(colorRampPalette(brewer.pal(12, "Paired"))(23)[c(1,3,5,7,8,9,20,21, 16,18,  11,19,22,23,13,21)],"black","green"))

#options(repr.plot.width=20, repr.plot.height=8)
#DimPlot(object = pbmc,label =F,label.size =8,raster=FALSE,repel=F,reduction='umap')+  theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.35,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20))  +
#scale_colour_manual(values=c(colorRampPalette(brewer.pal(12, "Paired"))(23)[c(1,3,5,7,8,9,20,21, 16,18,  11,19,22,23,13,21)],"black","green"))

pbmc1$type='SN'



pbmc2=readRDS("/sdc/xxjiang/rapl-sp/rapl-sc-114572.rds");   pbmc2;  names(pbmc2@meta.data)

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

pbmcx=pbmc2
Idents(pbmcx)="CX2"
options(repr.plot.width=20, repr.plot.height=8)
DimPlot(pbmcx,label = T,label.size = 6, raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.35,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=y[as.numeric(names(table(Idents(pbmcx))))],labels =rownames(table(pbmcx$Type2))) +guides(col = guide_legend(override.aes=list(size=3),ncol =4)) 
Idents(pbmcx)="CX1"
options(repr.plot.width=12, repr.plot.height=6)
DimPlot(pbmcx,label = T,label.size = 6, ncol=1,raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=z[as.numeric(names(table(Idents(pbmcx))))],labels =rownames(table(pbmcx$Type1))) +guides(col = guide_legend(override.aes=list(size=3),ncol =2)) 




other=pbmc1;pbmc1
HS=pbmc2;pbmc2



anchor1 <- FindTransferAnchors(
  reference= HS,
  query = other,
  normalization.method = "LogNormalize",
  reference.reduction='pca',
  dims = 1:30)



##saveRDS(anchor1, file = "sn148060_to_sc114572_anchor.rds")
anchor1=readRDS("sn148060_to_sc114572_anchor.rds");anchor1



other2 <- MapQuery(
  anchorset = anchor1,
  query = other,
  reference = HS,
  refdata = list(
    celltype1 = "CX2",
    celltype2 = "Type2",
    celltype3 = "CX1",
    celltype4 = "Type1"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

table(other2$predicted.celltype1)

options(repr.plot.width=35, repr.plot.height=8)
p0=DimPlot(other2, reduction = "ref.umap", cols=y[as.numeric(names(table(other2$predicted.celltype1)))],raster=FALSE,group.by = "predicted.celltype2", label = TRUE, label.size = 8,repel=F)+ theme(axis.title=element_blank(),axis.line=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.4,"inch"),text=element_text(size=20,face="plain"))  #+ NoLegend()+labs(title="")
p1=DimPlot(other2, reduction = "ref.umap", cols=z[as.numeric(names(table(other2$predicted.celltype1)))],raster=FALSE,group.by = "predicted.celltype4", label = TRUE, label.size = 6,repel=F)+ theme(axis.title=element_blank(),axis.line=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.4,"inch"),text=element_text(size=20,face="plain"))  #+ NoLegend()+labs(title="")
p2=DimPlot(other2, reduction = "umap",                                   raster=FALSE,group.by ="predicted.celltype2",     label = T,    label.size = 8,repel=F)+ theme(axis.title=element_blank(),axis.line=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.4,"inch"),text=element_text(size=20,face="plain")) 
p0|p1|p2

Idents(other2)='predicted.celltype1'
table(Idents(other2))
other2=subset(other2,idents=names(table(Idents(other2)))[as.numeric(table(Idents(other2)))>=5]);other2
levels(other2)=1:78
table(Idents(other2))
other2$CX2=Idents(other2)
Idents(other2)='predicted.celltype3'
other2=subset(other2,idents=names(table(Idents(other2)))[as.numeric(table(Idents(other2)))>=5]);other2
levels(other2)=1:38
table(Idents(other2))
other2$CX1=Idents(other2)



Idents(other2)="predicted.celltype2"
levels(other2)=names(table(Idents(other2)))[order( as.numeric(substr(names(table(Idents(other2))), 1, 2)))]
names(table(Idents(other2)))
other2$Type2=Idents(other2)
Idents(other2)="predicted.celltype4"
levels(other2)=names(table(Idents(other2)))[order( as.numeric(substr(names(table(Idents(other2))), 1, 2)))]
names(table(Idents(other2)))
other2$Type1=Idents(other2)

pbmcx=pbmc2
Idents(pbmcx)="CX2"
options(repr.plot.width=20, repr.plot.height=8)
DimPlot(pbmcx,label = T,label.size = 6, raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.35,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=y[as.numeric(names(table(Idents(pbmcx))))],labels =rownames(table(pbmcx$Type2))) +guides(col = guide_legend(override.aes=list(size=3),ncol =4)) 
Idents(pbmcx)="CX1"
options(repr.plot.width=12, repr.plot.height=6)
DimPlot(pbmcx,label = T,label.size = 6, ncol=1,raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=z[as.numeric(names(table(Idents(pbmcx))))],labels =rownames(table(pbmcx$Type1))) +guides(col = guide_legend(override.aes=list(size=3),ncol =2)) 


Idents(other2)="CX2"
options(repr.plot.width=12.6, repr.plot.height=8)
DimPlot(other2,label = F,label.size = 6, ncol=1,raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=y[as.numeric(names(table(Idents(other2))))],labels =rownames(table(other2$Type2))) +guides(col = guide_legend(override.aes=list(size=3),ncol =2)) 


Idents(other2)="CX1"
options(repr.plot.width=9.6, repr.plot.height=8)
DimPlot(other2,label = F,label.size = 6, ncol=1,raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.28,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=z[as.numeric(names(table(Idents(other2))))],labels =rownames(table(other2$Type1))) +guides(col = guide_legend(override.aes=list(size=3),ncol =1)) 



Idents(other2)="stage"
options(repr.plot.width=9.6, repr.plot.height=8)
DimPlot(other2, label=T,label.size = 8 ,raster=FALSE)+ theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.35,"inch"),text=element_text(size=20,face="plain"))

Idents(other2)="CX1"
options(repr.plot.width=16, repr.plot.height=14)
DimPlot(other2,label = T,label.size = 6,  raster=FALSE,repel=F,split.by="sample",ncol=4)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.40,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=z[as.numeric(names(table(Idents(other2))))],labels =rownames(table(other2$Type1))) +guides(col = guide_legend(override.aes=list(size=3),ncol =1)) 


Idents(other2)="CX2"
options(repr.plot.width=16, repr.plot.height=14)
DimPlot(other2,label = T,label.size = 6,  raster=FALSE,repel=F,split.by="sample",ncol=4)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.36,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=y[as.numeric(names(table(Idents(other2))))],labels =rownames(table(other2$Type2))) +guides(col = guide_legend(override.aes=list(size=3),ncol =1)) 

pbmc=pbmc3=other2;other2;names(other2@meta.data)

Idents(pbmc)='CX2'
Idents(pbmc3)='C16H'

#Idents(pbmc3)="CX2"
options(repr.plot.width=14, repr.plot.height=6)
DimPlot(pbmc3,label = T,reduction="umap",label.size = 6, ncol=1,raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=y[as.numeric(names(table(Idents(pbmc3))))],labels =rownames(table(pbmc3$Type2))) +guides(col = guide_legend(override.aes=list(size=3),ncol =3)) 


#Idents(pbmc3)="CX2"
options(repr.plot.width=14, repr.plot.height=6)
DimPlot(pbmc,label = T,reduction="umap",label.size = 6, ncol=1,raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=y[as.numeric(names(table(Idents(pbmc))))],labels =rownames(table(pbmc$Type2))) +guides(col = guide_legend(override.aes=list(size=3),ncol =3)) 


Idents(pbmc, cells = WhichCells(pbmc3,idents =c(16,14)))<-"53"

pbmc$CX3=Idents(pbmc)
table(Idents(pbmc))

Idents(pbmc3)='C39H'
Idents(pbmc, cells = WhichCells(pbmc3,idents =c(9)))<-"54"

pbmc$CX3=Idents(pbmc)
table(Idents(pbmc))

#Idents(pbmc3)="CX2"
options(repr.plot.width=14, repr.plot.height=6)
DimPlot(pbmc,label = T,reduction="umap",label.size = 6, ncol=1,raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=y[as.numeric(names(table(Idents(pbmc))))],labels =rownames(table(pbmc$Type2))) +guides(col = guide_legend(override.aes=list(size=3),ncol =3)) 


Idents(pbmc3)='C16H'

Idents(pbmc)='Type2';
Idents(pbmc, cells = WhichCells(pbmc3,idents =c(16,14)))<-"53 JZ-C6"
pbmc$Type3=Idents(pbmc)
table(Idents(pbmc))

Idents(pbmc3)='C39H'
Idents(pbmc, cells = WhichCells(pbmc3,idents =c(9)))<-"54 Unk."

pbmc$Type3=Idents(pbmc)
table(Idents(pbmc))

Idents(pbmc)="CX3"
options(repr.plot.width=12.6, repr.plot.height=6)
DimPlot(pbmc,label = T,reduction="umap",label.size = 6, ncol=1,raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=y[as.numeric(names(table(Idents(pbmc))))],labels =rownames(table(pbmc$Type3))) +guides(col = guide_legend(override.aes=list(size=3),ncol =3)) 


other2=pbmc

Idents(other2)="CX3"
levels(other2)=names(table(Idents(other2)))[order( as.numeric(substr(names(table(Idents(other2))), 1, 2)))]
names(table(Idents(other2)))
other2$CX3=Idents(other2)
Idents(other2)="Type3"
levels(other2)=names(table(Idents(other2)))[order( as.numeric(substr(names(table(Idents(other2))), 1, 2)))]
names(table(Idents(other2)))
other2$Type3=Idents(other2)


Idents(other2)="CX1"
options(repr.plot.width=9, repr.plot.height=6)
DimPlot(other2,label = T,reduction="umap",label.size = 6, ncol=1,raster=FALSE,repel=T)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=z[as.numeric(names(table(Idents(other2))))],labels =rownames(table(other2$Type1))) +guides(col = guide_legend(override.aes=list(size=3),ncol =2)) 


Idents(other2)="CX3"
options(repr.plot.width=12.6, repr.plot.height=6)
DimPlot(other2,label = T,reduction="umap",label.size = 6, ncol=1,raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=y[as.numeric(names(table(Idents(other2))))],labels =rownames(table(other2$Type3))) +guides(col = guide_legend(override.aes=list(size=3),ncol =3)) 


Idents(other2)="CX1"
options(repr.plot.width=11, repr.plot.height=6)
DimPlot(other2,label = T,reduction="ref.umap",label.size = 6, ncol=1,raster=FALSE,repel=T)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=z[as.numeric(names(table(Idents(other2))))],labels =rownames(table(other2$Type1))) +guides(col = guide_legend(override.aes=list(size=3),ncol =2)) 


Idents(other2)="CX3"
options(repr.plot.width=14, repr.plot.height=6)
DimPlot(other2,label = T,reduction="ref.umap",label.size = 6, ncol=1,raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.30,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20)) + 
scale_colour_manual(values=y[as.numeric(names(table(Idents(other2))))],labels =rownames(table(other2$Type3))) +guides(col = guide_legend(override.aes=list(size=3),ncol =3)) 


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(other2,cells.highlight = CellsByIdentities(object = other2, idents =c(53)),raster=FALSE,reduction="umap", label=F,label.size = 8) & theme(legend.key.height=unit(0.3,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))

options(repr.plot.width=8, repr.plot.height=6)
DimPlot(other2,cells.highlight = CellsByIdentities(object = other2, idents =c(20)),raster=FALSE,reduction="umap", label=F,label.size = 8) & theme(legend.key.height=unit(0.3,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))

options(repr.plot.width=8, repr.plot.height=6)
DimPlot(other2,cells.highlight = CellsByIdentities(object = other2, idents =c(4)),raster=FALSE,reduction="umap", label=F,label.size = 8) & theme(legend.key.height=unit(0.3,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))

options(repr.plot.width=9.6, repr.plot.height=8)
DimPlot(other2,cells.highlight = CellsByIdentities(object = other2, idents =c(47:48)), raster=FALSE,cols.highlight =rev(y[47:48]),pt.size=0.00001,sizes.highlight=0.5,reduction="umap", label=F,label.size = 8) & theme(legend.key.height=unit(0.3,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))

##### options(repr.plot.width=15, repr.plot.height=8)
DimPlot(other2,cells.highlight = CellsByIdentities(object = other2, idents =c(49:50)), raster=FALSE,cols.highlight =rev(y[49:50]),pt.size=0.00001,sizes.highlight=0.5,reduction="umap", label=F,label.size = 8) & theme(legend.key.height=unit(0.3,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))

##### options(repr.plot.width=15, repr.plot.height=8)
DimPlot(other2,cells.highlight = CellsByIdentities(object = other2, idents =c(51:52)), raster=FALSE,cols.highlight =rev(y[51:52]),pt.size=0.00001,sizes.highlight=0.5,reduction="umap", label=F,label.size = 8) & theme(legend.key.height=unit(0.3,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))

#options(repr.plot.width=15, repr.plot.height=8)
DimPlot(other2,cells.highlight = CellsByIdentities(object = other2, idents =c(53:54)),raster=FALSE, cols.highlight =rev(y[53:54]),pt.size=0.00001,sizes.highlight=0.5,reduction="umap", label=F,label.size = 8) & theme(legend.key.height=unit(0.3,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))

other2

names(other2@meta.data)

saveRDS(other2, file = "sn148053-sc114572-2nd.rds")

getwd()

















other2=readRDS("sn148053-sc114572-2nd.rds");other2;names(other2@meta.data)

pbmc1=other2;other2
median(pbmc1$nFeature_RNA)
options(repr.plot.width=6, repr.plot.height=6)
 VlnPlot(pbmc1,pt.size=0.0001, group.by = "Section",features = c("nFeature_RNA"),raster=FALSE) + ggtitle("nFeature_RNA / Number of genes")+ theme(axis.title=element_blank(),axis.text=element_text(size=20),plot.title=element_text(size=20,face='bold'))+NoLegend()
options(repr.plot.width=6, repr.plot.height=6)
 VlnPlot(pbmc1,pt.size=0, group.by = "Section",features = c("nFeature_RNA"),raster=FALSE) + ggtitle("nFeature_RNA / Number of genes")+ theme(axis.title=element_blank(),axis.text=element_text(size=20),plot.title=element_text(size=20,face='bold'))+NoLegend()
options(repr.plot.width=6, repr.plot.height=6)
 RidgePlot(pbmc1, group.by = "Section",features = c("nFeature_RNA")) + ggtitle("nFeature_RNA / Number of genes")+ theme(axis.title=element_blank(),axis.text=element_text(size=20),plot.title=element_text(size=20,face='bold'))+NoLegend()


