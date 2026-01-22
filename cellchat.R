getwd();library(Seurat);library(Matrix);suppressMessages(library(dplyr));library(ggplot2);suppressMessages(library('pagoda2'));suppressMessages(library("velocyto.R"));suppressMessages(library(monocle));
suppressMessages(library(cowplot));suppressMessages(library("patchwork"));suppressMessages(library(harmony))
options(repr.plot.width=7, repr.plot.height=6); library("RColorBrewer");library(magrittr);library(scales)
setwd("/sdc/xxjiang/rapl")
getwd()
packageVersion('Seurat')
packageVersion('patchwork')

library(CellChat)
library(ggplot2)
library(ggalluvial)
library(NMF)
library(circlize)

pbmc=readRDS("data/rapl169703-CDB.rds");pbmc;names(pbmc@meta.data)

Normal.cellchat=pbmc

Normal.cellchat$active.ident <- Normal.cellchat$Type5

Normal.cellchat@active.ident <- factor(Normal.cellchat$Type5)

normal.data.input <- Normal.cellchat@assays$RNA@data

normal.meta <- Normal.cellchat@meta.data

dim(normal.data.input);head(normal.meta)[,1:4]

Normal.cellchat <- createCellChat(object = Normal.cellchat, meta = normal.meta, group.by = "active.ident")

Normal.cellchat <- addMeta(Normal.cellchat, meta = normal.meta)

Normal.cellchat

Normal.cellchat <- setIdent(Normal.cellchat, ident.use = "active.ident")

levels(Normal.cellchat@idents)

#set the ligand-receptor interaction database
Normal.cellchatDB <- CellChatDB.human

head(Normal.cellchatDB$interaction)

dplyr::glimpse(Normal.cellchatDB$interaction)

Normal.cellchatDB$interaction$annotation

Normal.cellchat@DB <- Normal.cellchatDB

#processing the expression data for cell-cell communication
Normal.cellchat <- subsetData(Normal.cellchat)

Normal.cellchat <- identifyOverExpressedGenes(Normal.cellchat)

Normal.cellchat <- identifyOverExpressedInteractions(Normal.cellchat)

Normal.cellchat <- computeCommunProb(Normal.cellchat)

Normal.cellchat <- computeCommunProbPathway(Normal.cellchat)

Normal.cellchat <- aggregateNet(Normal.cellchat)


Normal.cellchat@netP$pathways


groupSize <- as.numeric(table(Normal.cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(Normal.cellchat@net$count,
                 vertex.weight=groupSize,
                 weight.scale=T,
                 label.edge=F,                  
                 title.name="Number_of_interactions"
                )
netVisual_circle(Normal.cellchat@net$weight,                  
                 vertex.weight = groupSize,                  
                 weight.scale = TRUE, 
                  label.edge = FALSE,                  
                 title.name = "Interaction weights/strength"
                )


options(repr.plot.width=12, repr.plot.height=6)
netVisual_aggregate(
  Normal.cellchat,
  signaling = "WNT", #pathway_name
  vertex.receiver = c(1,2,3), 
  layout = "circle"
)

