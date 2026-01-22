getwd();suppressMessages(library(Seurat));library(Matrix);suppressMessages(library(dplyr));library(ggplot2);suppressMessages(library(monocle));suppressMessages(library(SeuratDisk))
suppressMessages(library(cowplot));suppressMessages(library('pagoda2'));suppressMessages(library("velocyto.R"));suppressMessages(library("patchwork"))
library("RColorBrewer");library(scales);library(SeuratWrappers);library(argparser);suppressMessages(library(harmony))
suppressMessages(library(monocle3));suppressMessages(library(clusterProfiler));suppressMessages(library(org.Hs.eg.db));getwd();.libPaths();packageVersion('Seurat')
options(repr.plot.width=7, repr.plot.height=6); 


pbmc=readRDS("tr79843m-inte.rds");pbmc;names(pbmc@meta.data)

DefaultAssay(pbmc)='RNA';pbmc

Idents(pbmc)="CX1"
options(repr.plot.width=7.9, repr.plot.height=6)
DimPlot(pbmc,label = T,label.size = 8, ncol=1,raster=FALSE,repel=F)+ 
theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=22,face="plain"),legend.key.height=unit(0.45,"inch"),text=element_text(size=22,face="plain"),strip.text=element_text(size=22)) + 
scale_colour_manual(values=x[as.numeric(names(table(Idents(pbmc))))],labels =rownames(table(pbmc$Typex1))) +guides(col = guide_legend(override.aes=list(size=3),ncol =1)) 






library("monocle3")
library("ggplot2")
library("dplyr")

Idents(pbmc)="CX1"
options(repr.plot.width=6.7, repr.plot.height=6)
DimPlot(object = pbmc,label = T,label.size =8,repel=T)+  theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),legend.key.height=unit(0.4,"inch"),text=element_text(size=20,face="plain"),strip.text=element_text(size=20))  +
scale_colour_manual(values=x[as.numeric(names(table(Idents(pbmc))))],labels =rownames(table(Idents(pbmc)))) 



cds2 <- as.cell_data_set(pbmc)

options(repr.plot.width=6, repr.plot.height=6)
cds <- cluster_cells(cds2,k=100,resolution=0.0001)
plot_cells(cds, color_cells_by = "partition", group_label_size = 8)

cds <- learn_graph(cds)

cds@clusters@listData$UMAP$clusters=pbmc$CX1

options(repr.plot.width=6.05, repr.plot.height=6)
plot_cells(cds,label_groups_by_cluster = FALSE, label_leaves = FALSE,group_label_size = 8,cell_size = 1, label_branch_points = FALSE)+ theme( axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))+
scale_colour_manual(values=x[as.numeric(names(table(Idents(pbmc))))],labels =rownames(table(Idents(pbmc)))) 

Idents(pbmc)="CX1"
options(repr.plot.width=6.05, repr.plot.height=6)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE,group_label_size = 8,cell_size = 0.6,  graph_label_size = 0,
           label_branch_points = FALSE)+ theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=22,face="plain"),text=element_text(size=22,face="plain"))+
scale_colour_manual(values=x[as.numeric(names(table(Idents(pbmc))))],labels =rownames(table(Idents(pbmc))))


options(repr.plot.width=18, repr.plot.height=12)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE,group_label_size = 4,graph_label_size = 6,cell_size = 1, label_principal_points = T,label_branch_points =F)+ theme( axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

cds <- order_cells(cds, root_pr_nodes="Y_5")

options(repr.plot.width=7.75, repr.plot.height=6)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, cell_size = 0.5,  graph_label_size = 4,
    label_branch_points = FALSE)+ theme(legend.key.height=unit(0.5,"inch"), axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

options(repr.plot.width=7.75, repr.plot.height=6)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, cell_size = 0.3,  graph_label_size = 5,
    label_branch_points = FALSE)+ theme(axis.line=element_blank(),axis.title=element_blank(),legend.key.height=unit(0.5,"inch"), axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=22,face="plain"),text=element_text(size=22,face="plain"))#+
#scale_color_gradient2(low = "red", mid = "white", high = "green", midpoint = median(pseudotime(cds)))


Idents(pbmc)="CX1"
options(repr.plot.width=6.05, repr.plot.height=6)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE,group_label_size = 8,cell_size = 1,  graph_label_size = 4,
           label_branch_points = FALSE)+ theme( axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))+
scale_colour_manual(values=x[as.numeric(names(table(Idents(pbmc))))],labels =rownames(table(Idents(pbmc))))

saveRDS(cds, file = "mo/tr79843-mo3-CTB-STB1-CX1.rds",compress=F)

cds

pbmc

modulated_genes = graph_test(cds,                             cores=20)

modulated_genes

genes <- row.names(subset(modulated_genes, q_value== 0 & morans_I >0.25) )


library(ClusterGVis);library(ComplexHeatmap)

mat <- pre_pseudotime_matrix(
  cds_obj = cds,
  gene_list = genes
)

head(mat[1:6,1:6])

#options(repr.plot.width=6, repr.plot.height=6)
getClusters(mat)

options(repr.plot.width=6, repr.plot.height=8)

set.seed(123)  
ck <- clusterData(
  mat,
  cluster.method = "kmeans",
  cluster.num = 2    
)

table(ck$long.res$cluster)

options(repr.plot.width=10, repr.plot.height=8)
visCluster(object = ck,
       plot.type = "line") 

marker_genes <- intersect(rownames(mat),TF)
marker_genes

options(repr.plot.width=4.5, repr.plot.height=6)
visCluster(ck,plot.type='heatmap',
           add.sampleanno=F,
           markGenes=marker_genes,
           ht.col.list = list(col_range = c(-2, 0, 2), col_color = c("#00008B", "white",
    "orange")),
           cluster.order = c(1,2),
           ctAnno.col = NULL,
           textbox.size = 5,  
           column_names_rot = 45, 
           annnoblock.text = F,
           border=F,
           #line.size=0,
        #   add.line = FALSE,
           panel.arg = c(2, 0, 4, NA, NA),
           show_row_dend = F,  
      #     show_row_names=F,
           genes.gp = c('italic',fontsize = 15,color='black'))
