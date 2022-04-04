##RA-ILD sample re-clustering

library(Seurat)
library(qs)
library(rDBEC)
library(dplyr)
library(ggplot2)
library(data.table)

source("/datadrive/Rhapsody_analysis/Rscripts/library_source_Seurat.R")
mBC = qread("./Seurat/RAILD_combined_Seurat.qs", nthreads = 24)

mBC1 = mBC

#Epithelial sub-clustering

mBC = subset(mBC1, idents=c(8))
mBC = FindVariableFeatures(mBC, selection.method = "mvp", 
                           mean.cutoff=c(0.1, Inf), dispersion.cutoff=c(0.5, Inf))
message("Scaling data...")
mBC = ScaleData(object = mBC, vars.to.regress = c("nReads"))
hvg.number = length(x = mBC@assays$RNA@var.features)
mBC = RunPCA(object = mBC, features = VariableFeatures(object = mBC), 
             npcs=50, verbose = FALSE)


mBC = JackStraw(mBC, num.replicate = 100, dims = 50, verbose = FALSE)
mBC = ScoreJackStraw(mBC, dims = 1:50, score.thresh = 0.05)
#return p-value of Jackstraw analysis
JS_res = mBC@reductions$pca@jackstraw$overall.p.values
tmp = JS_res[JS_res[,2]> 1e-05,,drop=F]
dims =c(1:(min(tmp[,1])-1))

mBC = RunUMAP(mBC, reduction = "pca", dims = dims)
mBC = FindNeighbors(mBC, dims=dims)
mBC = FindClusters(object = mBC, resolution = 2, random.seed=42, verbose = FALSE)

p = FeaturePlot(mBC, features=c("FOXJ1", "SCGB1A1", "MUC5B", "SFTPC", "AGER", "HOPX", "TP63", "KRT17"), pt.size = 0.7)

ggsave(filename="Epi_featureplot_UMAP.png", plot=p, device="png", width=8, height=6, units="in", limitsize=F)

write.table(table(mBC@active.ident, mBC@meta.data$orig.ident), file="Epi_Cellnumber.txt", sep="\t", quote=F, row.names=T, col.names=T)

#find out marker genes
mBC.markers = rDBEC::FindMarkers_parallel_lite(mBC, AllcellsIdent=mBC@active.ident,
                                               test.use="wilcox", only.pos=TRUE, min.pct=0.1,
                                               features.use = NULL, nthreads = 7,
                                               adj.p.val.threshold=0.05)


mBC.markers$cluster = as.numeric(mBC.markers$cluster)
mBC.markers = mBC.markers %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
mBC.markers = as.data.frame(mBC.markers)

quiet(gc())

options(digits=3)
mBC.markers = mBC.markers[,c(10, 9, 8, 2, 3,4, 6,7,1)]
tmp1 = exp(mBC.markers$avg_logFC)
tmp1 = log2(tmp1)
mBC.markers$avg_logFC = tmp1
mBC.markers$avg_logFC = round(mBC.markers$avg_logFC, digits = 3)
mBC.markers$within_avg_exp = round(mBC.markers$within_avg_exp, digits = 3)
mBC.markers$without_avg_exp = round(mBC.markers$without_avg_exp, digits = 3)
mBC.markers$p_val_adj = mBC.markers$p_val_adj
mBC.markers$p_val = mBC.markers$p_val

mBC.markers_write = mBC.markers
mBC.markers_write$cluster=paste("cluster", mBC.markers_write$cluster, sep="")

fwrite(mBC.markers_write, "Epi_ALLmarkers_minpct0.1_Adj_p0.05.txt", row.names=F, col.names=T, sep="\t", quote=F)

new.ident = c("Epi_AT1", "Epi_abberant_basaloid", "Epi_club", "Epi_goblet", "Epi_abberant_basaloid",
              "Epi_ciliated", "Epi_AT2")
names(new.ident)=levels(mBC@active.ident)
mBC = RenameIdents(mBC, new.ident)

#add celltype info to metadata
tmp = as.character(mBC@active.ident)
names(tmp)=rownames(mBC@meta.data)
mBC = AddMetaData(mBC, metadata = tmp, col.name="celltype")

p = DimPlot(mBC, pt.size=1, label=TRUE, label.size = 3) + theme_void() + theme(legend.position = 'none')
ggsave(filename="Epi_UMAP_dimplot.png", plot=p, device="png", width=4, height=3, units="in", limitsize=F)

save(mBC, file="Epi_recluster.rda")

#####################
#Endo sub-clustering

mBC = subset(mBC1, idents=c(6,12,21))
mBC = FindVariableFeatures(mBC, selection.method = "mvp", 
                           mean.cutoff=c(0.1, Inf), dispersion.cutoff=c(0.5, Inf))
message("Scaling data...")
mBC = ScaleData(object = mBC, vars.to.regress = c("nReads"))
hvg.number = length(x = mBC@assays$RNA@var.features)
mBC = RunPCA(object = mBC, features = VariableFeatures(object = mBC), 
             npcs=50, verbose = FALSE)


mBC = JackStraw(mBC, num.replicate = 100, dims = 50, verbose = FALSE)
mBC = ScoreJackStraw(mBC, dims = 1:50, score.thresh = 0.05)
#return p-value of Jackstraw analysis
JS_res = mBC@reductions$pca@jackstraw$overall.p.values
tmp = JS_res[JS_res[,2]> 1e-05,,drop=F]
dims =c(1:(min(tmp[,1])-1))
dims=1:20
mBC = RunUMAP(mBC, reduction = "pca", dims = dims)
mBC = FindNeighbors(mBC, dims=dims)
mBC = FindClusters(object = mBC, resolution = 1.5, random.seed=42, verbose = FALSE)

p = FeaturePlot(mBC, features=c("CA4","TBX2", "PRX", "EDNRB", "COL15A1", "BMX", "HEY1", "NR2F2", "SELP", "ACKR1", "PROX1"),
                pt.size=0.5)
                                
ggsave(filename="Endo_featureplot_UMAP.png", plot=p, device="png", width=10, height=7, units="in", limitsize=F)

write.table(table(mBC@active.ident, mBC@meta.data$orig.ident), file="Endo_Cellnumber.txt", sep="\t", quote=F, row.names=T, col.names=T)

#find out marker genes
mBC.markers = rDBEC::FindMarkers_parallel_lite(mBC, AllcellsIdent=mBC@active.ident,
                                               test.use="wilcox", only.pos=TRUE, min.pct=0.1,
                                               features.use = NULL, nthreads = 7,
                                               adj.p.val.threshold=0.05)


mBC.markers$cluster = as.numeric(mBC.markers$cluster)
mBC.markers = mBC.markers %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
mBC.markers = as.data.frame(mBC.markers)

quiet(gc())

options(digits=3)
mBC.markers = mBC.markers[,c(10, 9, 8, 2, 3,4, 6,7,1)]
tmp1 = exp(mBC.markers$avg_logFC)
tmp1 = log2(tmp1)
mBC.markers$avg_logFC = tmp1
mBC.markers$avg_logFC = round(mBC.markers$avg_logFC, digits = 3)
mBC.markers$within_avg_exp = round(mBC.markers$within_avg_exp, digits = 3)
mBC.markers$without_avg_exp = round(mBC.markers$without_avg_exp, digits = 3)
mBC.markers$p_val_adj = mBC.markers$p_val_adj
mBC.markers$p_val = mBC.markers$p_val

mBC.markers_write = mBC.markers
mBC.markers_write$cluster=paste("cluster", mBC.markers_write$cluster, sep="")

fwrite(mBC.markers_write, "Endo_ALLmarkers_minpct0.1_Adj_p0.05.txt", row.names=F, col.names=T, sep="\t", quote=F)

new.ident = c("Endo_vascular_COL15A1hi", rep("Endo_venous", 4), "Endo_capillary", "Endo_arterial", "Endo_LEC")
names(new.ident)=levels(mBC@active.ident)
mBC = RenameIdents(mBC, new.ident)

#add celltype info to metadata
tmp = as.character(mBC@active.ident)
names(tmp)=rownames(mBC@meta.data)
mBC = AddMetaData(mBC, metadata = tmp, col.name="celltype")

p = DimPlot(mBC, pt.size=1, label=TRUE, label.size = 3) + theme_void() + theme(legend.position = 'none')
ggsave(filename="Endo_UMAP_dimplot.png", plot=p, device="png", width=4, height=3, units="in", limitsize=F)

save(mBC, file="Endo_recluster.rda")


#####################
#Fibro sub-clustering

mBC = subset(mBC1, idents=c(3,20))
mBC = FindVariableFeatures(mBC, selection.method = "mvp", 
                           mean.cutoff=c(0.1, Inf), dispersion.cutoff=c(0.5, Inf))
message("Scaling data...")
mBC = ScaleData(object = mBC, vars.to.regress = c("nReads"))
hvg.number = length(x = mBC@assays$RNA@var.features)
mBC = RunPCA(object = mBC, features = VariableFeatures(object = mBC), 
             npcs=50, verbose = FALSE)


mBC = JackStraw(mBC, num.replicate = 100, dims = 50, verbose = FALSE)
mBC = ScoreJackStraw(mBC, dims = 1:50, score.thresh = 0.05)
#return p-value of Jackstraw analysis
JS_res = mBC@reductions$pca@jackstraw$overall.p.values
tmp = JS_res[JS_res[,2]> 1e-05,,drop=F]
dims =c(1:(min(tmp[,1])-1))

mBC = RunUMAP(mBC, reduction = "pca", dims = dims)
mBC = FindNeighbors(mBC, dims=dims)
mBC = FindClusters(object = mBC, resolution = 1.0, random.seed=42, verbose = FALSE)

p = FeaturePlot(mBC, features=c("NPNT","FGF18", "PI16", "PDGFRA", "CTHRC1", "TCF21"),
                pt.size=0.5)

ggsave(filename="Fibro_featureplot_UMAP.png", plot=p, device="png", width=6, height=6, units="in", limitsize=F)

write.table(table(mBC@active.ident, mBC@meta.data$orig.ident), file="Fibro_Cellnumber.txt", sep="\t", quote=F, row.names=T, col.names=T)

#find out marker genes
mBC.markers = rDBEC::FindMarkers_parallel_lite(mBC, AllcellsIdent=mBC@active.ident,
                                               test.use="wilcox", only.pos=TRUE, min.pct=0.1,
                                               features.use = NULL, nthreads = 7,
                                               adj.p.val.threshold=0.05)


mBC.markers$cluster = as.numeric(mBC.markers$cluster)
mBC.markers = mBC.markers %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
mBC.markers = as.data.frame(mBC.markers)

quiet(gc())

options(digits=3)
mBC.markers = mBC.markers[,c(10, 9, 8, 2, 3,4, 6,7,1)]
tmp1 = exp(mBC.markers$avg_logFC)
tmp1 = log2(tmp1)
mBC.markers$avg_logFC = tmp1
mBC.markers$avg_logFC = round(mBC.markers$avg_logFC, digits = 3)
mBC.markers$within_avg_exp = round(mBC.markers$within_avg_exp, digits = 3)
mBC.markers$without_avg_exp = round(mBC.markers$without_avg_exp, digits = 3)
mBC.markers$p_val_adj = mBC.markers$p_val_adj
mBC.markers$p_val = mBC.markers$p_val

mBC.markers_write = mBC.markers
mBC.markers_write$cluster=paste("cluster", mBC.markers_write$cluster, sep="")

fwrite(mBC.markers_write, "Fibro_ALLmarkers_minpct0.1_Adj_p0.05.txt", row.names=F, col.names=T, sep="\t", quote=F)

new.ident = c("FB_alveolar","FB_adventitial", "FB_misc", "FB_misc", "FB_CTHRC1hi", "FB_alveolar")
names(new.ident)=levels(mBC@active.ident)
mBC = RenameIdents(mBC, new.ident)

#add celltype info to metadata
tmp = as.character(mBC@active.ident)
names(tmp)=rownames(mBC@meta.data)
mBC = AddMetaData(mBC, metadata = tmp, col.name="celltype")

p = DimPlot(mBC, pt.size=1, label=TRUE, label.size = 3) + theme_void() + theme(legend.position = 'none')
ggsave(filename="FB_UMAP_dimplot.png", plot=p, device="png", width=4, height=3, units="in", limitsize=F)

save(mBC, file="FB_recluster.rda")


#######
#Mac DC sub-clustering

mBC = subset(mBC1, idents=c(11, 14, 16, 22, 25))
mBC = FindVariableFeatures(mBC, selection.method = "mvp", 
                           mean.cutoff=c(0.1, Inf), dispersion.cutoff=c(0.5, Inf))
message("Scaling data...")
mBC = ScaleData(object = mBC, vars.to.regress = c("nReads"))
hvg.number = length(x = mBC@assays$RNA@var.features)
mBC = RunPCA(object = mBC, features = VariableFeatures(object = mBC), 
             npcs=50, verbose = FALSE)


mBC = JackStraw(mBC, num.replicate = 100, dims = 50, verbose = FALSE)
mBC = ScoreJackStraw(mBC, dims = 1:50, score.thresh = 0.05)
#return p-value of Jackstraw analysis
JS_res = mBC@reductions$pca@jackstraw$overall.p.values
tmp = JS_res[JS_res[,2]> 1e-05,,drop=F]
dims =c(1:(min(tmp[,1])-1))

dims=1:15

mBC = RunUMAP(mBC, reduction = "pca", dims = dims)
mBC = FindNeighbors(mBC, dims=dims)
mBC = FindClusters(object = mBC, resolution = 1.5, random.seed=42, verbose = FALSE)

p = FeaturePlot(mBC, features=c("MARCO","MSR1", "MRC1", "SPP1", "CD1C", "FABP4", "FOLR2", "FCN1", "TREM2", "C1QA"),
                pt.size=0.5)

ggsave(filename="MacDC_featureplot_UMAP.png", plot=p, device="png", width=10, height=7, units="in", limitsize=F)

write.table(table(mBC@active.ident, mBC@meta.data$orig.ident), file="Mac_Cellnumber.txt", sep="\t", quote=F, row.names=T, col.names=T)

#find out marker genes
mBC.markers = rDBEC::FindMarkers_parallel_lite(mBC, AllcellsIdent=mBC@active.ident,
                                               test.use="wilcox", only.pos=TRUE, min.pct=0.1,
                                               features.use = NULL, nthreads = 7,
                                               adj.p.val.threshold=0.05)


mBC.markers$cluster = as.numeric(mBC.markers$cluster)
mBC.markers = mBC.markers %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
mBC.markers = as.data.frame(mBC.markers)

quiet(gc())

options(digits=3)
mBC.markers = mBC.markers[,c(10, 9, 8, 2, 3,4, 6,7,1)]
tmp1 = exp(mBC.markers$avg_logFC)
tmp1 = log2(tmp1)
mBC.markers$avg_logFC = tmp1
mBC.markers$avg_logFC = round(mBC.markers$avg_logFC, digits = 3)
mBC.markers$within_avg_exp = round(mBC.markers$within_avg_exp, digits = 3)
mBC.markers$without_avg_exp = round(mBC.markers$without_avg_exp, digits = 3)
mBC.markers$p_val_adj = mBC.markers$p_val_adj
mBC.markers$p_val = mBC.markers$p_val

mBC.markers_write = mBC.markers
mBC.markers_write$cluster=paste("cluster", mBC.markers_write$cluster, sep="")

fwrite(mBC.markers_write, "Mac_ALLmarkers_minpct0.1_Adj_p0.05.txt", row.names=F, col.names=T, sep="\t", quote=F)

new.ident = c("Mac_C1QAhi","Mac_C1QAhi","Mac_SPP1hi", "Mo_classical", "Mac_C1QAhi", "Mac_AM", "DC_cDC2")
names(new.ident)=levels(mBC@active.ident)
mBC = RenameIdents(mBC, new.ident)

#add celltype info to metadata
tmp = as.character(mBC@active.ident)
names(tmp)=rownames(mBC@meta.data)
mBC = AddMetaData(mBC, metadata = tmp, col.name="celltype")

p = DimPlot(mBC, pt.size=1, label=TRUE, label.size = 3) + theme_void() + theme(legend.position = 'none')
ggsave(filename="Mac_UMAP_dimplot.png", plot=p, device="png", width=4, height=3, units="in", limitsize=F)

save(mBC, file="Mac_recluster.rda")

