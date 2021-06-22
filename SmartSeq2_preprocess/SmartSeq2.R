#re-analysis of Smart-seq2 data of murine lung of Tabula Muris data
options(stringsAsFactors = FALSE)

#please install Seurat v2.3.4. and rDBEC package before use (could be downloaded from this repository)
suppressWarnings(suppressMessages(source("library_source_Seurat.R")))

tableread_fast_csv = function(i){
  tmp = fread(i, header=TRUE, sep=",")
  tmp = as.data.frame(tmp)
  rownames(tmp) = tmp[,1]
  tmp = tmp[,2:ncol(tmp)]
  return(tmp)
}

##normal mouse Lung Smart-seq2 analysis (Tabula Muris data)
fileUrl <- "https://ndownloader.figshare.com/files/10700143"
download.file(fileUrl, destfile = "./FACS.zip", method = "wget")
unzip("FACS.zip")

##normal mouse Lung Smart-seq2 analysis (Tabula Muris data)

fnames = dir(pattern = "./FACS/Lung-counts.csv")
Lung_SMARTSeq2 = tableread_fast_csv(fnames)
Lung_SMARTSeq2 = Lung_SMARTSeq2[rowSums(Lung_SMARTSeq2)>0,]

#filter out ERCC data
ERCC = grep(pattern = "^ERCC-", x = rownames(Lung_SMARTSeq2), value = TRUE)
Lung_SMARTSeq2 = Lung_SMARTSeq2[!rownames(Lung_SMARTSeq2) %in% ERCC,]
colnames(Lung_SMARTSeq2) = paste("Smart-seq2", "_", c(1:ncol(Lung_SMARTSeq2)), sep='')

#Seurat analysis
Lung_SMARTSeq2 = CreateSeuratObject(raw.data = Lung_SMARTSeq2, min.cells = 5, min.genes = 500)

colnames(Lung_SMARTSeq2@meta.data)[2]="nReads"

ribo.genes = grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = Lung_SMARTSeq2@data), value = TRUE)
percent.ribo = Matrix::colSums(Lung_SMARTSeq2@raw.data[ribo.genes, ])/Matrix::colSums(Lung_SMARTSeq2@raw.data)
Lung_SMARTSeq2 = AddMetaData(object = Lung_SMARTSeq2, metadata = percent.ribo, col.name = "percent.ribo")

percent.RiboRNA = Matrix::colSums(Lung_SMARTSeq2@raw.data[c('Rn45s'), ])/Matrix::colSums(Lung_SMARTSeq2@raw.data)
Lung_SMARTSeq2 = AddMetaData(object = Lung_SMARTSeq2, metadata = percent.RiboRNA, col.name = "percent.RiboRNA")

nReads_log = log10(Matrix::colSums(Lung_SMARTSeq2@raw.data))
Lung_SMARTSeq2 = AddMetaData(object = Lung_SMARTSeq2, metadata = nReads_log, col.name = "nReads.log")

# filter cells
Lung_SMARTSeq2 = FilterCells(object = Lung_SMARTSeq2, subset.names = c("nGene", "nReads"),
                             low.thresholds = c(500, 10000), high.thresholds = c(Inf, Inf))

# normalizing data
Lung_SMARTSeq2 = NormalizeData(object = Lung_SMARTSeq2, scale.factor=1000000)
Lung_SMARTSeq2 = ScaleData(object = Lung_SMARTSeq2, vars.to.regress = c("nReads"), do.par=TRUE, num.cores=16)
Lung_SMARTSeq2 = FindVariableGenes(object = Lung_SMARTSeq2, mean.function = ExpMean, dispersion.function = LogVMR, 
                                   x.low.cutoff = 0.1, x.high.cutoff = Inf, y.cutoff = 0.5, do.plot=FALSE)
length(x = Lung_SMARTSeq2@var.genes)

# perform PCA
Lung_SMARTSeq2 = RunPCA(object = Lung_SMARTSeq2, pc.genes = Lung_SMARTSeq2@var.genes, do.print = FALSE, pcs.compute = 100, overwrite=TRUE)
Lung_SMARTSeq2 = ProjectPCA(object = Lung_SMARTSeq2, do.print = FALSE)
Lung_SMARTSeq2 = JackStraw(object = Lung_SMARTSeq2, num.replicate = 100, num.pc = 100,do.par=TRUE, num.cores=16)
Lung_SMARTSeq2 = JackStrawPlot(object = Lung_SMARTSeq2, PCs = 1:100)

#clustering and detection of cellular subsets
dims=Lung_SMARTSeq2@dr$pca@jackstraw@overall.p.values[Lung_SMARTSeq2@dr$pca@jackstraw@overall.p.values[,2]>0.05,1]
dims=c(1:(min(dims)-1))

Lung_SMARTSeq2 = FindClusters(object = Lung_SMARTSeq2, reduction.type = "pca", dims.use = dims,
                              resolution = 2, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
py_set_seed(42)
Lung_SMARTSeq2 = DoFItSNE(Lung_SMARTSeq2 , reduction_use = "pca", dims_use = as.integer(dims),
                          perplexity=100, n_jobs=as.integer(16), df = 0.9)

p1 = DimPlot(object = Lung_SMARTSeq2, reduction.use = "FItSNE", group.by = "orig.ident", 
             do.return = TRUE, vector.friendly = TRUE, pt.size = 1.0) +
  theme(axis.title.x = element_text(size=10, family = "Arial"), 
        axis.title.y = element_text(size=10, family = "Arial"), 
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"), 
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial")) +
  theme(panel.border = element_rect(fill = NA, size = 1.0))  
p2 = DimPlot(object = Lung_SMARTSeq2, reduction.use = "FItSNE", do.label = TRUE, label.size = 10,
             do.return = TRUE, vector.friendly = TRUE, pt.size = 1.0) +
  theme(axis.title.x = element_text(size=10, family = "Arial"), 
        axis.title.y = element_text(size=10, family = "Arial"), 
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"), 
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial")) +
  theme(panel.border = element_rect(fill = NA, size = 1.0)) 

ggsave(file = "FItSNE_plot_Lung_SmartSeq2.png", plot = plot_grid(p1, p2, ncol = 2), 
       device="png", units="in", dpi = 300, width = 8, height = 4, limitsize=FALSE)

#detection of subset-marker genes

tmp = table(Lung_SMARTSeq2@ident)
tmp = as.data.frame(tmp)
write.table(tmp, "subsetNumber.txt", row.names=F, sep="\t", quote=F)
Lung_SMARTSeq2.markers = rDBEC::FindMarkers_parallel_lite(Lung_SMARTSeq2, AllcellsIdent=Lung_SMARTSeq2@ident,
                                               test.use="wilcox", only.pos=TRUE, min.pct=0.1,
                                               features.use = NULL, nthreads = 16,
                                               adj.p.val.threshold=0.05)

Lung_SMARTSeq2.markers$cluster = as.numeric(Lung_SMARTSeq2.markers$cluster)
Lung_SMARTSeq2.markers = Lung_SMARTSeq2.markers %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
Lung_SMARTSeq2.markers = as.data.frame(Lung_SMARTSeq2.markers)

quiet(gc())

options(digits=3)
Lung_SMARTSeq2.markers = Lung_SMARTSeq2.markers[,c(10, 9, 8, 2, 3,4, 6,7,1)]
Lung_SMARTSeq2.markers$avg_logFC = round(Lung_SMARTSeq2.markers$avg_logFC, digits = 3)
Lung_SMARTSeq2.markers$within_avg_exp = round(Lung_SMARTSeq2.markers$within_avg_exp, digits = 3)
Lung_SMARTSeq2.markers$without_avg_exp = round(Lung_SMARTSeq2.markers$without_avg_exp, digits = 3)
Lung_SMARTSeq2.markers$p_val_adj = Lung_SMARTSeq2.markers$p_val_adj
Lung_SMARTSeq2.markers$p_val = Lung_SMARTSeq2.markers$p_val

Lung_SMARTSeq2.markers_write = Lung_SMARTSeq2.markers
Lung_SMARTSeq2.markers_write$cluster=paste("cluster", Lung_SMARTSeq2.markers_write$cluster, sep="")
fwrite(Lung_SMARTSeq2.markers_write, "ALLmarkers_minpct0.1_Adj_p0.05.txt", row.names=F, col.names=T, sep="\t", quote=F)


top20 = Lung_SMARTSeq2.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
top20 = as.data.frame(top20)
top20  = top20 [!duplicated(top20$gene),]
top20 = top20 %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
top20 = as.data.frame(top20)
p = DoHeatmap2(object = Lung_SMARTSeq2, genes.use = top20$gene, genes.ident = top20$cluster,
               slim.col.label = TRUE, remove.key = FALSE, cex.row=4, disp.min = -2.5, disp.max = 2.5)
ggsave(file = "marker.res2.png", plot = p, device="png", units="in", dpi = 300, width = 15, height = 30, limitsize=FALSE)

#export the number of each subsets
tmp = table(Lung_SMARTSeq2@ident, Lung_SMARTSeq2@meta.data$orig.ident)
write.table(tmp, "./Lung_SMARTSeq2_subsetNumber.txt", row.names=T, sep="\t", quote=F)
save.pigz(Lung_SMARTSeq2, file = "LungSmartseq2.rda", n.cores=16)


##########################

#DC analysis
sample.name=("/SMARTSeq2_DC_")
mBC = Lung_SMARTSeq2@raw.data[, WhichCells(object = Lung_SMARTSeq2, ident = c(8))]

mBC = CreateSeuratObject(raw.data = mBC, min.cells = 5, min.genes = 500)

# normalizing data

colnames(mBC@meta.data)[2]="nReads"
mBC = NormalizeData(object = mBC, scale.factor=1000000)
mBC = ScaleData(object = mBC, vars.to.regress = c("nReads"), do.par=TRUE, num.cores=16)
mBC = FindVariableGenes(object = mBC, mean.function = ExpMean, dispersion.function = LogVMR, 
                                   x.low.cutoff = 0.1, x.high.cutoff = Inf, y.cutoff = 0.5, do.plot=FALSE)
length(x = mBC@var.genes)
closeAllConnections()
gc()

# perform PCA
mBC = RunPCA(object = mBC, pc.genes = mBC@var.genes, do.print = FALSE, pcs.compute = 30, overwrite=TRUE)
mBC = ProjectPCA(object = mBC, do.print = FALSE)
mBC = JackStraw(object = mBC, num.replicate = 100, num.pc = 30,do.par=TRUE, num.cores=16)
mBC = JackStrawPlot(object = mBC, PCs = 1:30)

#clustering and detection of cellular subsets
dims=mBC@dr$pca@jackstraw@overall.p.values[mBC@dr$pca@jackstraw@overall.p.values[,2]>0.05,1]
dims=c(1:(min(dims)-1))

mBC = FindClusters(object = mBC, reduction.type = "pca", dims.use = dims,
                              resolution = 1.1, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
py_set_seed(42)
mBC = DoFItSNE(mBC , reduction_use = "pca", dims_use = as.integer(dims),
                          perplexity=10, n_jobs=as.integer(16), df = 0.9)

p1 = DimPlot(object = mBC, reduction.use = "FItSNE", group.by = "orig.ident", 
             do.return = TRUE, vector.friendly = TRUE, pt.size = 1.0) +
  theme(axis.title.x = element_text(size=10, family = "Arial"), 
        axis.title.y = element_text(size=10, family = "Arial"), 
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"), 
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial")) +
  theme(panel.border = element_rect(fill = NA, size = 1.0))  
p2 = DimPlot(object = mBC, reduction.use = "FItSNE", do.label = TRUE, label.size = 10,
             do.return = TRUE, vector.friendly = TRUE, pt.size = 1.0) +
  theme(axis.title.x = element_text(size=10, family = "Arial"), 
        axis.title.y = element_text(size=10, family = "Arial"), 
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"), 
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial")) +
  theme(panel.border = element_rect(fill = NA, size = 1.0)) 

ggsave(file = "FItSNE_plot_DC_SmartSeq2.png", plot = plot_grid(p1, p2, ncol = 2), 
       device="png", units="in", dpi = 300, width = 8, height = 4, limitsize=FALSE)

#detection of subset-marker genes

tmp = table(mBC@ident)
tmp = as.data.frame(tmp)
write.table(tmp, "DC_subsetNumber.txt", row.names=F, sep="\t", quote=F)
mBC.markers = rDBEC::FindMarkers_parallel_lite(mBC, AllcellsIdent=mBC@ident,
                                                          test.use="wilcox", only.pos=TRUE, min.pct=0.1,
                                                          features.use = NULL, nthreads = 16,
                                                          adj.p.val.threshold=0.05)

mBC.markers$cluster = as.numeric(mBC.markers$cluster)
mBC.markers = mBC.markers %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
mBC.markers = as.data.frame(mBC.markers)

quiet(gc())

options(digits=3)
mBC.markers = mBC.markers[,c(10, 9, 8, 2, 3,4, 6,7,1)]
mBC.markers$avg_logFC = round(mBC.markers$avg_logFC, digits = 3)
mBC.markers$within_avg_exp = round(mBC.markers$within_avg_exp, digits = 3)
mBC.markers$without_avg_exp = round(mBC.markers$without_avg_exp, digits = 3)
mBC.markers$p_val_adj = mBC.markers$p_val_adj
mBC.markers$p_val = mBC.markers$p_val

mBC.markers_write = mBC.markers
mBC.markers_write$cluster=paste("cluster", mBC.markers_write$cluster, sep="")
fwrite(mBC.markers_write, "DC_ALLmarkers_minpct0.1_Adj_p0.05.txt", row.names=F, col.names=T, sep="\t", quote=F)


top20 = mBC.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
top20 = as.data.frame(top20)
top20  = top20 [!duplicated(top20$gene),]
top20 = top20 %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
top20 = as.data.frame(top20)
p = DoHeatmap2(object = mBC, genes.use = top20$gene, genes.ident = top20$cluster,
               slim.col.label = TRUE, remove.key = FALSE, cex.row=4, disp.min = -2.5, disp.max = 2.5)
ggsave(file = "DC_marker.res1.1.png", plot = p, device="png", units="in", dpi = 300, width = 15, height = 30, limitsize=FALSE)

#export the number of each subsets
tmp = table(mBC@ident, mBC@meta.data$orig.ident)
write.table(tmp, "./LungDC_SMARTSeq2_DC_subsetNumber.txt", row.names=T, sep="\t", quote=F)
save.pigz(mBC, file = "LungSmartseq2_DC.rda", n.cores=8)

##########################

#Tcell analysis
sample.name=("/SMARTSeq2_Tcell_")
mBC = Lung_SMARTSeq2@raw.data[, WhichCells(object = Lung_SMARTSeq2, ident = c(13))]

mBC = CreateSeuratObject(raw.data = mBC, min.cells = 5, min.genes = 500)

# normalizing data

colnames(mBC@meta.data)[2]="nReads"
mBC = NormalizeData(object = mBC, scale.factor=1000000)
mBC = ScaleData(object = mBC, vars.to.regress = c("nReads"), do.par=TRUE, num.cores=16)
mBC = FindVariableGenes(object = mBC, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.1, x.high.cutoff = Inf, y.cutoff = 0.5, do.plot=FALSE)
length(x = mBC@var.genes)
closeAllConnections()
gc()

# perform PCA
mBC = RunPCA(object = mBC, pc.genes = mBC@var.genes, do.print = FALSE, pcs.compute = 30, overwrite=TRUE)
mBC = ProjectPCA(object = mBC, do.print = FALSE)
mBC = JackStraw(object = mBC, num.replicate = 100, num.pc = 30,do.par=TRUE, num.cores=16)
mBC = JackStrawPlot(object = mBC, PCs = 1:30)

#clustering and detection of cellular subsets
dims=mBC@dr$pca@jackstraw@overall.p.values[mBC@dr$pca@jackstraw@overall.p.values[,2]>0.05,1]
dims=c(1:(min(dims)-1))

mBC = FindClusters(object = mBC, reduction.type = "pca", dims.use = dims,
                   resolution = 0.8, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
py_set_seed(42)
mBC = DoFItSNE(mBC , reduction_use = "pca", dims_use = as.integer(dims),
               perplexity=10, n_jobs=as.integer(16), df = 0.9)

p1 = DimPlot(object = mBC, reduction.use = "FItSNE", group.by = "orig.ident", 
             do.return = TRUE, vector.friendly = TRUE, pt.size = 1.0) +
  theme(axis.title.x = element_text(size=10, family = "Arial"), 
        axis.title.y = element_text(size=10, family = "Arial"), 
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"), 
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial")) +
  theme(panel.border = element_rect(fill = NA, size = 1.0))  
p2 = DimPlot(object = mBC, reduction.use = "FItSNE", do.label = TRUE, label.size = 10,
             do.return = TRUE, vector.friendly = TRUE, pt.size = 1.0) +
  theme(axis.title.x = element_text(size=10, family = "Arial"), 
        axis.title.y = element_text(size=10, family = "Arial"), 
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"), 
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial")) +
  theme(panel.border = element_rect(fill = NA, size = 1.0)) 

ggsave(file = "FItSNE_plot_Tcell_SmartSeq2.png", plot = plot_grid(p1, p2, ncol = 2), 
       device="png", units="in", dpi = 300, width = 8, height = 4, limitsize=FALSE)

#because CD4/CD8 T cells could not separate clearly, we annotate CD4/CD8 T cells by Cd4/Cd8a/Cd8b1 expression
hoge = t(mBC@raw.data[c("Cd4", "Cd8a", "Cd8b1"),])
Cd4 = rownames(hoge[hoge[,1]>0,])
hoge2 = hoge[,c(2,3)]
Cd8 = rownames(hoge[rowSums(hoge2)>0,])
nonCd4Cd8 = rownames(hoge[!rownames(hoge)%in%c(Cd4, Cd8),])

Cd4_1 =rep("Tcell_CD4T", length(Cd4))
names(Cd4_1)=Cd4
Cd4_1 = Cd4_1[!names(Cd4_1) %in% c("Smart-seq2_1345", "Smart-seq2_124")]
Cd8_1 =rep("Tcell_CD8T", length(Cd8))
names(Cd8_1)=Cd8
Cd8_1 = Cd8_1[!names(Cd8_1) %in% c("Smart-seq2_1345", "Smart-seq2_124")]
nonCd4Cd8_1 = rep("misc", length(nonCd4Cd8)+2)
names(nonCd4Cd8_1)=c(nonCd4Cd8, "Smart-seq2_1345", "Smart-seq2_124")

hoge = c(Cd4_1, Cd8_1, nonCd4Cd8_1)

mBC = AddMetaData(mBC, metadata=c(Cd4_1, Cd8_1, nonCd4Cd8_1), col.name = "celltype")

#export the number of each subsets
tmp = table(mBC@meta.data$celltype, mBC@meta.data$orig.ident)
write.table(tmp, "./LungDC_SMARTSeq2_Tcell_subsetNumber.txt", row.names=T, sep="\t", quote=F)
save.pigz(mBC, file = "LungSmartseq2_Tcell.rda", n.cores=8)

#add annotations to original Seurat metadata and save Seurat object
res=NULL
load("LungSmartseq2_DC.rda")
hoge = mBC@meta.data$res.1.1
hoge[hoge %in% c("0")]="DC_cDC2"
hoge[hoge %in% c("1")]="DC_cDC1"
hoge[hoge %in% c("2")]="DC_migratory"
names(hoge)=rownames(mBC@meta.data)
res=c(res, hoge)

load("LungSmartseq2_Tcell.rda")
hoge = mBC@meta.data$celltype
names(hoge)=rownames(mBC@meta.data)
res=c(res, hoge)

hoge2 = Lung_SMARTSeq2@meta.data[!rownames(Lung_SMARTSeq2@meta.data)%in% 
                                            names(res),"res.2"]
hoge2[hoge2 %in% c("0", "2", "14")]="Endo_VEC"
hoge2[hoge2 %in% c("15")]="Endo_venous"
hoge2[hoge2 %in% c("5")]="Endo_capillary"
hoge2[hoge2 %in% c("7")]="Endo_arterial"
hoge2[hoge2 %in% c("10")]="Endo_LEC"
hoge2[hoge2 %in% c("4")]="Epi_AT2"
hoge2[hoge2 %in% c("20")]="Epi_ciliated"
hoge2[hoge2 %in% c("19")]="Pericyte"
hoge2[hoge2 %in% c("16")]="SMC"
hoge2[hoge2 %in% c("17")]="NK"
hoge2[hoge2 %in% c("18")]="Neutrophil"
hoge2[hoge2 %in% c("12")]="Bcell"
hoge2[hoge2 %in% c("3", "11")]="FB_Inmthi"
hoge2[hoge2 %in% c("1")]="FB_Dcnhi"
hoge2[hoge2 %in% c("6")]="Mo_Ly6Chi"
hoge2[hoge2 %in% c("9")]="Mo_Ly6Clo"

names(hoge2)=rownames(Lung_SMARTSeq2@meta.data[!rownames(Lung_SMARTSeq2@meta.data)%in% 
                                            names(res),])
res = c(hoge2, res)

Lung_SMARTSeq2=AddMetaData(Lung_SMARTSeq2, metadata=res, col.name="celltype")
save.pigz(Lung_SMARTSeq2, file="Lung_SMARTSeq2.rda", n.cores=16)

#export the number of each subsets
tmp = table(Lung_SMARTSeq2@meta.data$celltype, Lung_SMARTSeq2@meta.data$orig.ident)
write.table(tmp, "./Lung_SMARTSeq2_cellSubset.txt", row.names=T, sep="\t", quote=F)


