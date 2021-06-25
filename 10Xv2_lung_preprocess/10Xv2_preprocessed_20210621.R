#re-analysis of 10X Chromium v2 data of murine lung of Tabula Muris data
options(stringsAsFactors = FALSE)

#please install Seurat v2.3.4. and rDBEC package before use (could be downloaded from this repository)
suppressWarnings(suppressMessages(source("library_source_Seurat.R")))

##normal mouse Lung 10X analysis (Tabula Muris data)
fileUrl <- "https://ndownloader.figshare.com/files/10700167"
download.file(fileUrl, destfile = "./droplet.zip", method = "wget")
unzip("droplet.zip")

Lung10X_1 = Read10X(data.dir = "./droplet/Lung-10X_P7_8/")
Lung10X_2 = Read10X(data.dir = "./droplet/Lung-10X_P7_9/")
Lung10X_3 = Read10X(data.dir = "./droplet/Lung-10X_P8_12/")
Lung10X_4 = Read10X(data.dir = "./droplet/Lung-10X_P8_13/")

colnames(x = Lung10X_1) = paste('Lung10X_1', colnames(x = Lung10X_1), sep = '_')
colnames(x = Lung10X_2) = paste('Lung10X_2', colnames(x = Lung10X_2), sep = '_')
colnames(x = Lung10X_3) = paste('Lung10X_3', colnames(x = Lung10X_3), sep = '_')
colnames(x = Lung10X_4) = paste('Lung10X_4', colnames(x = Lung10X_4), sep = '_')

Lung10X_1 = CreateSeuratObject(raw.data = Lung10X_1, min.cells = 1, min.genes = 500)
Lung10X_2 = CreateSeuratObject(raw.data = Lung10X_2, min.cells = 1, min.genes = 500)
Lung10X_3 = CreateSeuratObject(raw.data = Lung10X_3, min.cells = 1, min.genes = 500)
Lung10X_4 = CreateSeuratObject(raw.data = Lung10X_4, min.cells = 1, min.genes = 500)


Lung10X = MergeSeurat(object1 = Lung10X_1, object2=Lung10X_2)
Lung10X = MergeSeurat(object1 = Lung10X, object2=Lung10X_3)
Lung10X = MergeSeurat(object1 = Lung10X, object2=Lung10X_4, min.cells = 5, min.genes = 500)
Lung10X = Lung10X@raw.data
Lung10X = as.matrix(Lung10X)
Lung10X = Lung10X[rowSums(Lung10X)>0,]
colnames(x = Lung10X) = paste('Lung10X', c(1:ncol(Lung10X)), sep = '_')

Lung10X = CreateSeuratObject(raw.data = Lung10X, min.cells = 5, min.genes = 500)

orig.ident = rep("10Xv2", nrow(Lung10X@meta.data))
names(orig.ident)=rownames(Lung10X@meta.data)
Lung10X = AddMetaData(object = Lung10X, metadata = orig.ident, col.name = "orig.ident")

ribo.genes = grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = Lung10X@data), value = TRUE)
percent.ribo = Matrix::colSums(Lung10X@raw.data[ribo.genes, ])/Matrix::colSums(Lung10X@raw.data)
Lung10X = AddMetaData(object = Lung10X, metadata = percent.ribo, col.name = "percent.ribo")

percent.RiboRNA = Lung10X@raw.data[c('Rn45s'), ]/Matrix::colSums(Lung10X@raw.data)
Lung10X = AddMetaData(object = Lung10X, metadata = percent.RiboRNA, col.name = "percent.RiboRNA")

#------------------------------------------------------------------

Lung10X = FilterCells(object = Lung10X, subset.names = c("nGene", "nUMI"),
                      low.thresholds = c(500, 1000), high.thresholds = c(Inf, Inf))

# normalizing data
Lung10X = NormalizeData(object = Lung10X, scale.factor=100000)
Lung10X = ScaleData(object = Lung10X, vars.to.regress = c("nUMI"), do.par=TRUE, num.cores=16)
Lung10X = FindVariableGenes(object = Lung10X, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.1, x.high.cutoff = Inf, y.cutoff = 0.5, do.plot=FALSE)
length(x = Lung10X@var.genes)

# perform PCA and select significant PCs
Lung10X = RunPCA(object = Lung10X, pc.genes = Lung10X@var.genes, do.print = FALSE, pcs.compute = 100, overwrite=TRUE)
Lung10X = ProjectPCA(object = Lung10X, do.print = FALSE)
Lung10X = JackStraw(object = Lung10X, num.replicate = 100, num.pc = 100,do.par=TRUE, num.cores=16)
Lung10X = JackStrawPlot(object = Lung10X, PCs = 1:100)
dims=Lung10X@dr$pca@jackstraw@overall.p.values[Lung10X@dr$pca@jackstraw@overall.p.values[,2]>0.05,1]
dims=c(1:(min(dims)-1))

#clustering and detection of cellular subsets
Lung10X = FindClusters(object = Lung10X, reduction.type = "pca", dims.use = dims,
                       resolution = 2.0, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

py_set_seed(42)
Lung10X = DoFItSNE(Lung10X, reduction_use = "pca", dims_use = as.integer(dims),
               perplexity=100, n_jobs=as.integer(16), df = 0.9)

p1 = DimPlot(object = Lung10X, reduction.use = "FItSNE", group.by = "orig.ident", 
             do.return = TRUE, vector.friendly = TRUE, pt.size = 1.0) +
  theme(axis.title.x = element_text(size=10, family = "Arial"), 
        axis.title.y = element_text(size=10, family = "Arial"), 
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"), 
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial")) +
  theme(panel.border = element_rect(fill = NA, size = 1.0))  
p2 = DimPlot(object = Lung10X, reduction.use = "FItSNE", do.label = TRUE, label.size = 10,
             do.return = TRUE, vector.friendly = TRUE, pt.size = 1.0) +
  theme(axis.title.x = element_text(size=10, family = "Arial"), 
        axis.title.y = element_text(size=10, family = "Arial"), 
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"), 
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial")) +
  theme(panel.border = element_rect(fill = NA, size = 1.0)) 

ggsave(file = "FItSNE_plot_Lung_10X.png", plot = plot_grid(p1, p2, ncol = 2), 
       device="png", units="in", dpi = 300, width = 8, height = 4, limitsize=FALSE)

#detection of subset-marker genes
mBC.markers = rDBEC::FindMarkers_parallel_lite(Lung10X, AllcellsIdent=Lung10X@ident,
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
fwrite(mBC.markers_write, "ALLmarkers_minpct0.1_Adj_p0.05.txt", row.names=F, col.names=T, sep="\t", quote=F)

top20 = mBC.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
top20 = as.data.frame(top20)
top20  = top20 [!duplicated(top20$gene),]
top20 = top20 %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
top20 = as.data.frame(top20)
p = DoHeatmap2(object = Lung10X, genes.use = top20$gene, genes.ident = top20$cluster,
               slim.col.label = TRUE, remove.key = FALSE, cex.row=4, disp.min = -2.5, disp.max = 2.5)
ggsave(file = "marker.res2.png", plot = p, device="png", units="in",
       dpi = 300, width = 15, height = 30, limitsize=FALSE)


#export the number of each subsets
tmp = table(Lung10X@ident, Lung10X@meta.data$orig.ident)
write.table(tmp, "./Lung10X_subsetNumber.txt", row.names=T, sep="\t", quote=F)

save.pigz(Lung10X, file = "Lung10X.rda", n.cores=8)

##########################

#DC analysis
sample.name=("/10X_DC_")
mBC = Lung10X@raw.data[, WhichCells(object = Lung10X, ident = c(18))]
mBC = CreateSeuratObject(raw.data = mBC, min.cells = 5, min.genes = 500)

# normalizing data
mBC = NormalizeData(object = mBC, scale.factor=100000)
mBC = ScaleData(object = mBC, vars.to.regress = c("nUMI"),
                do.par=TRUE, num.cores=16)
mBC = FindVariableGenes(object = mBC, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.1, x.high.cutoff = Inf, y.cutoff = 0.5, do.plot=FALSE)
length(x = mBC@var.genes)

closeAllConnections()
gc()
# perform PCA
mBC = RunPCA(object = mBC, pc.genes = mBC@var.genes, do.print = FALSE, pcs.compute = 30, overwrite=TRUE)
mBC = ProjectPCA(object = mBC, do.print = FALSE)

#load("Seurat_fibrosis_reads.rda")
mBC = JackStraw(object = mBC, num.replicate = 100, num.pc = 30,do.par=TRUE, num.cores=16)
mBC = JackStrawPlot(object = mBC, PCs = 1:30)
dev.off()
closeAllConnections()
gc()

#plot each principal component-associated genes
tmp=mBC@dr$pca@jackstraw@emperical.p.value.full
tmp = as.data.frame(mBC@dr$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05,1]
dims= c(1:(min(tmp1)-1))

#because significant PCs were less than 5, we use PC 1:5 for capture diversity of cells
dims=1:5

mBC = FindClusters(object = mBC, reduction.type = "pca", dims.use = dims,
                   resolution = 1.2, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

py_set_seed(42)
mBC = DoFItSNE(mBC, reduction_use = "pca", dims_use = as.integer(dims),
                   perplexity=10, n_jobs=as.integer(16), df = 0.9)

p1 = DimPlot(object = mBC, reduction.use = "FItSNE", group.by = "orig.ident", 
             do.return = TRUE, vector.friendly = TRUE, pt.size = 3.0, no.legend=FALSE) +
  theme(axis.title.x = element_text(size=10, family = "Arial"), 
        axis.title.y = element_text(size=10, family = "Arial"), 
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"), 
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial")) +
  theme(panel.border = element_rect(fill = NA, size = 1.0))  
p2 = DimPlot(object = mBC, reduction.use = "FItSNE", do.label = TRUE, label.size = 10,
             do.return = TRUE, vector.friendly = TRUE, pt.size = 3.0, no.legend=TRUE) +
  theme(axis.title.x = element_text(size=10, family = "Arial"), 
        axis.title.y = element_text(size=10, family = "Arial"), 
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"), 
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial")) +
  theme(panel.border = element_rect(fill = NA, size = 1.0)) 

dir.name="."
file.name=paste(dir.name, sample.name, "FItSNE_perplex10_reso1.2.png", sep='')
ggsave(file = file.name, plot = plot_grid(p1, p2), dpi = 100, width = 10, height = 5)


#detection of subset-marker genes
mBC.markers = rDBEC::FindMarkers_parallel_lite(mBC, AllcellsIdent=mBC@ident,
                                               test.use="wilcox", only.pos=TRUE, min.pct=0.1,
                                               features.use = NULL, nthreads = 4,
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
fwrite(mBC.markers_write, "ALLmarkers_minpct0.1_Adj_p0.05.txt", row.names=F, col.names=T, sep="\t", quote=F)

top20 = mBC.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
top20 = as.data.frame(top20)
top20  = top20 [!duplicated(top20$gene),]
top20 = top20 %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
top20 = as.data.frame(top20)
p = DoHeatmap2(object = mBC, genes.use = top20$gene, genes.ident = top20$cluster,
               slim.col.label = TRUE, remove.key = FALSE, cex.row=4, disp.min = -2.5, disp.max = 2.5)
ggsave(file = "DC_marker.res1.2.png", plot = p, device="png", units="in",
       dpi = 300, width = 5, height = 7, limitsize=FALSE)


#export the number of each subsets
tmp = table(mBC@ident, mBC@meta.data$orig.ident)
write.table(tmp, "./Lung10X_DC_subsetNumber.txt", row.names=T, sep="\t", quote=F)
save.pigz(mBC, file = "Lung10X_DC.rda", n.cores=8)


#add annotations to original Seurat metadata and save Seurat object
res=NULL
hoge = mBC@meta.data$res.1.2
hoge[hoge=="0"]="DC_cDC2"
hoge[hoge=="3"]="DC_pDC"
hoge[hoge=="2"]="Mac_IM"
hoge[hoge=="1"]="DC_cDC1"
names(hoge)=rownames(mBC@meta.data)
res=c(res, hoge)
save.pigz(mBC, file = "Lung10X_DC.rda", n.cores=8)

hoge1 = Lung10X@meta.data[!rownames(Lung10X@meta.data)%in% names(res),]
hoge2 = as.character(hoge1$res.2)

hoge2[hoge2 %in% c("0", "3", "4", "9", "19")]="FB_Inmthi"
hoge2[hoge2 %in% c("2")]="FB_Dcnhi"
hoge2[hoge2 %in% c("25")]="Mesothelial"
hoge2[hoge2 %in% c("8", "16")]="Endo_VEC"
hoge2[hoge2 %in% c("20")]="Endo_capillary"
hoge2[hoge2 %in% c("22")]="Epi_ciliated"
hoge2[hoge2 %in% c("17")]="Epi_AT2"
hoge2[hoge2 %in% c("24")]="Endo_venous"
hoge2[hoge2 %in% c("1", "7")]="NK"
hoge2[hoge2 %in% c("5")]="SMC"
hoge2[hoge2 %in% c("6", "23")]="Mac_AM"
hoge2[hoge2 %in% c("12")]="Tcell_CD4T"
hoge2[hoge2 %in% c("21")]="Tcell_CD8T"
hoge2[hoge2 %in% c("14")]="Neutrophil"
hoge2[hoge2 %in% c("13")]="Mo_Ly6Chi"
hoge2[hoge2 %in% c("10")]="Mo_Ly6Clo"
hoge2[hoge2 %in% c("15")]="Pericyte"
hoge2[hoge2 %in% c("11")]="Bcell"
hoge2[hoge2 %in% c("34")]="Endo_LEC"
hoge2[hoge2 %in% c("26")]="Mast_cell"

names(hoge2)=rownames(Lung10X@meta.data[!rownames(Lung10X@meta.data)%in% 
                                          names(res),])
res = c(hoge2, res)

Lung10X=AddMetaData(Lung10X, metadata=res, col.name="celltype")
save.pigz(Lung10X, file="Lung10X.rda", n.cores=16)

#export the number of each subsets
tmp = table(Lung10X@meta.data$celltype, Lung10X@meta.data$orig.ident)
write.table(tmp, "./Lung10Xv2_cellSubset.txt", row.names=T, sep="\t", quote=F)
