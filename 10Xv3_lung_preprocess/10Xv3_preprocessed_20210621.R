#re-analysis of 10X Chromium v3 data of murine lung (GSE145998)
options(stringsAsFactors = FALSE)
suppressWarnings(suppressMessages(source("library_source_Seurat.R")))

##normal mouse Lung 10X v3 analysis (GSE145998)
fileUrl <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE145nnn/GSE145998/suppl/GSE145998%5FCtrl%5Fmouse%5FScell1%2Edge%2Etxt%2Egz"
download.file(fileUrl, destfile = "GSE145998_Ctrl_mouse_Scell1.dge.txt.gz", method = "wget")

hoge = fread("./GSE145998_Ctrl_mouse_Scell1.dge.txt.gz")
hoge = as.data.frame(hoge)
rownames(hoge)=hoge[,1]
hoge = hoge[,2:ncol(hoge)]
hoge = as.matrix(hoge)
hoge = Matrix::Matrix(hoge, sparse = T)
Lung10Xv3 = CreateSeuratObject(hoge, min.cells = 5, min.genes=500)

#downsample cells to the cell number of 10X v2 dataset
Lung10Xv3 = SubsetData(Lung10Xv3, max.cells.per.ident = 5507)
orig.ident = rep("10Xv3", nrow(Lung10Xv3@meta.data))
names(orig.ident)=rownames(Lung10Xv3@meta.data)
Lung10Xv3 = AddMetaData(object = Lung10Xv3, metadata = orig.ident, col.name = "orig.ident")

mito.genes = grep(pattern = "^mt.", x = rownames(x = Lung10Xv3@data), value = TRUE)
percent.mito = Matrix::colSums(Lung10Xv3@raw.data[mito.genes, ])/Matrix::colSums(Lung10Xv3@raw.data)
Lung10Xv3 = AddMetaData(object = Lung10Xv3, metadata = percent.mito, col.name = "percent.mito")

ribo.genes = grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = Lung10Xv3@data), value = TRUE)
percent.ribo = Matrix::colSums(Lung10Xv3@raw.data[ribo.genes, ])/Matrix::colSums(Lung10Xv3@raw.data)
Lung10Xv3 = AddMetaData(object = Lung10Xv3, metadata = percent.ribo, col.name = "percent.ribo")

riboRNA.genes = grep(pattern = "^Rn[[:digit:]]", x = rownames(x = Lung10Xv3@data), value = TRUE)
percent.riboRNA = Matrix::colSums(Lung10Xv3@raw.data[riboRNA.genes, ])/Matrix::colSums(Lung10Xv3@raw.data)
Lung10Xv3 = AddMetaData(object = Lung10Xv3, metadata = percent.riboRNA, col.name = "percent.riboRNA")


Lung10Xv3 = FilterCells(object = Lung10Xv3, subset.names = c("nGene", "nUMI"),
                      low.thresholds = c(500, 1000), high.thresholds = c(Inf, Inf))

# normalizing data
Lung10Xv3 = NormalizeData(object = Lung10Xv3, scale.factor=100000)
Lung10Xv3 = ScaleData(object = Lung10Xv3, vars.to.regress = c("nUMI"), do.par=TRUE, num.cores=16)
Lung10Xv3 = FindVariableGenes(object = Lung10Xv3, mean.function = ExpMean, dispersion.function = LogVMR, 
                              x.low.cutoff = 0.1, x.high.cutoff = Inf, y.cutoff = 0.5, do.plot=FALSE)
length(x = Lung10Xv3@var.genes)


# perform PCA and select significant PCs
Lung10Xv3 = RunPCA(object = Lung10Xv3, pc.genes = Lung10Xv3@var.genes, do.print = FALSE, pcs.compute = 100, overwrite=TRUE)
Lung10Xv3 = ProjectPCA(object = Lung10Xv3, do.print = FALSE)
PCElbowPlot(object = Lung10Xv3, num.pc = 100) 
Lung10Xv3 = JackStraw(object = Lung10Xv3, num.replicate = 100, num.pc = 100,do.par=TRUE, num.cores=16)
Lung10Xv3 = JackStrawPlot(object = Lung10Xv3, PCs = 1:100)
dims=Lung10Xv3@dr$pca@jackstraw@overall.p.values[Lung10Xv3@dr$pca@jackstraw@overall.p.values[,2]>0.05,1]
dims=c(1:(min(dims)-1))


#clustering and detection of cellular subsets
Lung10Xv3 = FindClusters(object = Lung10Xv3, reduction.type = "pca", dims.use = dims,
                       resolution = 4.0, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

py_set_seed(42)
Lung10Xv3 = DoFItSNE(Lung10Xv3, reduction_use = "pca", dims_use = as.integer(dims),
               perplexity=100, n_jobs=as.integer(16), df = 0.9)

p1 = DimPlot(object = Lung10Xv3, reduction.use = "FItSNE", group.by = "orig.ident", 
             do.return = TRUE, vector.friendly = TRUE, pt.size = 1.0) +
  theme(axis.title.x = element_text(size=10, family = "Arial"), 
        axis.title.y = element_text(size=10, family = "Arial"), 
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"), 
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial")) +
  theme(panel.border = element_rect(fill = NA, size = 1.0))  
p2 = DimPlot(object = Lung10Xv3, reduction.use = "FItSNE", do.label = TRUE, label.size = 10,
             do.return = TRUE, vector.friendly = TRUE, pt.size = 1.0) +
  theme(axis.title.x = element_text(size=10, family = "Arial"), 
        axis.title.y = element_text(size=10, family = "Arial"), 
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"), 
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial")) +
  theme(panel.border = element_rect(fill = NA, size = 1.0)) 

ggsave(file = "FItSNE_plot_Lung_10Xv3.png", plot = plot_grid(p1, p2, ncol = 2), 
       device="png", units="in", dpi = 300, width = 8, height = 4, limitsize=FALSE)

#detection of subset-marker genes
mBC.markers = rDBEC::FindMarkers_parallel_lite(Lung10Xv3, AllcellsIdent=Lung10Xv3@ident,
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
p = DoHeatmap2(object = Lung10Xv3, genes.use = top20$gene, genes.ident = top20$cluster,
               slim.col.label = TRUE, remove.key = FALSE, cex.row=4, disp.min = -2.5, disp.max = 2.5)
ggsave(file = "marker.res4.png", plot = p, device="png", units="in",
       dpi = 300, width = 15, height = 30, limitsize=FALSE)


#export the number of each subsets
tmp = table(Lung10Xv3@ident, Lung10Xv3@meta.data$orig.ident)
write.table(tmp, "./Lung10Xv3_subsetNumber.txt", row.names=T, sep="\t", quote=F)

save.pigz(Lung10Xv3, file = "Lung10Xv3.rda", n.cores=16)

##########################

#lymphoid cells subclustering
sample.name=("/10Xv3_T_ILC_")
mBC = Lung10Xv3@raw.data[, WhichCells(object = Lung10Xv3, ident = c(4, 24))]
mBC = CreateSeuratObject(raw.data = mBC, min.cells = 5, min.genes = 500)

mBC = NormalizeData(object = mBC, scale.factor=100000)
mBC = ScaleData(object = mBC, vars.to.regress = c("nUMI"), do.par=TRUE, num.cores=16)
mBC = FindVariableGenes(object = mBC, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.1, x.high.cutoff = Inf, y.cutoff = 0.5, do.plot=FALSE)
length(x = mBC@var.genes)

# perform PCA and select significant PCs
mBC = RunPCA(object = mBC, pc.genes = mBC@var.genes, do.print = FALSE, pcs.compute = 50, overwrite=TRUE)
mBC = ProjectPCA(object = mBC, do.print = FALSE)

mBC = JackStraw(object = mBC, num.replicate = 100, num.pc = 50,do.par=TRUE, num.cores=16)
mBC = JackStrawPlot(object = mBC, PCs = 1:50)
dev.off()
closeAllConnections()
gc()

#plot each principal component-associated genes
tmp=mBC@dr$pca@jackstraw@emperical.p.value.full
tmp = as.data.frame(mBC@dr$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05,1]
dims= c(1:(min(tmp1)-1))


mBC = FindClusters(object = mBC, reduction.type = "pca", dims.use = dims,
                   resolution = 1.5, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

py_set_seed(42)
mBC = DoFItSNE(mBC, reduction_use = "pca", dims_use = as.integer(dims),
                   perplexity=50, n_jobs=as.integer(16), df = 0.9)

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
file.name=paste(dir.name, sample.name, "FItSNE_perplex50_reso1.5.png", sep='')
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
fwrite(mBC.markers_write, "T_ILC_ALLmarkers_minpct0.1_Adj_p0.05.txt", row.names=F, col.names=T, sep="\t", quote=F)

top30 = mBC.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
top30 = as.data.frame(top30)
top30  = top30 [!duplicated(top30$gene),]
top30 = top30 %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
top30 = as.data.frame(top30)
p = DoHeatmap2(object = mBC, genes.use = top30$gene, genes.ident = top30$cluster,
               slim.col.label = TRUE, remove.key = FALSE, cex.row=4, disp.min = -2.5, disp.max = 2.5)
ggsave(file = "T_ILC_marker.res1.5.png", plot = p, device="png", units="in",
       dpi = 300, width = 5, height = 7, limitsize=FALSE)


#export the number of each subsets
tmp = table(mBC@ident, mBC@meta.data$orig.ident)
write.table(tmp, "./Lung10Xv3_T_ILC_subsetNumber.txt", row.names=T, sep="\t", quote=F)
save.pigz(mBC, file = "Lung10Xv3_T_ILC.rda", n.cores=8)

#add annotations to original Seurat metadata and save Seurat object
res=NULL
hoge = mBC@meta.data$res.1.5
hoge[hoge %in% c("0", "6")]="ILC2"
hoge[hoge %in% c("1", "3")]="Tcell_CD8T"
hoge[hoge %in% c("2", "5")]="doublet"
hoge[hoge=="4"]="Tcell_CD4T"
names(hoge)=rownames(mBC@meta.data)
res=c(res, hoge)

hoge = Lung10Xv3@meta.data[!rownames(Lung10Xv3@meta.data)%in% names(res),]
hoge1 = hoge$res.4
hoge1[hoge1 %in% c("0", "14", "18")]="Endo_VEC"
hoge1[hoge1 %in% c("12", "17", "26")]="Endo_capillary"
hoge1[hoge1 %in% c("1", "2", "5", "11", "36")]="Mac_AM"
hoge1[hoge1 %in% c("13", "19", "30")]="FB_Inmthi"
hoge1[hoge1=="6"]="Pericyte"
hoge1[hoge1=="10"]="Epi_AT2"
hoge1[hoge1=="15"]="DC_cDC2"
hoge1[hoge1 %in% c("16", "25")]="Bcell"
hoge1[hoge1=="21"]="FB_Dcnhi"
hoge1[hoge1=="23"]="DC_cDC1"
hoge1[hoge1=="28"]="Neutrophil"
hoge1[hoge1=="27"]="Epi_AT1"
hoge1[hoge1=="31"]="SMC"
hoge1[hoge1=="9"]="doublet"
hoge1[hoge1=="32"]="doublet"
hoge1[hoge1=="35"]="DC_migratory"
hoge1[hoge1=="38"]="Mesothelial"
hoge1[hoge1=="37"]="DC_pDC"
hoge1[hoge1=="33"]="Endo_arterial"
hoge1[hoge1=="39"]="Endo_LEC"
hoge1[hoge1=="3"]="Mo_Ly6Chi"
hoge1[hoge1=="7"]="Mo_Ly6Clo"
hoge1[hoge1=="8"]="Endo_venous"
hoge1[hoge1=="20"]="NK"
hoge1[hoge1=="22"]="Bcell"
hoge1[hoge1=="29"]="Mac_IM"
hoge1[hoge1=="34"]="doublet"

names(hoge1)=rownames(Lung10Xv3@meta.data[!rownames(Lung10Xv3@meta.data)%in% 
                                          names(res),])
res = c(hoge1, res)

Lung10Xv3=AddMetaData(Lung10Xv3, metadata=res, col.name="celltype")
save.pigz(Lung10Xv3, file="Lung10Xv3.rda", n.cores=16)

#export the number of each subsets
tmp = table(Lung10Xv3@meta.data$celltype, Lung10Xv3@meta.data$orig.ident)
write.table(tmp, "./Lung10Xv3_cellSubset.txt", row.names=T, sep="\t", quote=F)
