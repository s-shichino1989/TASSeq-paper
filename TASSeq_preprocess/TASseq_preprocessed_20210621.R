#re-analysis of TAs-Seq data of murine lung
options(stringsAsFactors = FALSE)
suppressWarnings(suppressMessages(source("/datadrive/Rhapsody_analysis/Rscripts/library_source_Seurat.R")))

hoge = tableread_fast_sparse("./matrix_inflection_UTLung.txt.gz")
tablelist=list()
tablelist[[1]]=hoge
names(tablelist)[1]="TASSeq"
DBEC_filter = background_subtraction_Biex(tablelist, min.event=100, minimum.max.expr=8,
                                          species="mmu", min.ave=6, AutoThreshold = FALSE,
                                          min.diff=5.5, modelnames="E",
                                          uncert.thre=1, nthreads=24, sample.name="TASSeq")

names(DBEC_filter) = "TASSeq"
DBEC_res = apply_DBEC_filter(tablelist, DBEC_filter=DBEC_filter, nthreads=16, sample.name = "TASSeq")
names(DBEC_res) = "TASSeq"
colnames(DBEC_res[[1]])=paste("TAS-Seq", colnames(DBEC_res[[1]]), "_")
##create Seurat object (Seurat v2 workflow) and annotate by hashtag
Lung_BDRhapsody = CreateSeuratObject(DBEC_res[[1]], min.cells = 5, min.genes = 500)

colnames(Lung_BDRhapsody@meta.data)[2]="nReads"

mito.genes = grep(pattern = "^mt.", x = rownames(x = Lung_BDRhapsody@data), value = TRUE)
percent.mito = Matrix::colSums(Lung_BDRhapsody@raw.data[mito.genes, ])/Matrix::colSums(Lung_BDRhapsody@raw.data)
Lung_BDRhapsody = AddMetaData(object = Lung_BDRhapsody, metadata = percent.mito, col.name = "percent.mito")

ribo.genes = grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = Lung_BDRhapsody@data), value = TRUE)
percent.ribo = Matrix::colSums(Lung_BDRhapsody@raw.data[ribo.genes, ])/Matrix::colSums(Lung_BDRhapsody@raw.data)
Lung_BDRhapsody = AddMetaData(object = Lung_BDRhapsody, metadata = percent.ribo, col.name = "percent.ribo")

riboRNA.genes = grep(pattern = "^Rn[[:digit:]]", x = rownames(x = Lung_BDRhapsody@data), value = TRUE)
percent.riboRNA = Matrix::colSums(Lung_BDRhapsody@raw.data[riboRNA.genes, ])/Matrix::colSums(Lung_BDRhapsody@raw.data)
Lung_BDRhapsody = AddMetaData(object = Lung_BDRhapsody, metadata = percent.riboRNA, col.name = "percent.riboRNA")


nReads_log = log10(Matrix::colSums(Lung_BDRhapsody@raw.data))
Lung_BDRhapsody = AddMetaData(object = Lung_BDRhapsody, metadata = nReads_log, col.name = "nReads.log")


#filter out outliers
Lung_BDRhapsody = FilterCells(object = Lung_BDRhapsody, subset.names = c("percent.mito", "nGene"),
                  low.thresholds = c(-Inf, 500), high.thresholds = c(0.4, Inf))

message("Normalizing data...")
Lung_BDRhapsody = NormalizeData(object = Lung_BDRhapsody, scale.factor=1000000,display.progress = FALSE)
message("Scaling data...")
Lung_BDRhapsody = ScaleData(object = Lung_BDRhapsody, vars.to.regress = c("nReads", "percent.ribo", "percent.riboRNA"),
                do.par=TRUE, num.cores=16, display.progress = FALSE)
Lung_BDRhapsody = FindVariableGenes(object = Lung_BDRhapsody, mean.function = ExpMean, dispersion.function = LogVMR,
                        x.low.cutoff = 0.1, x.high.cutoff = Inf, y.cutoff = 0.5, do.plot=FALSE)
hvg.number = length(x = Lung_BDRhapsody@var.genes)

Lung_BDRhapsody = RunPCA(object = Lung_BDRhapsody, pc.genes = Lung_BDRhapsody@var.genes, do.print = FALSE, pcs.compute = 100, overwrite=TRUE)
Lung_BDRhapsody = ProjectPCA(object = Lung_BDRhapsody, do.print = FALSE)
Lung_BDRhapsody = JackStraw(object = Lung_BDRhapsody, num.replicate = 100, num.pc = 100,do.par=TRUE, num.cores=24, display.progress = FALSE)
Lung_BDRhapsody = JackStrawPlot(object = Lung_BDRhapsody, PCs = 1:100)
quiet(dev.off())
quiet(gc())

#clustering and detection of cellular subsets
tmp = as.data.frame(Lung_BDRhapsody@dr$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05,1]
if(length(tmp1)<1){
dims=1:100
} else {
dims= c(1:(min(tmp1)-1))
}

Lung_BDRhapsody = FindClusters(object = Lung_BDRhapsody, reduction.type = "pca", dims.use = dims,
                   resolution = 5.0, print.output = 0, save.SNN = FALSE, force.recalc = TRUE)
py_set_seed(42)
Lung_BDRhapsody = DoFItSNE(Lung_BDRhapsody, reduction_use = "pca", dims_use = as.integer(dims),
               perplexity=100, n_jobs=as.integer(16), df = 0.9, random.seed=42)

p1 = DimPlot(object = Lung_BDRhapsody, reduction.use = "FItSNE", do.label = TRUE,label.size = 6,
             do.return = TRUE, vector.friendly = TRUE, pt.size = 1, group.by="res.5",
             no.legend=FALSE, cols.use =custom_colors$discrete ) +
  theme(axis.title.x = element_text(size=10, family = "Arial"),
        axis.title.y = element_text(size=10, family = "Arial"),
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"),
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial"))
p1 = p1 + theme_void() + theme(legend.position='none')

p2 = DimPlot(object = Lung_BDRhapsody, reduction.use = "FItSNE", group.by = "orig.ident",
             do.label = FALSE,label.size = 5,
             do.return = TRUE, vector.friendly = TRUE, pt.size = 1,
             no.legend=FALSE) +
  theme(axis.title.x = element_text(size=10, family = "Arial"),
        axis.title.y = element_text(size=10, family = "Arial"),
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"),
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial"))
p2 = p2 + theme_void() + theme(legend.position='none')

ggsave(file = "./TASSeq_ALL_DimPlot.png",
       plot =plot_grid(p2, p1), device="png",
       units="in", dpi = 300,
       width = 6, height = 3, limitsize=FALSE)


tmp = table(Lung_BDRhapsody@ident)
tmp = as.data.frame(tmp)
write.table(tmp, "TASSeq_subsetNumber.txt", row.names=F, sep="\t", quote=F)

Lung_BDRhapsody.markers = rDBEC::FindMarkers_parallel_lite(Lung_BDRhapsody, AllcellsIdent=Lung_BDRhapsody@ident,
                                               test.use="wilcox", only.pos=TRUE, min.pct=0.1,
                                               features.use = NULL, nthreads = 28,
                                               adj.p.val.threshold=0.05)

Lung_BDRhapsody.markers$cluster = as.numeric(Lung_BDRhapsody.markers$cluster)
Lung_BDRhapsody.markers = Lung_BDRhapsody.markers %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
Lung_BDRhapsody.markers = as.data.frame(Lung_BDRhapsody.markers)

quiet(gc())

options(digits=3)
Lung_BDRhapsody.markers = Lung_BDRhapsody.markers[,c(10, 9, 8, 2, 3,4, 6,7,1)]
Lung_BDRhapsody.markers$avg_logFC = round(Lung_BDRhapsody.markers$avg_logFC, digits = 3)
Lung_BDRhapsody.markers$within_avg_exp = round(Lung_BDRhapsody.markers$within_avg_exp, digits = 3)
Lung_BDRhapsody.markers$without_avg_exp = round(Lung_BDRhapsody.markers$without_avg_exp, digits = 3)
Lung_BDRhapsody.markers$p_val_adj = Lung_BDRhapsody.markers$p_val_adj
Lung_BDRhapsody.markers$p_val = Lung_BDRhapsody.markers$p_val

Lung_BDRhapsody.markers_write = Lung_BDRhapsody.markers
Lung_BDRhapsody.markers_write$cluster=paste("cluster", Lung_BDRhapsody.markers_write$cluster, sep="")
fwrite(Lung_BDRhapsody.markers_write, "ALLmarkers_minpct0.1_Adj_p0.05.txt", row.names=F, col.names=T, sep="\t", quote=F)

top20 = Lung_BDRhapsody.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
top20 = as.data.frame(top20)
top20  = top20 [!duplicated(top20$gene),]
top20 = top20 %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
top20 = as.data.frame(top20)
file.name=paste(dir.name, sample.name, "marker_res5.png", sep='')
p = DoHeatmap2(object = Lung_BDRhapsody, genes.use = top20$gene, genes.ident = top20$cluster,
               slim.col.label = TRUE, remove.key = FALSE, cex.row=4, disp.min = -2.5, disp.max = 2.5)
ggsave(file = "marker.res5.png", plot = p, device="png", units="in", dpi = 300, width = 20, height = 40, limitsize=FALSE)


#export the number of each subsets and save file
tmp = table(Lung_BDRhapsody@ident, Lung_BDRhapsody@meta.data$orig.ident)
write.table(tmp, "./Lung_day00_Rhapsody_subsetNumber.txt", row.names=T, sep="\t", quote=F)
save.pigz(Lung_BDRhapsody, file = "Lung_day00_Rhapsody.rda", n.cores=8)

###################
#gdT_ILC_analysis

mBC = SubsetData(Lung_BDRhapsody, ident.use = 35, subset.raw=TRUE)

mBC = NormalizeData(object = mBC, scale.factor=1000000)
mBC = ScaleData(object = mBC, vars.to.regress = c("nReads"), do.par=TRUE, num.cores=16)
mBC = FindVariableGenes(object = mBC, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.1, x.high.cutoff = Inf, y.cutoff = 0.5, do.plot=FALSE)
length(x = mBC@var.genes)

# perform PCA and select significant PCs
mBC = RunPCA(object = mBC, pc.genes = mBC@var.genes, do.print = FALSE, pcs.compute = 15, overwrite=TRUE)
mBC = ProjectPCA(object = mBC, do.print = FALSE)

mBC = JackStraw(object = mBC, num.replicate = 100, num.pc = 15,do.par=TRUE, num.cores=16)
mBC = JackStrawPlot(object = mBC, PCs = 1:15)
dev.off()
closeAllConnections()
gc()

#plot each principal component-associated genes
tmp=mBC@dr$pca@jackstraw@emperical.p.value.full
tmp = as.data.frame(mBC@dr$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05,1]
dims= c(1:(min(tmp1)-1))

#set dims 1:5 to reflect diversity
dims=1:5

mBC = FindClusters(object = mBC, reduction.type = "pca", dims.use = dims,
                   resolution = 1.0, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

py_set_seed(42)
mBC = DoFItSNE(mBC, reduction_use = "pca", dims_use = as.integer(dims),
               perplexity=10, n_jobs=as.integer(16), df = 1)

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

file.name="gdT_ILC_FItSNE_perplex10_reso1.png"
ggsave(file = file.name, plot = plot_grid(p1, p2), dpi = 100, width = 10, height = 5)

#divide gdT and ILCs by the expression of Il1rl1 (ST2)
FeaturePlot(mBC, features.plot=c("Il1rl1"), reduction.use="FItSNE", pt.size=2)

hoge = colnames(mBC@raw.data[,mBC@raw.data["Il1rl1",]>0])
hoge1 = rownames(mBC@meta.data)[!rownames(mBC@meta.data) %in% hoge]
hoge2 = rep("ILC2", length(hoge))
hoge3 = rep("T_gdT", length(hoge1))
names(hoge2)=hoge
names(hoge3)=hoge1

mBC = AddMetaData(mBC, metadata = c(hoge2, hoge3), col.name="celltype")
save.pigz(mBC, file = "ILC_gdT_cluster35.rda", n.cores=8)

###################
#DC_IM analysis

mBC = SubsetData(Lung_BDRhapsody, ident.use = 20, subset.raw=TRUE)

mBC = NormalizeData(object = mBC, scale.factor=1000000)
mBC = ScaleData(object = mBC, vars.to.regress = c("nReads"), do.par=TRUE, num.cores=16)
mBC = FindVariableGenes(object = mBC, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.1, x.high.cutoff = Inf, y.cutoff = 0.5, do.plot=FALSE)
length(x = mBC@var.genes)

# perform PCA and select significant PCs
mBC = RunPCA(object = mBC, pc.genes = mBC@var.genes, do.print = FALSE, pcs.compute = 15, overwrite=TRUE)
mBC = ProjectPCA(object = mBC, do.print = FALSE)

mBC = JackStraw(object = mBC, num.replicate = 100, num.pc = 15,do.par=TRUE, num.cores=16)
mBC = JackStrawPlot(object = mBC, PCs = 1:15)
dev.off()
closeAllConnections()
gc()

#plot each principal component-associated genes
tmp=mBC@dr$pca@jackstraw@emperical.p.value.full
tmp = as.data.frame(mBC@dr$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05,1]
dims= c(1:(min(tmp1)-1))

#set dims 1:5 to reflect diversity
dims=1:5

mBC = FindClusters(object = mBC, reduction.type = "pca", dims.use = dims,
                   resolution = 2.0, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

py_set_seed(42)
mBC = DoFItSNE(mBC, reduction_use = "pca", dims_use = as.integer(dims),
               perplexity=10, n_jobs=as.integer(16), df = 1)

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

file.name="DC_IM_FItSNE_perplex10_reso1.png"
ggsave(file = file.name, plot = plot_grid(p1, p2), dpi = 100, width = 10, height = 5)


mBC = AddMetaData(mBC, metadata = c(hoge2, hoge3), col.name="celltype")
save.pigz(mBC, file = "DC_IM_cluster20.rda", n.cores=8)


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
fwrite(mBC.markers_write, "DC_IM_ALLmarkers_minpct0.1_Adj_p0.05.txt", row.names=F, col.names=T, sep="\t", quote=F)

top20 = mBC.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
top20 = as.data.frame(top20)
top20  = top20 [!duplicated(top20$gene),]
top20 = top20 %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
top20 = as.data.frame(top20)
file.name=paste(dir.name, sample.name, "marker_res2.png", sep='')
p = DoHeatmap2(object = mBC, genes.use = top20$gene, genes.ident = top20$cluster,
               slim.col.label = TRUE, remove.key = FALSE, cex.row=4, disp.min = -2.5, disp.max = 2.5)
ggsave(file = "DC_IM_marker.res2.png", plot = p, device="png", units="in", dpi = 300, width = 5, height = 7, limitsize=FALSE)


#export the number of each subsets and save file
tmp = table(mBC@ident, mBC@meta.data$orig.ident)
write.table(tmp, "./Lung_day00_Rhapsody_subsetNumber_DC_IM.txt", row.names=T, sep="\t", quote=F)
save.pigz(mBC, file = "DC_IM_cluster20.rda", n.cores=8)

#######
#annotate each cell subsets 
res=NULL
load("DC_IM_cluster20.rda")
hoge = mBC@meta.data$res.2
hoge[hoge=="1"]="DC_cDC2"
hoge[hoge=="0"]="Mac_IM"
hoge[hoge=="2"]="DC_cDC1"
hoge[hoge=="3"]="DC_cDC1"
names(hoge)=rownames(mBC@meta.data)
res = c(res, hoge)

load("ILC_gdT_cluster35.rda")

hoge = as.character(mBC@meta.data$celltype)
names(hoge)=rownames(mBC@meta.data)
res = c(res, hoge)


hoge1 = Lung_BDRhapsody@meta.data[!rownames(Lung_BDRhapsody@meta.data)%in% names(res),]
hoge2 = as.character(hoge1$res.5)

hoge2[hoge2 %in% c("0", "4")]="Epi_AT2"
hoge2[hoge2 %in% c("1", "11","16", "19", "22")]="FB_Inmthi"
hoge2[hoge2 %in% c("23")]="FB_Dcnhi"
hoge2[hoge2 %in% c("6", "36")]="doublet"
hoge2[hoge2 %in% c("2", "3", "5", "10")]="Endo_VEC"
hoge2[hoge2 %in% c("30")]="Endo_arterial"
hoge2[hoge2 %in% c("28")]="Endo_venous"
hoge2[hoge2 %in% c("34")]="Endo_LEC"
hoge2[hoge2 %in% c("18", "21")]="Endo_capillary"
hoge2[hoge2 %in% c("27")]="Epi_ciliated"
hoge2[hoge2 %in% c("9")]="Epi_AT1"
hoge2[hoge2 %in% c("8", "12", "24")]="Bcell"
hoge2[hoge2 %in% c("7")]="Mac_AM"
hoge2[hoge2 %in% c("14")]="Pericyte"
hoge2[hoge2 %in% c("13")]="SMC"
hoge2[hoge2 %in% c("17")]="Tcell_CD4T"
hoge2[hoge2 %in% c("25")]="Tcell_CD8T"
hoge2[hoge2 %in% c("33")]="NK"
hoge2[hoge2 %in% c("29")]="Mo_Ly6Chi"
hoge2[hoge2 %in% c("26")]="Mo_Ly6Clo"
hoge2[hoge2 %in% c("31")]="Neutrophil"
hoge2[hoge2 %in% c("32")]="Mesothelial"
hoge2[hoge2 %in% c("15")]="Epi_clara"
names(hoge2)=rownames(Lung_BDRhapsody@meta.data[!rownames(Lung_BDRhapsody@meta.data)%in% names(res),])

res = c(hoge2, res)

Lung_BDRhapsody=AddMetaData(Lung_BDRhapsody, metadata=res, col.name = "celltype")
save.pigz(Lung_BDRhapsody, file="Lung_day00_Rhapsody.rda", n.cores=16)




