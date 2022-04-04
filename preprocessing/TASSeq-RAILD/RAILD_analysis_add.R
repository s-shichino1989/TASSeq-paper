mBC = FindClusters(object = mBC, resolution = 2.5, verbose=FALSE, random.seed=42)

mBC.markers = rDBEC::FindMarkers_parallel_lite(mBC, AllcellsIdent=mBC@active.ident,
                                               test.use="wilcox", only.pos=TRUE, min.pct=0.1,
                                               features.use = NULL, nthreads = 26,
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

fwrite(mBC.markers_write, "ALLmarkers_minpct0.1_Adj_p0.05.txt", row.names=F, col.names=T, sep="\t", quote=F)


tmp = table(mBC@active.ident, mBC@meta.data$orig.ident)
tmp = as.data.frame(tmp)

tmp1 = tmp %>% tidyr::spread(key = Var2, value=Freq) 
colnames(tmp1)[1]="Seurat_Clusters"
hoge = c("doublet", "not-detected")
tmp1 = tmp1[,!colnames(tmp1) %in% hoge]

write.table(tmp1, "subsetNumber_RAILD.txt", row.names=F, sep="\t", quote=F)

