##TAS-Seq lung dataset2 subsetting

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 


mBC = qread("./Seurat/day00UT_results_Seurat.qs")
mBC = FindClusters(mBC, resolution=3.5, random.seed=42)

mBC1 = mBC

#DC

hoge = mBC1@assays$RNA@counts[,rownames(mBC1@meta.data[mBC1@meta.data$RNA_snn_res.3.5 %in% c(22),])]
mBC = CreateSeuratObject(hoge)
mBC = NormalizeData(mBC, scale.factor=1000000)
mBC = FindVariableFeatures(mBC, selection.method = "mvp", 
                           mean.cutoff=c(0.1, Inf), dispersion.cutoff=c(0.5, Inf))
message("Scaling data...")

colnames(mBC@meta.data) = c("orig.ident", "nReads", "nGene")
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


if(nrow(mBC@meta.data)/200 >=500){
  perplexity = 500
} else if (nrow(mBC@meta.data)/200 >=300){
  perplexity = 400 
} else if (nrow(mBC@meta.data)/200 >=200){
  perplexity = 300 
} else if (nrow(mBC@meta.data)/200 >=100){
  perplexity = 200 
} else if (nrow(mBC@meta.data)/200 >=50){
  perplexity = 100   
}else if (nrow(mBC@meta.data)/200 >=25){
  perplexity = 50   
}else if (nrow(mBC@meta.data)/200 >=12.5){
  perplexity = 25   
}else {
  perplexity = 10 
}

#perform FItSNE dimReduction
quiet(py_set_seed(42))
mBC = DoFItSNE(mBC, reduction_use = "pca", dims_use = as.integer(dims),
               perplexity=perplexity, n_jobs=as.integer(detectCores()/2), df = 0.9)
quiet(py_set_seed(42))
mBC = RunUMAP(mBC, reduction = "pca", dims = dims)

mBC = FindNeighbors(mBC, dims=dims)
mBC = FindClusters(object = mBC, resolution = 1.0, random.seed=42, verbose = FALSE)
DimPlot(mBC, reduction = "umap", pt.size=1, label = TRUE)
FeaturePlot(mBC, features=c("Xcr1", "Cd209a", "Siglech", "Fcgr1", "Cx3cr1", "Ccr7"), reduction="umap")

write.table(table(mBC@active.ident, mBC@meta.data$orig.ident), file="cDC1_cDC2_gate1_Cellnumber.txt", sep="\t", quote=F, row.names=T, col.names=T)
save(mBC, file="cDC1_cDC2_gate.rda")

