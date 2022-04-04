##smart-seq2 subsetting

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 


mBC = qread("./Seurat/UTlung_shallow_Seurat.qs")
mBC1 = mBC

#DC

mBC = subset(mBC1, idents=c(24))
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

resolution = c(0.6, 1.0, 1.5, 2.0, 2.5, 3.0)
mBC = FindNeighbors(mBC, dims=dims)
mBC = FindClusters(object = mBC, resolution = 1.8, random.seed=42, verbose = FALSE)
DimPlot(mBC, reduction = "FItSNE", pt.size=1, label = TRUE)
FeaturePlot(mBC, features=c("C1qa", "Cd209a", "Itgae", "Siglech", "Fscn1"), reduction="FItSNE")

write.table(table(mBC@active.ident, mBC@meta.data$orig.ident), file="DC_gate1_Cellnumber.txt", sep="\t", quote=F, row.names=T, col.names=T)
save(mBC, file="DC_IM_gate.rda")

