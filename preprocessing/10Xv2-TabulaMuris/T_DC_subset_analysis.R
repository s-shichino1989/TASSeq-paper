##smart-seq2 subsetting

mBC = qread("./Seurat/UTlung_10Xv2_TabulaMuris_Seurat.qs")

mBC1 = mBC

#Tcell

mBC = subset(mBC1, idents=c(1))
mBC = FindVariableFeatures(mBC, selection.method = "mvp", 
                           mean.cutoff=c(0.1, Inf), dispersion.cutoff=c(0.5, Inf))
message("Scaling data...")
mBC = ScaleData(object = mBC, vars.to.regress = c("nReads"))
hvg.number = length(x = mBC@assays$RNA@var.features)
mBC = RunPCA(object = mBC, features = VariableFeatures(object = mBC), 
             npcs=30, verbose = FALSE)


mBC = JackStraw(mBC, num.replicate = 100, dims = 30, verbose = FALSE)
mBC = ScoreJackStraw(mBC, dims = 1:30, score.thresh = 0.05)
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

dims=1:10 #manually set dimension for better CD4/CD8 separation
#perform FItSNE dimReduction
quiet(py_set_seed(42))
mBC = DoFItSNE(mBC, reduction_use = "pca", dims_use = as.integer(dims),
               perplexity=perplexity, n_jobs=as.integer(detectCores()/2), df = 0.9)
quiet(py_set_seed(42))
mBC = RunUMAP(mBC, reduction = "pca", dims = dims)
mBC = FindNeighbors(mBC, dims=dims)
mBC = FindClusters(object = mBC, resolution = 2, random.seed=42, verbose = FALSE)

DimPlot(mBC, pt.size=1, label=TRUE)
FeaturePlot(mBC, features=c("Cd4", "Cd8a", "Foxp3"))


write.table(table(mBC@active.ident, mBC@meta.data$orig.ident), file="Tcell_gate1_Cellnumber.txt", sep="\t", quote=F, row.names=T, col.names=T)

save(mBC, file="Tcell_gate.rda")


#DC
mBC = subset(mBC1, idents=c(8))
mBC = FindVariableFeatures(mBC, selection.method = "mvp", 
                           mean.cutoff=c(0.1, Inf), dispersion.cutoff=c(0.5, Inf))
message("Scaling data...")
mBC = ScaleData(object = mBC, vars.to.regress = c("nReads"))
hvg.number = length(x = mBC@assays$RNA@var.features)
mBC = RunPCA(object = mBC, features = VariableFeatures(object = mBC), 
             npcs=30, verbose = FALSE)


mBC = JackStraw(mBC, num.replicate = 100, dims = 30, verbose = FALSE)
mBC = ScoreJackStraw(mBC, dims = 1:30, score.thresh = 0.05)
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

dims=1:10
#perform FItSNE dimReduction
quiet(py_set_seed(42))
mBC = DoFItSNE(mBC, reduction_use = "pca", dims_use = as.integer(dims),
               perplexity=perplexity, n_jobs=as.integer(detectCores()/2), df = 0.9)
quiet(py_set_seed(42))
mBC = RunUMAP(mBC, reduction = "pca", dims = dims)

mBC = FindNeighbors(mBC, dims=dims)
mBC = FindClusters(object = mBC, resolution = 1.2, random.seed=42, verbose = FALSE)

DimPlot(mBC, pt.size=1, label=TRUE)
FeaturePlot(mBC, features=c("Itgae", "Cd209a", "Cx3cr1", "C1qa", "Ly6c2", "Ccr7", "Siglech"), pt.size = 0.5)


write.table(table(mBC@active.ident, mBC@meta.data$orig.ident), file="DC_IM_gate1_Cellnumber.txt", sep="\t", quote=F, row.names=T, col.names=T)
save(mBC, file="DC_IM_gate1.rda")




#DC gate2
mBC = subset(mBC, idents=c(1))
mBC = FindVariableFeatures(mBC, selection.method = "mvp", 
                           mean.cutoff=c(0.1, Inf), dispersion.cutoff=c(0.5, Inf))
message("Scaling data...")
mBC = ScaleData(object = mBC, vars.to.regress = c("nReads"))
hvg.number = length(x = mBC@assays$RNA@var.features)
mBC = RunPCA(object = mBC, features = VariableFeatures(object = mBC), 
               npcs=30, verbose = FALSE)


mBC = JackStraw(mBC, num.replicate = 100, dims = 30, verbose = FALSE)
mBC = ScoreJackStraw(mBC, dims = 1:30, score.thresh = 0.05)
#return p-value of Jackstraw analysis
JS_res = mBC@reductions$pca@jackstraw$overall.p.values
tmp = JS_res[JS_res[,2]> 1e-05,,drop=F]
dims =c(1:(min(tmp[,1])-1))


  perplexity = 10 


#perform FItSNE dimReduction
quiet(py_set_seed(42))
mBC = DoFItSNE(mBC, reduction_use = "pca", dims_use = as.integer(dims),
               perplexity=perplexity, n_jobs=as.integer(detectCores()/2), df = 0.9)
quiet(py_set_seed(42))
mBC = RunUMAP(mBC, reduction = "pca", dims = dims)

mBC = FindNeighbors(mBC, dims=dims)
mBC = FindClusters(object = mBC, resolution = 1, random.seed=42, verbose = FALSE)

DimPlot(mBC, pt.size=2, label=TRUE)
FeaturePlot(mBC, features=c("Itgae", "Cd209a", "Ccr7", "Siglech"), pt.size = 0.5)


write.table(table(mBC@active.ident, mBC@meta.data$orig.ident), file="DC_gate2_Cellnumber.txt", sep="\t", quote=F, row.names=T, col.names=T)
save(mBC, file="DC_gate2.rda")

