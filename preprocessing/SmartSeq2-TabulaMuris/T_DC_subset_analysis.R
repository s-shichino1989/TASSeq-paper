##smart-seq2 subsetting

load("./Seurat/UTlung_SmartSeq2_Seurat.rda")

mBC1 = mBC

#Tcell

mBC = subset(mBC1, idents=c(14))

mBC = FindVariableFeatures(mBC, selection.method = "vst",  
                           mean.cutoff=c(0.1, Inf), dispersion.cutoff=c(0.5, Inf))
mBC = ScaleData(object = mBC, vars.to.regress = c("nReads"))
hvg.number = length(x = mBC@assays$RNA@var.features)
mBC = RunPCA(object = mBC, features = VariableFeatures(object = mBC), 
             npcs=20, verbose = FALSE)

mBC = JackStraw(mBC, num.replicate = 100, dims = 20, verbose = FALSE)
mBC = ScoreJackStraw(mBC, dims = 1:20, score.thresh = 0.05)
JS_res = mBC@reductions$pca@jackstraw$overall.p.values
tmp = JS_res[JS_res[,2]> 1e-05,,drop=F]
dims =c(1:(min(tmp[,1])-1))
dims=1:4
mBC = RunUMAP(mBC, reduction = "pca", dims = dims)
mBC = FindNeighbors(mBC, dims=dims)
mBC = FindClusters(object = mBC, resolution = 1.2, random.seed=42, verbose = FALSE)
DimPlot(mBC, pt.size=1, label=TRUE)
FeaturePlot(mBC, features=c("Cd4", "Cd8a", "Cd8b1"), pt.size = 1)

write.table(table(mBC@active.ident, mBC@meta.data$orig.ident), file="T_Cellnumber.txt", sep="\t", quote=F, row.names=T, col.names=T)
save(mBC, file="T_gate.rda")

#IM/DC

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

perplexity = 10 
dims=1:10

#perform FItSNE dimReduction
quiet(py_set_seed(42))
mBC = DoFItSNE(mBC, reduction_use = "pca", dims_use = as.integer(dims),
               perplexity=perplexity, n_jobs=as.integer(detectCores()/2), df = 0.9)
quiet(py_set_seed(42))
mBC = RunUMAP(mBC, reduction = "pca", dims = dims)
mBC = FindNeighbors(mBC, dims=dims)
mBC = FindClusters(object = mBC, resolution = 2, random.seed=42, verbose = FALSE)
DimPlot(mBC, pt.size=1, label=TRUE)
FeaturePlot(mBC, features=c("Cd209a", "Ccr7", "Itgae", "Siglech"), pt.size = 1)

write.table(table(mBC@active.ident, mBC@meta.data$orig.ident), file="DC_Cellnumber.txt", sep="\t", quote=F, row.names=T, col.names=T)
save(mBC, file="DC_gate.rda")

FeaturePlot(mBC, features=c("Cd209a", "Ccr7", "Itgae", "Siglech"), pt.size = 1)

