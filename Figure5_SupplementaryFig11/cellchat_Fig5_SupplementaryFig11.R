library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(tibble)
library(tidyr)
library(dplyr)
library(pheatmap)
library(future)
options(stringsAsFactors = FALSE)
suppressWarnings(suppressMessages(source("./library_source_Seurat.R")))
species="mmu"
set.seed(42)

###custom functions
computeCommunProb2 <- function(object, type = c("triMean", "truncatedMean", "median", "mean"),
                               trim = NULL, LR.use = NULL, raw.use = FALSE,
                               population.size = FALSE, nboot = 100, seed.use = 1L, Kh = 0.5, n = 1,
                               nthreads=4, seed=42) {
  type <- match.arg(type)
  FunMean <- switch(type,
                    triMean = triMean,
                    truncatedMean = function(x) mean(x, trim = trim, na.rm = TRUE),
                    median = function(x) median(x, na.rm = TRUE),
                    mean =  function(x) mean(x, na.rm = TRUE))
  if (raw.use) {
    data <- as.matrix(object@data.signaling)
  } else {
    data <- object@data.project
  }
  if (is.null(LR.use)) {
    pairLR.use <- object@LR$LRsig
  } else {
    pairLR.use <- LR.use
  }
  complex_input <- object@DB$complex
  cofactor_input <- object@DB$cofactor
  my.sapply <- sapply

  ptm = Sys.time()

  pairLRsig <- pairLR.use
  group <- object@idents
  geneL <- as.character(pairLRsig$ligand)
  geneR <- as.character(pairLRsig$receptor)
  nLR <- nrow(pairLRsig)
  numCluster <- nlevels(group)
  if (numCluster != length(unique(group))) {
    stop("Please check `unique(object@idents)` and ensure that the factor levels are correct!
         You may need to drop unused levels using 'droplevels' function. e.g.,
         `meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))`")
  }
  # if (all(data[1:5, ] == floor(data[1:5, ]))) {
  #   stop("Please check your input data matrix and ensure that you use the normalized data instead of count data!")
  # }


  data.use <- data/max(data)
  nC <- ncol(data.use)

  # compute the expression of ligand and receptor
  dataL <- computeExpr_LR(geneL, data.use, complex_input)
  dataR <- computeExpr_LR(geneR, data.use, complex_input)
  # take account into the effect of co-activation and co-inhibition receptors
  dataR.co.A.receptor <- computeExpr_coreceptor(cofactor_input, data.use, pairLRsig, type = "A")
  dataR.co.I.receptor <- computeExpr_coreceptor(cofactor_input, data.use, pairLRsig, type = "I")
  dataR <- dataR * dataR.co.A.receptor/dataR.co.I.receptor
  # compute the average expression in each cell group
  dataLavg <- aggregate(t(dataL), list(group), FUN = FunMean)
  dataLavg <- t(dataLavg[,-1])
  rownames(dataLavg) <- geneL
  dataRavg <- aggregate(t(dataR), list(group), FUN = FunMean)
  dataRavg <- t(dataRavg[,-1])
  rownames(dataRavg) <- geneR

  dataL.binary = (dataL > 0)*1 ;dataR.binary = (dataR > 0)*1
  dataLavg2 <- aggregate(t(dataL.binary), list(group), FUN = sum)
  dataLavg2 <- t(dataLavg2[,-1])/nC
  dataRavg2 <- aggregate(t(dataR.binary), list(group), FUN = sum)
  dataRavg2 <- t(dataRavg2[,-1])/nC

  # compute the expression of agonist and antagonist
  index.agonist <- which(!is.na(pairLRsig$agonist) & pairLRsig$agonist != "")
  index.antagonist <- which(!is.na(pairLRsig$antagonist) & pairLRsig$antagonist != "")
  # quantify the communication probability
  set.seed(seed.use)
  permutation <- replicate(nboot, sample.int(nC, size = nC))
  Prob <- array(0, dim = c(numCluster,numCluster,nLR))
  Pval <- array(0, dim = c(numCluster,numCluster,nLR))

  cl = parallel::makeCluster(nthreads, type="PSOCK", useXDR=FALSE)
  doParallel::registerDoParallel(cl)
  doRNG::registerDoRNG(seed)
  res <- foreach::foreach(i = 1:nLR, .packages=c("future", "future.apply",
                                                 "Matrix", "pbapply",
                                                 "stats"),
                          .multicombine=TRUE,
                          .maxcombine=nLR,
                          .inorder=TRUE,
                          .noexport=setdiff(ls(),c("dataLavg", "dataRavg", "Kh", "n", "numCluster", "Prob", "Pval",
                                                   "index.agonist", "index.antagonist","data.use",
                                                   "pairLRsig","cofactor_input", "group","FunMean","dataLavg2",
                                                   "dataRavg2","population.size","dataL","dataR",
                                                   "dataL.binary","dataR.binary","my.sapply", "computeCommunProb_internal",
                                                   "nboot","permutation","nC", "trim")),
                          .options.RNG=seed
  ) %dopar% {
    computeCommunProb_internal(dataLavg = dataLavg,
                               dataRavg = dataRavg,
                               i,
                               Kh=Kh,
                               n=n,
                               numCluster=numCluster,
                               index.agonist = index.agonist,
                               index.antagonist = index.antagonist,
                               data.use = data.use,
                               pairLRsig = pairLRsig,
                               cofactor_input = cofactor_input,
                               group = group,
                               FunMean = FunMean,
                               dataLavg2 = dataLavg2,
                               dataRavg2 = dataRavg2,
                               population.size = population.size,
                               dataL = dataL,
                               dataR = dataR,
                               dataL.binary = dataL.binary,
                               dataR.binary = dataR.binary,
                               my.sapply = my.sapply,
                               nboot = nboot,
                               permutation = permutation,
                               nC = nC,
                               trim = trim,
                               Prob = Prob,
                               Pval = Pval)
  }
  parallel::stopCluster(cl)
  for (i in c(1:nLR)){
    Pval[, , i] = res[[i]][[1]]
    Prob[ , , i] = res[[i]][[2]]
  }

  Pval[Prob == 0] <- 1
  dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
  dimnames(Pval) <- dimnames(Prob)
  net <- list("prob" = Prob, "pval" = Pval)
  execution.time = Sys.time() - ptm
  object@options$run.time <- as.numeric(execution.time, units = "secs")
  object@options$parameter <- list(type.mean = type, trim = trim, raw.use = raw.use, population.size = population.size,  nboot = nboot, seed.use = seed.use, Kh = Kh, n = n)
  object@net <- net
  return(object)
}



#' internal function for computing the communication probability/strength between any interacting cell groups

computeCommunProb_internal <- function(dataLavg,
                                       dataRavg,
                                       i,
                                       Kh=0.5,
                                       n=1,
                                       numCluster,
                                       index.agonist,
                                       index.antagonist,
                                       data.use,
                                       pairLRsig,
                                       cofactor_input,
                                       group,
                                       FunMean,
                                       dataLavg2,
                                       dataRavg2,
                                       population.size = FALSE,
                                       dataL,
                                       dataR,
                                       dataL.binary,
                                       dataR.binary,
                                       my.sapply,
                                       nboot = 100,
                                       permutation,
                                       nC,
                                       trim,
                                       Prob,
                                       Pval) {
  dataLR <- Matrix::crossprod(matrix(dataLavg[i,], nrow = 1), matrix(dataRavg[i,], nrow = 1))
  P1 <- dataLR^n/(Kh^n + dataLR^n)
  if (sum(P1) == 0) {
    Pnull = P1
    Prob_res <- Pnull
    p = 1
    Pval_res <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
  } else {
    # agonist and antagonist
    if (is.element(i, index.agonist)) {
      data.agonist <- computeExprGroup_agonist(data.use = data.use, pairLRsig, cofactor_input, group = group,index.agonist = i, Kh = Kh, FunMean = FunMean, n = n)
      P2 <- Matrix::crossprod(matrix(data.agonist, nrow = 1))
    } else {
      P2 <- matrix(1, nrow = numCluster, ncol = numCluster)
    }
    if (is.element(i, index.antagonist)) {
      data.antagonist <- computeExprGroup_antagonist(data.use = data.use, pairLRsig, cofactor_input, group = group, index.antagonist = i, Kh = Kh, FunMean = FunMean, n = n)
      P3 <- Matrix::crossprod(matrix(data.antagonist, nrow = 1))
    } else {
      P3 <- matrix(1, nrow = numCluster, ncol = numCluster)
    }
    # number of cells
    if (population.size) {
      P4 <- Matrix::crossprod(matrix(dataLavg2[i,], nrow = 1), matrix(dataRavg2[i,], nrow = 1))
    } else {
      P4 <- matrix(1, nrow = numCluster, ncol = numCluster)
    }

    Pnull = P1*P2*P3*P4
    Prob_res <- Pnull

    Pnull <- as.vector(Pnull)
    dataL.i <- dataL[i,]; dataR.i <- dataR[i,];
    dataL2.i <- dataL.binary[i,]; dataR2.i <- dataR.binary[i,];
    #Pboot <- foreach(nE = 1:nboot) %dopar% {
    Pboot <- my.sapply(
      X = 1:nboot,
      FUN = function(nE) {
        groupboot <- group[permutation[, nE]]
        dataLavgB <- aggregate(matrix(dataL.i, ncol = 1), list(groupboot), FUN = FunMean)
        dataLavgB <- t(dataLavgB[,-1])
        dataLavgB <- matrix(dataLavgB, nrow = 1)

        dataRavgB <- aggregate(matrix(dataR.i, ncol = 1), list(groupboot), FUN = FunMean)
        dataRavgB <- t(dataRavgB[,-1])
        dataRavgB <- matrix(dataRavgB, nrow = 1)
        dataLRB = Matrix::crossprod(dataLavgB, dataRavgB)
        P1.boot <- dataLRB^n/(Kh^n + dataLRB^n)
        # agonist and antagonist
        if (is.element(i, index.agonist)) {
          data.agonist <- computeExprGroup_agonist(data.use = data.use, pairLRsig, cofactor_input, group = groupboot, index.agonist = i, Kh = Kh, FunMean = FunMean, n = n)
          P2.boot <- Matrix::crossprod(matrix(data.agonist, nrow = 1))
        } else {
          P2.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
        }
        if (is.element(i, index.antagonist)) {
          data.antagonist <- computeExprGroup_antagonist(data.use = data.use, pairLRsig, cofactor_input, group = groupboot,index.antagonist = i, Kh = Kh, FunMean = FunMean, n= n)
          P3.boot <- Matrix::crossprod(matrix(data.antagonist, nrow = 1))
        } else {
          P3.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
        }
        dataLavg2B <- by(matrix(dataL2.i, ncol = 1), groupboot, sum)/nC
        dataLavg2B <- matrix(dataLavg2B, nrow = 1)

        dataRavg2B <- by(matrix(dataR2.i, ncol = 1), groupboot, sum)/nC
        dataRavg2B <- matrix(dataRavg2B, nrow = 1)
        if (population.size) {
          P4.boot = Matrix::crossprod(dataLavg2B, dataRavg2B)
        } else {
          P4.boot = matrix(1, nrow = numCluster, ncol = numCluster)
        }

        Pboot = P1.boot*P2.boot*P3.boot*P4.boot
        return(as.vector(Pboot))
      }
    )
    Pboot <- matrix(unlist(Pboot), nrow=length(Pnull), ncol = nboot, byrow = FALSE)
    nReject <- rowSums(Pboot - Pnull >= 0)
    p = nReject/nboot
    Pval_res <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
  }
  res = list()
  res[[1]] = Pval_res
  res[[2]] = Prob_res
  return(res)
}


#outgoing incoming signaling contribution
netAnalysis_signalingRole_heatmap1 <- function(object, signaling = NULL,
                                               pattern = c("outgoing", "incoming","all"),
                                               slot.name = "netP"
){
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[,i] <- centr[[i]]$outdeg
    incoming[,i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  } else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  } else if (pattern == "all") {
    mat <- t(outgoing+ incoming)
    legend.name <- "Overall"
  }

  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  return(mat)
}



CellChatDB.use <- CellChatDB.mouse


##download annotated Seurat objects for mouse lung data
system("wget -L -O matrix_inflection_10Xv2-GSM3926540.txt.gz 
       https://tus.box.com/shared/static/wf4whgfdj2bx0kpbe05r66ntq1p7zltl.gz")
system("wget -L -O UTlung_deep_Seurat_annot.qs
       https://tus.box.com/shared/static/ohoohhgf3t12tcmupueajbbapu4xdxm6.qs")
system("wget -L -O UTlung_shallow_Seurat_annot.qs
       https://tus.box.com/shared/static/rhlcprxv2l3mtmulphsotwl2serl6u7o.qs")
system("wget -L -O UTlung_SmartSeq2_Seurat_annot.qs
       https://tus.box.com/shared/static/ej8tff4k7szdfqsw6niljnzrronnx910.qs")
system("wget -L -O UTlung_10Xv2_TabulaMuris_Seurat_annot.qs
       https://tus.box.com/shared/static/fqiabrfwxd7y15jgcwlovelhoz2tgz65.qs")


fnames = dir(pattern="UTlung_")
seu_list = lapply(fnames, qread, nthreads=24)

#rename Smart-seq2 orig.ident
seu_list[[5]]@meta.data$orig.ident = rep("Smart-seq2", nrow(seu_list[[5]]@meta.data))

#rename 10X v2 Tabula Muris data
seu_list[[2]]@meta.data$orig.ident = rep("10Xv2_TabulaMuris", nrow(seu_list[[2]]@meta.data))


#remove doublets
for(i in c(1:5)){
  doublet_key = seu_list[[i]]@meta.data$celltype %in% c("doublet", "misc")
  seu_list[[i]]=subset(seu_list[[i]], cells = rownames(seu_list[[i]]@meta.data[!doublet_key,]))
}

#merge data
for(i in c(1:5)){
  if(i>1){
    seu = merge(x=seu, y=seu_list[[i]])
  } else {
    seu = seu_list[[1]]
  }
}


#reorder identities
idents = unique(seu@meta.data$orig.ident)
idents = idents[c(1:2,9,6:8,3:5)]


Ncells = NULL
for (i in c(1:length(idents))){
  Ncells = c(Ncells, nrow(seu@meta.data[seu@meta.data$orig.ident == idents[i],]))
}


tmp = seu@meta.data$orig.ident
names(tmp)=rownames(seu@meta.data)
tmp = factor(tmp, levels = sort(unique(as.character(tmp))))
seu@active.ident = tmp

norm.data.input=list()
raw.data.input=list()
meta=list()

# simply use the default CellChatDB
#CellChatDB.use <- CellChatDB
CellChatDB.use.symbols = unique(CellChatDB.use$geneInfo$Symbol)

for (i in c(1:length(idents))){
  seu1 = subset(seu, idents = idents[i], downsample= min(Ncells))
  tmp = as.character(seu1@meta.data$celltype)
  names(tmp)=rownames(seu1@meta.data)
  tmp = factor(tmp)
  seu1@active.ident = tmp
  seu_list[[i]]=seu1
  norm.data.input[[i]] = seu1@assays$RNA@data[intersect(rownames(seu1@assays$RNA@data), CellChatDB.use.symbols),]
  norm.data.input[[i]] = norm.data.input[[i]][rowSums(norm.data.input[[i]])>0,]
  raw.data.input[[i]] = seu1@assays$RNA@counts[intersect(rownames(seu1@assays$RNA@data), CellChatDB.use.symbols),]
  raw.data.input[[i]] = raw.data.input[[i]][rowSums(raw.data.input[[i]])>0,]
  meta[[i]] = seu1@meta.data
  # a dataframe with rownames containing cell mata data
  meta[[i]]$celltype = as.character(meta[[i]]$celltype)
}

cellchat.list=list()
marker.list=list()
species="mmu"


#calculate thresh.pc = 0.05, 0.15, 0.25, 0.50, 0.75 for Gene filtering

for (j in c(0.05, 0.15, 0.25, 0.50, 0.75)){
  for (i in c(1:length(idents))){
    cellchat = createCellChat(object = norm.data.input[[i]], meta = meta[[i]],
                              group.by = "celltype") #specify Celltype column name
    cellchat@data.raw = raw.data.input[[i]]
    groupSize = as.numeric(table(cellchat@idents)) # number of cells in each cell group
    cellchat@DB = CellChatDB.use
    cellchat = CellChat::subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
    future::plan("multisession", workers = 32) # do parallel
    cellchat = identifyOverExpressedGenes(cellchat, thresh.p=0.05, thresh.pc = j)
    cellchat = identifyOverExpressedInteractions(cellchat)
    if(species=="mmu"){
      cellchat <- projectData(cellchat, PPI.mouse)
    } else if(species=="hsa"){
      cellchat <- projectData(cellchat, PPI.human)
    }
    future::plan("sequential") # paralell = OFF
    gc()
    message("DE analysis done")
    future::plan("multisession", workers = 4) # do parallel
    cellchat <- computeCommunProb2(cellchat, raw.use = TRUE, type = "mean",population.size = TRUE, nthreads=16)
    gc()
    closeAllConnections()
    message("commune probability analysis done")
    future::plan("multisession", workers = 16) # do parallel
    cellchat = filterCommunication(cellchat, min.cells = 5)
    cellchat = computeCommunProbPathway(cellchat, thresh = 0.05)
    cellchat = aggregateNet(cellchat, thresh = 0.05)
    cellchat = netAnalysis_computeCentrality(cellchat, slot.name="netP")
    future::plan("sequential") # paralell = OFF
    gc()
    cellchat.list[[i]]=cellchat
  }

  #naming cellchat objects
  names(cellchat.list) = idents

  #save Cellchat object list
  file.name=sprintf("%s.qs", j)
  file.name=paste("cellchat_list_thre", file.name, sep="")
  qsave(cellchat.list, file=file.name, nthreads=24)


}

percent = c(0.05, 0.15, 0.25, 0.50, 0.75)
res_mat_list=list()

#extract significant signaling pathways of cell-cell interactions of each dataset (Supplementary Figure 11a)
for (k in c(1:5)){
  j = percent[k]
  cellchat.list = qread(paste0("./cellchat_list_thre", j, ".qs"))
  mat=list()

  for (i in c(1:length(idents))){
    ht1 <- netAnalysis_signalingRole_heatmap1(cellchat.list[[i]], pattern = "all")
    ht1 <- sweep(ht1, 1L, apply(ht1, 1, max), '/', check.margin = FALSE)
    ht1 = as.data.frame(ht1)
    ht1 = ht1[rowSums(ht1)>0,]
    mat[[i]] = ht1
  }

  #visualize overlapped signaling pathways

  for (i in c(1:length(idents))){
    tmp = rownames(mat[[i]])
    tmp = as.data.frame(tmp)
    rownames(tmp)=tmp[,1]
    colnames(tmp)=idents[i]
    if(i>1){
      res_mat = dplyr::full_join(rownames_to_column(res_mat), rownames_to_column(tmp), by ="rowname")
    } else {
      res_mat = rownames_to_column(tmp)
    }
    rownames(res_mat)=res_mat$rowname
    res_mat = res_mat[,!colnames(res_mat) %in% "rowname", drop=F]
  }
  res_mat[is.na(res_mat)]=0
  hoge1=res_mat
  hoge1[hoge1!="0"]=1
  hoge1[hoge1=="0"]=0

  for(i in c(1:ncol(hoge1))){
    hoge1[,i]=as.numeric(hoge1[,i])
  }

  write.table(hoge1, file=paste0("./cellchat_pathways_ALL_thre", j, ".txt"), sep="\t", quote=F, row.names=T, col.names=T)

  hoge2 = hoge1[rowMin(as.matrix(hoge1))==0,]
  #hoge2 = hoge1[order(rownames(hoge1)),]

  m2b2g = colorRampPalette(c("grey","magenta"))
  p = pheatmap(as.data.frame(hoge2),
               scale="none",
               cluster_rows=TRUE,
               cluster_cols=FALSE,
               show_colnames=TRUE,
               show_rownames=TRUE,
               clustering_distance_rows = "euclidean",
               clustering_method = "ward.D",
               silent=TRUE, main=paste0("pct", j),
               color = m2b2g(1024)
  )

  hoge2 = hoge2[p$tree_row$order,]
  colnames(hoge2)=idents
  p = pheatmap(as.data.frame(hoge2),
               scale="none",
               cluster_rows=FALSE,
               cluster_cols=FALSE,
               show_colnames=TRUE,
               show_rownames=TRUE,
               silent=TRUE, main=paste0("pct", j),
               color = m2b2g(1024)
  )


  p = as.ggplot(p[[4]])

  ggsave(file = paste0("pct", j, "_pathways.png"),
         plot =p, device="png",
         units="in", dpi = 300,
         width = 3.5, height = round(max(8.5 * nrow(hoge1)/80, 5.5), digits=1), limitsize=FALSE)

  res_mat_list[[j]]=res_mat
}

qsave(res_mat_list, file="pathways_cellchat_list.qs", nthreads=4)


#create table for interaction numbers (Figure 5a)
######
percent = c(0.05, 0.15, 0.25, 0.50, 0.75)
res_mat_list=list()

for (k in c(1:5)){
  j = percent[k]
  cellchat.list = qread(paste0("./cellchat_list_thre", j, ".qs"))
  group.new=NULL
  for (i in c(1:length(cellchat.list))){
    group.new = union(group.new, levels(cellchat.list[[i]]@idents))
  }

  for (i in c(1:length(cellchat.list))){
    cellchat.list[[i]] = suppressMessages(suppressWarnings(liftCellChat(cellchat.list[[i]],
                                                                        group.new)))
  }

  cellchat = mergeCellChat(cellchat.list, add.names=idents)

  #extract total interaction number
  df <- as.data.frame(sapply(cellchat@net, function(x) sum(x$count)))
  colnames(df)=paste0("threshold_", j)
  if(k>1){
    interaction_res_df=cbind(interaction_res_df, df)
  } else {
    interaction_res_df = df
  }
}

#export result
write.table(interaction_res_df, file="./interaction_number.txt", sep="\t", quote=F, row.names=T, col.names=T)

#export the number of interaction pathways
for (k in c(1:5)){
  j = percent[k]
  cellchat.list = qread(paste0("./cellchat_list_thre", j, ".qs"))

  Npathway = NULL

  for (i in c(1:length(cellchat.list))){
    ht1 <- netAnalysis_signalingRole_heatmap1(cellchat.list[[i]], pattern = "all")
    ht1 = ht1[rowSums(ht1)>0,]
    Npathway = c(Npathway, nrow(ht1))
  }

  df <- as.data.frame(Npathway)
  rownames(df)=names(cellchat.list)
  colnames(df)=paste0("threshold_", j)
  if(k>1){
    interaction_res_df=cbind(interaction_res_df, df)
  } else {
    interaction_res_df = df
  }
}
#export result
write.table(interaction_res_df, file="./interaction_pathway_number.txt", sep="\t", quote=F, row.names=T, col.names=T)



#create outgoing/incoming signaling strength scatter plot at thresh.pc=0.15 (Fig. 5b)
######
cellchat.list = qread("./cellchat_list_thre0.15.qs", nthreads=12)

names(cellchat.list)=gsub("TASSeq", "TAS-Seq", names(cellchat.list))
names(cellchat.list)=gsub("10Xv2_", "10Xv2-", names(cellchat.list))


num.link <- sapply(cellchat.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(unlist(num.link)), round(max(unlist(num.link)), digits = -3))

num.link = NULL
for (i in c(1:9)){
  num.link1 = aggregateNet(cellchat.list[[i]], signaling = NULL, return.object = FALSE, remove.isolate = FALSE)$count
  num.link1 = rowSums(num.link1) + colSums(num.link1)-diag(num.link1)
  num.link = c(num.link, num.link1)
}
weight.MinMax <- c(min(num.link), max(num.link))

#set colors of cell subsets
hoge = NULL

for(i in c(1:9)){
  hoge = c(hoge, unique(as.character(cellchat.list[[i]]@meta$celltype)))
  hoge = unique(hoge)
}

hoge = sort(hoge)

#set custom color palette
custom_colors <- list()
colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51'
)
colors_spanish <- c(
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
)
colors_custom = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',
                  '#49beaa', '#611c35', '#2708a0')
custom_colors$discrete <- unique(c(colors_dutch, colors_spanish, colors_custom))


hoge1 = custom_colors$discrete[1:length(hoge)]
names(hoge1)=hoge

#extract legend
p = suppressMessages(netAnalysis_signalingRole_scatter(cellchat.list[[9]], color.use = hoge1,
                                                       title = names(cellchat.list[[9]]),
                                                       label.size = 4, font.size = 12,
                                                       font.size.title = 12,
                                                       weight.MinMax = weight.MinMax))
suppressWarnings(ggsave(file="Interaction_strength_Legend.png",
                        plot=p,
                        device="png", units="in",
                        dpi=300, width=5, height=5, limitsize = FALSE))

#plot interaction strength scatter plot
xbreaks = c(0, 0.3, 0.6, 0.9, 1.2)
ybreaks = c(0, 0.2, 0.4, 0.6, 0.8)

gg_list=list()
for(i in c(1:9)){
  gg_list[[i]] <- suppressMessages(netAnalysis_signalingRole_scatter(cellchat.list[[i]], color.use = hoge1,
                                                                     title = names(cellchat.list)[i],
                                                                     label.size = 4, font.size = 8,
                                                                     font.size.title = 10, show.legend = F,
                                                                     weight.MinMax = weight.MinMax)) +
    scale_x_continuous(breaks=xbreaks, labels=xbreaks, limits = c(0, 1.2)) +
    scale_y_continuous(breaks=ybreaks, labels=ybreaks, limits = c(0, 0.8)) +
    xlab("") + ylab("") + theme(axis.text=element_text(size=12),
                                plot.margin = unit(c(2,2,2,2), "mm"))
}

suppressWarnings(ggsave(file="Interaction_strength_thre0.15.png",
                        plot=plot_grid(plotlist=gg_list, ncol=3, nrow=3),
                        device="png", units="in",
                        dpi=300, width=8, height=8, limitsize = FALSE))

######


#create outgoing/incoming network circle plot of major contributing stromal cells
#at thresh.pc=0.15 (Figure 5c)

res=NULL
for (i in c(1:9)){
  tmp = max(cellchat.list[[i]]@net$weight)
  res = c(res, tmp)
}

cell_populations = c("FB_alveolar", "FB_adventitial", "Epi_AT2",
                     "Endo_capillary", "Endo_aerocyte", "Mo_Ly6Chi")

for (i in c(1:9)){
  groupSize <- as.numeric(table(cellchat.list[[i]]@idents))
  names(groupSize)=names(table(cellchat.list[[i]]@idents))
  groupSize = groupSize[cell_populations]
  hoge2 = cellchat.list[[i]]@net$weight
  hoge2 = hoge2[cell_populations, cell_populations]
  filename1=paste("interactions_weight_", names(cellchat.list)[i], ".png", sep="")
  png(filename = filename1, width = 2.8, height = 2.8, units = "in", res = 300, bg = "transparent")
  netVisual_circle(hoge2, vertex.weight = groupSize, top=1,
                   color.use = hoge1[cell_populations],
                   edge.weight.max = max(res),
                   vertex.label.cex = 0.5,
                   weight.scale = T, label.edge= F)
  dev.off()
}

#create outgoing/incoming network circle plot of ALL cells. Only show top 10% interactions (Supplementary Figure 11b)
for (i in c(1:9)){
  groupSize <- as.numeric(table(cellchat.list[[i]]@idents))
  names(groupSize)=names(table(cellchat.list[[i]]@idents))
  hoge2 = cellchat.list[[i]]@net$weight
  filename1=paste("interactions_weight_ALL", names(cellchat.list)[i], ".png", sep="")
  png(filename = filename1, width = 5.5, height = 5.5, units = "in", res = 300, bg = "transparent")
  netVisual_circle(hoge2, vertex.weight = groupSize, top=0.1,
                   color.use = hoge1[rownames(hoge2)],
                   edge.weight.max = max(res),
                   vertex.label.cex = 0.5,
                   weight.scale = T, label.edge= F)
  dev.off()
}
