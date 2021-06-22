#Cellchat analysis for make figure 2k and Supplementary Figure 6
library(CellChat)
library(patchwork)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE)
#please install rDBEC package before use (could be downloaded from this repository)
suppressWarnings(suppressMessages(source("library_source_Seurat.R")))
suppressMessages(library(stringr))
species="mmu"

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

#Please Place preprocessed Seurat v2.3.4 objects into working directory.

#load data
load("./Lung10X.rda")
load("Lung_day00_Rhapsody.rda")
load("Lung_SMARTSeq2.rda")
load("Lung10Xv3.rda")
seu_list=list()


#remove doublets and make seurat object list
mBC = Lung_BDRhapsody
doublets = mBC@meta.data$celltype %in% c("doublet", "misc")
mBC = SubsetData(object = mBC, cells.use=rownames(mBC@meta.data[!doublets,]),
                 subset.raw = TRUE)
seu_list[[1]] = mBC

mBC = Lung_SMARTSeq2
doublets = mBC@meta.data$celltype %in% c("doublet", "misc")
mBC = SubsetData(object = mBC, cells.use=rownames(mBC@meta.data[!doublets,]),
                 subset.raw = TRUE)
seu_list[[2]] = mBC

mBC = Lung10X
doublets = mBC@meta.data$celltype %in% c("doublet", "misc")
mBC = SubsetData(object = mBC, cells.use=rownames(mBC@meta.data[!doublets,]),
                 subset.raw = TRUE)
seu_list[[3]] = mBC

mBC = Lung10Xv3
doublets = mBC@meta.data$celltype %in% c("doublet", "misc")
mBC = SubsetData(object = mBC, cells.use=rownames(mBC@meta.data[!doublets,]),
                 subset.raw = TRUE)
seu_list[[4]] = mBC

idents = c("TAS-Seq","Smart-seq2", "10X_v2", "10X_v3")
names(seu_list) = idents

# set cell cluster colors
hoge = unique(c(as.character(unique(seu_list[[1]]@meta.data$celltype)),
                as.character(unique(seu_list[[2]]@meta.data$celltype)),
                as.character(unique(seu_list[[3]]@meta.data$celltype)),
                as.character(unique(seu_list[[4]]@meta.data$celltype))))
hoge1 = custom_colors$discrete[1:length(hoge)]
names(hoge1)=hoge

#make datasets for CellChat imput (downsample to the cell number pf Smart-seq2 dataset)

norm.data.input=list()
raw.data.input=list()
meta=list()
species = "mmu"

if(species=="mmu"){
  CellChatDB <- CellChatDB.mouse
} else if (species=="hsa"){
  CellChatDB <- CellChatDB.human
}
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

CellChatDB.use.symbols = unique(CellChatDB.use$geneInfo$Symbol)
Ncells = c(nrow(seu_list[[1]]@meta.data),
           nrow(seu_list[[2]]@meta.data),
           nrow(seu_list[[3]]@meta.data),
           nrow(seu_list[[4]]@meta.data)
)


for (i in c(1:length(seu_list))){
  seu = seu_list[[i]]

  tmp = as.character(seu@meta.data$orig.ident)
  names(tmp)=rownames(seu@meta.data)
  seu@ident = as.factor(tmp)
  seu = SubsetData(seu, ident.use = as.character(unique(seu@meta.data$orig.ident)),
                   max.cells.per.ident = min(Ncells), subset.raw = TRUE)
  tmp = as.character(seu@meta.data$celltype)
  names(tmp)=rownames(seu@meta.data)
  tmp = factor(tmp)
  seu@ident = tmp
  seu = NormalizeData(seu, scale.factor=1000000) #scale factor=1000000 for unifying global scaling factor between datasets
  norm.data.input[[i]] = seu@data[intersect(rownames(seu@data), CellChatDB.use.symbols),]
  raw.data.input[[i]] = seu@raw.data[intersect(rownames(seu@data), CellChatDB.use.symbols),]
  meta[[i]] = seu@meta.data
  # a dataframe with rownames containing cell mata data
  meta[[i]]$orig.ident = as.character(meta[[i]]$orig.ident)
  meta[[i]]$celltype = as.character(meta[[i]]$celltype)
}

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

#calculate thresh.pc = 0.05, 0.15, 0.25, 0.50, 0.75 for Gene filtering


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

for (j in c(0.05, 0.15, 0.25, 0.50, 0.75)){
  cellchat.list=list()
  for (i in c(1:length(seu_list))){
    cellchat = createCellChat(object = norm.data.input[[i]], meta = meta[[i]],
                              group.by = "celltype") #specify Celltype column name
    cellchat@data.raw = as.matrix(raw.data.input[[i]])
    groupSize = as.numeric(table(cellchat@idents)) # number of cells in each cell group
    cellchat@DB = CellChatDB.use
    cellchat = CellChat::subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
    future::plan("multisession", workers = 16) # do parallel
    cellchat = identifyOverExpressedGenes(cellchat, thresh.p=0.05, thresh.pc = j)
    cellchat = identifyOverExpressedInteractions(cellchat)
    if(species=="mmu"){
      cellchat <- projectData(cellchat, PPI.mouse)
    } else if(species=="hsa"){
      cellchat <- projectData(cellchat, PPI.human)
    }
    future::plan("multisession", workers = 4) # paralell = OFF
    cellchat <- computeCommunProb2(cellchat, raw.use = TRUE, type = "mean",
                                   population.size = TRUE, nthreads = 12)
    future::plan("multisession", workers = 16) # do parallel
    cellchat <- filterCommunication(cellchat, min.cells = 5)
    cellchat <- computeCommunProbPathway(cellchat, thresh = 0.05)
    cellchat <- aggregateNet(cellchat, thresh = 0.05)
    cellchat = netAnalysis_computeCentrality(cellchat, slot.name="netP")
    cellchat.list[[i]]=cellchat
  }


  idents = c("TAS-Seq","Smart-seq2", "10X_v2", "10X_v3")
  names(cellchat.list) = idents
  
  #extract significant signaling pathways of cell-cell interactions of each dataset
  mat=list()
  
  for(i in c(1:4)){
    ht1 <- netAnalysis_signalingRole_heatmap1(cellchat.list[[i]], pattern = "all")
    ht1 <- sweep(ht1, 1L, apply(ht1, 1, max), '/', check.margin = FALSE)
    mat[[i]] = as.data.frame(ht1)
  }
  
  #visualize overlapped outgoing signaling pathways
  TASseq = rownames(mat[[1]])
  SmartSeq2 = rownames(mat[[2]])
  TenXv2 = rownames(mat[[3]])
  TenXv3 = rownames(mat[[4]])
  
  
  fuga1 = as.data.frame(TASseq)
  rownames(fuga1)=fuga1[,1]
  fuga2 = as.data.frame(SmartSeq2)
  rownames(fuga2)=fuga2[,1]
  fuga3 = as.data.frame(TenXv2)
  rownames(fuga3)=fuga3[,1]
  fuga4 = as.data.frame(TenXv3)
  rownames(fuga4)=fuga4[,1]
  
  fuga5 = dplyr::full_join(rownames_to_column(fuga1), rownames_to_column(fuga2), by ="rowname")
  fuga5 = dplyr::full_join(fuga5, rownames_to_column(fuga3), by ="rowname")
  fuga5 = dplyr::full_join(fuga5, rownames_to_column(fuga4), by ="rowname")
  fuga5[is.na(fuga5)]=0
  
  file.name=sprintf("%s.txt", j)
  file.name=paste("pathways_", file.name, sep="")
  write.table(fuga5, file=file.name, sep="\t", quote=F)
  
  #create interaction number plot (for Figure 2k)
  ######
  group.new=NULL
  for (i in c(1:length(cellchat.list))){
    group.new = union(group.new, levels(cellchat.list[[i]]@idents))
  }
  cellchat.list1=list()
  for (i in c(1:length(cellchat.list))){
    cellchat.list1[[i]] = suppressMessages(suppressWarnings(liftCellChat(cellchat.list[[i]],
                                                                         group.new)))
  }
  hoge = cellchat.list1[[1]]@meta
  hoge = hoge[,c(3,9)]
  cellchat.list1[[1]]@meta = hoge
  hoge = cellchat.list1[[2]]@meta
  hoge = hoge[,c(3,8)]
  cellchat.list1[[2]]@meta = hoge
  hoge = cellchat.list1[[3]]@meta
  hoge = hoge[,c(3,7)]
  cellchat.list1[[3]]@meta = hoge
  hoge = cellchat.list1[[4]]@meta
  hoge = hoge[,c(3,8)]
  cellchat.list1[[4]]@meta = hoge
  
  cellchat = mergeCellChat(cellchat.list1, add.names=idents)
  gg1 <- compareInteractions(cellchat, show.legend = F)
  gg2 <- compareInteractions(cellchat, show.legend = F, measure = "weight")
  
  gg1 = as.ggplot(gg1)
  gg2 = as.ggplot(gg2)
  file.name=sprintf("%s.png", j)
  file.name=paste("compareInteractions_", file.name, sep="")
  ggsave(file= file.name, plot=plot_grid(gg1, gg2), device="png", units="in",
         dpi=300, width=4, height=3, limitsize = FALSE)
  ######
  
  #create outgoing/incoming signaling strength scatter plot at thresh.pc=0.15 (Supplementary Figure 7a)
  ######
  if(j==0.15){
    num.link <- sapply(cellchat.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
    weight.MinMax <- c(min(unlist(num.link)), round(max(unlist(num.link)), digits = -3))
    
    num.link = NULL
    for (i in c(1:4)){
      num.link1 = aggregateNet(cellchat.list[[i]], signaling = NULL, return.object = FALSE, remove.isolate = FALSE)$count
      num.link1 = rowSums(num.link1) + colSums(num.link1)-diag(num.link1)
      num.link = c(num.link, num.link1)
    }
    weight.MinMax <- c(min(num.link), max(num.link))
    
    gg_list=list()
    gg_list[[1]] <- suppressMessages(netAnalysis_signalingRole_scatter(cellchat.list[[1]], color.use = hoge1,
                                                                       title = names(cellchat.list[[1]]),
                                                                       label.size = 4, font.size = 12,
                                                                       font.size.title = 12,
                                                                       weight.MinMax = weight.MinMax)) +
      xlim(0, 1.5) + ylim(0,1) + theme(axis.text=element_text(size=12))
    gg_list[[2]] <- suppressMessages(netAnalysis_signalingRole_scatter(cellchat.list[[2]], color.use = hoge1,
                                                                       title = names(cellchat.list[[2]]),
                                                                       label.size = 4, font.size = 12,
                                                                       font.size.title = 12,
                                                                       weight.MinMax = weight.MinMax)) +
      xlim(0, 0.8) + ylim(0,0.4) + theme(axis.text=element_text(size=12))
    gg_list[[3]] <- suppressMessages(netAnalysis_signalingRole_scatter(cellchat.list[[3]], color.use = hoge1,
                                                                       title = names(cellchat.list[[3]]),
                                                                       label.size = 4, font.size = 12,
                                                                       font.size.title = 12,
                                                                       weight.MinMax = weight.MinMax)) +
      xlim(0, 0.6) + ylim(0,0.3) + theme(axis.text=element_text(size=12))
    gg_list[[4]] <- suppressMessages(netAnalysis_signalingRole_scatter(cellchat.list[[4]], color.use = hoge1,
                                                                       title = names(cellchat.list[[4]]),
                                                                       label.size = 4, font.size = 12,
                                                                       font.size.title = 12,
                                                                       weight.MinMax = weight.MinMax)) +
      xlim(0, 0.6) + ylim(0,0.3) + theme(axis.text=element_text(size=12))
    
    suppressWarnings(ggsave(file="Interaction_strength.png",
                            plot=plot_grid(plotlist=gg_list, ncol=2, nrow=2),
                            device="png", units="in",
                            dpi=300, width=5, height=5, limitsize = FALSE))
    
    
    gg_list=list()
    gg_list[[1]] <- suppressMessages(netAnalysis_signalingRole_scatter(cellchat.list[[1]], color.use = hoge1,
                                                                       title = names(cellchat.list[[1]]),
                                                                       label.size = 3.8, font.size = 12,
                                                                       font.size.title = 12,show.legend = FALSE,
                                                                       weight.MinMax = weight.MinMax)) +
      xlim(0, 1.5) + ylim(0,1) + ggtitle(idents[1]) +
      theme(axis.text=element_text(size=12), axis.title = element_blank(), plot.title=element_text(size=12))
    gg_list[[2]] <- suppressMessages(netAnalysis_signalingRole_scatter(cellchat.list[[2]], color.use = hoge1,
                                                                       title = names(cellchat.list[[2]]),
                                                                       label.size = 3.8, font.size = 12,
                                                                       font.size.title = 12,show.legend = FALSE,
                                                                       weight.MinMax = weight.MinMax)) +
      xlim(0, 0.8) + ylim(0,0.4) + ggtitle(idents[2]) +
      theme(axis.text=element_text(size=12), axis.title = element_blank(), plot.title=element_text(size=12))
    gg_list[[3]] <- suppressMessages(netAnalysis_signalingRole_scatter(cellchat.list[[3]], color.use = hoge1,
                                                                       title = names(cellchat.list[[3]]),
                                                                       label.size = 3.8, font.size = 12,
                                                                       font.size.title = 12,show.legend = FALSE,
                                                                       weight.MinMax = weight.MinMax)) +
      xlim(0, 0.6) + ylim(0,0.3)+ ggtitle(idents[3]) +
      theme(axis.text=element_text(size=12), axis.title = element_blank(), plot.title=element_text(size=12))
    gg_list[[4]] <- suppressMessages(netAnalysis_signalingRole_scatter(cellchat.list[[4]], color.use = hoge1,
                                                                       title = names(cellchat.list[[4]]),
                                                                       label.size = 3.8, font.size = 12,
                                                                       font.size.title = 12,show.legend = FALSE,
                                                                       weight.MinMax = weight.MinMax)) +
      xlim(0, 0.6) + ylim(0,0.3) + ggtitle(idents[4]) +
      theme(axis.text=element_text(size=12), axis.title = element_blank(), plot.title=element_text(size=12))
    
    suppressWarnings(ggsave(file="Interaction_strength1.png",
                            plot=plot_grid(plotlist=gg_list, ncol=2, nrow=2),
                            device="png", units="in",
                            dpi=300, width=5, height=5, limitsize = FALSE))
  }
  ######
  
  
  #create outgoing/incoming network circle plot of major contributing stromal cells
  #at thresh.pc=0.15 (Figure 2l)
  
  if(j==0.15){
  res=NULL
  for (i in c(1:4)){
    tmp = max(cellchat.list[[i]]@net$weight)
    res = c(res, tmp)
  }

  cell_populations = c("FB_Inmthi", "FB_Dcnhi", "Epi_AT2", 
                       "Endo_VEC", "Endo_capillary", "Mo_Ly6Chi")
  
  for (i in c(1:4)){
    groupSize <- as.numeric(table(cellchat.list[[i]]@idents))
    names(groupSize)=names(table(cellchat.list[[i]]@idents))
    groupSize = groupSize[cell_populations]
    hoge2 = cellchat.list[[i]]@net$weight
    hoge2 = hoge2[cell_populations, cell_populations]
    filename1=paste("interactions_weight_", names(cellchat.list)[i], ".png", sep="")
    png(filename = filename1, width = 2.8, height = 2.8, units = "in", res = 300)
    netVisual_circle(hoge2, vertex.weight = groupSize, top=1,
                     color.use = hoge1[cell_populations],
                     edge.weight.max = max(res),
                     vertex.label.cex = 0.5,
                     weight.scale = T, label.edge= F)
    dev.off()
  }
  }

  #save cellchat object
  file.name=sprintf("%s.rda", j)
  file.name=paste("cellchat_list_", file.name, sep="")
  save.pigz(cellchat.list, file = file.name, n.cores=16)
}
