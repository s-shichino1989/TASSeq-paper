##draw Figure 6

library(Seurat)
library(ggplot2)
library(cowplot)
library(org.Mm.eg.db) 
library(AnnotationDbi)
library(GO.db)
library(qs)

set.seed(42)
dir.name="./"

DotPlot2 = function (object, assay = NULL, features, cols = c("lightgrey",
                                                              "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,

                     idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE,
                     scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA)
{


  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in%
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = scale_size,
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features),
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE,
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  data.features <- FetchData(object = object, vars = features,
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits,
                              sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)),
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident,
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove,
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot,
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  if (length(x = levels(x = data.plot$id)) == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled",
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot ==
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min,
                                                  max = col.max)
                             }
                             else {
                               data.use <- log1p(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled,
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot,
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(X = as.character(x = data.plot$id),
                         FUN = gsub, FUN.VALUE = character(length = 1L),
                         pattern = paste0("^((", paste(sort(x = levels(x = object),
                                                            decreasing = TRUE), collapse = "|"), ")_)"),
                         replacement = "", USE.NAMES = FALSE)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors",
                     no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot],
                                       levels = unique(x = feature.groups))
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot",
                                                        y = "id")) +
    geom_point(mapping = aes_string(size = "pct.exp",
                                    color = color.by)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    scale_size(breaks=c(5, 25, 50, 75, 100), limits = c(5, 100)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = "Percent\nexpressed")) +
    labs(x = "Features", y = ifelse(test = is.null(x = split.by),
                                    yes = "Identity", no = "Split Identity")) + theme_cowplot()
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(facets = ~feature.groups,
                              scales = "free_x", space = "free_x", switch = "y") +
      theme(panel.spacing = unit(x = 1, units = "lines"),
            strip.background = element_blank())
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  }
  else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = "Average\nExpression"))
  }
  return(plot)
}

################

##download annotated Seurat objects for mouse lung data
system("wget -L -O UTlung_deep_Seurat_annot.qs https://tus.box.com/shared/static/ohoohhgf3t12tcmupueajbbapu4xdxm6.qs")
system("wget -L -O UTlung_SmartSeq2_Seurat_annot.qs https://tus.box.com/shared/static/ej8tff4k7szdfqsw6niljnzrronnx910.qs")
system("wget -L -O UTlung_10Xv2_TabulaMuris_Seurat_annot.qs https://tus.box.com/shared/static/fqiabrfwxd7y15jgcwlovelhoz2tgz65.qs")


fnames = dir(pattern="UTlung_")
seu_list = lapply(fnames, qread, nthreads=24)


fnames = dir(pattern="UTlung_")
seu_list = lapply(fnames, qread, nthreads=24)

#rename Smart-seq2 orig.ident
seu_list[[3]]@meta.data$orig.ident = rep("Smart-seq2", nrow(seu_list[[5]]@meta.data))

#rename 10X v2 Tabula Muris data
seu_list[[1]]@meta.data$orig.ident = rep("10Xv2_TabulaMuris", nrow(seu_list[[2]]@meta.data))


#remove doublets
for(i in c(1:3)){
  doublet_key = seu_list[[i]]@meta.data$celltype %in% c("doublet", "misc")
  seu_list[[i]]=subset(seu_list[[i]], cells = rownames(seu_list[[i]]@meta.data[!doublet_key,]))
}

#merge data
for(i in c(1:3)){
  if(i>1){
    seu = merge(x=seu, y=seu_list[[i]])
  } else {
    seu = seu_list[[1]]
  }
}

idents = c("10Xv2-TabulaMuris","TAS-Seq.deep", "Smart-seq2")

#plot interleukin genes
for(i in c(1:3)){
interleukins = grep(pattern="Il[1-9].*", x = rownames(seu_list[[i]]@assays$RNA@counts), value=TRUE)
interleukins = interleukins[!interleukins %in% c(grep(pattern="r", x = interleukins, value=TRUE), "IL18bp", "Il1bos", "Il4i1", "Il6st")]
interleukins=sort(interleukins)

p = DotPlot2(seu_list[[i]], cols = c("grey", "red"), features=interleukins, scale = FALSE, dot.min = 0.05,
            group.by="celltype") +
   theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), ) +
   xlab("") + ylab("") +
   ggtitle(idents[i])+
       theme(plot.title=element_text(hjust = 0.5),
          text=element_text(size=10),
          axis.text.x = element_text(size=10, family="italic"))

ggsave(filename=paste0(idents[i], "_IL_dotplot.png"), plot = p, device = "png", width = 7, height=8, units = "in", dpi = 300,
       limitsize = FALSE, bg = "white")
}


#extract growth factor genes
data2 = rownames(seu@assays$RNA@counts)
res = select(org.Mm.eg.db, keys = data2, keytype = "SYMBOL",columns = c("SYMBOL", "GOALL"))
data3 = res[res$ONTOLOGYALL %in% "MF",]
key1 = as.character(as.vector(data3[,2]))
key2 = Term(key1)
data4 = transform(data3,Target=key2)

#specify GO term, you can search by using QuickGO
annotation = c("growth factor activity")
data5 = data4[as.character(data4[,5]) %in% annotation,]
annotation_res = unique(data5[,1])

#remove interleukin genes
annotation_res = annotation_res[!annotation_res %in% interleukins]
annotation_res = sort(annotation_res)

#plot GF genes
for(i in c(1:3)){
  GF_genes = annotation_res[annotation_res %in% rownames(seu_list[[i]]@assays$RNA@counts)]

  p = DotPlot2(seu_list[[i]], cols = c("grey", "red"), features=GF_genes, scale = FALSE, dot.min = 0.05,
              group.by="celltype") +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), ) +
    xlab("") + ylab("") +
    ggtitle(idents[i])+
        theme(plot.title=element_text(hjust = 0.5),
          text=element_text(size=10),
          axis.text.x = element_text(size=10, family="italic"))

  ggsave(filename=paste0(idents[i], "_GF_dotplot.png"), plot = p, device = "png", width = 9, height=8, units = "in", dpi = 300,
         limitsize = FALSE, bg = "white")
}

