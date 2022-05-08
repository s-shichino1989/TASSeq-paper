#source file of DBEC preprocessing
#background subtration of BD Rhapsody WTA data based on bimodal distribution of lognormal raw read counts
#of each gene by using mclust package
#Written by Shigeyuki Shichino, ver0.2 20190603


reticulate::use_python(
  python = "/usr/bin/python3",
  required = TRUE
)
library(reticulate)
py_set_seed(42, disable_hash_randomization = TRUE)
suppressMessages(library(BiocParallel))
suppressMessages(library(rDBEC))
suppressMessages(library(tradeSeq))
suppressMessages(library(slingshot))
suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(Matrix))
suppressMessages(library(mclust))
suppressMessages(library(Seurat))
suppressMessages(library(recommenderlab))
suppressMessages(library(MASS))
suppressMessages(library(future.apply))
suppressMessages(library(future))
suppressMessages(library(doFuture))
suppressMessages(library(fastSave))
suppressMessages(library(stats))
suppressMessages(library(factoextra))
suppressMessages(library(grDevices))
suppressMessages(library(ggplotify))
suppressMessages(library(ggpubr))
suppressMessages(library(cowplot))
suppressMessages(library(pheatmap))
suppressMessages(library(knitr))
suppressMessages(library(rmdformats))
suppressMessages(library(stringr))
suppressMessages(library(DT))
suppressMessages(library(cluster))
suppressMessages(library(parallelDist))
suppressMessages(library(WGCNA))
suppressMessages(library(flashClust))
suppressMessages(library(SingleR))
suppressMessages(library(ggsci))
suppressMessages(library(qs))
set.seed(42)
options(future.globals.maxSize= 100000*1024^3)


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
colors_temp = pal_igv(palette = "default")
colors_temp = colors_temp(51)
custom_colors$discrete <- unique(c(colors_dutch, colors_spanish, colors_custom, colors_temp))


quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

tableread_fast_hashtag = function(x, sep="\t", header=TRUE){
  tmp = data.table::fread(x, header=header, sep=sep, quote="")
  tmp = as.data.frame(tmp)
  return(tmp)
}
