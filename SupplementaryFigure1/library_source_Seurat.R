#source file of DBEC preprocessing
#background subtration of BD Rhapsody WTA data based on bimodal distribution of lognormal raw read counts
#of each gene by using mclust package
#Written by Shigeyuki Shichino, ver0.2 20190603

library(reticulate)
reticulate::use_python(
  python = "/usr/local/Python3-3.7.2/bin/python3.7",
  required = TRUE
)
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
suppressMessages(library(tidyr))
suppressMessages(library(RColorBrewer))
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
suppressMessages(library(ggsci))

set.seed(42)
options(future.globals.maxSize= 10000*1024^3)

