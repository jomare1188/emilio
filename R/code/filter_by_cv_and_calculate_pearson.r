# load libraries

library(dplyr)
library(tidyr)
library(DESeq2)
library(WGCNA)
library(ggrepel)
library(wesanderson)
library(tibble)
library(HiClimR)
library(reshape2)
library(igraph)
# set up things
HOME_DIR = "/Storage/data1/jorge.munoz/emilio/"
setwd(HOME_DIR)
# read vst matrix
datExpr <- as.matrix(read.table("../data/vst_all.csv", sep = ",", header = T, row.names = 1))

# calculate cv coefficient to filter
row_means <- rowMeans(datExpr)
row_sds <- rowSds(datExpr)
cv <- (row_sds / row_means) * 100

# Set a threshold for the minimum coefficient of variation you want to keep
min_cv_threshold <- 18

# Filter rows (genes) based on the coefficient of variation threshold
filtered_genes <- datExpr[cv >= min_cv_threshold, ]
dim(filtered_genes)

# make histrogram of vst
#hist <- hist(filtered_vst, xlim=c(2,20), breaks = 18)
# calculate pearson and format correlation table in triples
pcor <- melt(fastCor(t(filtered_genes), nSplit = 40, upperTri = TRUE, optBLAS = TRUE, verbose = TRUE))
dim(pcor)

# filter strong correlations only and erase diagonal 
person_fil <- pcor[abs(pcor$value) > 0.90,] %>% filter(Var1 != Var2)
dim(person_fil)

# you need positive  values to run mcl
person_fil$value <- abs(person_fil$value)

# we import to igraph to simplify i.e. remove loops
colnames(person_fil) <- c("Var1", "Var2", "weight")
graph_raw <- graph_from_data_frame(person_fil, directed = FALSE)
graph_simple <- simplify(graph_raw, remove.multiple = T)
edgelist_simple <-as.data.frame(get.edgelist(graph_simple))
## we add the weight as an atrribute to the edge liste (triples)
edgelist_simple$weight <- abs(E(graph_simple)$weight)
## write triples
write.table(edgelist_simple, file = "../results/full_p90_cv18.triples", row.names = F, col.names = F , quote = F, sep = " ")

