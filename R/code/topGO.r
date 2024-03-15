load_cluster_libraries <-function(){
library(rlang, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(graph, lib.loc = "/Storage/data1/jorge.munoz/TOPGO/data")
library(GO.db, lib.loc = "/Storage/data1/jorge.munoz/TOPGO/data")
library(SparseM, lib.loc = "/Storage/data1/jorge.munoz/TOPGO/data")
library(topGO, lib.loc = "/Storage/data1/jorge.munoz/TOPGO/data")
library(dplyr, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(Rgraphviz, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(ggplot2, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(scales, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
}
results_path_cluster <-"/Storage/data1/jorge.munoz/ANGELICA/redes_angelica/enrichment"
dir.create(results_path_cluster)

load_cluster_files <- function(){
## read clusters size
number.clusters <<- read.table("/Storage/data1/jorge.munoz/ANGELICA/redes_angelica/formated_mclI2.number", header = F)
## read list of genes for all the network 
geneID2GO <<- readMappings(file = "/Storage/data1/jorge.munoz/ANGELICA/redes_angelica/GOterms_ok.tsv")
# read genes in graph 
genes_in_graph <<- as.data.frame(rownames(read.table("/Storage/data1/jorge.munoz/ANGELICA/redes_angelica/mcl_I2.formated.csv", header =F, row.names=1)))
#genes_in_graph <<- as.data.frame(rownames(read.table("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/modules_specc_diff_expr.txt")))
colnames(genes_in_graph) <<- "gene"
# filter background genes to genes in the network
geneID2GO_filtered <<- geneID2GO[genes_in_graph[,1]]
#geneNames <<- names(geneID2GO_filtered)
geneNames <<- genes_in_graph$gene
# read modules
#dynamicMods <<- read.table(file = "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/modules_specc_diff_expr.txt")
dynamicMods <<- read.table(file = "/Storage/data1/jorge.munoz/ANGELICA/redes_angelica/mcl_I2.formated.csv", header = F, row.names=1)
colnames(dynamicMods) <<- "module_No"
}

# load things in cluster
load_cluster_libraries()
load_cluster_files()

# make annot fun
anot_modules <- function(Module_no, results_path){
paste0("MODULO:", Module_no,"      ","NODOS:",number.clusters[Module_no,1])
myInterestingGenes <- genes_in_graph[dynamicMods$module_No==Module_no,]
geneList <<- as.factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
GOdata <<- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO_filtered)
allGO=usedGO(GOdata)
Classic <<- runTest(GOdata, algorithm = "classic", statistic = "fisher")
# Make results  table
#table <- GenTable(GOdata, Classic = resultClassic, Weight01 = resultsWeight01, topNodes = length(allGO), orderBy = 'Classic')
table <- GenTable(GOdata, Classic = Classic, topNodes = length(allGO), orderBy = 'Classic')
# Filter not significant values for classic algorithm
table1 <- filter(table, Classic < 0.05 )
# Performing BH correction on our p values FDR
p.adj <- round(p.adjust(table1$Classic,method="BH"),digits = 4)
# Create the file with all the statistics from GO analysis
all_res_final <<- cbind(table1,p.adj)
all_res_final <<- all_res_final[order(all_res_final$p.adj),]
# Get list of significant GO before multiple testing correction
results.table.p = all_res_final[which(all_res_final$Classic <=0.05),]
# Get list of significant GO after multiple testing correction
results.table.bh = all_res_final[which(all_res_final$p.adj<=0.05),]
# Save first top 50 ontolgies sorted by adjusted pvalues
write.table(all_res_final[1:50,], file = paste0(results_path, "module_", Module_no, ".csv"), quote=FALSE, row.names=FALSE, sep = ",")
}

for (i in 1:dim(table(dynamicMods))){ 
	if(number.clusters[i,1]>=4){
		print(i)
		anot_modules(i,results_path_cluster)
		
}
}


