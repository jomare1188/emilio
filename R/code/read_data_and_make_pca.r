# load libraries
library(dplyr)
library(tidyr)
library(DESeq2)
library(WGCNA)
library(ggrepel)
library(wesanderson)
library(tibble)

# read raw matrix
HOME_DIR = "/Storage/data1/jorge.munoz/emilio"
setwd(HOME_DIR)
raw <- read.table("../data/GSE153345_TIS_counts.txt", head = T, row.names = 1)
metadata <- read.table("../data/run_selector_moni.txt", sep = ",", head = T)
row_data <- read.table("../data/metadata.txt")
colnames(raw) <- row_data$V2

# make pca with ALL the data
# We make r-base PCA it worth to filter by cv before making PCA to only left the most variat genes or simply using the plotPCA function of DESEQ2
raw.good <- raw[goodGenes(t(raw)),]
raw.vst <- vst(as.matrix(raw.good))
raw.pca <- prcomp(x = t(raw.vst), scale = TRUE)
var_explained <- raw.pca$sdev^2/sum(raw.pca$sdev^2)

# set fancy colors
colors = wesanderson::wes_palettes$GrandBudapest1

# from metadata we break the columns to assing groups (e.g. wild-type_mock_T12h_rep1 turns wild-type, T12h_rep1, rep1)
# this way we can indentify the points in the PCA 

pca_data <- raw.pca$x %>% 
  as.data.frame %>%
  rownames_to_column("group") %>%
  separate(group,c("type", "status", "time"), sep = "_") %>%
  unite(reps, c("type", "status", "time"), remove = F)

p <-  ggplot(data=, pca_data, aes(x=PC1,y=PC2, label = time, fill = type )) + 
  geom_text(check_overlap = TRUE, vjust = 1.5) +
  geom_point(aes(color= type, shape = status ), size = 4 ) + 
  scale_color_manual(values = colors) +
  theme_bw(base_size=20) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top",
         text = element_text(family = "Times", size=20), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text.x=element_text(colour="black", angle = 8, vjust = 0.7, hjust=0.5),
         axis.text.y=element_text(colour="black"),
         axis.line = element_line(colour = "black"))

ggsave(p,filename = "../results/PCA_raw.png" ,units = "cm", width = 25, height = 25,dpi = 320)
# now we agreggate replicates (this for sure can be done in a better way)
# make a vector for each condition anb mean values (replicates),
wild.type_mock_T12h <- raw %>% select(1:5) %>% rowMeans() %>% as.data.frame() 
colnames(wild.type_mock_T12h) <- "wild.type_mock_T12h"

wild.type_infected_T12h <- raw %>% select(6:10) %>% rowMeans() %>% as.data.frame() 
colnames(wild.type_infected_T12h) <- "wild.type_infected_T12h"

wild.type_mock_T24h <- raw %>% select(11:15) %>% rowMeans() %>% as.data.frame() 
colnames(wild.type_mock_T24h) <- "wild.type_mock_T24h"

wild.type_infected_T24h <- raw %>% select(16:20) %>% rowMeans() %>% as.data.frame() 
colnames(wild.type_infected_T24h) <- "wild.type_infected_T24h"

wild.type_mock_T48h <- raw %>% select(21:25) %>% rowMeans() %>% as.data.frame() 
colnames(wild.type_mock_T48h) <- "wild.type_mock_T48h"

wild.type_infected_T48h <- raw %>% select(26:30) %>% rowMeans() %>% as.data.frame() 
colnames(wild.type_infected_T48h) <- "wild.type_infected_T48h"

wild.type_mock_T5d <- raw %>% select(31:35) %>% rowMeans() %>% as.data.frame() 
colnames(wild.type_mock_T5d) <- "wild.type_mock_T5d"

wild.type_infected_T5d <- raw %>% select(36:40) %>% rowMeans() %>% as.data.frame() 
colnames(wild.type_infected_T5d) <- "wild.type_infected_T5d"

wild.type_mock_T10d <- raw %>% select(41:45) %>% rowMeans() %>% as.data.frame() 
colnames(wild.type_mock_T10d) <- "wild.type_mock_T10d"

wild.type_infected_T10d <- raw %>% select(46:50) %>% rowMeans() %>% as.data.frame() 
colnames(wild.type_infected_T10d) <- "wild.type_infected_T10d"

wild.type_mock_T20d <- raw %>% select(51:55) %>% rowMeans() %>% as.data.frame() 
colnames(wild.type_mock_T20d) <- "wild.type_mock_T20d"

wild.type_infected_T20d <- raw %>% select(56:60) %>% rowMeans() %>% as.data.frame() 
colnames(wild.type_infected_T20d) <- "wild.type_infected_T20d"

wild.type_mock_T30d <- raw %>% select(61:65) %>% rowMeans() %>% as.data.frame() 
colnames(wild.type_mock_T30d) <- "wild.type_mock_T30d"

wild.type_infected_T30d <- raw %>% select(66:70) %>% rowMeans() %>% as.data.frame() 
colnames(wild.type_infected_T30d) <- "wild.type_infected_T30d"

AtCKX2_mock_T12h <- raw %>% select(71:75) %>% rowMeans() %>% as.data.frame() 
colnames(AtCKX2_mock_T12h) <- "AtCKX2_mock_T12h"

AtCKX2_infected_T12h <- raw %>% select(76:80) %>% rowMeans() %>% as.data.frame() 
colnames(AtCKX2_infected_T12h) <- "AtCKX2_infected_T12h"

AtCKX2_mock_T24h <- raw %>% select(81:85) %>% rowMeans() %>% as.data.frame() 
colnames(AtCKX2_mock_T24h) <- "AtCKX2_mock_T24h"

AtCKX2_infected_T24h <- raw %>% select(86:90) %>% rowMeans() %>% as.data.frame() 
colnames(AtCKX2_infected_T24h) <- "AtCKX2_infected_T24h"

AtCKX2_mock_T48h <- raw %>% select(91:95) %>% rowMeans() %>% as.data.frame() 
colnames(AtCKX2_mock_T48h) <- "AtCKX2_mock_T48h"

AtCKX2_infected_T48h <- raw %>% select(96:100) %>% rowMeans() %>% as.data.frame() 
colnames(AtCKX2_infected_T48h) <- "AtCKX2_infected_T48h"

AtCKX2_mock_T5d <- raw %>% select(101:105) %>% rowMeans() %>% as.data.frame() 
colnames(AtCKX2_mock_T5d) <- "AtCKX2_mock_T5d"

AtCKX2_infected_T5d <- raw %>% select(106:110) %>% rowMeans() %>% as.data.frame() 
colnames(AtCKX2_infected_T5d) <- "AtCKX2_infected_T5d"

AtCKX2_mock_T10d <- raw %>% select(111:115) %>% rowMeans() %>% as.data.frame() 
colnames(AtCKX2_mock_T10d) <- "AtCKX2_mock_T10d"

AtCKX2_infected_T10d <- raw %>% select(116:120) %>% rowMeans() %>% as.data.frame() 
colnames(AtCKX2_infected_T10d) <- "AtCKX2_infected_T10d"
 
AtCKX2_mock_T30d <- raw %>% select(121:125) %>% rowMeans() %>% as.data.frame() 
colnames(AtCKX2_mock_T30d) <- "AtCKX2_mock_T30d"

AtCKX2_infected_T30d <- raw %>% select(126:130) %>% rowMeans() %>% as.data.frame() 
colnames(AtCKX2_infected_T30d) <- "AtCKX2_infected_T30d"

# AtCKX2_infected
# 12h, 24h, 48h, 5d, 10d, 30d
AtCKX2_infected <- cbind(AtCKX2_infected_T12h, AtCKX2_infected_T24h, AtCKX2_infected_T48h, AtCKX2_infected_T5d, AtCKX2_infected_T10d, AtCKX2_infected_T30d)

# AtCKX2_mock
# 12h, 24h, 48h, 5d, 10d, 30d
AtCKX2_mock <- cbind(AtCKX2_mock_T12h, AtCKX2_mock_T24h, AtCKX2_mock_T48h, AtCKX2_mock_T5d, AtCKX2_mock_T10d, AtCKX2_mock_T30d)

# wild_type_infected
# 12h, 24h, 48h, 5d, 10d, 20d, 30d 
wild.type_infected <- cbind(wild.type_infected_T12h, wild.type_infected_T24h, wild.type_infected_T48h, wild.type_infected_T5d, wild.type_infected_T10d, wild.type_infected_T20d, wild.type_infected_T30d)

# wild_type_mock
# 12h, 24h, 48h, 5d, 10d, 20d, 30d
wild.type_mock <- cbind(wild.type_mock_T12h, wild.type_mock_T24h, wild.type_mock_T48h, wild.type_mock_T5d, wild.type_mock_T10d, wild.type_mock_T20d, wild.type_mock_T30d)
###################################################################################3

## and we make a big matrix with all conditions 
aggregated <- cbind(AtCKX2_infected, AtCKX2_mock, wild.type_infected, wild.type_mock)

# PCA with aggregated replicates
aggregated.good <- round(aggregated[goodGenes(t(aggregated)),])
aggregated.vst <- vst(as.matrix(aggregated.good))
write.table(aggregated.vst, "vst_aggregated.csv", col.names =T , sep=",", quote = F)


aggregated.pca <- prcomp(x = t(aggregated.vst), scale = F)
# no necessary to scale???????
#aggregated.pca <- prcomp(x = vst, scale = T)

var_explained <- aggregated.pca$sdev^2/sum(aggregated.pca$sdev^2)

colors = wesanderson::wes_palettes$GrandBudapest1

pca_data2 <- aggregated.pca$x %>% 
  as.data.frame %>%
  rownames_to_column("group") %>%
  separate(group,c("type", "status", "time"), sep = "_") %>%
  unite(reps, c("type", "status", "time"), remove = F)

q <- ggplot(data = pca_data2, aes(x=PC1,y=PC2, label = time, fill = type )) + 
  geom_text(check_overlap = TRUE, vjust = 1.5) +
  geom_point(aes(color= type, shape = status ), size = 4 ) + 
  scale_color_manual(values = colors) +
  theme_bw(base_size=30) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top",
        text = element_text(family = "Times New Roman", size=20), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(colour="black", angle = 8, vjust = 0.7, hjust=0.5),
        axis.text.y=element_text(colour="black"),
        axis.line = element_line(colour = "black"))

ggsave(q,filename = "..results/PCA_agreggated.png" ,units = "cm", width = 25, height = 25,dpi = 320)

# We write the matrix 
write.table(file = "../results/vst_all.csv", x = aggregated.vst, sep = ",", row.names = T, col.names = T, quote = F)
write.table(file = "../results/counts_all.csv", x = aggregated.good, sep = ",", row.names = T, col.names = T, quote = F)
