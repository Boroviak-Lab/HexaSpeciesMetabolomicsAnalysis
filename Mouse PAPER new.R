library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)
library("pheatmap")
library(data.table)
library(tidyverse)
library(readr)
library(stringr)

#######################################################################################
 
deng_counts <- read.table("PAPER/Datasets/Mouse/newDeng/DengData.csv",
                      sep=",", header=T, row.names=1)

deng_key <- read.table("PAPER/Datasets/Mouse/Deng/DengIDGG.csv", 
                       header=T, sep=",", row.names=1)

isinvit <- deng_key$Isin

labs <- deng_key$Type
labs <- labs[which(isinvit>0)]

deng_data <-  CreateSeuratObject(counts = deng_counts[,which(isinvit>0)], assay = "RNA",min.cells = 0, min.features = 0)
deng_data$CT <- labs
deng_data$dataset <- "Deng"

#Change the annotation of each cell to e.g. its type Amnnion etc
Idents(deng_data) <- deng_data@meta.data$CT

deng_data <- subset(deng_data, idents = c("Zy", "4C", "8C"))

#Normalise the data
deng_data <- NormalizeData(deng_data, verbose = FALSE)

#Get variable genes (reduce number)
deng_data<- FindVariableFeatures(deng_data, selection.method = "vst", nfeatures = 30000)

#scaling 
deng_data <- ScaleData(deng_data, verbose = FALSE)

deng_data <- RunPCA(deng_data, npcs = 30, verbose = FALSE) #,features = rownames(Markers))

DimPlot(deng_data, dims=c(1,2), reduction="pca", label= TRUE, repel=TRUE)

#deng_data <- FindNeighbors(deng_data, dims = 1:20, k.param = 10)
#deng_data <- FindClusters(deng_data, resolution=2)

#####################################################################################
# Mouse 3D 

mouse_counts2<-read.table("PAPER/Datasets/Mouse/3D/featurecountsAllMouse3D_all.csv", sep = ",", header = TRUE,row.names = 1)

mouse_key2 <-read.table("PAPER/Datasets/Mouse/3D/mouse3D_ALLKEY.csv", sep = ",", header = TRUE, row.names = 1)

isinvit <- mouse_key2$Isin
labs <-mouse_key2$Ident
labs <- labs[which(isinvit>0)]

mouse_data2 <-  CreateSeuratObject(counts = mouse_counts2[,which(isinvit>0)], assay = "RNA",min.cells = 0, min.features = 0)
mouse_data2$CT <- labs
mouse_data2$dataset <- "3D"

#Change the annotation of each cell to e.g. its type Amnnion etc
Idents(mouse_data2) <- mouse_data2@meta.data$CT

#Normalise the data
mouse_data2 <- NormalizeData(mouse_data2, verbose = FALSE)
#Get variable genes (reduce number)
mouse_data2<- FindVariableFeatures(mouse_data2, selection.method = "vst", nfeatures = 30000)

gene_list_mouse2 <- mouse_data2@assays[["RNA"]]@var.features

mouse_data2 <- FindVariableFeatures(mouse_data2, selection.method = "vst", nfeatures = 20000)

#scaling 
mouse_data2 <- ScaleData(mouse_data2, verbose = FALSE)

mouse_data2 <- RunPCA(mouse_data2, npcs = 30, verbose = FALSE) #,features = rownames(Markers))

DimPlot(mouse_data2, dims=c(1,2), reduction="pca", label= TRUE, repel=TRUE)
ggsave("PAPER/PCAs/3D_PCA_mouse.pdf", width=6, height=4)


oxphos <- read.csv("GeneOnt/oxphos.csv")
oxphos <-  pull(oxphos,Column2)
oxphos <- str_to_title(oxphos)
oxphos <-unique(oxphos)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_mouse2) %>% 
  unique()

glycolysis <- read.csv("GeneOnt/glycolysis.csv")
glycolysis <-  pull(glycolysis,Column2)
glycolysis <- str_to_title(glycolysis)
glycolysis <-unique(glycolysis)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_mouse2)%>% 
  unique()

##### Scores 

mouse_data2 <- AddModuleScore(mouse_data2, features = glycolysis, 
                              name = "glycolysis")
mouse_data2 <- AddModuleScore(mouse_data2, features = oxphos, 
                              name = "oxphos")

scores_mouse2 <- data.frame(mouse_data2@meta.data$oxphos1,mouse_data2@meta.data$glycolysis1, 
                            mouse_data2@active.ident)
scores_mouse2 <- rownames_to_column(scores_mouse2,var="cell")

scores_mouse2 <- reshape2::melt(scores_mouse2, measure.vars=c(2:3), variable.name="score")

colnames(scores_mouse2) <- c("cell","CT", "score", "value")

write.csv(scores_mouse2, "PAPER/Scores/mouse_3D_scores.csv")

ggplot(scores_mouse2, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.5, aes(fill=score))+
  scale_fill_manual(values = c("#6280e3", "#b83939"),
                    name="", labels = c("OXPHOS", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("Mouse 3D")+
  ylim(-0.5,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()
ggsave("PAPER/Scores/Mouse_3D_scores.pdf", width=6, height=4)


#####################################################################################
#################################################################################
# Chen et al
#### Load data 
#D1<-read.table("Mouse files/",
#               sep="\t",header = T, row.names=1)

library(purrr)
library(stringr)

D1 <- read.table ("Mouse files/a/GSE74155.txt") 

library(readr)

#### change gene ids
names <- rownames(D1) 
names <- as.data.frame(names)
gene_id <- read_tsv("Mouse files/mart_export-27.txt")

colnames(gene_id) <- c("stable_id","gene")

gene_id <- 
  filter(gene_id, gene_id$stable_id %in% names$names)

D1$stable_id <- rownames(D1)

D1 <- merge(D1, gene_id, by="stable_id")

cols <- colnames(D1)
View(cols)

select(D1, stable_id, gene)
D1 <- select(D1, -stable_id)

D1$gene <- make.unique(D1$gene)
rownames(D1) <- D1$gene
D1 <- select(D1, -gene)

### Seurat
mouse_data <- CreateSeuratObject(counts = D1,
                                 assay = "RNA",min.cells = 0, min.features = 0)

cluster2 <- WhichCells(object = mouse_data, ident = str_subset(mouse_data@active.ident, pattern="M"))

mouse_data <- SetIdent(
  object = mouse_data,
  cells = cluster2,
  value = 'N')


data <- rownames(mouse_data@meta.data)

cluster3 <- WhichCells(object = mouse_data, 
                       cells = str_subset(data, pattern="ES_2i"))

mouse_data <- SetIdent(
  object = mouse_data,
  cells = cluster3,
  value = 'ES2i')


###Idents(mouse_data) <- c("Epi", "ES")
#mouse_data <- subset(mouse_data, idents = c("Epi", "ES"), invert = FALSE)
#
mouse_data <- subset(mouse_data, subset = nFeature_RNA > 0)

#log normalise
mouse_data <- NormalizeData(mouse_data, verbose = FALSE)
#Get variable genes (reduce number)
mouse_data <- FindVariableFeatures(mouse_data, selection.method = "vst", nfeatures = 20000)



#Rename to keep variable convention for the join species modelling
mouse_data <- ScaleData(mouse_data, verbose = FALSE)

#Dimesionality reduction and clustering
mouse_data <- RunPCA(mouse_data, npcs = 30, verbose = FALSE) #,features = rownames(Markers))

DimPlot(mouse_data, dims=c(1,2), reduction="pca", label= TRUE, repel=TRUE)
ggsave("Mouse files/PCA_mouse.pdf", width=6 , height=4)

#####Subclustering 
mouse_data_Epi <- subset(mouse_data, idents = "Epi")
mouse_data_Epi <- RunPCA(mouse_data_Epi,npcs = 30, verbose = FALSE)
mouse_data_Epi <- FindNeighbors(mouse_data_Epi, reduction = "pca", dims = 1:20, k.param = 10) #k.param reduce for later applications? 
mouse_data_Epi <- FindClusters(mouse_data_Epi, resolution = 0.01)
DimPlot (mouse_data_Epi, dims=c(1,2), reduction= "pca")

#markers <- FindMarkers(mouse_data_Epi, ident.1=0, ident.2=1, min.pct = 0.25, logfc.threshold = 0.25)

cells <- WhichCells(mouse_data_Epi, idents = 0)

mouse_data <- SetIdent(
  object = mouse_data,
  cells = cells ,
  value = 'Epi_delay')

mouse_data <- subset(mouse_data, idents=c("ES2i","Epi", "ES", "N"))

mouse_data<- RunPCA(mouse_data, npcs=30, verbose=FALSE)
DimPlot(mouse_data, dims=c(1,2), label=TRUE, repel=TRUE, cols = col)
ggsave(filename="PAPER/PCAs/chen_dimplot.pdf", width=6, height=4)

pct <- mouse_data@reductions[["pca"]]@stdev /sum(mouse_data@reductions[["pca"]]@stdev) *100


levels(mouse_data) <- c ("ES2i", "ES", "Epi", "N")
col <- c("#f9cbcc", "#f4a57d", "#e41f85", "#c75fb9")

meta <- mouse_data@meta.data
mouse_data <- FindVariableFeatures(mouse_data, selection.method = "vst", nfeatures = 30000)
mouse_data <- ScaleData(mouse_data, verbose = FALSE)

###### Mammalian scores 
mouse_data <- FindVariableFeatures(mouse_data, selection.method = "vst", nfeatures = 40000)
gene_list <- mouse_data@assays[["RNA"]]@var.features

oxphos <- read.csv("GeneOnt/oxphos.csv")
oxphos <-  pull(oxphos,Column2)
oxphos <- str_to_title(oxphos)
oxphos <-unique(oxphos)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list) %>% 
  unique()

glycolysis <- read.csv("GeneOnt/glycolysis.csv")
glycolysis <-  pull(glycolysis,Column2)
glycolysis <- str_to_title(glycolysis)
glycolysis <-unique(glycolysis)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list)%>% 
  unique()

glycolysis <- c("HK1", "HK2", "HK3", "GCK", "GPI1", "PFKL", 
                "PFKM", "PFKP", "ALDOA", "ALDOB", "ALDOC", 
                "TPI1", "GAPDH", "GAPDHS", "PGK1", "ENO1", "ENO2", "ENO3", 
                "PKLR", "PKM", "ADPGK", "BPGM", "FOXK1", "FOXK2", 
                "HIF1A", "LDHA", "LIN28A", "LIN28B", "MYC", "PDK1",
                "PDK3", "PGAM1", "PGAM2", "PGM1", "PGM2L1", "PYGL")
glycolysis <- str_to_title(glycolysis)


##### Scores 

mouse_data <- AddModuleScore(mouse_data, features = glycolysis, 
                             name = "glycolysis")
mouse_data <- AddModuleScore(mouse_data, features = oxphos, 
                             name = "oxphos")


######Plotting 
scores <- data.frame(mouse_data@meta.data$oxphos1,mouse_data@meta.data$glycolysis1, 
                     mouse_data@active.ident)

scores <- reshape2::melt(scores, measure.vars=c(1:2), variable.name="score")

colnames(scores) <- c("CT", "score", "value")

scores$CT <- factor(scores$CT , levels=c("ES2i", "ES", "Epi", "N"))

write.csv(scores, "PAPER/Scores/Chen_mouse_scores.csv")

ggplot(scores, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.25, aes(fill=score))+
  scale_fill_manual(values = c("#289aeb", "#eb2838"),
                    name="", labels = c("OxPhos", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("")+
  ylim(-0.32,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()

ggsave("PAPER/Scores/chen_scores.pdf", width=6, height=4)

mark_m <- FindMarkers(mouse_data, ident.1="M", ident.2=c("ES2i","ES","Epi_delay","Epi"))
mark_m

write.csv(glycolysis, "Mouse files/Chen_glycolysis.csv")
write.csv(oxphos, "Mouse files/Chen_oxphos.csv")
write.csv(glycolysis_to_lactate, "Mouse files/Chen_glycolysis_to_lactate.csv")

mouse_data <- FindVariableFeatures(mouse_data, selection.method = "vst", nfeatures = 30000)

levels(mouse_data) <- c("ES2i", "ES", "Epi", "N")

levels(mouse_data2) <- c("E3.5", "E4.5_Epiblast", "E5.5_Em", "E6.5_Em", "E4.5_PrE", "E6.5_VE")

####################################################################################
oxphos <- sort(oxphos$.)
glycolysis <- sort(glycolysis$.)

avgexp_mouse_oxphos <- AverageExpression(mouse_data, features=oxphos, use.scale=FALSE)
avgexp_mouse_glycolysis <- AverageExpression(mouse_data, features=glycolysis, use.scale=FALSE)

avgexp_mouse_glycolysis$glycolysis <- "TRUE"
avgexp_mouse_glycolysis$oxphos <- "FALSE"

avgexp_mouse_oxphos$oxphos <- "TRUE"
avgexp_mouse_oxphos$glycolysis <- "FALSE"

avgexp_mouse_glycolysis <- as.data.frame(avgexp_mouse_glycolysis)
avgexp_mouse_oxphos <- as.data.frame(avgexp_mouse_oxphos)

avgexp <- rbind(avgexp_mouse_glycolysis, avgexp_mouse_oxphos)

write.csv(avgexp, "PAPER/avgexp_mouse_invitro.csv")

glycolysis <- str_to_title(glycolysis)


DoHeatmap(mouse_data, features=oxphos, group.colors = col)
ggsave("PAPER/Heatmaps/Chen_heatmap_oxphos.pdf", width = 8, height=12)
DoHeatmap(mouse_data, features=glycolysis, group.colors = col)
ggsave("PAPER/Heatmaps/Chen_heatmap_gly.pdf", width = 8, height=6)

meta<- deng_data@meta.data 
meta<- mouse_data@meta.data
meta <-mouse_moh@meta.data 

###################################################################################
#Mouse data from the 2017 Cell Reports papers 

mouse_counts <- read.table("Mouse files Cell 2017/GSE100597/GSE100597_count_table_QC_filtered.txt")

#### change gene ids from gene-chr-location to just gene
names <- rownames(mouse_counts) 
names <- str_split_fixed(names, "_", n=3)
names <- as.data.frame(names)
names <- pull(names, V1)
mouse_counts$gene <- names
mouse_counts <- mouse_counts %>%
  select(gene, everything())

mouse_counts$gene <- make.unique(mouse_counts$gene)
rownames(mouse_counts) <- mouse_counts$gene
mouse_counts <- select(mouse_counts, -gene)

#### making a seurat object
mouse_moh <- CreateSeuratObject(counts = mouse_counts,
                                  assay = "RNA",min.cells = 0, min.features = 0)

key_mohammed <- read.table("PAPER/Datasets/Mouse/Mohammed_key.csv", sep=",", header=T, row.names=1)
labs <- key_mohammed
mouse_moh$CT <- labs

Idents(mouse_moh) <- mouse_moh@meta.data$CT

#mouse_moh@active.ident

################################################################################
#re-clustering of Mohammed data
#log normalise
mouse_moh <- NormalizeData(mouse_moh, verbose = FALSE)
#Get variable genes (reduce number)
mouse_moh <- FindVariableFeatures(mouse_moh, selection.method = "vst", nfeatures = 30000)
gene_list <- mouse_moh@assays[["RNA"]]@var.features

mouse_moh <- FindVariableFeatures(mouse_moh, selection.method = "vst", nfeatures = 20000)

#Rename to keep variable convention for the join species modelling
mouse_moh <- ScaleData(mouse_moh, verbose = FALSE)

#Dimesionality reduction and clustering
mouse_moh <- RunPCA(mouse_moh, npcs = 30, verbose = FALSE) #,features = rownames(Markers))

#plot pca
DimPlot(mouse_moh, dims=c(1,2), reduction="pca", label=TRUE, repel=TRUE)

#FeaturePlot(mouse_moh, features="Nanog")

mouse_moh <- FindNeighbors(mouse_moh, dims = 1:20, k.param = 10)
mouse_moh <- FindClusters(mouse_moh, resolution=0.2)

#mouse_moh_mark <- FindAllMarkers(mouse_moh, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#Get top 100
#top20 <- mouse_moh_mark %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
#Plot heatmap
#DoHeatmap(mouse_moh, features = top20$gene) + NoLegend()
#ggsave(filename=paste("Markers_mouse_moh",".pdf",sep=""),width = 40, height = 100,limitsize = FALSE)
#FeaturePlot(mouse_moh, features=c("Myo6"))

#DimPlot(mouse_moh, dims=c(1,2), reduction="pca", label=TRUE, repel=TRUE)
#ggsave("Mouse files Cell 2017/mouse_moh_PCA.pdf", width=6, height=4)

#DotPlot(mouse_moh, features=c("T", "Snai1", "Lef1", "Evx1","Mesp1", 
#                                "Greb1", "Ncoa5", "Eid1", "H2afy2"))

cells <- WhichCells(mouse_moh, idents = c(2,8))
mouse_moh <- SetIdent(
  object = mouse_moh,
  cells = cells ,
  value = 'E3.5')

cells <- WhichCells(mouse_moh, idents = c(4))
mouse_moh <- SetIdent(
  object = mouse_moh,
  cells = cells ,
  value = 'E4.5_PrE')

cells <- WhichCells(mouse_moh, idents = c(6))
mouse_moh <- SetIdent(
  object = mouse_moh,
  cells = cells ,
  value = 'E4.5_Epiblast')

cells <- WhichCells(mouse_moh, idents = c(0))
mouse_moh <- SetIdent(
  object = mouse_moh,
  cells = cells ,
  value = 'E5.5')

cells <- WhichCells(mouse_moh, idents = c(5))
mouse_moh <- SetIdent(
  object = mouse_moh,
  cells = cells ,
  value = 'E6.5_VE')

cells <- WhichCells(mouse_moh, idents = c(7,3,1))
mouse_moh <- SetIdent(
  object = mouse_moh,
  cells = cells ,
  value = 'E6.5_Em')

key_mohammed <- mouse_moh@active.ident
write.csv(key_mohammed,"PAPER/Datasets/Mouse/Mohammed_key.csv")


DimPlot(mouse_moh)

mouse_moh$dataset <- "Mohammed"

####################################################################################################
#merging 
#mouse_1 <- merge(mouse_moh, mouse_data)
mouse_combined <- merge(mouse_moh, deng_data)

mouse_combined <- NormalizeData(mouse_combined, verbose = FALSE)
#Get variable genes (reduce number)
mouse_combined <- FindVariableFeatures(mouse_combined, selection.method = "vst", nfeatures = 20000)

#scaling 
mouse_combined <- ScaleData(mouse_combined, verbose = FALSE)

cType <- c("Zy","4C","8C","E3.5","E4.5_Epiblast", "E5.5", "E6.5_Em","E4.5_PrE","E6.5_VE")

BaseCol <- c("#82F5F1", "#A6D9E7", "#5AC0CE",
             "#5B859E", "#2AB6B9", "#3F94D1", "#3064AD", "#E5B405", "#E94E12")

VlnPlot(mouse_combined, features= "Hif1a")  

#cType <- c("E3.5","E4.5_Epiblast", "E5.5", "E6.5_Em","E4.5_PrE","E6.5_VE", 
#           "E5.5_emb", "E6.0_emb", "E6.5_emb", "E7.0_emb", "E7.5_emb", "E7.0_end",
#           "ES2i", "ES", "Epi", "N")

#BaseCol <- c("#BCDFD6", "#5C859F", "#2DB6B8", "#3165AE", "#E5B40C", "#E95017",
#             "#2DB6B8", "#3F94D1", "#3165AE","#197e80" , "#284597", "#5554A0",
#             "#F9CBCC", "#D3A5B4", "#F193AF","#E42386")

mouse_combined@active.ident


colind <- integer(length(levels(Idents(mouse_combined))))
for (i in 1:length( levels(Idents(mouse_combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mouse_combined))[i])
}
coluse <- BaseCol[colind]

#Dimesionality reduction and clustering
mouse_combined <- RunPCA(mouse_combined, npcs = 30, verbose = FALSE) #,features = rownames(Markers))
levels(mouse_combined) <- c("Zy", "4C", "8C",
                            "E3.5","E4.5_Epiblast", "E5.5", "E6.5_Em","E4.5_PrE","E6.5_VE")

DimPlot(mouse_combined, dims=c(1,2), pt.size = 2.5, label=TRUE, repel=TRUE, shape.by = "dataset", cols=coluse)
ggsave("PAPER/PCAs/mouse_deng_mohammed_new.pdf", width=6, height=4)


pct <- mouse_combined@reductions[["pca"]]@stdev /sum(mouse_combined@reductions[["pca"]]@stdev) *100
  
  


###### Mammalian scores 
mouse_combined <- FindVariableFeatures(mouse_combined, selection.method = "vst", nfeatures = 30000)
gene_list <- mouse_combined@assays[["RNA"]]@var.features

oxphos <- read.csv("GeneOnt/oxphos.csv")
oxphos <-  pull(oxphos,Column2)
oxphos <- str_to_title(oxphos)
oxphos <-unique(oxphos)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list) %>% 
  unique()

glycolysis <- read.csv("GeneOnt/glycolysis.csv")
glycolysis <-  pull(glycolysis,Column2)
glycolysis <- str_to_title(glycolysis)
glycolysis <-unique(glycolysis)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list)%>% 
  unique()

mouse_combined <- AddModuleScore(mouse_combined, features = glycolysis, 
                              name = "glycolysis")
mouse_combined <- AddModuleScore(mouse_combined, features = oxphos, 
                              name = "oxphos")

scores_mouse <- data.frame(mouse_combined@meta.data$oxphos1,mouse_combined@meta.data$glycolysis1, 
                            mouse_combined@active.ident)
scores_mouse <- rownames_to_column(scores_mouse,var="cell")

scores_mouse <- reshape2::melt(scores_mouse, measure.vars=c(2:3), variable.name="score")

colnames(scores_mouse) <- c("cell","CT", "score", "value")

write.csv(scores_mouse, "PAPER/Scores/mouse_deng_mohammed_scores.csv")

scores_mouse <- filter(scores_mouse, scores_mouse$CT %in% c("Zy", "4C", "8C",
                                "E3.5","E4.5_Epiblast", "E5.5", "E6.5_Em"))

ggplot(scores_mouse, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.25, aes(fill=score))+
  scale_fill_manual(values = c("#289aeb", "#eb2838"),
                    name="", labels = c("OxPhos", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("")+
  ylim(-0.32,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()
ggsave("PAPER/Scores/Good_scores_deng_mohammed_all.pdf", width=6, height=4)


#################################################################################
# Heatmap
oxphos <- sort(oxphos$.)
glycolysis <- sort(glycolysis$.)

oxphos <- sort(oxphos$.)
glycolysis <- sort(glycolysis$.)

glycolysis_order <- read.csv("PAPER/Scores/glycolysis_order.csv")
glycolysis <- glycolysis_order

avgexp_mouse_oxphos <- AverageExpression(mouse_combined, features=oxphos, use.scale=FALSE)
avgexp_mouse_glycolysis <- AverageExpression(mouse_combined, features=glycolysis, use.scale=FALSE)

avgexp_mouse_glycolysis$glycolysis <- "TRUE"
avgexp_mouse_glycolysis$oxphos <- "FALSE"

avgexp_mouse_oxphos$oxphos <- "TRUE"
avgexp_mouse_oxphos$glycolysis <- "FALSE"

avgexp_mouse_glycolysis <- as.data.frame(avgexp_mouse_glycolysis)
avgexp_mouse_oxphos <- as.data.frame(avgexp_mouse_oxphos)

avgexp <- rbind(avgexp_mouse_glycolysis, avgexp_mouse_oxphos)

write.csv(avgexp, "PAPER/avgexp_mouse_invivo.csv")

########################################################################################
mouse_combined <- ScaleData(mouse_combined, verbose = FALSE)

levels (mouse_combined) <- c("Zy", "4C", "8C",
                             "E3.5","E4.5_Epiblast", "E5.5", "E6.5_Em","E4.5_PrE","E6.5_VE")

glycolysis <- glycolysis$.

DoHeatmap(mouse_combined, features = oxphos, group.colors = coluse, raster=FALSE)
ggsave("PAPER/Heatmaps/Mouse_final_final_heatmap_oxphos_rasterfalse.pdf", width=6, height=12)

DoHeatmap(mouse_combined, features = glycolysis, group.colors = coluse)
ggsave("PAPER/Heatmaps/mouse_final_heatmap_glyco_AA.pdf", width=6, height=6)

#############################################################################################################################################################
#pheatmap test
oxphos <- sort(oxphos$.)
glycolysis <- sort(glycolysis$.)

heatmap_df<- mouse_combined[["RNA"]]@scale.data
heatmap_df <- as.data.frame(heatmap_df)
heatmap_df <- rownames_to_column(heatmap_df, "genes") 
heatmap_df <- filter(heatmap_df, heatmap_df$genes %in% glycolysis)
heatmap_df <- column_to_rownames(heatmap_df, var="genes")
heatmap_df <- as.matrix(heatmap_df)

annotation_df <- as.data.frame(mouse_combined$CT)
annotation_t <- t(annotation_df)


annotation_df <- factor(annotation_df, levels = c("Zy", "4C", "8C",
                           "E3.5","E4.5_Epiblast", "E5.5", "E6.5_Em","E4.5_PrE","E6.5_VE"))

heatmap_df$
heatmap_df <- factor(annotation_df, levels = c("Zy", "4C", "8C",
                                                  "E3.5","E4.5_Epiblast", "E5.5", "E6.5_Em","E4.5_PrE","E6.5_VE"))



pheatmap(mat = heatmap_df, show_colnames= FALSE, 
         cluster_rows = FALSE, cluster_cols = FALSE, 
         annotation_col = annotation_df, kmeans_k = NA, 
         color= redblue1(20), breaks= mat_breaks, border_color = NA)

color =  redblue1(20),gaps_col=c(4,9,14,20,24,28),
gaps_row=c(3,6,9,17,23,34,39),breaks = mat_breaks, border_color = NA,
annotation_col = annotation_col, annotation_colors = anno_colors, 
cluster_rows=FALSE,cluster_cols=FALSE,  filename = "~/Desktop/HM1.pdf"



mat_breaks <- seq(0, 2, length.out = 20)



#quantile_breaks <- function(xs, n = 10) {
#  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
#  breaks[!duplicated(breaks)]
#}

#mat_breaks <- quantile_breaks(heatmap_df, n = 10)

library(RColorBrewer)
pal <- brewer.pal(n=8,name="RdBu")