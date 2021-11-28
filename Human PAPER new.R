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

###### 
#loading data

key_petro <- read.table("PAPER/Datasets/Human/Stirparo_petro_key.csv", row.names=1, sep=",", header=T)

key_blake <- read.table("PAPER/Datasets/Human/Stirparo_blake_key.csv", sep=",", header=T, row.names=1)

key_yan <- read.table("PAPER/Datasets/Human/Stirparo_yan_key.csv", sep=",", header=T, row.names=1)


petro_human <-read.csv("PAPER/Datasets/Human/petrocounts.csv", sep=",", header=T, row.names=1)

blake_human <- read.table("PAPER/Datasets/Human/featurecountsNiakan.txt", sep="", header=T, row.names=1)

yan_human <- read.table("PAPER/Datasets/Human/featurecountsYan.txt", sep="", header=T, row.names=1)

#####################################################################################
#Yan data
cells <- rownames(key_yan)
yan_human <- select(yan_human, c(all_of(cells)))

labs_yan <- key_yan$Ident_2

yan_data <- CreateSeuratObject(counts = yan_human, assay = "RNA",min.cells = 0, min.features = 0)

yan_data$CT <- labs_yan

yan_data$dataset <- "Yan"

Idents(yan_data) <- yan_data@meta.data$CT #Change the annotation of each cell to e.g. its type Amnnion etc

#Blakeley data
cells <- rownames(key_blake)

blake_human <- select(blake_human, c(all_of(cells)))

labs_blake <- key_blake$Ident_2

blake_data <- CreateSeuratObject(counts = blake_human, assay = "RNA",min.cells = 0, min.features = 0)

blake_data$CT <- labs_blake

blake_data$dataset <- "Blakeley"

Idents(blake_data) <- blake_data@meta.data$CT #Change the annotation of each cell to e.g. its type Amnnion etc

#Petropoulus data
cells <- rownames(key_petro)

petro_human <- select(petro_human, c(all_of(cells)))

labs_petro <- key_petro$Ident_2

petro_data <- CreateSeuratObject(counts = petro_human, assay = "RNA",min.cells = 0, min.features = 0)
petro_data$CT <- labs_petro

petro_data$dataset <- "Petropoulus"


Idents(petro_data) <- petro_data@meta.data$CT #Change the annotation of each cell to e.g. its type Amnnion etc

##########################################################################################
#CS7 human data
CS7_counts <- read_rds("PAPER/Datasets/Human/Human_CS7.rds")
CS7_subset <- read_rds("PAPER/Datasets/Human/HumanCS7_subs.rds")

CS7_counts <- t(CS7_counts)
CS7_counts <- as.data.frame(CS7_counts)

#write.csv(CS7_counts, "PAPER/Datasets/Human/Human_CS7_counts.csv")
#CS7_counts <- read.csv("PAPER/Datasets/Human/Human_CS7_counts.csv") 


rownames(CS7_counts)[rownames(CS7_counts) == "ATP5F1A"] <- "ATP5A1"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5F1B"] <- "ATP5B"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5F1C"] <- "ATP5C1"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5F1D"] <- "ATP5D"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5MC1"] <- "ATP5G1"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5MC3"] <- "ATP5G3"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5PO"] <- "ATP5O"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5PB"] <- "ATP5F1"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5PD"] <- "ATP5H"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5MF"] <- "ATP5J2"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5MG"] <- "ATP5L"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5PF"] <- "ATP5J"

CS7_data <- CreateSeuratObject(counts = CS7_counts, min.cells=1 , assay = "RNA", min.features = 0)

CS7_data$dataset <- "Tyser"

key <- CS7_subset@meta.data 

cells <- rownames(key)
#cells <-str_sub(cells, 2, str_length(cells))

CS7_data <- subset(CS7_data, cells=cells)

labs <- CS7_subset@active.ident

CS7_data$CT <- labs 

Idents(CS7_data) <- CS7_data@meta.data$CT

########################################################################################
#Merging
blake_yan <- merge(x=blake_data, y=yan_data)

human_combined1 <- merge(x=blake_yan, y=CS7_data)

human_combined <- merge(x=human_combined1, y=petro_data)

remove(blake_yan)
remove(human_combined1)

human_combined <- RenameIdents(object = human_combined, 'early_Tb_CS3' = 'Tb_CS3')
human_combined <- RenameIdents(object = human_combined, 'late_Tb_CS3' = 'Tb_CS3')
human_combined <- RenameIdents(object = human_combined, 'late_PrE_CS3' = 'Hyp_CS3')
human_combined <- RenameIdents(object = human_combined, 'Epi' = 'Epi_CS7')


human_combined <- subset(human_combined, 
                         idents = c("Zy_CS1", "4cell_CS2", "8cell_CS2", 
                                    "cMor_CS3", "early_ICM_CS3", "late_Epi_CS3",
                                    "Epi_CS7", "Tb_CS3", "Hyp_CS3"))
                                                  

levels(human_combined) <- c("Zy_CS1", "4cell_CS2", "8cell_CS2", 
                            "cMor_CS3", "early_ICM_CS3", "late_Epi_CS3",
                            "Epi_CS7", "Tb_CS3", "Hyp_CS3")


cType <- c("Zy_CS1", "4cell_CS2", "8cell_CS2", 
           "cMor_CS3", "early_ICM_CS3", "late_Epi_CS3",
           "Epi_CS7", "Tb_CS3", "Hyp_CS3")

BaseCol <-c("#7dfff8", "#A6D9E7", '#BBDFD5', 
            '#5AC0CE', '#5B859E', '#2AB6B9', 
            '#284597', "#F7EF8E", '#E5B405')

colind <- integer( length( levels(Idents(human_combined)) )  )
for (i in 1:length( levels(Idents(human_combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(human_combined))[i])
}
coluse <- BaseCol[colind]

#Normalise the data
human_combined <- NormalizeData(human_combined, verbose = FALSE)
#Get variable genes (reduce number)
human_combined <- FindVariableFeatures(human_combined, selection.method = "vst", nfeatures = 20000)

gene_list_human <- human_combined@assays[["RNA"]]@var.features

#scaling 
human_combined <- ScaleData(human_combined, verbose = FALSE)

#Idents(human_combined) <- human_combined@meta.data$CT #Change the annotation of each cell to e.g. its type Amnnion etc

#Dimesionality reduction and clustering
human_combined <- RunPCA(human_combined, npcs = 30, verbose = FALSE) #,features = rownames(Markers))
DimPlot(human_combined, reduction = "pca", cols=coluse, pt.size= 2, dims=c(1,2), shape.by="dataset", repel = TRUE, label = FALSE)
ggsave("PAPER/PCAs/Human_final_10.pdf", width=6, height=4)


pct <- human_combined@reductions[["pca"]]@stdev/sum(human_combined@reductions[["pca"]]@stdev) *100

human_combined <- FindVariableFeatures(human_combined, selection.method = "vst", nfeatures = 40000)

gene_list_human <- human_combined@assays[["RNA"]]@var.features

oxphos <- read.csv("GeneOnt/oxphos.csv")
oxphos <-  pull(oxphos,Column2)
oxphos <- str_to_upper(oxphos)
oxphos <-unique(oxphos)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_human)%>% 
  unique()

glycolysis <- read.csv("GeneOnt/glycolysis.csv")
glycolysis <-  pull(glycolysis,Column2)
glycolysis <- str_to_upper(glycolysis)
glycolysis <-unique(glycolysis)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_human)%>% 
  unique()

write.csv(glycolysis, "PAPER/Gene modules/glycolysis_human_Yan_blakeley_tyser.csv")
write.csv(oxphos,"PAPER/Gene modules/oxphos_human_yan_blakeley_tyser.csv")

human_combined <- AddModuleScore(human_combined, features = oxphos, 
                           name = "oxphos")
human_combined <- AddModuleScore(human_combined, features = glycolysis, 
                             name = "glycolysis")

human_scores <- data.frame(human_combined@meta.data$oxphos1,human_combined@meta.data$glycolysis1,
                           human_combined@active.ident)

colnames(human_scores) <- c("oxphos","glycolysis", "CT")

human_scores <- rownames_to_column(human_scores,var="cell")

human_scores <- reshape2::melt(human_scores, measure.vars=c(2:3), variable.name="score")

colnames(human_scores) <- c("cell","CT", "score", "value")


human_scores_extraemb <- filter(human_scores, human_scores$CT %in% extraemb)

emb<- c("Zy_CS1", "4cell_CS2", "8cell_CS2", "cMor_CS3", "early_ICM_CS3",
        "late_Epi_CS3", "Epi_CS7")

extraemb <- c("Tb_CS3", "Hyp_CS3")
write.csv(human_scores, "PAPER/Scores/human_scores_final_all.csv")

ggplot(human_scores, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.25, aes(fill=score))+
  scale_fill_manual(values = c("#289aeb", "#eb2838"),
                    name="", labels = c("OxPhos", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("")+
  ylim(-0.32,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()

ggsave("PAPER/Scores/human_all_final_extra_new.pdf", width=6, height=4)
cells <- human_combined@meta.data
write.csv(human_scores,"PAPER/Scores/human_all_figure.csv")


#### avg exp
oxphos <- sort(oxphos$.)
glycolysis <- sort(glycolysis$.)

avgexp_human_oxphos <- AverageExpression(human_combined, features=oxphos, use.scale=FALSE)
avgexp_human_glycolysis <- AverageExpression(human_combined, features=glycolysis, use.scale=FALSE)

avgexp_human_glycolysis$glycolysis <- "TRUE"
avgexp_human_glycolysis$oxphos <- "FALSE"

avgexp_human_oxphos$oxphos <- "TRUE"
avgexp_human_oxphos$glycolysis <- "FALSE"

avgexp_human_glycolysis <- as.data.frame(avgexp_human_glycolysis)
avgexp_human_oxphos <- as.data.frame(avgexp_human_oxphos)

avgexp <- rbind(avgexp_human_glycolysis, avgexp_human_oxphos)

write.csv(avgexp, "PAPER/avg_exp_human_revised.csv")

##########################################################################################
#CS7 human data
CS7_counts <- read_rds("PAPER/Datasets/Human/Human_CS7.rds")
CS7_subset <- read_rds("PAPER/Datasets/Human/HumanCS7_subs.rds")

CS7_counts <- t(CS7_counts)
CS7_counts <- as.data.frame(CS7_counts)

write.csv(CS7_counts, "PAPER/Datasets/Human/Human_CS7_counts.csv")

rownames(CS7_counts)[rownames(CS7_counts) == "ATP5F1A"] <- "ATP5A1"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5F1B"] <- "ATP5B"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5F1C"] <- "ATP5C1"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5F1D"] <- "ATP5D"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5MC1"] <- "ATP5G1"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5MC3"] <- "ATP5G3"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5PO"] <- "ATP5O"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5PB"] <- "ATP5F1"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5PD"] <- "ATP5H"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5MF"] <- "ATP5J2"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5MG"] <- "ATP5L"
rownames(CS7_counts)[rownames(CS7_counts) == "ATP5PF"] <- "ATP5J"

CS7_data <- CreateSeuratObject(counts = CS7_counts, min.cells=1 , assay = "RNA", min.features = 0)

CS7_data$dataset <- "Tyser"

key <- CS7_subset@meta.data 

cells <- rownames(key)
#cells <-str_sub(cells, 2, str_length(cells))

CS7_data <- subset(CS7_data, cells=cells)

labs <- CS7_subset@active.ident

CS7_data$CT <- labs 

Idents(CS7_data) <- CS7_data@meta.data$CT

CS7_data <- NormalizeData(CS7_data, verbose = FALSE)
#Get variable genes (reduce number)

CS7_data <- subset (CS7_data, idents=("Epi"))

CS7_data <- FindVariableFeatures(CS7_data, selection.method = "vst", nfeatures = 30000)

gene_list_human <- CS7_data@assays[["RNA"]]@var.features

#scaling 
CS7_subset <- ScaleData(CS7_subset, verbose = FALSE)

CS7_subset <- FindVariableFeatures(CS7_subset, selection.method = "vst", nfeatures = 40000)

gene_list_human <- CS7_subset@assays[["RNA"]]@var.features


#Idents(CS7_subset) <- CS7_subset@meta.data$CT #Change the annotation of each cell to e.g. its type Amnnion etc

#Dimesionality reduction and clustering
CS7_subset <- RunPCA(CS7_subset, npcs = 30, verbose = FALSE) #,features = rownames(Markers))
DimPlot(CS7_subset, reduction = "pca", dims=c(1,2), repel = TRUE, label = FALSE)
ggsave("PAPER/PCAs/CS7_human_full_pca.pdf", width=6, height=4)

oxphos <- read.csv("GeneOnt/oxphos.csv")
oxphos <-  pull(oxphos,Column2)
oxphos <- str_to_upper(oxphos)
oxphos <-unique(oxphos)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_human)%>% 
  unique()

glycolysis <- read.csv("GeneOnt/glycolysis.csv")
glycolysis <-  pull(glycolysis,Column2)
glycolysis <- str_to_upper(glycolysis)
glycolysis <-unique(glycolysis)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_human)%>% 
  unique()

write.csv(glycolysis, "PAPER/Gene modules/glycolysis_human_Cs7.csv")
write.csv(oxphos,"PAPER/Gene modules/oxphos_human_Cs7.csv")

CS7_data <- AddModuleScore(CS7_data, features = oxphos, 
                           name = "oxphos")
CS7_data <- AddModuleScore(CS7_data, features = glycolysis, 
                           name = "glycolysis")

CS7 <- data.frame(CS7_data@meta.data$oxphos1,CS7_data@meta.data$glycolysis1,
                           CS7_data@active.ident)
colnames(CS7) <- c("oxphos","glycolysis", "CT")


CS7 <- rownames_to_column(CS7,var="cell")

CS7 <- reshape2::melt(CS7, measure.vars=c(2:3), variable.name="score")

colnames(CS7) <- c("cell","CT", "score", "value")
#
#########################################################################################
# Xiang et al 2020 data

xiang_key<-read.table("Datasets/Li/humanIDs.csv",sep=",",header = T, row.names=1)
xiang_counts<- read.table("Datasets/Li/featurecountsLi.csv", sep=",",header = T, row.names=1)

xiang_counts <- select(xiang_counts, -Length)
xiang_data <- CreateSeuratObject(counts = xiang_counts, assay = "RNA",min.cells = 0, min.features = 0)

xiang_key <- xiang_key$Type3

xiang_data$CT <- xiang_key
Idents(xiang_data) <- xiang_data@meta.data$CT 

xiang_data$dataset <- "Xiang"

xiang_data <- NormalizeData(xiang_data, verbose = FALSE)
#Get variable genes (reduce number)

xiang_data <- FindVariableFeatures(xiang_data, selection.method = "vst", nfeatures = 20000)

#scaling 
xiang_data <- ScaleData(xiang_data, verbose = FALSE)

gene_list_human <- xiang_data@assays[["RNA"]]@var.features

levels(xiang_data@active.ident)


levels(xiang_data) <- c("ICM_CS3", "EPI_CS3", "EPI_CS4", "EmDisc_CS5", "EmDisc_CS6", 
                        "Tb_CS4",  "Tb_CS5", "Tb_CS6", "HYP_CS3", "VE_CS4", 
                        "VE_CS5", "SYS_CS6")

xiang_data <- RunPCA(xiang_data, npcs = 30, verbose = FALSE) #,features = rownames(Markers))
DimPlot(xiang_data, reduction = "pca", cols=col, dims=c(1,2), repel = TRUE, label = FALSE)

pct <- xiang_data@reductions[["pca"]]@stdev/sum(xiang_data@reductions[["pca"]]@stdev) *100

col <- c("#5AC0CE", "#5B859E", "#2AB6B9", "#3F94D1",
        "#284597", 
        "#F7EF8E", "#f5d950", "#E3DF1E", 
        "#E5B405", "#E94E12", "#B43C0E", 
        "#E58603")
ggsave("PAPER/PCAs/xiang.pdf", width=6, height=4)


oxphos <- read.csv("GeneOnt/oxphos.csv")
oxphos <-  pull(oxphos,Column2)
oxphos <- str_to_upper(oxphos)
oxphos <-unique(oxphos)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_human)%>% 
  unique()

glycolysis <- read.csv("GeneOnt/glycolysis.csv")
glycolysis <-  pull(glycolysis,Column2)
glycolysis <- str_to_upper(glycolysis)
glycolysis <-unique(glycolysis)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_human)%>% 
  unique()

oxphos <- sort(oxphos$.)
glycolysis <- sort(glycolysis$.)


avgexp_xiang_oxphos <- AverageExpression(xiang_data, features=oxphos, use.scale=FALSE)
avgexp_xiang_glycolysis <- AverageExpression(xiang_data, features=glycolysis, use.scale=FALSE)

avgexp_xiang_glycolysis$glycolysis <- "TRUE"
avgexp_xiang_glycolysis$oxphos <- "FALSE"

avgexp_xiang_oxphos$oxphos <- "TRUE"
avgexp_xiang_oxphos$glycolysis <- "FALSE"

avgexp_xiang_glycolysis <- as.data.frame(avgexp_xiang_glycolysis)
avgexp_xiang_oxphos <- as.data.frame(avgexp_xiang_oxphos)

avgexp_xiang <- rbind(avgexp_xiang_glycolysis, avgexp_xiang_oxphos)

write.csv(avgexp, "PAPER/avgexp_human.csv")
write.csv(avgexp_xiang, "PAPER/avgexp_human_xiang.csv")

meta <- xiang_data@meta.data 

write.csv(meta, "PAPER/meta_xiang.csv")

xiang_data <- AddModuleScore(xiang_data, features = oxphos, 
                           name = "oxphos")
xiang_data <- AddModuleScore(xiang_data, features = glycolysis, 
                           name = "glycolysis")

human_scores <- data.frame(xiang_data@meta.data$oxphos1,xiang_data@meta.data$glycolysis1,
                           xiang_data@active.ident)

human_scores <- rownames_to_column(human_scores,var="cell")

human_scores <- reshape2::melt(human_scores, measure.vars=c(2:3), variable.name="score")

colnames(human_scores) <- c("cell","CT", "score", "value")

write.csv(human_scores, "PAPER/Scores/Scores_human_xiang.csv")

ggplot(human_scores, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.25, aes(fill=score))+
  scale_fill_manual(values = c("#289aeb", "#eb2838"),
                    name="", labels = c("OxPhos", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("")+
  ylim(-0.32,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()
ggsave("PAPER/Scores/xiang_scores_final.pdf", width=6, height=4)

DoHeatmap(xiang_data, features=glycolysis, group.colors=coluse)
ggsave("PAPER/Heatmaps/xiang_glyco.pdf", width = 8, height=6)

DoHeatmap(xiang_data, features=oxphos, group.colors = coluse)
ggsave("PAPER/Heatmaps/xiang_ox.pdf", width = 8, height = 12)


##########################################################################################
human_combined1 <- merge(x=human_combined, y=xiang_data)

#Normalise the data
human_combined1 <- NormalizeData(human_combined1, verbose = FALSE)
#Get variable genes (reduce number)
human_combined1 <- FindVariableFeatures(human_combined1, selection.method = "vst", nfeatures = 30000)

gene_list_human <- human_combined@assays[["RNA"]]@var.features

human_combined1 <- FindVariableFeatures(human_combined1, selection.method = "vst", nfeatures = 20000)

#scaling 
human_combined1 <- ScaleData(human_combined1, verbose = FALSE)

#Idents(human_combined) <- human_combined@meta.data$CT #Change the annotation of each cell to e.g. its type Amnnion etc

oxphos <- read.csv("GeneOnt/oxphos.csv")
oxphos <-  pull(oxphos,Column2)
oxphos <- str_to_upper(oxphos)
oxphos <-unique(oxphos)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_human1)%>% 
  unique()

glycolysis <- read.csv("GeneOnt/glycolysis.csv")
glycolysis <-  pull(glycolysis,Column2)
glycolysis <- str_to_upper(glycolysis)
glycolysis <-unique(glycolysis)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_human)%>% 
  unique()
View(glycolysis)

#Dimesionality reduction and clustering
#human_combined <- RunPCA(human_combined, npcs = 30, verbose = FALSE) #,features = rownames(Markers))

#levels(human_combined1)<- c("Zy_CS1", "4cell_CS2", "8cell_CS2", 
#                                     "cMor_CS3", 
#                                     "early_ICM_CS3", "late_Epi_CS3", 
#                                     "Tb_CS3", 
#                             "late_PrE_CS3", "ICM_CS3", "EPI_CS3", "EPI_CS4", "EmDisc_CS5", 
#                                     "EmDisc_CS6", "Tb_CS4", "Tb_CS5", "Tb_CS6",
#                                     "HYP_CS3", "VE_CS4", "VE_CS5",
#                                     "SYS_CS6")

#DimPlot(human_combined, reduction = "pca", cols=coluse, dims=c(1,2), shape.by="dataset", repel = TRUE, label = FALSE)

#ggsave("PAPER/PCAs/Human_test_7.pdf", width=8, height=6)

#######################################################################################
#colours

cType <- c("Zy_CS1", "4cell_CS2", "8cell_CS2", 
           "cMor_CS3", 
           "early_ICM_CS3", "late_Epi_CS3", 
           "early_Tb_CS3", "late_Tb_CS3", 
           "late_PrE_CS3", "ICM_CS3", "EPI_CS3", "EPI_CS4", "EmDisc_CS5", 
           "EmDisc_CS6", "Tb_CS4", "Tb_CS5", "Tb_CS6",
           "HYP_CS3", "VE_CS4", "VE_CS5",
           "SYS_CS6")

BaseCol <-c("#DAEDEC","#A6D9E7","#BBDFD5",
            "#5AC0CE",
            "#5B859E", "#2AB6B9", 
            "#FFF5C4","#F7EF8E",
            "#E5B405","#82F5F1", "#197E80","#3F94D1", "#3064AD",
            "#284597", "#F5D950", "#E3DF1E", "#E6C800",
            "#E94E12", "#B43C0E", "#DA3843", 
            "#D17612")

colind <- integer( length( levels(Idents(human_combined)) )  )
for (i in 1:length( levels(Idents(human_combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(human_combined))[i])
}
coluse <- BaseCol[colind]

#######################################################################################
#PCA 3D

# Export PC loadings into a dataframe organized by type
embed_mat <- as.data.frame(human_combined@reductions$pca@cell.embeddings)
meta_data <- Idents(human_combined)
stage <- human_combined$dataset

#check order
check_meta = NULL; check_stage = NULL;
for (i in 1:nrow(embed_mat))  {
  check_meta <- append(check_meta,rownames(embed_mat)[i] == labels(meta_data)[i])
  check_stage<- append(check_stage,rownames(embed_mat)[i] == labels(stage)[i])
}
nrow(embed_mat) == sum(check_meta)
nrow(embed_mat) == sum(check_stage)

# add metadata into counts_matrix
pca3d_df <- cbind(meta_data,stage, embed_mat)
colnames(pca3d_df)[1:2] <- c("type", "stage")


names(coluse) = cType[colind]
pca3dcolors = rep(NA,nrow(pca3d_df))
for(j in 1:length(pca3dcolors)){
  temp <- which(rownames(as.data.frame(coluse))  ==  as.character(pca3d_df$type)[j] )
  pca3dcolors[j] <- as.character(as.data.frame(coluse)$coluse)[temp]
}

pca3dpch <- rep(NA,nrow(pca3d_df))
for(j in 1:length(pca3dpch)){
  if(as.character(pca3d_df$stage)[j]   == 'Blakeley') {pca3dpch[j] <- 19}
  if(as.character(pca3d_df$stage)[j]   == 'Petropoulus') {pca3dpch[j] <- 17}
  if(as.character(pca3d_df$stage)[j]   == 'Yan') {pca3dpch[j] <- 3}
  if(as.character(pca3d_df$stage)[j]   == 'Xiang') {pca3dpch[j] <- 15}
  
}

library("plot3D")

# Add small dots on basal plane and on the depth plane
# Can be adapted to plot on any plane
#scatter3D_fancy <- function(x, y, z,..., colvar = z)
#{
#  panelfirst <- function(pmat) {
#    XY <- trans3D(x, y, z = rep(min(z), length(z)), pmat = pmat)
#    scatter2D(XY$x, XY$y, col= coluse, colvar = colvar, pch = ".", 
#              cex = 2, add = TRUE, colkey = FALSE)

# XY <- trans3D(x, y= rep(max(y), length(y)), z, pmat = pmat)
# scatter2D(XY$x, XY$y, col=coluse,  colvar = colvar, pch = ".",
# cex = 2, add = TRUE, colkey = FALSE)
#  }
#  scatter3D(x, y, z, ..., colvar = colvar, panel.first=panelfirst,
#            colkey = FALSE)
#}

#phi = 30
#theta = 270

#print a pdf with lots of versions!
pdf('PAPER/PCAs/3D_human_test_3.pdf', height=10,width=10, useDingbats = FALSE)
#scatter3D_fancy(pca3d_df$PC_1, pca3d_df$PC_2, pca3d_df$PC_3, col=coluse, colvar = as.integer(pca3d_df$type), pch = c(24,21)[pca3d_df$stage],
#                phi = 35, theta = 270, type = 'h',cex = 1.3,  #colkey = FALSE,
#                xlab = "PC1", ylab="PC2", zlab = "PC3", main = "PCA",
#                ticktype = "detailed", nticks = 4)
scatter3D(pca3d_df$PC_1, pca3d_df$PC_2, pca3d_df$PC_3, col=coluse, bty = "b2",
          colvar = as.integer(pca3d_df$type), 
          phi = 40, theta = 215, colkey = FALSE, cex = 0.6 , pch = pca3dpch,
          xlab = "PC1", ylab="PC2", zlab = "PC3", main = "PCA",
          ticktype = "detailed", nticks =4)
dev.off()

###################################################################################
# Scores 

write.csv(glycolysis, "PAPER/Gene modules/glycolysis_human_pgc.csv")
write.csv(oxphos,"PAPER/Gene modules/oxphos_human_pgc.csv")

pgc_data <- AddModuleScore(pgc_data, features = oxphos, 
                             name = "oxphos")
pgc_data <- AddModuleScore(pgc_data, features = glycolysis, 
                             name = "glycolysis")

human_scores <- data.frame(pgc_data@meta.data$oxphos1,pgc_data@meta.data$glycolysis1,
                          pgc_data@active.ident)

human_scores <- rownames_to_column(human_scores,var="cell")

human_scores <- reshape2::melt(human_scores, measure.vars=c(2:3), variable.name="score")

colnames(human_scores) <- c("cell","CT", "score", "value")

write.csv(human_scores, "PAPER/Scores/Scores_human_pgc.csv")

#############
# Plot 

emb <- c("Zy_CS1","4cell_CS2","8cell_CS2", "cMor_CS3", "early_ICM_CS3", "late_Epi_CS3", "ICM_CS3", 
         "EPI_CS3", "EPI_CS4", "EmDisc_CS5", "EmDisc_CS6")

extraemb <- c("early_Tb_CS3","late_Tb_CS3", "Tb_CS5", "Tb_CS6","late_PrE_CS3", 
              "VE_CS4","VE_CS5", "SYS_CS6")

scores_human_pgc <- c ("PGC_4W_M","PGC_4W_F", "PGC_7W_M",
"PGC_8W_F", "PGC_10W_M", "PGC_10W_F", 
"PGC_11W_M", "PGC_11W_F", "PGC_17W_F", 
"PGC_19W_M", "Soma_4W_F", "Soma_7W_M", 
"Soma_8W_F","Soma_10W_M", "Soma_10W_F",
"Soma_11W_M", "Soma_11W_F", "Soma_17W_F", 
"Soma_19W_M")

pgc_fem <- c("PGC_4W_F","PGC_8W_F","PGC_10W_F","PGC_11W_F", "PGC_17W_F", "Soma_8W_F","Soma_10W_F", "Soma_11W_F","Soma_17W_F")
pgc_male <- c("PGC_4W_M", "PGC_7W_M", "PGC_10W_M", "PGC_11W_M", "PGC_19W_M", "Soma_7W_M", "Soma_10W_M","Soma_11W_M", "Soma_19W_M")


scores_human_emb <- human_scores %>% 
  filter(CT %in% emb)

scores_human_emb$CT <- factor(scores_human_emb$CT , levels=c("Zy_CS1","4cell_CS2","8cell_CS2", "cMor_CS3", "early_ICM_CS3", "late_Epi_CS3", "ICM_CS3", 
                                                             "EPI_CS3", "EPI_CS4", "EmDisc_CS5", "EmDisc_CS6"))

scores_human_extraemb <- human_scores %>% 
  filter(CT %in% extraemb)

scores_human_extraemb$CT <- factor(scores_human_extraemb$CT , levels=c("early_Tb_CS3","late_Tb_CS3", "Tb_CS5", "Tb_CS6","late_PrE_CS3", 
                                                                  "VE_CS4","VE_CS5", "SYS_CS6"))


scores_human_fem <- human_scores %>% 
  filter(CT %in% pgc_fem)

scores_human_male <- human_scores %>% 
  filter(CT %in% pgc_male)


ggplot(scores_human_male, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.5, aes(fill=score))+
  scale_fill_manual(values = c("#6280e3", "#b83939"),
                    name="", labels = c("OXPHOS", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("Extraembryonic")+
  ylim(-0.5,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()
ggsave("PAPER/Scores/Human_scores_pgc_only_male.pdf", width=8,height=4)

#########################################################################################
# Heatmaps 

oxphos <- sort(oxphos$.)
glycolysis <- sort(glycolysis$.)

human_combined <- FindVariableFeatures(human_combined, selection.method = "vst", nfeatures = 40000)

#scaling 
human_combined <- ScaleData(human_combined, verbose = FALSE)

DoHeatmap(human_combined, features = oxphos, group.colors = coluse)
ggsave("PAPER/Heatmaps/Human_oxphos_2.pdf", width=8, height=12)

DoHeatmap(human_combined, features = glycolysis, group.colors = coluse)
ggsave("PAPER/Heatmaps/Human_glycolysis_2.pdf", width=8, height=6)

####################################################################################
heatmap_df<- human_combined[["RNA"]]@scale.data
heatmap_df <- as.data.frame(heatmap_df)
heatmap_df <- rownames_to_column(heatmap_df, "genes") 
heatmap_df <- filter(heatmap_df, heatmap_df$genes %in% glycolysis)
setcolorder(heatmap_df,cells$cell)


heatmap_df <- column_to_rownames(heatmap_df, var="genes")
heatmap_df <- as.matrix(heatmap_df)

annotation_df <- as.data.frame(human_combined$CT)
annotation_df$`human_combined$CT` <- factor(annotation_df$`human_combined$CT`, levels=order)

cells <- human_combined@meta.data
cells <- select(cells, CT)
cells <- rownames_to_column(cells, "cell")
cells$CT <- factor(cells$CT,levels=order)

heatmap_df <- rownames_to_column(heatmap_df, "genes")
heatmap_df <- as.data.table(heatmap_df)

setcolorder(heatmap_df,cells$cell)

levels(annotation_df) <- c("Zy_CS1", "4cell_CS2", "8cell_CS2", 
                           "cMor_CS3", 
                           "early_ICM_CS3", "late_Epi_CS3", 
                           "early_Tb_CS3", "late_Tb_CS3", 
                           "late_PrE_CS3", "ICM_CS3", "EPI_CS3", "EPI_CS4", "EmDisc_CS5", 
                           "EmDisc_CS6", "Tb_CS4", "Tb_CS5", "Tb_CS6",
                           "HYP_CS3", "VE_CS4", "VE_CS5",
                           "SYS_CS6")
order <- c("Zy_CS1", "4cell_CS2", "8cell_CS2", 
           "cMor_CS3", 
           "early_ICM_CS3", "late_Epi_CS3", 
           "early_Tb_CS3", "late_Tb_CS3", 
           "late_PrE_CS3", "ICM_CS3", "EPI_CS3", "EPI_CS4", "EmDisc_CS5", 
           "EmDisc_CS6", "Tb_CS4", "Tb_CS5", "Tb_CS6",
           "HYP_CS3", "VE_CS4", "VE_CS5",
           "SYS_CS6")


pheatmap(mat = heatmap_df, show_colnames= FALSE, 
         cluster_rows = FALSE, cluster_cols = FALSE, 
         annotation_col = annotation_df, kmeans_k = NA, 
        color = pal)

library(RColorBrewer)
library(viridis)

pal <- brewer.pal("RdBu", 10)

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(heatmap_df, n = 10)


###
D1 <- as.data.frame(Subs1[["umap"]]@cell.embeddings)
D2 <- as.data.frame(Subs1[["pca"]]@cell.embeddings)
###
##################################################################################
cType <-  c("PGC_4W_M","PGC_4W_F", "PGC_7W_M",
              "PGC_8W_F", "PGC_10W_M", "PGC_10W_F", 
              "PGC_11W_M", "PGC_11W_F", "PGC_17W_F", 
              "PGC_19W_M", "Soma_4W_F", "Soma_7W_M", 
              "Soma_8W_F","Soma_10W_M", "Soma_10W_F",
              "Soma_11W_M", "Soma_11W_F", "Soma_17W_F", 
              "Soma_19W_M")

BaseCol <-c("#b3ffd9"," #ffb3ec", "#80ffbf", "#ff80df", "#1aff8c", '#ff4dd2', '#00cc66', 
            "")

colind <- integer( length( levels(Idents(human_combined)) )  )
for (i in 1:length( levels(Idents(human_combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(human_combined))[i])
}
coluse <- BaseCol[colind]

#####################################################################################
#human naive and primed

naive <- read.csv("PAPER/Datasets/Human/counts_naive.csv", header=T, row.names=1, sep=",")
primed <- read.csv("PAPER/Datasets/Human/counts_primed.csv", header=T, row.names=1, sep=",")

primed <- select(primed, -Length)
naive <- select(naive, -Length)

naive <- CreateSeuratObject(counts = naive, assay = "RNA",min.cells = 0, min.features = 0)

naive$CT <- "ESCs_naive"

primed <- CreateSeuratObject(counts = primed, assay = "RNA",min.cells = 0, min.features = 0)

primed$CT <- "ESCs_primed"

Idents(primed) <- primed@meta.data$CT
Idents(naive) <-naive@meta.data$CT

human_ESCs <- merge(naive, primed)

human_ESCs <- NormalizeData(human_ESCs, verbose = FALSE)
#Get variable genes (reduce number)
human_ESCs <- FindVariableFeatures(human_ESCs, selection.method = "vst", nfeatures = 20000)

#scaling 
human_ESCs <- ScaleData(human_ESCs, verbose = FALSE)

#Dimesionality reduction and clustering
human_ESCs <- RunPCA(human_ESCs, npcs = 30, verbose = FALSE) #,features = rownames(Markers))
DimPlot(human_ESCs, reduction = "pca", cols= coluse, pt.size= 2, dims=c(1,2), repel = TRUE, label = FALSE)

pct <- human_ESCs@reductions[["pca"]]@stdev/sum(human_ESCs@reductions[["pca"]]@stdev) *100


coluse <- c("#F9CBCC", "#F192AE")
ggsave("PAPER/PCAs/human_cells_pca.pdf", width=6, height=4)

human_ESCs <- FindVariableFeatures(human_ESCs, selection.method = "vst", nfeatures = 40000)

gene_list_human <- human_ESCs@assays[["RNA"]]@var.features

oxphos <- read.csv("GeneOnt/oxphos.csv")
oxphos <-  pull(oxphos,Column2)
oxphos <- str_to_upper(oxphos)
oxphos <-unique(oxphos)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_human)%>% 
  unique()

glycolysis <- read.csv("GeneOnt/glycolysis.csv")
glycolysis <-  pull(glycolysis,Column2)
glycolysis <- str_to_upper(glycolysis)
glycolysis <-unique(glycolysis)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_human)%>% 
  unique()

human_ESCs <- AddModuleScore(human_ESCs, features = oxphos, 
                                 name = "oxphos")
human_ESCs <- AddModuleScore(human_ESCs, features = glycolysis, 
                                 name = "glycolysis")

human_scores <- data.frame(human_ESCs@meta.data$oxphos1,human_ESCs@meta.data$glycolysis1,
                           human_ESCs@active.ident)

colnames(human_scores) <- c("oxphos","glycolysis", "CT")

human_scores <- rownames_to_column(human_scores,var="cell")

human_scores <- reshape2::melt(human_scores, measure.vars=c(2:3), variable.name="score")

colnames(human_scores) <- c("cell","CT", "score", "value")

write.csv(human_scores, "PAPER/Scores/scores_human_ESCs.csv")

ggplot(human_scores, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.25, aes(fill=score))+
  scale_fill_manual(values = c("#289aeb", "#eb2838"),
                    name="", labels = c("OxPhos", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("")+
  ylim(-0.32,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()

ggsave("PAPER/Scores/human_cells_scores.pdf", width=4, height=4)

oxphos <- sort(oxphos$.)
glycolysis <- sort(glycolysis$.)

avgexp_human_oxphos <- AverageExpression(human_ESCs, features=oxphos, use.scale=FALSE)
avgexp_human_glycolysis <- AverageExpression(human_ESCs, features=glycolysis, use.scale=FALSE)

avgexp_human_glycolysis$glycolysis <- "TRUE"
avgexp_human_glycolysis$oxphos <- "FALSE"

avgexp_human_oxphos$oxphos <- "TRUE"
avgexp_human_oxphos$glycolysis <- "FALSE"

avgexp_human_glycolysis <- as.data.frame(avgexp_human_glycolysis)
avgexp_human_oxphos <- as.data.frame(avgexp_human_oxphos)

avgexp <- rbind(avgexp_human_glycolysis, avgexp_human_oxphos)

write.csv(avgexp, "PAPER/avgexp_humanESCs.csv")

oxphos <- sort(oxphos$.)

glycolysis <- c("HK1", "HK2", "HK3", "GCK", "GPI", "PFKL", 
                "PFKM", "PFKP", "ALDOA", "ALDOB", "ALDOC", 
                "TPI1", "GAPDH", "GAPDHS", "PGK1", "ENO1", "ENO2", "ENO3", 
                "PKLR", "PKM", "ADPGK", "BPGM", "FOXK1", "FOXK2", 
                "HIF1A", "LDHA", "LIN28A", "LIN28B", "MYC", "PDK1",
                "PDK3", "PGAM1", "PGAM2", "PGM1", "PGM2L1", "PYGL")


DoHeatmap(human_ESCs, features=oxphos, group.colors=coluse)
ggsave("PAPER/Heatmaps/human_cells_2.pdf", width=8, height=12)

DoHeatmap(human_ESCs, features=glycolysis, group.colors=coluse)
ggsave("PAPER/Heatmaps/human_cells.pdf", width =8 , height = 6)
