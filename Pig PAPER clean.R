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

##################################################################################
#I ran addmodulescore separately for Ramos_Ibeas and Kong (from GEO) and 
#the PCA using the merged dataset from the dataset folder (one
#of the last sections of the script)

# Ramos_Ibeas 
key_pig <- read.table("PAPER/Datasets/Pig/Porcine_labels.csv",sep=",",header = T, row.names=1)

#counts_pig1 <-read.table("PAPER/Datasets/Pig/Porcine_TPMs_All.csv", sep=",",header = T, row.names=1)
counts_pig <- read.table("PAPER/Datasets/Pig/GSE112380_Pig_RawCounts.txt/GSE112380_Pig_RawCounts.txt", sep="", header=T, row.names=1)
gene_names <- read.csv("PAPER/Datasets/Pig/biomart.csv", sep=",", header=T)

counts_pig$ï..Gene_ID <- rownames(counts_pig)
#counts_pig$ï..Gene_ID
#gene_names$ï..Gene_ID

counts_pig <- left_join(counts_pig, gene_names, by="ï..Gene_ID")

#counts_pig$Gene_name
#write.csv(counts_pig, "PAPER/counts_pig.csv")
counts_pig <- read.csv("PAPER/counts_pig.csv", header=1, row.names=1)

counts_pig <- select(counts_pig, -"ï..Gene_ID")

#counts_pig$Gene_name <- make.unique(counts_pig$Gene_name)
#View(counts_pig)

#counts_pig$Gene_name

#counts_pig <- counts_pig[complete.cases(counts_pig[ ,221]),]

counts_pig <- subset(counts_pig, !duplicated(Gene_name))
rownames(counts_pig) <- counts_pig$Gene_name
counts_pig <- select(counts_pig, -Gene_name)


labs <- key_pig$label

###### Seurat 
#Create a Seurat object
pig_data <- CreateSeuratObject(counts = counts_pig, assay = "RNA",min.cells = 0, min.features = 0)

#annotate
pig_data$CT <- labs
Idents(pig_data) <- pig_data@meta.data$CT #Change the annotation of each cell to e.g. its type Amnnion etc

pig_data$dataset <- "Ramos-Ibeas"

#Normalise the data
pig_data <- NormalizeData(pig_data, verbose = FALSE)

#Get variable genes (reduce number)
pig_data <- FindVariableFeatures(pig_data, selection.method = "vst", nfeatures = 40000)

gene_list_pig <- pig_data@assays[["RNA"]]@var.features

pig_data <- ScaleData(pig_data, verbose = FALSE)

pig_data <- RunPCA(pig_data, npcs = 30, verbose = FALSE) #,features = rownames(Markers))

levels(pig_data) <- c("cMor_CS3","ICM_CS3","Epi_CS3",
                          "EmDisc_CS5","Tb_CS3","Hyp_CS3","VE_CS5")

cType <- c("cMor_CS3","ICM_CS3","Epi_CS3",
           "EmDisc_CS5","Tb_CS3","Hyp_CS3","VE_CS5")

BaseCol <-c("#5AC0CE","#5B859E","#2AB6B9", 
            "#3F94D1","#F7EF8E", "#E5B405","#E94E12")

colind <- integer( length( levels(Idents(pig_data)) )  )
for (i in 1:length( levels(Idents(pig_data)) ) ) {
  colind[i] <- which(cType==levels(Idents(pig_data))[i])
}
coluse <- BaseCol[colind]

DimPlot(pig_data, dims=c(1,2), reduction="pca",  pt.size= 2,label= TRUE,repel=TRUE, shape.by = "dataset")

#D1 <- Embeddings(pig_data, reduction = "pca")
#D1 <- as.data.frame(pig_data[["pca"]]@cell.embeddings)

gene_list_pig <- pig_data@assays[["RNA"]]@var.features

oxphos <- read.csv("GeneOnt/oxphos.csv")
oxphos <-  pull(oxphos,Column2)
oxphos <- str_to_upper(oxphos)
oxphos <-unique(oxphos)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_pig) %>% 
  unique()

glycolysis <- read.csv("GeneOnt/glycolysis.csv")
glycolysis <-  pull(glycolysis,Column2)
glycolysis <- str_to_upper(glycolysis)
glycolysis <-unique(glycolysis)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_pig)%>% 
  unique()

#glycolysis <- c("HK1", "HK2", "HK3", "GCK", "GPI", "PFKL", 
                "PFKM", "PFKP", "ALDOA", "ALDOB", "ALDOC", 
                "TPI1", "GAPDH", "GAPDHS", "PGK1", "ENO1", "ENO2", "ENO3", 
                "PKLR", "PKM", "ADPGK", "BPGM", "FOXK1", "FOXK2", 
                "HIF1A", "LDHA", "LIN28A", "LIN28B", "MYC", "PDK1",
                "PDK3", "PGAM1", "PGAM2", "PGM1", "PGM2L1", "PYGL")

oxphos <- sort(oxphos$.)

DoHeatmap(pig_data, features = oxphos, group.colors = coluse)
ggsave("PAPER/Heatmaps/pig_ramos_oxphos.pdf", width=8, height=12)

DoHeatmap(pig_data, features = glycolysis, group.colors = coluse)
ggsave("PAPER/Heatmaps/pig_combined_glycolysis.pdf", width=8, height=6)

####################################################################################
#### zygote to 8cell, Kong et al 2020, FASEB
#Original dataset 
pig_counts2 <- read.table("PAPER/Datasets/Pig/fsb220038-sup-0002-tables1.csv", sep=",", header=T, row.names=1)
pig_counts2 <- subset(pig_counts2, !duplicated(NAME))
rownames(pig_counts2) <- NULL
pig_counts2<- column_to_rownames(pig_counts2, var="NAME")

View(pig_counts2)

pig_key2 <- read.table("PAPER/Datasets/Pig/Kong_key.csv", sep=",", header=T, row.names=1)
labs <- pig_key2$CT

#create a Seurat object
pig_data2 <-  CreateSeuratObject(counts = pig_counts2, assay = "RNA",min.cells = 0, min.features = 0)
pig_data2$CT <- labs
pig_data2$dataset <- "Kong"
Idents(pig_data2) <- pig_data2@meta.data$CT 

pig_data2<-subset(pig_data2, idents = c("Zy_CS1", "4cell_CS2", "8cell_CS2"))

pig_data2 <- NormalizeData(pig_data2, verbose = FALSE)

pig_data2 <- FindVariableFeatures(pig_data2, selection.method = "vst", nfeatures = 40000)

pig_data2 <- ScaleData(pig_data2, verbose = FALSE)

levels(pig_data2) <- c("Zy_CS1","4cell_CS2","8cell_CS2")

cType <- c("Zy_CS1","4cell_CS2","8cell_CS2")

BaseCol <-c("#DAEDEC","#A6D9E7","#BBDFD5")

colind <- integer( length( levels(Idents(pig_data2)) )  )
for (i in 1:length( levels(Idents(pig_data2)) ) ) {
  colind[i] <- which(cType==levels(Idents(pig_data2))[i])
}
coluse <- BaseCol[colind]

gene_list_pig <- pig_data2@assays[["RNA"]]@var.features

oxphos <- sort(oxphos$.)

oxphos <- read.csv("GeneOnt/oxphos.csv")
oxphos <-  pull(oxphos,Column2)
oxphos <- str_to_upper(oxphos)
oxphos <-unique(oxphos)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_pig) %>% 
  unique()

glycolysis <- read.csv("GeneOnt/glycolysis.csv")
glycolysis <-  pull(glycolysis,Column2)
glycolysis <- str_to_upper(glycolysis)
glycolysis <-unique(glycolysis)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_pig)%>% 
  unique()

DoHeatmap(pig_data2, features=glycolysis, group.colors = coluse)
ggsave("PAPER/Heatmaps/kong_glyco.pdf", width=8, height=6)

DoHeatmap(pig_data2, features=oxphos, group.colors = coluse)
ggsave("PAPER/Heatmaps/kong_oxphos.pdf", width=8, height=6)



#pig_data2 <- RunPCA(pig_data2, npcs = 30, verbose = FALSE) #,features = rownames(Markers))

#pig_data2 <- subset(pig_data2, cells= colnames(x = pbmc_small)[3:43])

#DimPlot(pig_data2, dims=c(1,2), reduction="pca",  pt.size= 2,label= TRUE,repel=TRUE, shape.by = "dataset" , cols=coluse)

#D2 <- as.data.frame(pig_data2[["pca"]]@cell.embeddings)
#D2<- D2[row.names(D2) != "Zygote.2",]

#PCAemb <-rbind(D1, D2)
#plot(PCAemb, x=PC_1, y=PC_2, main="Scatterplot Example",
#     xlab="", ylab="", pch=19)

#PCAemb <- cbind(PCAemb, idents)

#idents<- idents[row.names(idents) != "Zygote.2",]

#filt <- c("Zy_CS1","4cell_CS2","8cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3", 
  "EmDisc_CS5","Hyp_CS3","VE_CS5","Tb_CS3")

#idents_fil <- PCAemb$ident.idents

#ggplot(PCAemb, aes(x=PC_1, y=PC_2, color = idents))+
#  geom_point()+
#  scale_colour_manual()



#########################################################################################
### merge 
pig_combined <- merge(pig_data2, pig_data, merge.data = TRUE)
idents <- pig_combined@active.ident
idents <- as.data.frame(idents)

#Normalise the data
pig_combined <- NormalizeData(pig_combined, verbose = FALSE)

#Get variable genes (reduce number)
pig_combined <- FindVariableFeatures(pig_combined, selection.method = "vst", nfeatures = 20000)

#gene_list_pig <- pig_combined@assays[["RNA"]]@var.features

pig_combined <- ScaleData(pig_combined, verbose = FALSE)

pig_combined <- RunPCA(pig_combined, npcs = 30, verbose = FALSE) #,features = rownames(Markers))

D2 <- as.data.frame(pig_combined[["pca"]]@cell.embeddings)


###################################################################################
### colours 

levels(pig_combined) <- c("Zy_CS1","4cell_CS2","8cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3",
                          "EmDisc_CS5","Hyp_CS3","VE_CS5","Tb_CS3")

cType <- c("Zy_CS1","4cell_CS2","8cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3",
           "EmDisc_CS5","Hyp_CS3","VE_CS5","Tb_CS3")

BaseCol <-c("#71eac4","#A6D9E7","#BBDFD5","#5AC0CE","#5B859E","#2AB6B9", 
            "#3F94D1","#E5B405","#E94E12","#F7EF8E")

colind <- integer( length( levels(Idents(pig_combined)) )  )
for (i in 1:length( levels(Idents(pig_combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(pig_combined))[i])
  }
coluse <- BaseCol[colind]

plot <- DimPlot(pig_combined, dims=c(1,2), reduction="pca",  pt.size= 2,label= TRUE,repel=TRUE, shape.by = "dataset" , cols=coluse)
plot[[1]] + xlim(c(-150,50)) +ylim(c(-100, 70))

ggsave("PAPER/PCAs/Dimplot_pig_legend_ptsize15.pdf",width =8, height=6)

####################################################################################
#Scores for

pig_combined <- pig_data2

pig_combined <- FindVariableFeatures(pig_combined, selection.method = "vst", nfeatures = 40000)

gene_list_pig <- pig_combined@assays[["RNA"]]@var.features

sort(oxphos$.)

oxphos <- read.csv("GeneOnt/oxphos.csv")
oxphos <-  pull(oxphos,Column2)
oxphos <- str_to_upper(oxphos)
oxphos <-unique(oxphos)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_pig) %>% 
  unique()

glycolysis <- read.csv("GeneOnt/glycolysis.csv")
glycolysis <-  pull(glycolysis,Column2)
glycolysis <- str_to_upper(glycolysis)
glycolysis <-unique(glycolysis)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_pig)%>% 
  unique()

write.csv(oxphos, "PAPER/Gene modules/oxphos_pig.csv")
write.csv(glycolysis, "PAPER/Gene modules/glycolysis_pig.csv")

##### Add Module Scores 
#
pig_data <- AddModuleScore(pig_data, features = glycolysis, 
                               name = "glycolysis")
pig_data <- AddModuleScore(pig_data, features = oxphos, 
                               name = "oxphos")

levels(pig_combined)

scores_pig <- data.frame(pig_data@meta.data$oxphos1,pig_data@meta.data$glycolysis1, 
                         pig_data@active.ident)

scores_pig <- rownames_to_column(scores_pig,var="cell")

scores_pig <- reshape2::melt(scores_pig, measure.vars=c(2:3), variable.name="score")

colnames(scores_pig) <- c("cell","CT", "score", "value")

#write.csv(scores_pig, "PAPER/Scores/Pig_scores_combined.csv")

scores_pig <- filter(scores_pig,CT %in% c("Zy_CS1", "4cell_CS2", "8cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3",
                                                 "EmDisc_CS5"))

scores_pig$CT <- factor(scores_pig$CT , levels=c("Zy_CS1", "4cell_CS2", "8cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3",
                                                 "EmDisc_CS5"))

scores_pig <- filter(scores_pig,CT %in% c("Tb_CS3", "Hyp_CS3", "VE_CS5"))

scores_pig$CT <- factor(scores_pig$CT , levels=c("Tb_CS3", "Hyp_CS3", "VE_CS5"))


scores_pig <- read.csv("PAPER/Scores/pig_scores_late_separate.csv")

levels(scores_pig)<- c("cMor_CS3", "ICM_CS3","Epi_CS3", "EmDisc_CS5","Tb_CS3", "Hyp_CS3", "VE_CS5")

ggplot(scores_pig, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.25, aes(fill=score))+
  scale_fill_manual(values = c("#289aeb", "#eb2838"),
                    name="", labels = c("OxPhos", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("")+
  ylim(-0.32,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()

ggsave("PAPER/Scores/Pig_new_extra.pdf", width=6,height=4)
write.csv(scores_pig,"PAPER/Scores/Pig_FINAL_RAMOS_Scores.csv")


meta<- pig_data@meta.data
###################################################################################
#Average Expression for Supplementary 
oxphos <- sort(oxphos$.)
glycolysis <- sort(glycolysis$.)

levels(pig_data) <- c("cMor_CS3","ICM_CS3","Epi_CS3",
                      "EmDisc_CS5","Tb_CS3", "Hyp_CS3", "VE_CS5")

avgexp_pig_oxphos <- AverageExpression(pig_combined, features=oxphos, use.scale=FALSE)
avgexp_pig_glycolysis <- AverageExpression(pig_combined, features=glycolysis, use.scale=FALSE)

avgexp_pig_glycolysis$glycolysis <- "TRUE"
avgexp_pig_glycolysis$oxphos <- "FALSE"

avgexp_pig_oxphos$oxphos <- "TRUE"
avgexp_pig_oxphos$glycolysis <- "FALSE"

avgexp_pig_glycolysis <- as.data.frame(avgexp_pig_glycolysis)
avgexp_pig_oxphos <- as.data.frame(avgexp_pig_oxphos)

avgexp <- rbind(avgexp_pig_glycolysis, avgexp_pig_oxphos)

write.csv(avgexp, "PAPER/avgexp_pig_kong.csv")

meta <- pig_combined@meta.data


###################################################################################
####
#Heatmap Seurat
oxphos <- sort(oxphos$.)
glycolysis <- sort(glycolysis$.)

cType <- c("Zy_CS1","4cell_CS2","8cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3",
           "EmDisc_CS5","Hyp_CS3","VE_CS5","Tb_CS3")

BaseCol <-c("#DAEDEC","#A6D9E7","#BBDFD5","#5AC0CE","#5B859E","#17bc80", 
            "#3F94D1","#F7EF8E", "#E5B405", "#E94E12")

levels(pig_data) <- c("cMor_CS3","ICM_CS3","Epi_CS3",
                      "EmDisc_CS5","Hyp_CS3","VE_CS5","Tb_CS3")

cols <-c("#5AC0CE","#5B859E","#2AB6B9", 
         "#3F94D1","#E5B405","#E94E12","#F7EF8E")

DoHeatmap(pig_combined, features = oxphos, group.colors = coluse)
ggsave("PAPER/Heatmaps/pig_combined_oxphos.pdf", width=8, height=10)

DoHeatmap(pig_combined, features = glycolysis, group.colors = coluse)
ggsave("PAPER/Heatmaps/pig_combined_glycolysis.pdf", width=8, height=6)


####################################################################################
#dataset from the 5 species folder - Ramos-Ibeas and Kong merged 
#gives weird scores

key_pig <- read.table ("PAPER/Datasets/Pig/new_pig/PigKey.csv", header=T, row.names=1, sep=",")
counts_pig <- read.table ("PAPER/Datasets/Pig/new_pig/PigData.csv", header=T, row.names=1, sep=",")
labs <- key_pig$Type

#counts_pig <- select(counts_pig, -Chr, -Start, -End, -Strand, -Length)
pig_data <- CreateSeuratObject(counts = counts_pig, 
                               assay = "RNA",min.cells = 0, min.features = 0)
#annotate
pig_data$CT <- labs
Idents(pig_data) <- pig_data@meta.data$CT #Change the annotation of each cell to e.g. its type Amnnion etc
pig_data <- subset(pig_data, idents = c("Zy_CS1", "4-cell_CS2", "8-cell_CS2", "cMor_CS3", 
                                        "ICM_CS3", "Epi_CS3", "EmDisc_CS5", "Tb_CS3", "Hyp_CS3", "Hyp_CS5"))


#Normalise the data
pig_data <- NormalizeData(pig_data, verbose = FALSE)

pig_data <- FindVariableFeatures(pig_data, selection.method = "vst", nfeatures = 20000)

gene_list_pig <- pig_data@assays[["RNA"]]@var.features

#scaling 
pig_data <- ScaleData(pig_data, verbose = FALSE)

# Scores
oxphos <- read.csv("GeneOnt/oxphos.csv")
oxphos <-  pull(oxphos,Column2)
oxphos <- str_to_upper(oxphos)
oxphos <-unique(oxphos)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_pig) %>% 
  unique()

glycolysis <- read.csv("GeneOnt/glycolysis.csv")
glycolysis <-  pull(glycolysis,Column2)
glycolysis <- str_to_upper(glycolysis)
glycolysis <-unique(glycolysis)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_pig)%>% 
  unique()

##### Scores 
#
pig_data <- AddModuleScore(pig_data, features = glycolysis, 
                           name = "glycolysis")
pig_data <- AddModuleScore(pig_data, features = oxphos, 
                           name = "oxphos")

scores_pig <- data.frame(pig_data@meta.data$oxphos1,pig_data@meta.data$glycolysis1, 
                         pig_data@active.ident)
scores_pig <- rownames_to_column(scores_pig,var="cell")

scores_pig <- reshape2::melt(scores_pig, measure.vars=c(2:3), variable.name="score")

colnames(scores_pig) <- c("cell","CT", "score", "value")

#
pig_data <- subset(pig_data, idents= c("Zy_CS1", "4-cell_CS2", "8-cell_CS2", "cMor_CS3", 
                                             "ICM_CS3", "Epi_CS3", "EmDisc_CS5", "Tb_CS3", "Hyp_CS3"))

###
scores_pig$CT <- factor(scores_pig$CT, 
                        levels=c("Zy_CS1", "4-cell_CS2", "8-cell_CS2", "cMor_CS3", 
                                 "ICM_CS3", "Epi_CS3", "EmDisc_CS5", "Tb_CS3", "Hyp_CS3", "VE_CS5"))

ggplot(scores_pig, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.5, aes(fill=score))+
  scale_fill_manual(values = c("#6280e3", "#b83939"),
                    name="", labels = c("OXPHOS", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("Embryonic")+
  ylim(-0.5,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()
ggsave("PAPER/Scores/scores_wrong_combined_dataset.pdf", width=6,height=4)

###############################################################################################
pig_data <- FindVariableFeatures(pig_data, selection.method = "vst", nfeatures = 20000)
pig_data <- ScaleData(pig_data, verbose = FALSE)

meta_data <- read.table("PAPER/Datasets/Pig/meta_data.csv",sep=",",header = T, row.names=1)

pig_data$dataset <- meta_data$dataset

#meta_data <- pig_data@meta.data
#write.csv(meta_data, "PAPER/Datasets/Pig/meta_data.csv")

pig_data <- RunPCA(pig_data, npcs = 30, verbose = FALSE) #,features = rownames(Markers))
DimPlot(pig_data, dims=c(1,2), reduction="pca", pt.size= 3, label= TRUE, repel=TRUE, cols=BaseCol, shape.by = "dataset")

pct <- pig_data@reductions[["pca"]]@stdev/sum(pig_data@reductions[["pca"]]@stdev) *100



ggsave("PAPER/PCAs/Pig_final.pdf", width=6, height=4)


pig_data <- subset(pig_data, idents = c("Zy_CS1", "4-cell_CS2", "8-cell_CS2", "cMor_CS3", 
                                        "ICM_CS3", "Epi_CS3", "EmDisc_CS5", "Tb_CS3", "Hyp_CS3", "Hyp_CS5"))

levels(pig_data) <-c("Zy_CS1", "4-cell_CS2", "8-cell_CS2", "cMor_CS3", 
                     "ICM_CS3", "Epi_CS3", "EmDisc_CS5", "Tb_CS3", "Hyp_CS3", "Hyp_CS5")
