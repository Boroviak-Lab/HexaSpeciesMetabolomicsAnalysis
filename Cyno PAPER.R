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

Key_cyn<- read.table("PAPER/Datasets/Cyno/cyKey.csv",sep=",",header = T, row.names=1)
counts_cyn <-read.table("PAPER/Datasets/Cyno/cyData.csv", sep=",",header = T, row.names=1)

isinvit <- Key_cyn$All 
labs <- Key_cyn$LABEL_4
labs <- labs[which(isinvit>0)]


# subset(cyno_data, idents = c("), invert=FALSE) 

#Create a Seurat object
cyno_data <- CreateSeuratObject(counts = counts_cyn[which(isinvit>0)], assay = "RNA",min.cells = 0, min.features = 0)
cyno_data$CT <- labs

Idents(cyno_data) <- cyno_data@meta.data$CT #Change the annotation of each cell to e.g. its type Amnnion etc

cyno_data <- subset(cyno_data, idents= c("ICM_CS3","Epi_CS3","EmDisc_CS5", 
                                         "EmDisc_CS6/7",
                                         "Hyp_CS3","VE_CS5",
                                         "SYS_CS6/7",
                                         "Tb_CS3",
                                         "Tb_CS6/7",
                                         "ExMes_CS5","ExMes_CS6/7",
                                         "PGC_CS5", "PGC_CS6/7",
                                         "ESC_primed"))

#Normalise the data
cyno_data <- NormalizeData(cyno_data, verbose = FALSE)
#Get variable genes (reduce number)
cyno_data <- FindVariableFeatures(cyno_data, selection.method = "vst", nfeatures = 20000)

gene_list_cyn <- cyno_data@assays[["RNA"]]@var.features

#scaling 
cyno_data <- ScaleData(cyno_data, verbose = FALSE)

#Idents(cyno_data) <- cyno_data@meta.data$CT #Change the annotation of each cell to e.g. its type Amnnion etc

#Dimesionality reduction and clustering
cyno_data <- RunPCA(cyno_data, npcs = 30, verbose = FALSE) #,features = rownames(Markers))

pct <- cyno_data@reductions[["pca"]]@stdev/sum(cyno_data@reductions[["pca"]]@stdev) *100

levels(cyno_data) <- c("ICM_CS3","Epi_CS3","EmDisc_CS5","EmDisc_CS6/7",
                                                                "Hyp_CS3","VE_CS5",
                                                                "SYS_CS6/7",
                                                                "Tb_CS3",
                                                                "Tb_CS6/7",
                                                                "ExMes_CS5","ExMes_CS6/7",
                                                                "PGC_CS5", "PGC_CS6/7",
                                                                "ESC_primed")

cType <- c("ICM_CS3","Epi_CS3","EmDisc_CS5", 
           "EmDisc_CS6/7")

BaseCol <-c("#5AC0CE","#5B859E","#2AB6B9",
            "#284597"
          )

colind <- integer( length( levels(Idents(cyno_data)) )  )
for (i in 1:length( levels(Idents(cyno_data)) ) ) {
  colind[i] <- which(cType==levels(Idents(cyno_data))[i])
}
coluse <- BaseCol[colind]

cyno_data$dataset <- "Nakamura"

DimPlot(cyno_data, reduction = "pca", dims=c(1,2), cols=coluse, pt.size= 2, shape.by="dataset", repel = TRUE, label = TRUE )
ggsave("PAPER/PCAs/Dimcyno.pdf", width=6, height=4)


##### Mammalian lists 
#gene_list <- cyno_data@assays[["RNA"]]@var.features

cyno_data <- FindVariableFeatures(cyno_data, selection.method = "vst", nfeatures = 30000)

gene_list_cyn <- cyno_data@assays[["RNA"]]@var.features

oxphos <- read.csv("GeneOnt/oxphos.csv")
oxphos <-  pull(oxphos,Column2)
oxphos <- str_to_upper(oxphos)
oxphos <-unique(oxphos)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_cyn)%>% 
  unique()

glycolysis <- read.csv("GeneOnt/glycolysis.csv")
glycolysis <-  pull(glycolysis,Column2)
glycolysis <- str_to_upper(glycolysis)
glycolysis <-unique(glycolysis)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_cyn)%>% 
  unique()

write.csv(oxphos,"PAPER/Scores/oxphos_nakamura.csv")
write.csv(glycolysis, "PAPER/Scores/glycolysis_nakamura.csv")

cyno_data <- AddModuleScore(cyno_data, features = oxphos, 
                            name = "oxphos")
cyno_data <- AddModuleScore(cyno_data, features = glycolysis, 
                            name = "glycolysis")

scores_cyn <- data.frame(cyno_data@meta.data$oxphos1,cyno_data@meta.data$glycolysis1,
                     cyno_data@active.ident)

scores_cyn <- rownames_to_column(scores_cyn,var="cell")

scores_cyn <- reshape2::melt(scores_cyn, measure.vars=c(2:3), variable.name="score")

colnames(scores_cyn) <- c("cell","CT", "score", "value")

emb <- c("ICM_CS3", "Epi_CS3", "EmDisc_CS5", "EmDisc_CS6/7")

scores_cyn_emb <- filter(scores_cyn, CT %in% emb)


ggplot(scores_cyn, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.25, aes(fill=score))+
  scale_fill_manual(values = c("#289aeb", "#eb2838"),
                    name="", labels = c("OxPhos", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("")+
  ylim(-0.32,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()

ggsave("PAPER/Scores/scores_cyn_new_emb_nakamura_all.pdf", width=6, height=4)

write.csv(scores_cyn, "PAPER/Scores/scores_cyn_invivo_new_emb.csv")

#################################################################################
oxphos <- sort(oxphos$.)
glycolysis <- sort(glycolysis$.)

avgexp_cyn_oxphos <- AverageExpression(cyno_data, features=oxphos, use.scale=FALSE)
avgexp_cyn_glycolysis <- AverageExpression(cyno_data, features=glycolysis, use.scale=FALSE)

avgexp_cyn_glycolysis$glycolysis <- "TRUE"
avgexp_cyn_glycolysis$oxphos <- "FALSE"

avgexp_cyn_oxphos$oxphos <- "TRUE"
avgexp_cyn_oxphos$glycolysis <- "FALSE"

avgexp_cyn_glycolysis <- as.data.frame(avgexp_cyn_glycolysis)
avgexp_cyn_oxphos <- as.data.frame(avgexp_cyn_oxphos)

avgexp <- rbind(avgexp_cyn_glycolysis, avgexp_cyn_oxphos)

write.csv(avgexp, "PAPER/avgexp_cyno.csv")

###
emb <- c("ICM_CS3","Epi_CS3","EmDisc_CS5","EmDisc_CS6/7")

extraemb <- c("Hyp_CS3","VE_CS5","SYS_CS6/7", "Tb_CS3",
              "Tb_CS6/7","ExMes_CS5","ExMes_CS6/7")
pgc <- c("PGC_CS5","PGC_CS6/7")

cells <- c("ESC_primed")

scores_cyn_emb <- scores_cyn %>% 
  filter(CT %in% emb)

scores_cyn_emb$CT <- factor(scores_cyn_emb$CT , levels=c("ICM_CS3", "Epi_CS3", "EmDisc_CS5", "EmDisc_CS6/7"))

ggplot(scores_cyn_emb, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.5, aes(fill=score))+
  scale_fill_manual(values = c("#6280e3", "#b83939"),
                    name="", labels = c("OXPHOS", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("Embryonic")+
  ylim(-0.5,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()
ggsave("PAPER/Scores/cyn_embryo.pdf", width=6,height=4)

#### extraembryonic 
scores_cyn_extraemb <- scores_cyn %>% 
  filter(CT %in% extraemb)

scores_cyn_extraemb$CT <- factor(scores_cyn_extraemb$CT , levels=c("Hyp_CS3", "VE_CS5", "Tb_CS3","Tb_CS6/7", "SYS_CS6/7", "ExMes_CS5","ExMes_CS6/7"))

ggplot(scores_cyn_extraemb, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.5, aes(fill=score))+
  scale_fill_manual(values = c("#6280e3", "#b83939"),
                    name="", labels = c("OXPHOS", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("Extraembryonic")+
  ylim(-0.5,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()
ggsave("PAPER/Scores/cyn_extraembryo.pdf", width=6,height=4)

### pgc
scores_cyn_pgc <- scores_cyn %>% 
  filter(CT %in% pgc)

scores_cyn_pgc$CT <- factor(scores_cyn_pgc$CT , levels=c("PGC_CS5", "PGC_CS6/7"))

ggplot(scores_cyn_pgc, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.5, aes(fill=score))+
  scale_fill_manual(values = c("#6280e3", "#b83939"),
                    name="", labels = c("OXPHOS", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("PGCs")+
  ylim(-0.5,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()
ggsave("PAPER/Scores/cyn_pgc.pdf", width=4,height=4)

### cells 
scores_cyn_cells <- scores_cyn %>% 
  filter(CT %in% cells)

ggplot(scores_cyn_cells, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.5, aes(fill=score))+
  scale_fill_manual(values = c("#6280e3", "#b83939"),
                    name="", labels = c("OXPHOS", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("Cell culture")+
  ylim(-0.5,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()
ggsave("PAPER/Scores/cyn_cells.pdf", width=4,height=4)

#################################################################################################
###CYNO DATA IN VITRO 
##

cyn_counts2 <- read.table("PAPER/Datasets/Cyno/cyInVitData.csv", sep=",", header=TRUE, row.names=1)

cyn_key2 <- read.table("PAPER/Datasets/Cyno/cyinvitLab.csv", sep= ",", header=TRUE, row.names=1)

labs <- select(cyn_key2, Type)
  
  
# subset(cyno_data, idents = c("), invert=FALSE) 

#Create a Seurat object
cyno_data2 <- CreateSeuratObject(counts = cyn_counts2, assay = "RNA",min.cells = 0, min.features = 0)
cyno_data2$CT <- labs

Idents(cyno_data2) <- cyno_data2@meta.data$CT #Change the annotation of each cell to e.g. its type Amnnion etc

#Normalise the data
cyno_data2 <- NormalizeData(cyno_data2, verbose = FALSE)
#Get variable genes (reduce number)
cyno_data2 <- FindVariableFeatures(cyno_data2, selection.method = "vst", nfeatures = 20000)

gene_list_cyn2 <- cyno_data2@assays[["RNA"]]@var.features

#scaling 
cyno_data2 <- ScaleData(cyno_data2, verbose = FALSE)

#Idents(cyno_data2) <- cyno_data2@meta.data$CT #Change the annotation of each cell to e.g. its type Amnnion etc

#Dimesionality reduction and clustering
cyno_data2 <- RunPCA(cyno_data2, npcs = 30, verbose = FALSE) #,features = rownames(Markers))
DimPlot(cyno_data2, reduction = "pca", dims=c(1,2), repel = TRUE, label = FALSE, cols=BaseCol)
ggsave("PAPER/Ma_cyno.pdf", width=6, height=4)


pct <- cyno_data2@reductions[["pca"]]@stdev/sum(cyno_data2@reductions[["pca"]]@stdev) *100



##### Mammalian lists 
#gene_list <- cyno_data2@assays[["RNA"]]@var.features

oxphos <- read.csv("GeneOnt/oxphos.csv")
oxphos <-  pull(oxphos,Column2)
oxphos <- str_to_upper(oxphos)
oxphos <-unique(oxphos)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_cyn2)%>% 
  unique()

glycolysis <- read.csv("GeneOnt/glycolysis.csv")
glycolysis <-  pull(glycolysis,Column2)
glycolysis <- str_to_upper(glycolysis)
glycolysis <-unique(glycolysis)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_cyn2)%>% 
  unique()

cyno_data2 <- AddModuleScore(cyno_data2, features = oxphos, 
                            name = "oxphos")
cyno_data2 <- AddModuleScore(cyno_data2, features = glycolysis, 
                            name = "glycolysis")

scores_cyn2 <- data.frame(cyno_data2@meta.data$oxphos1,cyno_data2@meta.data$glycolysis1,
                         cyno_data2@active.ident)

scores_cyn2 <- rownames_to_column(scores_cyn2,var="cell")

scores_cyn2 <- reshape2::melt(scores_cyn2, measure.vars=c(2:3), variable.name="score")

colnames(scores_cyn2) <- c("cell","CT", "score", "value")

write.csv(scores_cyn2, "PAPER/Scores/scores_cyn_invitro.csv")

#######

oxphos <- sort(oxphos$.)
glycolysis <- sort(glycolysis$.)

avgexp_cyn_oxphos <- AverageExpression(cyno_data2, features=oxphos, use.scale=FALSE)
avgexp_cyn_glycolysis <- AverageExpression(cyno_data2, features=glycolysis, use.scale=FALSE)

avgexp_cyn_glycolysis$glycolysis <- "TRUE"
avgexp_cyn_glycolysis$oxphos <- "FALSE"

avgexp_cyn_oxphos$oxphos <- "TRUE"
avgexp_cyn_oxphos$glycolysis <- "FALSE"

avgexp_cyn_glycolysis <- as.data.frame(avgexp_cyn_glycolysis)
avgexp_cyn_oxphos <- as.data.frame(avgexp_cyn_oxphos)

avgexp <- rbind(avgexp_cyn_glycolysis, avgexp_cyn_oxphos)

write.csv(avgexp, "PAPER/avgexp_cyno_invitro.csv")
#Plotting

emb <- c("EmDisc_CS5","EmDisc_CS6","EmDisc_CS7")

extraemb <- c("VE_CS5","Tb_CS5", "Tb_CS6", "Tb_CS7","Am_CS5", 
              "Am_CS6","Am_CS7", "SYS_CS6","SYS_CS7",
              "ExMes_CS5","ExMes_CS6", "ExMes_CS7")
pgc <- c("PGC_CS6")

scores_cyn2_emb <- scores_cyn2 %>% 
  filter(CT %in% emb)

scores_cyn2_emb$CT <- factor(scores_cyn2_emb$CT , levels=c("EmDisc_CS5","EmDisc_CS6","EmDisc_CS7"))


ggsave("PAPER/Scores/cyn_invitro_embryo.pdf", width=4,height=4)

#### extraembryonic 
scores_cyn2 <- scores_cyn2 %>% 
  filter(CT %in% c("EmDisc_CS5", "EmDisc_CS6", "EmDisc_CS7", "Tb_CS5", "Tb_CS6", "Tb_CS7", 
                   "VE_CS5"))

scores_cyn2$CT <- factor(scores_cyn2$CT , levels=c("EmDisc_CS5", "EmDisc_CS6", "EmDisc_CS7", "VE_CS5","Tb_CS5", "Tb_CS6", "Tb_CS7","Am_CS5", 
                                                                     "Am_CS6","Am_CS7", "SYS_CS6","SYS_CS7",
                                                                     "ExMes_CS5","ExMes_CS6", "ExMes_CS7", "PGC_CS6"))

ggplot(scores_cyn2, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.25, aes(fill=score))+
  scale_fill_manual(values = c("#289aeb", "#eb2838"),
                    name="", labels = c("OxPhos", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("")+
  ylim(-0.32,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()

ggsave("PAPER/Scores/cyn_invitro_all_new.pdf", width=6,height=4)


View(cyno_data@meta.data)

### pgc
scores_cyn2_pgc <- scores_cyn2 %>% 
  filter(CT %in% pgc)

ggplot(scores_cyn2_pgc, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.5, aes(fill=score))+
  scale_fill_manual(values = c("#6280e3", "#b83939"),
                    name="", labels = c("OXPHOS", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("PGCs")+
  ylim(-0.5,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()
ggsave("PAPER/Scores/cyn_invitro_pgc.pdf", width=3,height=4)

#######################################################################################
#### MERGE for PCA
labs <- select(cyn_key2, Type)

cyno_data$CT <- labs
Idents(cyno_data) <- cyno_data2@meta.data$CT

cyno_combined <- merge(cyno_data, y = cyno_data2, add.cell.ids = c("1", "2"), project = "cyno_combined")

meta <- cyno_combined@meta.data

write.csv(meta, "PAPER/Datasets/Cyno/key_merged.csv")
#added a column with dataset name based on cell id to the key_merged.csv file
key_merged <- read.table("PAPER/Datasets/Cyno/key_merged.csv", sep=",", header = TRUE, row.names=1)

cyno_combined$dataset <- key_merged$data

cyno_combined <- FindVariableFeatures(cyno_combined, selection.method = "vst", nfeatures = 20000)

#scaling 
cyno_combined <- ScaleData(cyno_combined, verbose = FALSE)

#Idents(cyno_combined) <- cyno_combined@meta.data$CT #Change the annotation of each cell to e.g. its type Amnnion etc

meta<- select(key_merged,data)

#Dimesionality reduction and clustering
cyno_combined <- RunPCA(cyno_combined, npcs = 30, verbose = FALSE) #,features = rownames(Markers))

levels(cyno_combined) <- c("ICM_CS3","Epi_CS3","EmDisc_CS5", "EmDisc_CS6", "EmDisc_CS6/7", "EmDisc_CS7",
                          "Am_CS5","Am_CS6","Am_CS7",
                           "Hyp_CS3","VE_CS5","SYS_CS6","SYS_CS6/7","SYS_CS7", 
                          "Tb_CS3",
                           "Tb_CS5","Tb_CS6","Tb_CS6/7", "Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS6/7","ExMes_CS7",
                           "PGC_CS5","PGC_CS6", "PGC_CS6/7",
                          "ESC_primed")

#### 
colours 
cType <- c("EmDisc_CS5", 
           "EmDisc_CS6", "EmDisc_CS7",
           "Am_CS5","Am_CS6","Am_CS7",
           "VE_CS5",
           "SYS_CS6","SYS_CS7", 
           "Tb_CS5",
           "Tb_CS6", "Tb_CS7",
           "ExMes_CS5","ExMes_CS6","ExMes_CS7",
           "PGC_CS6")
#switched to invitro
BaseCol <-c("#2AB6B9",
            "#3F94D1", "#284597",
            "#F7EF8E", "#F5D950", "#E3DF1E",
            "#E94E12",
            "#8379B7", "#5554A0",
            "#2A2A64",
            "#E58603", "#D17612",
            "#D79E72", "#A35B2F","#521602",
            "#800080")

colind <- integer( length( levels(Idents(cyno_combined)) )  )
for (i in 1:length( levels(Idents(cyno_combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(cyno_combined))[i])
}
coluse <- BaseCol[colind]


DimPlot(cyno_combined, reduction = "pca", dims=c(1,2), repel = TRUE, 
        label = FALSE, shape.by = 'dataset', cols=coluse)

ggsave("PAPER/PCAs/Dimplot_cyno.pdf", width=10, height=6)

DoHeatmap(cyno_combined, group.colors = coluse, features=glycolysis)
ggsave("PAPER/Heatmaps/cyno_combined_colour_scheme.pdf", width = 12, height=6)

##################################################################################
#PCA 3D

# Export PC loadings into a dataframe organized by type
embed_mat <- as.data.frame(cyno_combined@reductions$pca@cell.embeddings)
meta_data <- Idents(cyno_combined)
stage <- cyno_combined$dataset

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
  if(as.character(pca3d_df$stage)[j]   == 'Nakamura') {pca3dpch[j] <- 19}
  if(as.character(pca3d_df$stage)[j]   == 'Ma') {pca3dpch[j] <- 17}
}

library("plot3D")

# Add small dots on basal plane and on the depth plane
# Can be adapted to plot on any plane
scatter3D_fancy <- function(x, y, z,..., colvar = z)
{
  panelfirst <- function(pmat) {
    XY <- trans3D(x, y, z = rep(min(z), length(z)), pmat = pmat)
    scatter2D(XY$x, XY$y, col= coluse, colvar = colvar, pch = ".", 
              cex = 2, add = TRUE, colkey = FALSE)
    
    # XY <- trans3D(x, y= rep(max(y), length(y)), z, pmat = pmat)
    # scatter2D(XY$x, XY$y, col=coluse,  colvar = colvar, pch = ".",
    # cex = 2, add = TRUE, colkey = FALSE)
  }
  scatter3D(x, y, z, ..., colvar = colvar, panel.first=panelfirst,
            colkey = FALSE)
}

phi = 30
theta = 270

#print a pdf with lots of versions!
pdf('PAPER/PCAs/3D_cyno.pdf', height=10,width=10, useDingbats = FALSE)
scatter3D_fancy(pca3d_df$PC_1, pca3d_df$PC_2, pca3d_df$PC_3, col=coluse, colvar = as.integer(pca3d_df$type), pch = c(24,21)[pca3d_df$stage],
                phi = 35, theta = 270, type = 'h',cex = 1.3,  #colkey = FALSE,
                xlab = "PC1", ylab="PC2", zlab = "PC3", main = "PCA",
                ticktype = "detailed", nticks = 4)
scatter3D(pca3d_df$PC_1, pca3d_df$PC_2, pca3d_df$PC_3, col=coluse, bty = "b2",
          colvar = as.integer(pca3d_df$type), 
          phi = 30, theta = 180, colkey = FALSE, cex = 0.6 , pch = pca3dpch,
          xlab = "PC1", ylab="PC2", zlab = "PC3", main = "PCA",
          ticktype = "detailed", nticks =4)
dev.off()

####################################################################################
# Heatmaps

## Nakamura
# run the oxphos and glyco script for nakamura
oxphos <- sort(oxphos$.)
glycolysis <- sort(glycolysis$.)

levels(cyno_data) <- c("ICM_CS3","Epi_CS3","EmDisc_CS5","EmDisc_CS6/7",
                       "Hyp_CS3","VE_CS5","SYS_CS6/7", "Tb_CS3",
  "Tb_CS6/7","ExMes_CS5","ExMes_CS6/7", "PGC_CS5","PGC_CS6/7","ESC_primed")

cols <-c("#5AC0CE","#5B859E","#2AB6B9",
            "#3F94D1", 
            "#E5B405", "#E94E12",
            "#E58603", 
            "#fff5c4",
            "#E3DF1E", 
            "#D79E72", "#7C2C10", 
            "#BA55D3","#4f064f",
            "#F192AE")

DoHeatmap(cyno_data, features = oxphos, group.colors = cols)
ggsave("PAPER/Heatmaps/cyno_nakamura_oxphos.pdf", width=8, height=12)

DoHeatmap(cyno_data, features = glycolysis, group.colors = cols)
ggsave("PAPER/Heatmaps/cyno_nakamura_glycolysis.pdf", width=6, height=6)

### MA
#run the oxphos and glyco script for MA
oxphos <- sort(oxphos$.)
glycolysis <- sort(glycolysis$.)

cType <- c("EmDisc_CS5", 
           "EmDisc_CS6", "EmDisc_CS7",
           "Tb_CS5",
           "Tb_CS6", "Tb_CS7",
           "VE_CS5",
           "Am_CS5","Am_CS6","Am_CS7",
           "SYS_CS6","SYS_CS7", 
           "ExMes_CS5","ExMes_CS6","ExMes_CS7",
           "PGC_CS6")

BaseCol <-c("#2AB6B9",
            "#3F94D1", "#284597",
            "#E94E12",
            "#E58603", "#BE7114",
            "#F7EF8E",
            "#8379B7", "#5554A0", "#2A2A64",
            "#f5d950", "#E6C800",
            "#D79E72", "#A35B2F","#521602",
            "#800080")


levels(cyno_data2) <- c("EmDisc_CS5","EmDisc_CS6","EmDisc_CS7",
                        "Tb_CS5", "Tb_CS6", "Tb_CS7",
                        "VE_CS5",
                        "Am_CS5","Am_CS6","Am_CS7",
                        "SYS_CS6","SYS_CS7",
                        "ExMes_CS5","ExMes_CS6", "ExMes_CS7",
                        "PGC_CS6")
cols2 <-c("#2AB6B9","#3F94D1", "#284597",
          "#F7EF8E","#f5d950", "#E6C800",
            "#E94E12",
          "#8379B7", "#5554A0", "#2A2A64",
          "#E58603", "#BE7114",
            "#D79E72", "#A35B2F", "#521602",
            "#800080")


DoHeatmap(cyno_data2, features = oxphos, group.colors = cols2)
ggsave("PAPER/Heatmaps/cyno_ma_oxphos.pdf", width=8, height=12)

DoHeatmap(cyno_data2, features = glycolysis, group.colors = cols2)
ggsave("PAPER/Heatmaps/cyno_ma_glycolysis.pdf", width=8, height=6)
