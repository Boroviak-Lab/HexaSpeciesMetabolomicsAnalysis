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

###
#load data
BS<-read.table("Erins data/Keycorrect_CPAll2.csv",sep=",",header = T, row.names=1)

#Load the marmoset data ... to get ordering
raw_counts <- read.table("Erins data/featurecountsAll_CAPProcessed (1).csv", 
                         sep=",",header = T, row.names=1)

isinvit <- BS$QC #1 selected +QC
labs <- BS$Sample_type
labs1 <- BS$Carnegie.Stage
labs2 <- BS$Annotation2
labs <- labs[which(isinvit>0)]
labs1 <- labs1[which(isinvit>0)]

labs2 <- labs2[which(isinvit>0)]


# subset(marmoset_data, idents = c("), invert=FALSE) 

#Create a Seurat object
marmoset_data <- CreateSeuratObject(counts = raw_counts[,which(isinvit>0)], assay = "RNA",min.cells = 0, min.features = 0)
marmoset_data$CS <- labs1
marmoset_data$CT <- labs2
marmoset_data$sample_type <- labs
meta <- marmoset_data@meta.data


CS3 <- c("cMor_CS3","ICM_CS3","Tb_CS3", "Epi_CS3", "Hyp_CS3")
CS5 <- c("Am_CS5","EmDisc_CS5","ExMes_CS5", "PGC_CS5","VE_CS5","Tb_CS5","SYS_CS5")
CS6 <- c("EmDisc_CS6","Am_CS6","ExMes_CS6","PGC_CS6","SYS_CS6","Tb_CS6","VE_CS6")
CS1<- c("Zy_CS1", "4-cell_CS2", "8-cell_CS2")
CS7 <- c("EmDisc_CS7","Am_CS7","ExMes_CS7","PGC_CS7","SYS_CS7","Tb_CS7","VE_CS7")
Cell <- c("ESC_naive", "ESC_primed", "Neonate_neural_cells")

Idents(marmoset_data) <- marmoset_data@meta.data$CT #Change the annotation of each cell to e.g. its type Amnnion etc

marmoset_data <- subset(marmoset_data, idents = c("Zy_CS1", "4-cell_CS2", "8-cell_CS2", "cMor_CS3","ICM_CS3","Tb_CS3", "Epi_CS3", "Hyp_CS3", 
                                                  "Am_CS5","EmDisc_CS5","ExMes_CS5", "PGC_CS5","VE_CS5","Tb_CS5","SYS_CS5", 
                                                  "EmDisc_CS6","Am_CS6","ExMes_CS6","PGC_CS6","SYS_CS6","Tb_CS6","VE_CS6", 
                                                  "EmDisc_CS7","Am_CS7","ExMes_CS7","SYS_CS7","Tb_CS7", 
                                                  "ESC_naive", "ESC_primed", "Neonate_neural_cells"),invert=FALSE) 
###
#Normalise the data
marmoset_data <- NormalizeData(marmoset_data, verbose = FALSE)
#Get variable genes (reduce number)
marmoset_data <- FindVariableFeatures(marmoset_data, selection.method = "vst", nfeatures = 20000)

gene_list_marmoset <- marmoset_data@assays[["RNA"]]@var.features

#scaling 
marmoset_data <- ScaleData(marmoset_data, verbose = FALSE)

marmoset_data <- FindVariableFeatures(marmoset_data, selection.method = "vst", nfeatures = 20000)

#Idents(marmoset_data) <- marmoset_data@meta.data$CT #Change the annotation of each cell to e.g. its type Amnnion etc

#Dimesionality reduction and clustering
marmoset_data <- RunPCA(marmoset_data, npcs = 30, verbose = FALSE) #,features = rownames(Markers))

pct <- marmoset_data@reductions[["pca"]]@stdev/sum(marmoset_data@reductions[["pca"]]@stdev) *100



marmoset_data <- subset(marmoset_data, idents = c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3",
                                                  "EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7",
                                                  "Hyp_CS3","VE_CS5","VE_CS6","SYS_CS5","SYS_CS6","SYS_CS7", "Tb_CS3",
                                                  "Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7",
                                                  "PGC_CS5","PGC_CS6"))

meta <- marmoset_data@meta.data

levels(marmoset_data) <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3",
                           "EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7",
                           "Hyp_CS3","VE_CS5","VE_CS6","SYS_CS5","SYS_CS6","SYS_CS7", "Tb_CS3",
                           "Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7",
                           "PGC_CS5","PGC_CS6", "ESC_naive","ESC_primed","Neonate_neural_cells")

cType <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3",
           "EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7",
           "Hyp_CS3","VE_CS5","VE_CS6","SYS_CS5","SYS_CS6","SYS_CS7", "Tb_CS3",
           "Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7",
           "PGC_CS5","PGC_CS6",
           "ESC_naive","ESC_primed","Neonate_neural_cells")

BaseCol <-c("#DAEDEC","#A6D9E7","#BBDFD5","#5AC0CE","#5B859E","#2AB6B9", 
            "#3F94D1", "#3064AD", "#284597","#8379B7", "#5554A0", "#2A2A64", 
            "#E5B405", "#E94E12", "#B43C0E", "#E58603", "#D17612", "#BE7114",
            "#F7EF8E", "#f5d950", "#E3DF1E", "#E6C800", "#D79E72", "#A35B2F", "#7C2C10", "#BA55D3", 
            "#800080",  "#F9CBCC", "#F192AE", "#E41F85")

colind <- integer( length( levels(Idents(marmoset_data)) )  )
for (i in 1:length( levels(Idents(marmoset_data)) ) ) {
  colind[i] <- which(cType==levels(Idents(marmoset_data))[i])
}
coluse <- BaseCol[colind]

marmoset_data$dataset <- "Bergmann"

DimPlot(marmoset_data, reduction = "pca", dims=c(1, 2), pt.size =2, shape.by = "dataset", label = FALSE, repel=FALSE, cols=coluse)
ggsave("PAPER/PCAs/Dimplot_marmoset_new.pdf", width =8, height=4)

##### Mammalian lists 
gene_list_marmoset <- marmoset_data@assays[["RNA"]]@var.features

oxphos <- read.csv("GeneOnt/oxphos.csv")
oxphos <-  pull(oxphos,Column2)
oxphos <- str_to_upper(oxphos)
oxphos <-unique(oxphos)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_marmoset)%>% 
  unique()

glycolysis <- read.csv("GeneOnt/glycolysis.csv")
glycolysis <-  pull(glycolysis,Column2)
glycolysis <- str_to_upper(glycolysis)
glycolysis <-unique(glycolysis)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list_marmoset)%>% 
  unique()

oxphos <- sort(oxphos$.)
glycolysis <- sort(glycolysis$.)

avgexp_marmoset_oxphos <- AverageExpression(marmoset_data, features=oxphos, use.scale=FALSE)
avgexp_marmoset_glycolysis <- AverageExpression(marmoset_data, features=glycolysis, use.scale=FALSE)

avgexp_marmoset_glycolysis$glycolysis <- "TRUE"
avgexp_marmoset_glycolysis$oxphos <- "FALSE"

avgexp_marmoset_oxphos$oxphos <- "TRUE"
avgexp_marmoset_oxphos$glycolysis <- "FALSE"

avgexp_marmoset_glycolysis <- as.data.frame(avgexp_marmoset_glycolysis)
avgexp_marmoset_oxphos <- as.data.frame(avgexp_marmoset_oxphos)

avgexp <- rbind(avgexp_marmoset_glycolysis, avgexp_marmoset_oxphos)

write.csv(avgexp, "PAPER/avgexp_marmoset.csv")



write.csv(glycolysis, "PAPER/Gene modules/Marmoset_glycolysis.csv")
write.csv(oxphos, "PAPER/Gene modules/Marmoset_oxphos.csv")

marmoset_data <- FindVariableFeatures(marmoset_data, selection.method = "vst", nfeatures = 40000)

marmoset_data <- AddModuleScore(marmoset_data, features = oxphos, 
                                name = "oxphos", control=130)
marmoset_data <- AddModuleScore(marmoset_data, features = glycolysis, 
                                name = "glycolysis", control=35)

scores_marmoset <- data.frame(marmoset_data@meta.data$oxphos1,marmoset_data@meta.data$glycolysis1,
                     marmoset_data@active.ident)
scores_marmoset <- rownames_to_column(scores_marmoset,var="cell")

scores_marmoset <- reshape2::melt(scores_marmoset, measure.vars=c(2:3), variable.name="score")

colnames(scores_marmoset) <- c("cell","CT", "score", "value")

write.csv(scores_marmoset, "PAPER/Scores/marmoset_scores.csv")


emb <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3",
         "Epi_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7")

extraemb <- c("Tb_CS3",
              "Tb_CS5","Tb_CS6","Tb_CS7","Hyp_CS3","VE_CS5","VE_CS6", "VE_CS7")
pgc <- c("PGC_CS5","PGC_CS6")

cells <- c("ESC_naive","ESC_primed","Neonate_neural_cells")

scores_marmoset_emb <- scores_marmoset %>% 
  filter(CT %in% emb)

scores_marmoset_emb$CT <- factor(scores_marmoset_emb$CT , levels=c("Zy_CS1", "4-cell_CS2", "8-cell_CS2",
                                         "cMor_CS3", "ICM_CS3", "Epi_CS3", 
                                         "EmDisc_CS5", "EmDisc_CS6", "EmDisc_CS7"))

scores_marmoset_extraemb$CT <- factor(scores_marmoset_extraemb$CT , levels=c("Tb_CS3",
                                  "Tb_CS5","Tb_CS6","Tb_CS7","Hyp_CS3","VE_CS5","VE_CS6", "VE_CS7"))



ggplot(scores_marmoset_emb, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.25, aes(fill=score))+
  scale_fill_manual(values = c("#289aeb", "#eb2838"),
                    name="", labels = c("OxPhos", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("")+
  ylim(-0.32,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()

ggplot(scores_marmoset_emb, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.5, aes(fill=score))+
  scale_fill_manual(values = c("#6280e3", "#b83939"),
                    name="", labels = c("OXPHOS", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("Embryonic")+
  ylim(-0.5,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()
ggsave("PAPER/Scores/Marmoset_embryo_new.pdf", width=6,height=4)

#### extraembryonic 
scores_marmoset_extraemb <- scores_marmoset %>% 
  filter(CT %in% extraemb)

ggplot(scores_marmoset_extraemb, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.25, aes(fill=score))+
  scale_fill_manual(values = c("#289aeb", "#eb2838"),
                    name="", labels = c("OxPhos", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("")+
  ylim(-0.32,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()
ggsave("PAPER/Scores/Marmoset_extraembryo_new.pdf", width=6,height=4)

### pgc
scores_marmoset_pgc <- scores_marmoset %>% 
  filter(CT %in% pgc)

ggplot(scores_marmoset_pgc, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.5, aes(fill=score))+
  scale_fill_manual(values = c("#6280e3", "#b83939"),
                    name="", labels = c("OXPHOS", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("PGCs")+
  ylim(-0.5,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()
ggsave("PAPER/Scores/Marmoset_pgc.pdf", width=4,height=4)

### cells 
scores_marmoset_cells <- scores_marmoset %>% 
  filter(CT %in% cells)

ggplot(scores_marmoset_cells, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.5, aes(fill=score))+
  scale_fill_manual(values = c("#6280e3", "#b83939"),
                    name="", labels = c("OXPHOS", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("Cell culture")+
  ylim(-0.5,1)+
  geom_hline(yintercept=0, colour="black", size=1.2)+ RotatedAxis()
ggsave("PAPER/Scores/Marmoset_cells.pdf", width=4,height=4)

######################################################################################
#Heatmaps

oxphos <- sort(oxphos$.)
glycolysis <- sort(glycolysis$.)

DoHeatmap(marmoset_data, features = oxphos, group.colors = coluse)
ggsave("PAPER/Heatmaps/marmoset_oxphos_new.pdf", width=10, height=12)

DoHeatmap(marmoset_data, features = glycolysis, group.colors = coluse)
ggsave("PAPER/Heatmaps/marmoset_glycolysis.pdf", width=10, height=6)

##################################################################################
#PCA 3D

# Export PC loadings into a dataframe organized by type
embed_mat <- as.data.frame(marmoset_data@reductions$pca@cell.embeddings)
meta_data <- Idents(marmoset_data)
marmoset_data$stage <- "cell"
stage <- marmoset_data$stage

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
  if(as.character(pca3d_df$stage)[j]   == 'cell') {pca3dpch[j] <- 19}
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
pdf('PAPER/PCAs/3D_marmoset.pdf', height=10,width=10, useDingbats = FALSE)
scatter3D_fancy(pca3d_df$PC_1, pca3d_df$PC_2, pca3d_df$PC_3, col=coluse, colvar = as.integer(pca3d_df$type), pch = c(24,21)[pca3d_df$stage],
                phi = 35, theta = 270, type = 'h',cex = 1.3,  #colkey = FALSE,
                xlab = "PC1", ylab="PC2", zlab = "PC3", main = "PCA",
                ticktype = "detailed", nticks = 4)
scatter3D(pca3d_df$PC_1, pca3d_df$PC_2, pca3d_df$PC_3, col=coluse, bty = "b2",
          colvar = as.integer(pca3d_df$type), 
          phi = 20, theta = 160, colkey = FALSE, cex = 0.6 , pch = pca3dpch,
          xlab = "PC1", ylab="PC2", zlab = "PC3", main = "PCA",
          ticktype = "detailed", nticks =4)
dev.off()

############################################################################################
#Pheatmaps

heatmap_df<- marmoset_data[["RNA"]]@scale.data %>% 
  rownames_to_column(as.data.frame(heatmap_df), "genes") %>% 
  filter(heatmap_df, heatmap_df$genes %in% glycolysis)
heatmap_df <- column_to_rownames(heatmap_df, var="genes")
heatmap_df <- as.matrix(heatmap_df)

annotation_df <- as.data.frame(marmoset_data$CT)
annotation_t <- t(annotation_df)

levels(annotation_df) <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3",
  "EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7",
  "Hyp_CS3","VE_CS5","VE_CS6","SYS_CS5","SYS_CS6","SYS_CS7", "Tb_CS3",
  "Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7",
  "PGC_CS5","PGC_CS6",
  "ESC_naive","ESC_primed","Neonate_neural_cells")

heatmap_df <- heatmap_df[rownames(annotation_df), ]
left_join(heatmap_df, 

heatmap_df %>%
  as.data.table()

heatmap_df <- match(heatmap_df, annotation_df) 
left_join(data.frame(name=target),df,by="name")

  
pheatmap(mat = heatmap_df, show_colnames= FALSE, 
         cluster_rows = FALSE, cluster_cols = FALSE, 
         annotation_col = annotation_df, kmeans_k = NA, )

