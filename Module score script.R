library(Seurat)
library(stringr)
library(tidyverse)
library(dplyr)
library(reshape2)
library(data.table)

# I repeat this code for all datasets (just gave cynomolgus as an example),
# depending on the dataset there are differences when annotating
# seurat clusters and/or formatting of the raw counts file.
# Once everything is annotated properly the rest is the same.

#Set the seed (for repeatable analysis)
set.seed(1)

######## Loading data
Keycyn<- read.table("Datasets/Cyno/cyKey.csv",sep=",",header = T, row.names=1)
raw_counts_cyn <-read.table("Datasets/Cyno/Cynomologous.csv", sep=",",header = T, row.names=1)

isinvit <- Keycyn$All 
labs <- Keycyn$LABEL_4
labs <- labs[which(isinvit>0)]

#Create a Seurat object
cyno_data <- CreateSeuratObject(counts = raw_counts_cyn[which(isinvit>0)], assay = "RNA",min.cells = 0, min.features = 0)
cyno_data$CT <- labs

Idents(cyno_data) <- cyno_data@meta.data$CT #Change the annotation of each cell to e.g. its type Amnnion etc

#Normalise the data
cyno_data <- NormalizeData(cyno_data, verbose = FALSE)
#Get variable genes (reduce number)
cyno_data <- FindVariableFeatures(cyno_data, selection.method = "vst", nfeatures = 30000)
#increase nfeatures to catch all genes

gene_list <- cyno_data@assays[["RNA"]]@var.features

##### Gene lists 
#generating species-specific gene list  

oxphos <- read.csv("GeneOnt/oxphos.csv")
oxphos <-  pull(oxphos,Column2)
oxphos <- str_to_upper(oxphos)
oxphos <-unique(oxphos)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list)%>% 
  unique()

glycolysis <- read.csv("GeneOnt/glycolysis.csv")
glycolysis <-  pull(glycolysis,Column2)
glycolysis <- str_to_upper(glycolysis)
glycolysis <-unique(glycolysis)%>% 
  as.data.frame() %>% 
  filter(. %in% gene_list)%>% 
  unique()

#write.csv(glycolysis, "Datasets/Cyno_glycolysis.csv")
#write.csv(oxphos, "Datasets/Cyno_oxphos.csv")

#Module score
cyno_data <- AddModuleScore(cyno_data, features = oxphos, 
                           name = "oxphos")
cyno_data <- AddModuleScore(cyno_data, features = glycolysis, 
                           name = "glycolysis")

###### Plotting
scores <- data.frame(cyno_data@meta.data$oxphos1,cyno_data@meta.data$glycolysis1,
                     cyno_data@active.ident)

scores <- reshape2::melt(scores, measure.vars=c(1:2), variable.name="score")

colnames(scores) <- c("CT", "score", "value")

#filtering 

embryonic <- c("ICM_CS3","Epi_CS3", "EmDisc_CS5","EmDisc_CS6/7","ESC_primed")

scores <- scores %>% 
  filter(scores$CT %in% embryonic)

#arrange 
scores$CT <- factor(scores$CT , levels=c("ICM_CS3","Epi_CS3", "EmDisc_CS5",
                                         "EmDisc_CS6/7","ESC_primed"))

#plot 
ggplot(scores, aes(x=CT, y=value)) +
  geom_boxplot(alpha=1, size=0.8, aes(fill=score))+
  scale_fill_manual(values = c("#6280e3", "#b83939"),
                    name="", labels = c("OXPHOS", "Glycolysis"))+
  ylab("Score")+
  xlab("Cell type")+
  ggtitle("Cynomolgus")+
  ylim(-0.5,1)+
  geom_hline(yintercept=0, colour="black", size=1.2) + RotatedAxis()
ggsave("GfS/Cyno_metabolism.pdf", width=4, height=4)
