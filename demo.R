#data
data <- read.csv("data.csv",row.names= 1,header=TRUE)
celltype <- read.csv("celltype.csv",row.names= 1,header=TRUE)
batch <- read.csv("batch.csv",row.names= 1,header=TRUE)

#cKBET results
rejection <- cKEBT(data,k=20,celltype,batch)

#Visualization
#t-SNE
library(Rtsne)
batch.id=factor(as.matrix(batch))
cell.type=factor(as.matrix(celltype))
tsne <- Rtsne(data,is.distance=T)
plotdat <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
ggplot(plotdat,aes(x=TSNE1, y=TSNE2)) + geom_point(size=2,aes(color=cell.type, shape=batch.id)) + 
  scale_shape_manual(values = c(3, 20))+
   theme_classic()+theme(title = element_text(size = rel(1.5)))
#UMAP
library(umap)
iris.umap = umap::umap(data)
plotdat <- data.frame(umap1=iris.umap$layout[,1],umap2=iris.umap$layout[,2])
ggplot(plotdat,aes(x=umap1, y=umap2)) + geom_point(size=2,aes(color=cell.type, shape=batch.id)) + 
  scale_shape_manual(values = c(3, 20))+
   theme_classic()+theme(title = element_text(size = rel(1.5)))

