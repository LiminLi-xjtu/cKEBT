#cKBET
install.packages("kBET")
library(kBET)

##input:
#data£ºgene expression matrix, each row represents a cell with each column meaning one gene
#celltype: information on cell types
#batch: batch information of cells
#k: neighborhood size

cKEBT <- function(data, k, celltype, batch){
  celltype <- as.matrix(celltype)
  batch <- as.matrix(batch)
  cellinfor <- cbind(celltype, batch)
  rownames(cellinfor) <- rownames(data)
  colnames(cellinfor) <- c('celltype','batch')
  clust.name <- unique(celltype)
  num <- matrix(data = NA,nrow = nrow(clust.name),ncol = 1)
  cKBET <- matrix(data = NA,nrow = nrow(clust.name),ncol = 1)
  
  for (i in 1:nrow(clust.name)){
    #Selection based on cell type
    cellinfori <- subset(cellinfor,celltype==clust.name[i,])
    #Corresponding batch information
    batchi <- cellinfori[,2]
    batchi <- as.factor(batchi)
    #Corresponding gene expression matrix
    datai <- data[which(rownames(data)%in%rownames(cellinfori)),]
    num[i,]<-nrow(datai)
    ckBETi <- kBET(k0 = k, datai, batchi)
    cKBET[i,]<-mean(ckBETi$stats$kBET.observed)
  }
  
  rejection_rate <- t(num)%*%cKBET / nrow(data)
  
  return(rejection_rate)
}
