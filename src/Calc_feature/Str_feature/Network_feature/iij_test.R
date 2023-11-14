library(bio3d)
library(igraph)
setwd("E:/wendang/pycharm1/PTM_function/alphafold2/02calc_feature/str_feature/NACEN")
source('PSN_One.R')

out_path = "iij"

#if (!file.exists(out_path)){
 # dir.create(file.path(out_path))
#}
pdb_name = 'A1X283'
out_file <- paste('E:/wendang/pycharm1/PTM_function/alphafold2/02calc_feature/str_feature/NACEN/feature/',pdb_name,'.txt',sep = "")
pdb_file ='A1X283.pdb'

pdb <- read.pdb(pdb_file)
print(pdb_file)
Iij = PSNConstructor(pdb)
Iij <- (Iij - min(Iij)) / (max(Iij) - min(Iij))
rownames(Iij) <- c(1:dim(Iij)[1])
colnames(Iij) <- c(1:dim(Iij)[2])
dimnames(Iij)
g <- graph_from_adjacency_matrix(Iij, weighted=TRUE, mode='undirected')

clo <- closeness(g)
pr <- page_rank(g)
write.table('closeness', out_file, append=TRUE,row.names=FALSE,col.names=FALSE)
write.table(clo, out_file,append=TRUE,col.names=FALSE)
write.table('page_rank', out_file, append=TRUE,row.names=FALSE,col.names=FALSE)
write.table(pr$vector, out_file,append=TRUE,col.names=FALSE)


