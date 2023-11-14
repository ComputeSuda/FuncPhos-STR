#! /usr/bin/env Rscript

# Rscript NACEN.R -d /opt/anaconda3/bin/mkdssp -p protein.pdb

#library(argparse)

# Creating parameter resolution objects.
#parser <- ArgumentParser()

#parser$add_argument("-d", "--dssp_file", type = "character", default = "dssp-3.0.0.exe", help = "The directory of DSSP executive file.")
#parser$add_argument("-p", "--pdb_path", type = "character", help = "The directory of protein complex pdb files.")

#args <- parser$parse_args()

library(NACEN)

# pdb_files_path <- dir('E:/wendang/pycharm1/PTM_function/alphafold2/01str_and_seq/02del_low_pldd_acc')

# dssp_file <- args$dssp_file
# pdb_path <- args$pdb_path
pdb_name <-'A1X283.pdb'
pdb_file <- paste('E:/wendang/pycharm1/PTM_function/alphafold2/01str_and_seq/02del_low_pldd_acc/A1X283_new.pdb')
print(pdb_file)
out_file <- paste('feature/',substring(pdb_name, 1, 6),'_nacen.csv',sep = "")

if(file.exists(out_file))
{
  print(out_file)
  next
}
Net <- NACENConstructor(pdb_file, WeightType = "xxx", exefile = "dssp-3.0.0.exe")
print(" *** Network constructed! *** ")

Net_analyzed <- NACENAnalyzer(Net$AM, Net$NodeWeights)
print(" *** Network analyzed! *** ")

NetP <- data.frame(Net_analyzed$NetP)
# output <- data.frame(cbind(as.character(NetP$ID), as.character(NetP$B), as.character(NetP$chain), as.character(NetP$Resid), as.character(NetP$Res)))
output <- data.frame(cbind(as.character(NetP$ID), as.character(NetP$chain), as.character(NetP$Resid), as.character(NetP$Res), as.character(NetP$K), as.character(NetP$B), as.character(NetP$C), as.character(NetP$Kw), as.character(NetP$Bw), as.character(NetP$Cw)))
colnames(output) <- c("ID", "Chain", "ResID", "ResName", "Unweighted Degree", "Unweighted Betweenness", "Unweighted Closeness", "Node-weighted degree", "Node-weighted Betweenness", "Node-weighted Closeness")

print(out_file)
write.table(output, out_file, quote = F, row.names = F, sep = "\t")
  

