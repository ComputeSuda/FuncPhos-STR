PSNConstructor <- function(pdb, ResidueDis=4.5){
  ca.inds <-  bio3d::atom.select(pdb, "calpha")#select the C_alpha index
  print(ca.inds)
  CaAtom <- pdb$atom[ca.inds$atom,]
  inds <- bio3d::atom.select(pdb,elety=c('CA','C','O','N'),inverse=T) #select the side chain atom
  print(inds)
  SideChainXYZ <- pdb$xyz[inds$xyz]
  SideChainAtom <- pdb$atom[inds$atom,]
  
  print("Calculating the distance between side-chain atoms...")
  SideChainDM <- bio3d::dm.xyz(SideChainXYZ)#calculate the distance between side-chain atom pairs
  
  SideChainAM <- matrix(0,nrow = nrow(SideChainDM),ncol = nrow(SideChainDM))
  SideChainAM[which(SideChainDM<ResidueDis,arr.ind = T)] <- 1
  SideChainAM2 <- SideChainAM + t(SideChainAM)
  
  print("Calculating the Iij between residues...")
  nij <- matrix(0,nrow = nrow(CaAtom),ncol = nrow(CaAtom))#calculate the num of side-chain atom pairs
  for (i in 1:nrow(CaAtom)){
    for (j in 1:nrow(CaAtom)){
      if (i != j){
        #select the side-chain atoms of each residues
        tempN <- SideChainAM2[which(SideChainAtom$chain==CaAtom$chain[i] & SideChainAtom$resno==CaAtom$resno[i]),
                              which(SideChainAtom$chain==CaAtom$chain[j] & SideChainAtom$resno==CaAtom$resno[j])]
        nij[i,j] <- sum(tempN)
      }
    }
  }
  
  load("NormFactor.RData")#path need change
  
  tempNi <- data.frame(Residue=CaAtom$resid,No=1:nrow(CaAtom))
  tempNi2 <- merge(tempNi,NormFactor,by='Residue')
  Ni <- data.matrix(tempNi2[order(tempNi2$No),]$Norm)
  Iij <- nij * 100/((Ni%*%t(Ni))^0.5)

  
  NodeNames <- paste(CaAtom$chain,CaAtom$resno,CaAtom$resid,sep=":")
  colnames(Iij) <- NodeNames
  rownames(Iij) <- NodeNames
  return(Iij)
}




