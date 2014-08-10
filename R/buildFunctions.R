# Identifies the children alleles according to their parents (left allele inherited from the father and right allele inherited from the mother)
idAlleles=function(child,parents){    
  identified=0
  if (any(is.na(child))&!any(is.na(parents[1,]))&parents[1,1]==parents[1,2]) {child[1] <- parents[1,1]; identified=1}
  if (any(is.na(child))&!any(is.na(parents[2,]))&parents[2,1]==parents[2,2]) {child[2] <- parents[2,1]; identified=1}
  else if (!any(is.na(child))&!any(is.na(parents[1,]))&length(inP2<-which(!child%in%parents[1,]))>0) {child <- child[c(3-inP2,inP2)]; identified=1} 
  else if (!any(is.na(child))&!any(is.na(parents[2,]))&length(inP1<-which(!child%in%parents[2,]))>0) {child <- child[c(inP1,3-inP1)]; identified=1}
  else if (!any(is.na(child))&child[1]==child[2]) identified=1
  return(c(child,identified))  
}
# Sorts the parents alleles according to previous identification of children
idParents=function(child,parents){
  if (!any(is.na(parents[1,]))&parents[1,1]!=parents[1,2]&!is.na(child[1])&parents[1,1]!=child[1]) parents[1,] <- parents[1,2:1]
  if (!any(is.na(parents[2,]))&parents[2,1]!=parents[2,2]&!is.na(child[2])&parents[2,1]!=child[2]) parents[2,] <- parents[2,2:1]
  return(parents)
}
# Inserts the parents haplotypes (when it exists, at least, one completely identified child)
insParHaps=function(id,idOffsHap,otHapsP,haplos){
  idOffsHap <- unique(idOffsHap)
  if (length(idOffsHap)==1) {
    if (is.na(haplos[id,1])) haplos[id,1] <- idOffsHap else haplos[id,2] <- idOffsHap 
    if (!any(is.na(otHapsP))&any(is.na(haplos[id,2]))) {
      hapMatch1 <- otHapsP[which(otHapsP[,1]%in%idOffsHap),2]
      hapMatch2 <- otHapsP[which(otHapsP[,2]%in%idOffsHap),1]
      haplos[id,2] <- if (length(hapMatch1)>0) hapMatch1 else if (length(hapMatch2)>0) hapMatch2 
    }
  } 
  else if (length(idOffsHap)==2) haplos[id,] <- idOffsHap
  return(haplos[id,])
}
# Inserts the parents haplotypes (when it exists, at least, one partial identified child)
reinsParHaps=function(id,haplos){
  uniqHaps <- na.omit(unique(haplos[3:nrow(haplos),id]))
  if (length(uniqHaps)==2) haplos[id,] <- uniqHaps
  else if (length(uniqHaps)==1) {
    if (all(is.na(haplos[id,]))) haplos[id,1] <- uniqHaps[1]
    else if (uniqHaps!=haplos[id,1]) haplos[id,2] <- uniqHaps[1]
  }
  return(haplos[id,])
}
# Creates all the haplotypic combinations when one or more markers are unidentified
otherHaps=function(vec,posPs,indvMks){
  if (!all(is.na(indvMks))) {
    if (length(posPs)>0) {
      indvMks <- t(indvMks)
      m <- permutations(2,length(posPs),c(0,1),repeats.allowed=TRUE)
      omLeft <- indvMks[,seq(1,ncol(indvMks),2)]; omLeft2 <- omLeft[posPs]
      omRight <- indvMks[,seq(2,ncol(indvMks),2)]; omRight2 <- omRight[posPs] 
      omLeftRep <- do.call(rbind, replicate(nrow(m),omLeft,simplify=FALSE))        
      omRightRep <- do.call(rbind, replicate(nrow(m),omRight,simplify=FALSE))    
      omLeftRep[,posPs] <- t(apply(m,1,function(x) ifelse(x==0,omLeft2,omRight2)))
      omRightRep[,posPs] <- t(apply(m,1,function(x) ifelse(x==0,omRight2,omLeft2)))
      leftHaps <- as.matrix(apply(omLeftRep,1,function(x) paste(x,collapse='')))
      rightHaps <- as.matrix(apply(omRightRep,1,function(x) paste(x,collapse='')))
    } 
    else {
      leftHaps <- unidHaps(vec,indvMks)[1]
      rightHaps <- unidHaps(vec,indvMks)[2]
    } 
  } 
  else {leftHaps=NA;rightHaps=NA}
  return(unique(cbind(leftHaps,rightHaps)))
} 
# Solves partial unidentified markers 
newFamMks=function(ids2,parentId,childId,posPUM,posPs,famMks){
  conflict=FALSE
  hapsOffs <- list(father=famMks[childId,seq(1,ncol(famMks),by=2)],mother=famMks[childId,seq(2,ncol(famMks),by=2)])
  for (id in parentId) {
    hapsPar <- rbind(famMks[id,seq(1,ncol(famMks),by=2)],famMks[id,seq(2,ncol(famMks),by=2)])
    for (iter in 1:2) {                 # First iteration to solve parental haplotypes and second to solve offspring markers
      for (ch in 1:length(childId)) {
        chinPar <- if (!is.null(ncol(hapsOffs[[id]]))) hapsOffs[[id]][ch,] else t(hapsOffs[[id]])[ch,] # Extraemos los marcadores del hijo
        posRowMatch <- which(apply(hapsPar,1,function(x) all(chinPar[-posPs]==x[-posPs])))             # Se comparan los haps del hijo con el del padre (menos los marcadores no identificados). En caso de dos coincidencias se escoge la primera posicion
        posids0 <- which(ids2[ch,]==0)
        posFix <- if (iter==1) posPUM[!posPUM%in%posids0] else if (iter==2) posPUM[posPUM%in%posids0]
        if (length(posRowMatch)==1) {   # Checks if one children haplotype has been inherited from one parental haplotype
          for (k in posFix) {
            colsK <- c(k*2-1,k*2)
            if (hapsPar[posRowMatch,k]!=chinPar[k]) {
              if (iter==1) { hapsPar[,k] <- hapsPar[2:1,k]; famMks[id,colsK] <- famMks[id,colsK[2:1]] }
              else if (iter==2) {
                chK1 <- hapsOffs[[id]][ch,k]; chK2 <- hapsOffs[[3-id]][ch,k]
                hapsOffs[[id]][ch,k] <- chK2; hapsOffs[[3-id]][ch,k] <- chK1
                famMks[childId[ch],colsK] <- famMks[childId[ch],colsK[2:1]]
              }
            }
          }
        } else conflict=TRUE             # If conflict=TRUE is not possible to solve the familiar markers (one children haplotype inherited from two identical parental haplotypes)
      }
    }
  }
  return(list(famMks=famMks,conflict=conflict))
}
# Creates individual haplotypes according to the number of unidentified markers
unidHaps=function(idsR,indvMks){
  indvMks <- t(indvMks)
  offsMksOdd <- ifelse(idsR==0,"\u00B7",ifelse(is.na(indvMks[,seq(1,ncol(indvMks),2)]),"\u00B7",indvMks[,seq(1,ncol(indvMks),2)]))
  offsMksEven <- ifelse(idsR==0,"\u00B7",ifelse(is.na(indvMks[,seq(2,ncol(indvMks),2)]),"\u00B7",indvMks[,seq(2,ncol(indvMks),2)]))
  unidHaplo1 <- paste(offsMksOdd,collapse="")
  unidHaplo2 <- paste(offsMksEven,collapse="")
  return(c(unidHaplo1,unidHaplo2))
}
# Creates parental haplotypes according to the number of unidentified markers
unidParHaps=function(id,ids,famMks,haplos){
  hapsOffs <- famMks[3:nrow(famMks),seq(id,ncol(famMks),by=2)]
  colsIds1 <- which(colSums(unique(ids))>1)
  rowsNoMatch=NULL; rowsMatch=NULL
  for (j in colsIds1) {
    rowsIds1 <- which(ids[,j]==1)
    if (length(rowsIds1)>1&length(unique(hapsOffs[rowsIds1,j]))>1){
      levs <- levels(as.factor(hapsOffs[rowsIds1,j]))
      posRowMatch <- match(levs[1:2],hapsOffs[rowsIds1,j])
      namesOffs <- names(hapsOffs[rowsIds1[posRowMatch],j])
      rowsNoMatch <- match(namesOffs,names(hapsOffs[,j]))
    }
  } 
  if (length(rowsNoMatch)>1) haplos[id,] <- haplos[rowsNoMatch[1:2]+2,id] 
  else haplos[id,1] <- haplos[which.max(rowSums(ids))+2,id]
  return(haplos[id,])
}
# Formats the final output of the familiar haplotypes
hapsOut=function(nMkrs,posPs,famMks,haplos){
  haplos[1,] <- sort(haplos[1,],na.last=TRUE)
  haplos[2,] <- sort(haplos[2,],na.last=TRUE) 
  chRows <- which(is.na(haplos[3:nrow(haplos),1]))+2; vec=rep(1,nMkrs); vec[posPs] <- 0
  if (length(chRows)>0) for (ch in chRows) haplos[ch,] <- unidHaps(vec,famMks[ch,])
  familyRows1 <- which(is.na(haplos[,1])); familyRows2 <- which(is.na(haplos[,2]))
  for (fam1 in familyRows1) haplos[fam1,1] <- unidHaps(rep(1,nMkrs),famMks[fam1,])[1]
  for (fam2 in familyRows2) haplos[fam2,2] <- unidHaps(rep(1,nMkrs),famMks[fam2,])[2]  
  return(haplos)
}
