#' @title Generation of haplotypes from a dataset composed by families
#' @description By analyzing all possible combinations of a parent-offspring pedigree in which parents may be missing (missParProb>0), as long as one child was genotyped, it is possible an unequivocal reconstruction of many parental haplotypes. When neither parent was genotyped (missParProb==1), also it is possible to reconstruct at least two parental haplotypes in certain cases. Regarding offspring haplotype reconstruction, if both parents are completely genotyped (missParProb==0), in majority of cases partial offspring haplotypes may be successfully reconstructed (missOffProb>0).
#' @param datasetAlls Dataframe containing clinical and genetic information of families
#' @return familiar haplotypes dataframe
#' @import gtools
#' @export buildHaps
#' @references Medina-Rodriguez, N. Santana A. et al. (2014) alleHap: an efficient algorithm to reconstruct zero-recombinant haplotypes from parent-offspring pedigrees. BMC Bioinformatics, 15, A6 (S-3).
#' @examples
#' ## Haplotype Reconstruction of a dataset composed by 100 families and 12 markers 
#' reconstHaps <- buildHaps(simFams(100,12)[[1]])
#' 
#' ## Haplotype Reconstruction of a dataset composed by 100 families and 12 markers 
#' # with missing parental data (Reconstruction Rate may lead a value lower than one, RR<1)
#' reconstHaps <- buildHaps(simFams(100,12,missParProb=0.1)[[1]])
#' 
#' ## Haplotype Reconstruction of a dataset composed by 100 families and 12 markers 
#' # with missing offspring data (Reconstruction Rate may lead a value lower than one, RR<1)
#' reconstHaps <- buildHaps(simFams(100,12,missOffProb=0.05)[[1]])
buildHaps=function(datasetAlls){  
  idFams <- unique(datasetAlls[,1])                                # Number of identification of the families
  nMrks2 <- ncol(datasetAlls[,-(1:4)])                             # Number of markers * 2
  nMkrs <- nMrks2/2                                                # Number of Markers
  famsHaps <- NULL                                                 # Initialization of variable famHaps
  for (f in idFams) {                                              # Browses the familiar ids
    famMks <- as.matrix(datasetAlls[datasetAlls[,1]==f,-(1:4)])    # Extracts the markers of the family
    nOffs <- nrow(famMks)-2; offspring <- famMks[3:nrow(famMks),]  # Counts the number of offspring
    refRow <- if (!all(is.na(offspring))&length(offspring)>nMrks2) # Extracts the row of reference child (in current family)
                  which(apply(offspring,1,function(x) all(!is.na(x))))[1]+2 else 3  
    ids <- matrix(0,nrow=nOffs,ncol=nMkrs)                         # Initialization of variable ids
    for (j in seq(1,nMrks2,by=2)){                                 # Browses the markers (two by two)
      parents <- famMks[1:2,c(j,j+1)]                              # Extracts the parents alleles
      for (i in 3:nrow(famMks)) {                                  # Browses through all the individuals of the family
        child <- famMks[i,c(j,j+1)]                                # Extracts the child alleles
        idA <- idAlleles(child,parents)                            # Identifies the alleles from the parents in the offpring
        famMks[i,c(j,j+1)] <- idA[1:2]                             # Replaces the identified and sorted alleles of the child
        ids[i-2,(j+1)/2] <- as.integer(idA[3])                     # Saves the identification variable
      }
      refChild <- famMks[refRow,c(j,j+1)]                          # Extracts the alleles of reference child
      famMks[1:2,c(j,j+1)] <- idParents(refChild,parents)          # Identifies the alleles from the parents in the offpring and then replaces the sorted alleles in famMks
    }
    {
      NAoffs <- if (nOffs>1)                                       # Extracts the rows (positions) of children with missing alleles
        which(apply(famMks[3:nrow(famMks),],1,function(x) any(is.na(x)))) else 
          which(any(is.na(famMks[3,])))                           
      ids2 <- if (length(NAoffs)>0) ids[-NAoffs,] else ids         # Removes the rows with missing values from ids matrix
      nOffs2 <- nOffs-length(NAoffs)                               # Removes the number of children with missing values from nOffs vector
      idSum <- colSums(as.matrix(ids))                             # Sums the identificable variable of children by columns
      idRowSum <- rowSums(as.matrix(ids))                          # Sums the identificable variable of children by rows
      idChRows <- which(idRowSum==nMkrs)                           # Positions (rows) of identificable children
      unidChRows <- which(idRowSum>0&idRowSum<ncol(ids))           # Positions (rows) of unidentificable children
      posPs <- which(idSum>=0&idSum<nOffs)                         # Positions (columns) of unidentificable markers
      posPUM <- which(idSum>0&idSum<nOffs)                         # Positions (columns) of partial unidentificable markers
      posTUM <- which(idSum==0)                                    # Positions (columns) of totally unidentificable markers
      posPs2 <- which(idSum>=0&idSum<nOffs2)                       # Positions (columns) of unidentificable markers (excluding rows with missing values)
      posPUM2 <- which(idSum>0&idSum<nOffs2)                       # Positions (columns) of partial unidentificable markers (excluding rows with missing values)
      haplos <- cbind(hap_Parent1=rep(NA,nOffs+2),
                      hap_Parent2=rep(NA,nOffs+2))                 # Initialization of variable haplos
      parentId <- which(!is.na(famMks[1:2,1]))                     # Counts the identified parents
      childId <- if (length(NAoffs)>0)                             # Counts the identified children
                    (3:nrow(famMks))[-NAoffs] else 3:nrow(famMks)
      famRows <- c(parentId,childId)                               # Group the identified parents and children
      idOffsHaps <- NULL                                           # Initialization of variable idOffsHaps
      vec <- rep(1,nMkrs)                                          # Initialization of variable vec 
    }
    if (length(idChRows)>0) {                                       # Checks if there is completely identified children
      for (idCh in idChRows) {                                      # Browses the identified children
        idChHaps <- otherHaps(vec,NULL,famMks[idCh+2,])                 # Creates the identified child haplotypes
        haplos[idCh+2,] <- idChHaps                                 # Inserts the identified child haplotypes in the variable haplos
        idOffsHaps <- as.matrix(unique(rbind(idOffsHaps,idChHaps))) # Saves the identified children haplotypes in the variable idOffsHaps
      }
      haplos[1,] <- insParHaps(1,unique(idOffsHaps[,1]),otherHaps(vec,posPs2,famMks[1,]),haplos) # Inserts the identified father haplotypes in the variable haplos (using identified child haplotypes)
      haplos[2,] <- insParHaps(2,unique(idOffsHaps[,2]),otherHaps(vec,posPs2,famMks[2,]),haplos) # Inserts the identified mother haplotypes in the variable haplos (using identified child haplotypes)
    }    
    if (length(unidChRows)>0&length(childId)>1) {                          # Checks if there is completely identified together with unidentified children
      if (length(parentId)>0) {                                            # Checks if there is completely identified children
        lnfm <- newFamMks(ids2,parentId,childId,posPUM2,posPs2,famMks)     # Creates a list with partial unidentified markers solved and conflict=TRUE if not possible
        if (!lnfm$conflict) famMks <- lnfm$famMks                          # Replaces new familiar markers if there is not conflict (one children haplotype inherited from two identical parental haplotypes)
        else vec[posPs] <- 0                                               # Updates variable vec if there is conflict
      }
      if (length(posTUM)==0)                                               # Checks if there are not positions (markers) totally unidentified
        for (indv in famRows) haplos[indv,] <- unidHaps(vec,famMks[indv,]) # Inserts the identified haplotypes in the variable haplos
      if (any(is.na(haplos[1,]))) haplos[1,] <- reinsParHaps(1,haplos)     # Checks if father haplotypes still remain unidentified. If so, it tries to insert identified children haplotypes  
      if (any(is.na(haplos[2,]))) haplos[2,] <- reinsParHaps(2,haplos)     # Checks if mother haplotypes still remain unidentified. If so, it tries to insert identified children haplotypes  
    }
    if (length(posTUM)>0) {                                                # Checks if there are positions (markers) totally unidentified
      vec[posTUM] <- 0                                                     # Updates variable vec according to totally unidentified markers
      for (id in parentId) haplos[id,] <- unidHaps(vec,famMks[id,])        # Creates the parental haplotypes in the variable haplos (inserting the character "\u00B7" in unidentified positions)
      for (ch in childId) haplos[ch,] <- if (length(parentId)>0)           # Creates the offspring haplotypes in the variable haplos (inserting the character "\u00B7" in unidentified positions)
        unidHaps(vec,famMks[ch,]) else unidHaps(ids[ch-2,],famMks[ch,]) 
      parNoId <- (1:2)[!(1:2)%in%parentId]                                 # Counts the non-identified parents
      for (id in parNoId) haplos[id,] <- unidParHaps(id,ids,famMks,haplos) # Creates the parental haplotypes in the variable haplos according to the number of unidentified markers
    }  
    haplos <- hapsOut(nMkrs,posPs,famMks,haplos)                           # Formats the final output of the familiar haplotypes
    lPoints <- apply(haplos,2,function(x) round(1-ifelse(is.na(x),nMkrs,nchar(gsub("[^\u00B7]","",x)))/nMkrs,4)) # Counts the number of "\u00B7" characters in each haplotype and calculates the reconstruction rate (RR)
    RR <- apply(lPoints,2,prettyNum)                                       # Updates the format the reconstruction rates
    colnames(RR)=c("RR1","RR2")                                            # Re-labels the reconstruction rates columns
    familyHaplos <- cbind(datasetAlls[datasetAlls[,1]==f,1:4],haplos,RR)   # Binds by columns of clinical information together with created haplotypes and RRs  
    famsHaps <- rbind(famsHaps,familyHaplos,deparse.level=0)               # Binds by rows the variable familyHaplos with the previous familiar information
  }
  return(famsHaps)                                                 # Exports the previously generated data.frame
}
