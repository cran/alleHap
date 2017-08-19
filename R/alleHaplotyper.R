#' @title Haplotyping of a dataset composed by several families.
#' @description By analyzing all possible combinations of a parent-offspring pedigree in which parents may be missing (missParProb>0), as long as one child was genotyped, it is possible an unequivocal reconstruction of many parental haplotypes. When neither parent was genotyped (missParProb==1), also it is possible to reconstruct at least two parental haplotypes in certain cases. Regarding offspring haplotypes, if both parents are completely genotyped (missParProb==0), in majority of cases partial offspring haplotypes may be successfully obtained (missOffProb>0).
#' @param data Data containing non-genetic and genetic information of families (or PED file path).
#' @param NAsymbol Icon which will be placed in the NA values of the haplotypes.
#' @param alleSep Icon which will be used as separator of the haplotype alleles.
#' @param invisibleOutput Data are not shown by default.
#' @param dataSummary A summary of the data is shown by default.
#' @return Re-imputed alleles and haplotypes for each loaded family.
#' @import abind
#' @import stats
#' @import utils
#' @export alleHaplotyper
#' @references Medina-Rodriguez, N. Santana A. et al. (2014) alleHap: an efficient algorithm to reconstruct zero-recombinant haplotypes from parent-offspring pedigrees. BMC Bioinformatics, 15, A6 (S-3).
#' @examples
#' 
#' ## Haplotype reconstruction for 3 families without missing data.
#' simulatedFams <- alleSimulator(3,3,6)  
#' (famsAlls <- simulatedFams[[1]])      # Original data 
#' famsList <- alleHaplotyper(famsAlls)  # List containing families' alleles and haplotypes
#' famsList$reImputedAlls                # Re-imputed alleles
#' famsList$haplotypes                   # Reconstructed haplotypes
#' 
#' ## Haplotype reconstruction of a family containing missing data in a parent. 
#' infoFam <- data.frame(famID="FAM002",indID=1:6,patID=c(0,0,1,1,1,1),
#'                      matID=c(0,0,2,2,2,2),sex=c(1,2,1,2,1,2),phenot=c(2,1,1,2,1,2))
#' Mkrs <- rbind(c(1,4,2,5,3,6),rep(NA,6),c(1,7,2,3,3,2),
#'               c(4,7,5,3,6,2),c(1,1,2,2,3,3),c(1,4,2,5,3,6))
#' colnames(Mkrs) <- c("Mk1_1","Mk1_2","Mk2_1","Mk2_2","Mk3_1","Mk3_2")
#' (family <- cbind(infoFam,Mkrs))    # Original data 
#' famList <- alleHaplotyper(family)  # List containing family's alleles and haplotypes
#' famList$reImputedAlls              # Re-imputed alleles
#' famList$haplotypes                 # Reconstructed haplotypes
#' 
#' ## Haplotype reconstruction from a PED file
#' pedFamPath <- file.path(find.package("alleHap"), "examples", "example3.ped") # PED file path
#' pedFamAlls <- alleLoader(pedFamPath,dataSummary=FALSE) 
#' pedFamList <- alleHaplotyper(pedFamAlls)
#' pedFamAlls                # Original data 
#' pedFamList$reImputedAlls  # Re-imputed alleles 
#' pedFamList$haplotypes     # Reconstructed haplotypes
#' 
alleHaplotyper=function(data,NAsymbol="?",alleSep="",invisibleOutput=TRUE,dataSummary=TRUE){

  ############################################# I. INTERNAL FUNCTIONS ###########################################
  { 
    # Haplotyping per family (main function)
    famHaplotyper=function(family){  
      #=== INTERNAL FUNCTIONS ===#
      as.familyArray=function(Mkrs){
        famSize <- nrow(Mkrs)
        nOffs <- famSize-2
        nMarkers <- ncol(Mkrs)/2
        outArray <- array(t(Mkrs),dim=c(2,nMarkers,famSize), 
                          dimnames=list(haplo=c("hap1","hap2"),marker=paste("Mkr",1:nMarkers,sep=""),
                                        subject=c("father","mother",paste("child",1:nOffs,sep=""))))
        return(outArray)
      }
      as.familyMatrix=function(familyArray){
        outArray <- t(array(familyArray,dim=c(2*dim(familyArray)[2],dim(familyArray)[3])))
        return(outArray)
      }
      excludeCase5Mkrs=function(fam,case5){
        fam$nMarkers=fam$nMarkers-length(case5$mrk)
        fam$imputedAlls=fam$imputedAlls[,-case5$mrkCols]
        fam$allelesNumber=fam$allelesNumber[-case5$mrk]
        fam$reImpAlls=fam$reImpAlls[,-case5$mrk,,drop=FALSE]
        return(fam)
      }
      findHaps=function(fam){
        ### Internal Functions ###
        canBeEqual=function(hap1,hap2){
          compareAlls=function(v1,v2) all(v1==v2,na.rm=TRUE)|all(v1==v2[2:1],na.rm=TRUE)
          n=dim(hap1)[2]
          comps=vector("logical",n)
          for (j in 1:n) comps[j]=compareAlls(hap1[,j],hap2[,j])
          return(all(comps))
        }
        identifyChildHaps=function(childIndex){
          child=fam$reImpAlls[,,childIndex+2]
          compatible=which(sapply(possibleHapPairs,canBeEqual,child))
          if (length(compatible)==1) idHaps=hapPairOrigin[compatible,]
          else if (length(compatible)>1){
            idHaps=c(0,0)
            hpo=hapPairOrigin[compatible,]
            for (j in 1:2){
              ph=unique(hpo[,j])
              if (length(ph)==1) idHaps[j]=ph 
            }
          }
          else if (length(compatible)==0) idHaps=c(NA,NA) ## RECOMBINATION EVENT !!!!
          idHaps
        }
        updateGeno=function(genotype,genotypeIDS,markersToUpdate,identifiedAlleles,rowToUpdate){
          for (j in 1:length(markersToUpdate)){
            allele=identifiedAlleles[j]
            mk=markersToUpdate[j]
            if (is.na(genotype[rowToUpdate,mk])){
              if (is.na(genotype[3-rowToUpdate,mk])){ 
                genotype[rowToUpdate,mk]=allele
                genotypeIDS[rowToUpdate,mk]=1
              } # Case 10
              else if (genotype[3-rowToUpdate,mk]==allele){ 
                if (genotypeIDS[3-rowToUpdate,mk]==0){ 
                  genotype[,mk]=genotype[2:1,mk]
                  genotypeIDS[rowToUpdate,mk]=1
                } # Case 6
                else{ 
                  genotype[rowToUpdate,mk]=allele
                  genotypeIDS[rowToUpdate,mk]=1
                } # Case 7
              } # Cases 6,7
              else{ 
                if (genotypeIDS[3-rowToUpdate,mk]==0){
                  genotype[rowToUpdate,mk]=allele
                  genotypeIDS[1:2,mk]=1
                } # Case 8
                else{
                  genotype[rowToUpdate,mk]=allele
                  genotypeIDS[rowToUpdate,mk]=1
                } # Case 9
              } # Cases 8,9
            } # Cases 6,7,8,9,10
            else{ # when !is.na(genotype[rowToUpdate,mk])
              if (is.na(genotype[3-rowToUpdate,mk])){
                if (genotype[rowToUpdate,mk]==allele){
                  genotypeIDS[rowToUpdate,mk]=1
                } # Case 1
                else {
                  genotype[,mk]=genotype[2:1,mk]
                  genotype[rowToUpdate,mk]=allele
                  genotypeIDS[1:2,mk]=1
                } # Case 5
              } # Cases 1,5
              else if (genotype[rowToUpdate,mk]==allele){
                genotypeIDS[1:2,mk]=1
              } # Case 2
              else{
                if (genotype[3-rowToUpdate,mk]==allele){
                  genotype[,mk]=genotype[2:1,mk]
                  genotypeIDS[1:2,mk]=1
                } # Case 3
                else{
                  return(list(genotype=NA,genotypeIDS=NA, errorMkr=mk))
                } # Case 4
              } # Cases 3,4
            } # Cases 1,2,3,4,5
          }
          return(list(genotype=genotype,genotypeIDS=genotypeIDS))
        }
        updateFamHaps=function(childIndex,fam){
          child=fam$reImpAlls[,,childIndex+2]
          childIDS=fam$IDS[,,childIndex+2]
          parents=fam$reImpAlls[,,1:2]
          parentsIDS=fam$IDS[,,1:2]
          idh=fam$idHaps[childIndex,]
          for (i in 1:2){
            difids=childIDS[i,]-parentsIDS[idh[i],,i]
            par2ch=which(difids==-1)  # alleles identified in parent but not in child
            ch2par=which(difids==1)   # alleles identified in child but not in parent
            if (length(par2ch)>0){
              updatedChild=updateGeno(child,childIDS,par2ch,parents[idh[i],par2ch,i],rowToUpdate=i)
              if (exists("errorMkr",updatedChild)){
                fam$hapIncidences=paste("Genotyping error or recombination event in marker",
                                        updatedChild$mk,"in subject",childIndex+2)
                return(fam)
              }
              fam$reImpAlls[,,childIndex+2]=updatedChild$genotype
              fam$IDS[,,childIndex+2]=updatedChild$genotypeIDS
            }
            if (length(ch2par)>0){
              updatedParent=updateGeno(parents[,,i],parentsIDS[,,i],ch2par,child[i,ch2par],rowToUpdate=idh[i])
              if (exists("errorMkr",updatedParent)){
                fam$hapIncidences=paste("Genotyping error or recombination event in marker",updatedParent$mk,"in subject",i)
                return(fam)
              }
              fam$reImpAlls[,,i]=updatedParent$genotype
              fam$IDS[,,i]=updatedParent$genotypeIDS
              idhaps0=which(apply(fam$idHaps,1,function(x) any(x==0)))
              if (length(idhaps0)>0){
                parentsHaps=fam$reImpAlls[,,1:2]
                parentsHaps[fam$IDS[,,1:2]==0]=NA
                possibleHapPairs=list(t(parentsHaps[1,,]),rbind(parentsHaps[1,,1],parentsHaps[2,,2]),
                                      rbind(parentsHaps[2,,1],parentsHaps[1,,2]),t(parentsHaps[2,,]))
                for (k in idhaps0) fam$idHaps[k,]=identifyChildHaps(k)
              }
            }
          }
          return(fam)
        }
        
        ### Initialization ###
        hapPairOrigin=matrix(rbind(c(1,1),c(1,2),c(2,1),c(2,2)),ncol=2)
        parentsHaps=fam$reImpAlls[,,1:2]
        parentsHaps[fam$IDS[,,1:2]==0]=NA
        possibleHapPairs=list(t(parentsHaps[1,,]),rbind(parentsHaps[1,,1],parentsHaps[2,,2]),
                              rbind(parentsHaps[2,,1],parentsHaps[1,,2]),t(parentsHaps[2,,]))
        fam$idHaps=matrix(0,nrow=fam$famSize-2,ncol=2)
        for (k in 1:(fam$famSize-2)) fam$idHaps[k,]=identifyChildHaps(k)
        if (any(is.na(fam$idHaps))){
          idhapsNA=which(apply(fam$idHaps,1,function(x) any(is.na(x))))
          fam$hapIncidences=paste("Genotyping error or recombination event in subject(s)",idhapsNA,collapse=",")
          fam$haps=writeFamilyHaps(fam)
        } 
        else {
          ### Haplotyping ###
          sumIDS0=sum(fam$IDS)
          repeat{
            for (k in 1:(fam$famSize-2)) {
              fam=updateFamHaps(k,fam)
              if (exists("hapIncidences",fam)){
                fam$haps <- writeFamilyHaps(fam)
                return(fam)
              } 
            }
            sumIDS1=sum(fam$IDS)
            if (sumIDS1==sumIDS0) break else sumIDS0=sumIDS1
          }
          fam$haps=writeFamilyHaps(fam)
        }
        return(fam)
      }
      includeCase5Mkrs=function(fam,fam0,case5,NAsymbol="?"){
        fam0$imputedAlls[,-case5$mrkCols]=fam$imputedAlls
        fam0$reImpAlls[,-case5$mrk,]=fam$reImpAlls
        fam0$HMZ=matrix(0,nrow=fam0$famSize,ncol=fam0$nMarkers)
        fam0$HMZ[,-case5$mrk]=fam$HMZ
        fam0$HMZ[,case5$mrk]=t(apply(fam0$reImpAlls[,case5$mrk,,drop=FALSE],3,function(x) apply(x,2,function(a) as.integer(a[1]==a[2]))))
        fam0$IDS=array(0,dim=c(2,fam0$nMarkers,fam0$famSize))
        fam0$IDS[,-case5$mrk,]=fam$IDS
        fam0$idHaps=fam$idHaps
        fam0$haps=writeFamilyHaps(fam0)
        nf0=names(fam0)
        nf=names(fam)
        nv=nf[which(!nf%in%nf0)]
        for (k in nv) fam0[[k]]=fam[[k]]
        return(fam0)
      }
      inheritedHap=function(hap,parentHaps,hapIDS,pHapsIDS){
        # This function determines which one of parentHaps is (or can be) equal to 'hap'; it returns 0 
        # if both parentHaps are compatible with 'hap'
        if(is.null(dim(parentHaps))){ 
          dim(hap) <- dim(hapIDS) <- c(1,1)
          dim(parentHaps) <- dim(pHapsIDS) <- c(2,1)
        } 
        canBeSameHap=vector("logical",2)
        for (i in 1:2){
          pos=which(hapIDS+pHapsIDS[i,]==2)   
          if (length(pos)==0) canBeSameHap[i]=TRUE
          else canBeSameHap[i]=all(hap[pos]==parentHaps[i,pos])  
          if (length(pos)==ncol(pHapsIDS)&canBeSameHap[i]) return(i)
        }
        inhrtdHap=which(canBeSameHap)
        if (length(inhrtdHap)==2) inhrtdHap=0
        else if (length(inhrtdHap)==0) inhrtdHap=3
        return(inhrtdHap)
      }
      initializeIDS=function(fam){  
        reImpAlls <- as.familyArray(fam$imputedAlls)
        ### HMZ initialization
        HMZ <- t(apply(reImpAlls,3,function(x) apply(x,2,function(a) as.integer(a[1]==a[2]))))
        if (nrow(HMZ)==1) HMZ <- t(HMZ)
        ### IDS initialization
        IDS=with(fam,array(NA,dim=c(2,nMarkers,famSize), 
                           dimnames=list(haplo=c("hap1","hap2"),marker=paste("Mkr",1:nMarkers,sep=""),
                                         subject=c("father","mother",paste("child",1:(famSize-2),sep="")))))
        IDS[1,,] = IDS[2,,] = t(HMZ)
        IDS[is.na(IDS)]<-0
        ### IDS in children
        for (i in 3:fam$famSize) 
          for (j in 1:fam$nMarkers){
            NAchild <- which(is.na(reImpAlls[,j,i]))
            NAfather <- which(is.na(reImpAlls[,j,1]))
            NAmother <- which(is.na(reImpAlls[,j,2]))
            # A child has only one NA value if the other value has been imputed from an homozygous parent or
            # maybe after a pass by the three children algorithm, which may have re-imputed only one allele.
            if (length(NAchild)==1) {  
              hompar <- which(HMZ[1:2,j]==1)
              chAl <- reImpAlls[3-NAchild,j,i]
              proc <- which(c(chAl %in% reImpAlls[,j,1],chAl %in% reImpAlls[,j,2]))
              if (length(proc)==1){
                if (NAchild==proc){
                  reImpAlls[,j,i]=reImpAlls[2:1,j,i]
                  NAchild=3-NAchild
                } 
                IDS[proc,j,i]=1
              } else if (length(proc)==2){
                if (length(hompar)==1){
                  if (NAchild==hompar) {
                    reImpAlls[,j,i]=reImpAlls[2:1,j,i]
                    NAchild=3-NAchild
                  }
                  IDS[hompar,j,i]=1
                }
              }
              if (length(hompar)==1){
                parAl=reImpAlls[1,j,hompar]
                if (chAl!=parAl){
                  reImpAlls[NAchild,j,i]=parAl
                  IDS[NAchild,j,i]=1
                }
              } else if (length(hompar)==2){
                reImpAlls[,j,i]=reImpAlls[1,j,1:2]
                IDS[,j,i]=1
              }
            } 
            else if (length(NAchild)==0){  # Child without missing alleles
              if (length(NAfather)+length(NAmother)==0&HMZ[i,j]==0){ # Parents without missing alleles
                fromWhichParent<-list(fromFather=which(reImpAlls[,j,i]%in%reImpAlls[,j,1]),
                                      fromMother=which(reImpAlls[,j,i]%in%reImpAlls[,j,2]))
                ia=lapply(fromWhichParent,length)
                if (ia[[1]]+ia[[2]]<4){ 
                  ia1=which(ia==1)[1]
                  if (fromWhichParent[[ia1]]!=ia1) reImpAlls[,j,i]=reImpAlls[2:1,j,i]
                  IDS[,j,i]=1        
                }
              } else if (length(NAfather)>0&length(NAmother)==0){ # Missing alleles only in father
                fromFather=which(!reImpAlls[,j,i]%in%reImpAlls[,j,2])  
                if (length(fromFather)==1){ 
                  if (fromFather==2) reImpAlls[,j,i]=reImpAlls[2:1,j,i]
                  IDS[,j,i]=1
                }
              } else if (length(NAfather)==0&length(NAmother)>0){ # Missing alleles only in mother
                fromMother=which(!reImpAlls[,j,i]%in%reImpAlls[,j,1])
                if (length(fromMother)==1){ 
                  if (fromMother==1) reImpAlls[,j,i]=reImpAlls[2:1,j,i]
                  IDS[,j,i]=1
                }
              }
            }
          }
        ### IDS in parents 
        # IDS in parents is initialized from the child with more IDS==1
        for (i in 1:2){
          if (fam$famSize==3) keyCh<-1
          else{
            nrIDS1=if (fam$nMarkers==1) IDS[i,,-(1:2)] else colSums(IDS[i,,-(1:2)],na.rm=TRUE)  ## MODIFICADO 20/8/2015
            keyCh<-which.max(nrIDS1)   ## MODIFICADO 20/8/2015
          }
          keyMks=which(IDS[i,,keyCh+2]==1)
          for (j in keyMks){
            pap=which(match(reImpAlls[,j,i],reImpAlls[i,j,keyCh+2])==1)
            if (!is.na(pap[1])) if (pap[1]==2) reImpAlls[,j,i]=reImpAlls[2:1,j,i]
            IDS[which(!is.na(reImpAlls[,j,i])),j,i]=1
          }
        }
        # Which parents allele has been inherited by child 
        idHaps<-array(0,dim=c(fam$famSize,2))
        for (i in 3:fam$famSize)
          for (j in 1:2)
            idHaps[i,j]=inheritedHap(reImpAlls[j,,i],reImpAlls[,,j],IDS[j,,i],IDS[,,j])
        ### New variables attached
        fam$reImpAlls<-reImpAlls
        fam$HMZ <- HMZ
        fam$IDS <- IDS
        fam$idHaps <- idHaps
        return(fam)
      }
      makeHapsFromThreeChildren=function(fam){
        makeFamilyFromChildren=function(fam){ 
          findHapsfrom3or4Children=function(children,parentsInfo=NULL,NAsymbol="?"){
            ### Internal Functions 
            canBeEqual=function(hapPair1,hapPair2){
              # Compares two haplotype pairs
              for (j in 1:ncol(hapPair1)){
                equalMrk=(all(hapPair1[,j]==hapPair2[,j],na.rm=TRUE))|
                  (all(hapPair1[,j]==hapPair2[,j][2:1],na.rm=TRUE))
                if (!equalMrk) return(FALSE)
              }
              return(TRUE)
            }  
            commonHap=function(twoChildren){
              # Finds the possible common haplotypes between two children (stored in the array 'twoChildren'). 
              commonAllelles <- vector("list")  
              nMarkers <- dim(twoChildren)[2]
              for (j in 1:nMarkers) commonAllelles[[j]] <- intersect(twoChildren[,j,1],twoChildren[,j,2])          
              commonHap <- expand.grid(commonAllelles,KEEP.OUT.ATTRS=FALSE)
              names(commonHap) <- paste("Mkr",1:ncol(commonHap),sep="")
              return(commonHap)
            }
            complementHap=function(hap,child){
              # Determines the complementary haplotype for a given one. 
              # This function supose that "hap" is placed in "child", and it finds the other haplotype.
              comp <- child[1,]
              flip <- (comp==hap)
              comp[flip] <- child[2,flip]
              return(comp)
            }
            findParHaps=function(chldHaps){
              # Function to determine the parental haplotypes from the haplotypes found in children:
              foundHaps=t(apply(chldHaps,3,function(x) apply(x,1,paste,collapse="")))
              nHaps=length(unique(as.vector(foundHaps))) # Num de haplotipos distintos
              hmzCh=sum(foundHaps[,1]==foundHaps[,2])    # num de hijos homocigotos
              if(nHaps==3&hmzCh==0)  
                return(list(sortedParents=NULL, sortedChildren=chldHaps, 
                            hapIncidences <- "Multiple compatible haplotypes in parents",
                            alignedParents=FALSE))
              hapIncidences=""
              parentOrdered=all(apply(foundHaps,2,function(x) length(unique(x)))<3)
              if (!parentOrdered){
                i=1
                repeat{
                  foundHaps[i,]=foundHaps[i,2:1]
                  parentOrdered=all(apply(foundHaps,2,function(x) length(unique(x)))<3)
                  if (!parentOrdered){
                    foundHaps[i,]=foundHaps[i,2:1]
                    i=i+1
                    if (i==4) {
                      parents=array(NA,dim=c(2,nMrk,2))
                      hapIncidences <- "Genotyping error, recombination event or inheritance from non-declared parent"
                      break
                    }
                  } 
                  else{
                    chldHaps[,,i]=chldHaps[2:1,,i]
                    break
                  }
                }
              }
              if (hapIncidences==""){
                parents=abind(t(unique(chldHaps[1,,],MARGIN=2)), t(unique(chldHaps[2,,],MARGIN=2)),along=3)
                parents=abind(sortByRows(parents[,,1]),sortByRows(parents[,,2]),along=3)
                if (is.null(parentsInfo)) alignedParents=FALSE else {
                  eqPar=matrix(NA,2,2)
                  for (i in 1:2)
                    for (j in 1:2)
                      eqPar[i,j]=canBeEqual(parents[,,i],parentsInfo[,,j])
                  if (any(rowSums(eqPar)==0)|any(colSums(eqPar)==0)) {
                    alignedParents=FALSE
                    parents[,,1]<-parents[,,2]<-NA
                    hapIncidences <- "Parental information is not compatible with haplotypes found in children"
                  } else{
                    weq=which(!eqPar,arr.ind=TRUE)
                    if (length(weq)==0){
                      alignedParents=FALSE
                    } else{
                      alignedParents=TRUE
                      if (weq[1,1]==weq[1,2]){
                        parents=parents[,,2:1]
                        chldHaps=chldHaps[2:1,,]
                      } 
                    }
                  }
                }      
              }    
              return(list(sortedParents=parents, sortedChildren=chldHaps, 
                          hapIncidences=hapIncidences, alignedParents=alignedParents))
            }  
            identifyFourthChild=function(famHaps,fourthChild){
              # Function for identifying haplotypes in fourth child
              c4possibleHaps=list(t(famHaps$sortedParents[1,,1:2]),t(famHaps$sortedParents[2,,1:2]),
                                  rbind(famHaps$sortedParents[1,,1],famHaps$sortedParents[2,,2]),
                                  rbind(famHaps$sortedParents[2,,1],famHaps$sortedParents[1,,2]))
              c4possibleHapsSorted=lapply(c4possibleHaps,function(ph) apply(ph,2,sort))
              fourthChild=apply(fourthChild,2,sort)
              c4Haps=which(sapply(c4possibleHapsSorted,function(x) all(x==fourthChild | x==fourthChild[2:1,])))
              if (length(c4Haps)>0) famHaps$sortedChildren=abind(famHaps$sortedChildren, c4possibleHaps[[c4Haps]],along=3) 
              else{
                famHaps$sortedChildren=abind(famHaps$sortedChildren, fourthChild,along=3) 
                famHaps$hapIncidences="Haplotypes in one child are not compatible with the haplotypes found in the rest of the offspring"
              } 
              return(famHaps)
            }
            hapInChild=function(hap,child){
              # Determines if the haplotype 'hap' is located in a given child. 
              # In such a case returns the children haplotypes including 'hap' and its complementary haplotype.
              # If a haplotype is not found, the function returns "NULL".
              if (is.null(dim(child))) child=matrix(child,ncol=1)
              isin=(t(child)==hap)
              if (any(rowSums(isin)==0)) return(NULL)
              s=which(!isin[,1])
              child[,s]=child[2:1,s]
              return(matrix(child,2))
            }
            sortByRows=function(arr) {
              # Function for sorting a matrix row by row
              return(arr[do.call(order, lapply(1:NCOL(arr), function(i) arr[, i])),])
            }
            
            ### Initialization 
            nMrk <- dim(children)[2]
            nchld <- dim(children)[3]
            if (nchld==4) {
              fourthChild=children[,,4]
              children=children[,,1:3]
            }
            if (nchld<3){
              return(list(sortedParents=array(NA,dim=c(2,nMrk,2)),sortedChildren=children,alignedParents=FALSE,
                          hapIncidences="Less than three children. Haplotypes can not be generated at this stage",
                          childrenIDS=array(0,dim=c(2,nMrk,nchld)), parentsIDS=array(0,dim=c(2,nMrk,2))))
            }
            
            ### Haplotype finding from three unique children without missing values nor recombination.
            commonHaps <- vector("list")
            for (i in 1:3) commonHaps[[i]] <- commonHap(children[,,-i])
            nch <- sapply(commonHaps,nrow)
            refH <- which(nch>0)
            chHaps <- list()
            for (ref in refH) {    
              chld=(1:3)[-ref]
              # The complemmentary haplotypes are added in each child
              possibleHaps <- array(rbind(t(commonHaps[[ref]]), apply(commonHaps[[ref]],1,complementHap,children[,,chld[1]]),
                                          apply(commonHaps[[ref]],1,complementHap,children[,,chld[2]])),dim=c(nMrk,3,nch[ref]))    
              for (ih in 1:nch[ref]){
                identifiedHaps <- t(possibleHaps[,,ih])
                uHaps <- unique(identifiedHaps)
                for (j in 1:nrow(uHaps)) {
                  hapsChRef <- hapInChild(uHaps[j,],children[,,ref])
                  if (!is.null(hapsChRef)){
                    nHaps <- nrow(unique(rbind(uHaps,hapsChRef)))
                    if (nHaps<=4){              
                      chh <- rbind(identifiedHaps[c(1,2,1,3),],hapsChRef)
                      chh <- chh[order(2*c(sort(rep((1:3)[-ref],2)),ref,ref)-c(1,0)),]            
                      nh=nrow(unique(chh)) # Number of different haplotypes in the three children
                      ndh=nrow(unique(chh[duplicated(chh),])) # Num. of haplotypes which are repeated
                      if (!(nh==4&ndh==1)) {
                        chHaps[[length(chHaps)+1]]<-array(cbind(sortByRows(chh[1:2,]),sortByRows(chh[3:4,]),
                                                                sortByRows(chh [5:6,])),dim=c(2,nMrk,3))    
                      }
                    }
                  }
                }
              }
            }
            chHaps <- unique(chHaps)
            
            ### Identification and reconstruction of the posible parental haplotypes
            famHaps=lapply(chHaps,findParHaps)
            noIncidences=which(sapply(famHaps,function(x) nchar(x$hapIncidences))==0)
            # If there is only one posible haplotype structure in children
            famHaps <- if (length(noIncidences)==0) famHaps[1] else famHaps[noIncidences]
            # If there is a fourth child (diferent from the rest of the offspring)
            if (nchld>3){
              famHaps <- lapply(famHaps,identifyFourthChild,fourthChild)
              noIncidences=which(sapply(famHaps,function(x) nchar(x$hapIncidences))==0)
              famHaps <- if (length(noIncidences)==0) famHaps[1] else famHaps[noIncidences]
            }
            
            ### IDS calculation
            if (length(famHaps)>1) return(NULL) else famHaps=famHaps[[1]]
            ids <- if (famHaps$hapIncidences=="") 1 else 0
            famHaps$childrenIDS=array(ids,dim=c(2,nMrk,nchld))
            famHaps$parentsIDS=array(ids,dim=c(2,nMrk,2))
            
            ### Function output
            return(famHaps)
          }
          isCase5=function(mrk){
            # Indicates if a marker is in case 5
            globAlNr=length(unique(as.vector(mrk[!is.na(mrk)])))  # Global number of alleles in marker
            indivAlNr=apply(mrk,3,function(x) if(any(is.na(x))) 2 else length(unique(x)))
            return((globAlNr==2) & all(indivAlNr==2))
          }
          selChMk=function(ndx,children){
            alleMCh=alleMat[ndx,]
            noNAMks=which(colSums(alleMCh)==0)
            if (length(noNAMks)>0){
              c5=which(sapply(noNAMks,function(k) isCase5(children[,k,ndx,drop=FALSE])))
              if (length(c5)>0) noNAMks=noNAMks[-c5]
            }
            if (length(noNAMks)>1){
              sch=array(apply(children[,noNAMks,ndx,drop=FALSE],3, function(M) apply(M,2,sort)),dim=c(2,length(noNAMks),3))
              idupCh=duplicated(sch,MARGIN=3)
              uch=ndx[!idupCh]  # unique complete Children
              if (length(uch)>=3) return(list(children=uch,markers=noNAMks))
            }
          }
          updateParentsIDS=function(fam){
            ## IDS in parents is initialized from the child with more IDS==1
            for (i in 1:2){
              fam$IDS[,,i]=0
              nrIDS1 <- if (fam$nMarkers==1) fam$IDS[i,,-(1:2)] else colSums(fam$IDS[i,,-(1:2)],na.rm=TRUE)
              keyCh <- which.max(nrIDS1)   
              keyMks <- which(fam$IDS[i,,keyCh+2]==1)
              for (j in keyMks){
                chAl <- fam$reImpAlls[i,j,keyCh+2]
                newa <- which(!chAl%in%fam$reImpAlls[,j,i])
                if(length(newa)>0){
                  pos <- which(is.na(fam$reImpAlls[,j,i]))[1:length(newa)]
                  fam$reImpAlls[pos,j,i]=chAl[newa]
                }
                pap=which(match(fam$reImpAlls[,j,i],fam$reImpAlls[i,j,keyCh+2])==1)
                if (!is.na(pap[1])) if (pap[1]==2) fam$reImpAlls[,j,i]=fam$reImpAlls[2:1,j,i]
                fam$IDS[which(!is.na(fam$reImpAlls[,j,i])),j,i]=1
              }
            }
            # Which parents allele has been inherited by child 
            fam$idHaps<-array(0,dim=c(fam$famSize,2))
            for (i in 3:fam$famSize)
              for (j in 1:2)
                fam$idHaps[i,j]=inheritedHap(fam$reImpAlls[j,,i],fam$reImpAlls[,,j],fam$IDS[j,,i],fam$IDS[,,j])
            return(fam)
          }
          compareNA=function(v1,v2) {
            same <- (v1 == v2) | (is.na(v1) & is.na(v2))
            same[is.na(same)] <- FALSE
            return(same)
          }
          find2Haps=function(mat){
            canBeEqual=function(v1,v2) all(v1==v2,na.rm=TRUE)
            nr=nrow(mat)
            differents=NULL
            for (i in 1:(nr-1)){
              for(j in 2:nr)
                if (!canBeEqual(mat[i,],mat[j,])){
                  differents=c(i,j)
                  break
                }
              if (!is.null(differents)) break
            }
            if (is.null(differents)){
              nna=apply(mat,1,function(x) length(which(is.na(x))))
              return(rbind(mat[which.min(nna),],rep(NA,ncol(mat))))
            } else{
              hap=rbind(mat[i,],mat[j,])
              hta=(1:nr)[-c(i,j)]
              repeat{
                natot0=length(which(is.na(hap)))
                for (k in hta){
                  if(canBeEqual(mat[k,],hap[1,])&!canBeEqual(mat[k,],hap[2,])){
                    h1na=intersect(which(is.na(hap[1,])),which(!is.na(mat[k,])))
                    hap[1,h1na]=mat[k,h1na]
                    hta=hta[-which(hta==k)]
                  } 
                  else if (!canBeEqual(mat[k,],hap[1,])&canBeEqual(mat[k,],hap[2,])){
                    h2na=intersect(which(is.na(hap[2,])),which(!is.na(mat[k,])))
                    hap[2,h2na]=mat[k,h2na]
                    hta=hta[-which(hta==k)]
                  }
                }
                natot1=length(which(is.na(hap)))
                if (natot0==natot1) break
              }
              return(hap)
            }
          }
          tryParentAlign=function(fam,selectedMarkers,parents,idsPr,children,idsCh){
            changeMarkers=function(fam,selectedMarkers,parents,idsPr,children,idsCh,parentsOrder=1:2){
              fam$reImpAlls[,selectedMarkers,]=abind(parents[,,parentsOrder],children[parentsOrder,,],along=3)
              fam$IDS[,selectedMarkers,]=abind(idsPr[,,parentsOrder],idsCh[parentsOrder,,],along=3)
              return(fam)
            }
            f1=f2=fam
            f1=changeMarkers(fam,selectedMarkers,parents,idsPr,children,idsCh)
            f2=changeMarkers(fam,selectedMarkers,parents,idsPr,children,idsCh,2:1)
            f1=updateParentsIDS(f1)
            f2=updateParentsIDS(f2)
            compat1=!any(f1$idHaps==3)
            compat2=!any(f2$idHaps==3)
            if (compat1&!compat2){
              fam=findHaps(f1)
            } else if (compat2&!compat1){
              fam=findHaps(f2)
            } else if (compat1&compat2){
              fam1=findHaps(f1)
              fam2=findHaps(f2)
              for (k in 1:nChildren){
                nam=which(apply(fam$reImpAlls[,selectedMarkers,k+2],2,function(x) any(is.na(x))))
                if (length(nam)>0){
                  for (nn in nam){
                    nalle=which(is.na(fam$reImpAlls[,selectedMarkers[nn],k+2]))
                    nallefam1=fam1$reImpAlls[,selectedMarkers[nn],k+2]
                    nallefam2=fam2$reImpAlls[,selectedMarkers[nn],k+2]
                    commonAll=nallefam1[which(nallefam1%in%nallefam2)]
                    commonAll=commonAll[!is.na(commonAll)]
                    if (length(nalle)==2){
                      if (length(commonAll)==2) fam$reImpAlls[,selectedMarkers[nn],k+2]=commonAll
                      else if (length(commonAll)==1) fam$reImpAlls[1,selectedMarkers[nn],k+2]=commonAll
                    }
                    else{ #length(nalle)==1
                      if (length(commonAll)==1){
                        if (commonAll!=fam$reImpAlls[3-nalle,selectedMarkers[nn],k+2])
                          fam$reImpAlls[nalle,selectedMarkers[nn],k+2]=commonAll
                      } else if (length(commonAll)==2){
                        if (commonAll[1]==commonAll[2]) fam$reImpAlls[nalle,selectedMarkers[nn],k+2]=commonAll[1]
                        else{
                          wo=which(commonAll!=fam$reImpAlls[3-nalle,selectedMarkers[nn],k+2])
                          fam$reImpAlls[nalle,selectedMarkers[nn],k+2]=commonAll[wo]
                        }
                      }
                    }
                    if(!any(is.na(fam$reImpAlls[,selectedMarkers[nn],k+2]))) 
                      if (fam$reImpAlls[1,selectedMarkers[nn],k+2]==fam$reImpAlls[2,selectedMarkers[nn],k+2])
                        fam$IDS[,selectedMarkers[nn],k+2]=c(1,1)
                  }
                }
              }
              fam$alignedParents=FALSE
            } else {
              fam$hapIncidences <- "Genotyping error, recombination event or inheritance from non-declared parent"
            }
            return(fam)
          }
          
          if (!exists("alignedParents",fam)){
            if (all(is.na(fam$reImpAlls[,,1:2]))) fam$alignedParents=FALSE
            else{
              idsM1=ifelse(apply(fam$IDS[,,-(1:2)],2,sum)>0,TRUE,FALSE)  # Markers with ids=1 in some sib
              c5=apply(fam$reImpAlls[,,-(1:2)],2,function(x){
                x=array(x,dim=c(2,1,fam$famSize-2))
                isCase5(x)}
              )
              # Determination of which haps come from father and which from mother (which markers can be used for alignment of parents)
              alinMk=which(!apply(fam$reImpAlls[,,1:2],2,function(M) all(compareNA(M[,1],M[,2])))&idsM1&!c5)  
              fam$alignedParents=if (length(alinMk)>0) TRUE else FALSE
              fam$alignmentMks=if (fam$alignedParents) list(alinMk=alinMk,reImpAlls=fam$reImpAlls[,alinMk,,drop=FALSE],
                                                            IDS=fam$IDS[,alinMk,,drop=FALSE]) else NULL
            }
          }
          nChildren=fam$famSize-2
          alleMat=t(apply(fam$reImpAlls[,,-c(1,2)],3,function(x) apply(x,2,function(x) any(is.na(x))))+0)
          # The sum by rows indicates the number of missing markers on each child      
          combChildren=combn(nChildren,3) # Search for marker combinations for three full children (at least)
          scm=apply(combChildren,2,selChMk,fam$reImpAlls[,,-c(1,2)])  
          if (!is.null(scm)) scm=scm[-which(sapply(scm,is.null))]
          nmk=sapply(scm,function(cm){length(cm$markers)})
          if (all(nmk<=1)) return(fam)
          
          sl=order(nmk,decreasing = TRUE)
          sl=sl[nmk[sl]>1]
          selSet=NULL
          selMks=NULL
          for (k in sl){
            nwm=which(!scm[[k]]$markers%in%selMks)
            if (length(nwm)>0){
              selSet=c(selSet,k)
              selMks=c(selMks,scm[[k]]$markers[nwm])
            }
          }
          for (s in selSet) if(any(is.na(fam$reImpAlls[,scm[s][[1]]$markers,]))){
            selected=scm[s][[1]]
            idsCh=fam$IDS[,selected$markers,-(1:2)]
            children=fam$reImpAlls[,selected$markers,-c(1,2)]
            partialFam=findHapsfrom3or4Children(fam$reImpAlls[,selected$markers,selected$children+2],
                                                parentsInfo=fam$reImpAlls[,selected$markers,1:2])
            if (!is.null(partialFam$sortedParents)){
              children[,,selected$children]=partialFam$sortedChildren
              idsCh[,,selected$children]=partialFam$childrenIDS
              if (length(selected$children)<nChildren){
                unSelectedChildren=(1:(fam$famSize-2))[-selected$children]
                uHMZ=t(apply(children[,,unSelectedChildren],3,
                             function(x) apply(x,2,function(a) as.integer(a[1]==a[2]))))
                if (nrow(uHMZ)==1) uHMZ=t(uHMZ) 
                uHMZ[is.na(uHMZ)]=0
                for (u in 1:length(unSelectedChildren)) idsCh[,,unSelectedChildren[u]]=rbind(uHMZ[u,],uHMZ[u,])
                partFam=list(reImpAlls=abind(partialFam$sortedParents,children,along=3),
                             IDS=abind(partialFam$parentsIDS,idsCh,along=3))
                partFam$famSize=fam$famSize
                partFam=findHaps(partFam)
                children[,,-selected$children]=partFam$reImpAlls[,,unSelectedChildren+2]
                idsCh[,,-selected$children]=partFam$IDS[,,unSelectedChildren+2]
              }
              if (partialFam$alignedParents){
                fam$reImpAlls[,selected$markers,1:2]=partialFam$sortedParents
                fam$reImpAlls[,selected$markers,-(1:2)]=children
                fam$IDS[,selected$markers,-(1:2)]=idsCh
                fam=updateParentsIDS(fam)
              } else {
                if (!is.null(fam$alignmentMks)) 
                  fam=tryParentAlign(fam,selected$markers,partialFam$sortedParents,partialFam$parentsIDS,children,idsCh)
                if (is.null(fam$alignmentMks)|!fam$alignedParents){
                  oldCh=fam$reImpAlls[,,-(1:2)]
                  oldCh[fam$IDS[,,-(1:2)]==0]=NA
                  newCh=children
                  newCh[idsCh==0]=NA
                  newCh1<-newCh2<-oldCh
                  newCh1[,selected$markers,]=newCh
                  newCh2[,selected$markers,]=newCh[2:1,,]
                  # Checking of those markers with ids=1 which are common in oldCh and newCh
                  ids1old=t(apply(fam$IDS[,selected$markers,-(1:2)],3,function(M) ifelse(colSums(M)>0,1,0)))
                  ids1new=t(apply(idsCh,3,function(M) ifelse(colSums(M)>0,1,0)))
                  ids1comm=ifelse(ids1old+ids1new==2,1,0)
                  # If any is homozygous, it can be useful to the alignment
                  ids1comm=which(ids1comm==1,arr.ind=TRUE)
                  het=apply(ids1comm,1,function(v) children[1,v[2],v[1]]!=children[2,v[2],v[1]])
                  hetMk=which(het)
                  if (length(hetMk)>0){
                    hcm=ids1comm[hetMk[1],]
                    child0=fam$reImpAlls[,selected$markers,hcm[1]+2]
                    gmk0=child0[,hcm[2]]
                    gmk=children[,hcm[2],hcm[1]]
                    if (all(gmk==gmk0)) {
                      fam$reImpAlls[,selected$markers,-(1:2)]=children
                      fam$IDS[,selected$markers,-(1:2)]=idsCh
                      fam$reImpAlls[,selected$markers,(1:2)]=partialFam$sortedParents
                    } else{
                      fam$reImpAlls[,selected$markers,-(1:2)]=children[2:1,,]
                      fam$IDS[,selected$markers,-(1:2)]=idsCh[2:1,,]
                      fam$reImpAlls[,selected$markers,(1:2)]=partialFam$sortedParents[,,2:1]
                    }
                    fam=updateParentsIDS(fam)
                  } 
                  else if (all(compareNA(newCh1,newCh2[2:1,,]))){ 
                    # If there are homozygous and there are not common ones which lead to problems
                    fam$reImpAlls[,selected$markers,1:2]=partialFam$sortedParents
                    fam$reImpAlls[,selected$markers,-(1:2)]=children
                    fam$IDS[,selected$markers,-(1:2)]=idsCh
                    unSelectedMarkers=(1:ncol(fam$IDS))[-selected$markers]
                    for (u in unSelectedMarkers){
                      cma=intersect(fam$reImpAlls[,u,1],fam$reImpAlls[,u,2])
                      if (length(cma)==1) cma=c(cma,NA) else if (length(cma)==0) cma=c(NA,NA)
                      fam$reImpAlls[,u,1]<-fam$reImpAlls[,u,2]<-cma
                    }
                    fam$alignedParents=FALSE  # This means that parents will be assigned haplotypes 
                    # but it is not known which haps really correspond to father and which correspond to mother.
                    fam=updateParentsIDS(fam)
                  } 
                  else{
                    cpr=compareNA(newCh1,newCh2[2:1,,])
                    fam$reImpAlls[,selected$markers,-(1:2)]=children
                    fam$IDS[,selected$markers,-(1:2)]=idsCh
                    fam$IDS[,,-(1:2)]=(cpr+0)*fam$IDS[,,-(1:2)]
                    fch=fam$reImpAlls[,,-(1:2)]
                    ids0=fam$IDS[,,-(1:2)]==0
                    fch[ids0]=NA
                    posParHaps=list(unique(t(fch[1,,])),unique(t(fch[2,,]))) # Possible parental haplotypes
                    phaps=lapply(posParHaps,find2Haps)
                    for (i in 1:2) fam$reImpAlls[,,i]=phaps[[i]]
                    fam=updateParentsIDS(fam)
                  }
                }
              }
              fam=findHaps(fam)
            }
          }
          if (!fam$alignedParents){
            fam$unAlignedParentsHaps=fam$reImpAlls[,,1:2]
            dimnames(fam$unAlignedParentsHaps)[[3]]=c("Parent1","Parent2")
            cnf=list(fam$reImpAlls[,,1]==fam$reImpAlls[,,2],
                     fam$reImpAlls[,,1]==fam$reImpAlls[2:1,,2])
            fid=fam$IDS[,,1]*fam$IDS[,,2]
            cnf=lapply(cnf,function(id) {id[is.na(id)]=FALSE; id*fid})
            idsC=cnf[[which.max(sapply(cnf,sum))]]
            phParAlls=fam$reImpAlls[,,1]
            phParAlls[!idsC]=NA
            fam$reImpAlls[,,1:2]=phParAlls
            fam$IDS[,,1:2]=idsC
            fam$haps=writeFamilyHaps(fam)
          }
          return(fam)
        }
        na0=sum(fam$IDS)
        repeat{
          fam <- tryCatch(makeFamilyFromChildren(fam),error=function(e) {
            fam$incidences="Some error in processing data"
            return(fam)
          })
          na1=sum(fam$IDS)
          if (na1==na0) break
          na0=na1
        }
        return(fam)
      }
      makeHapsFromTwoChildren=function(fam){
        # Function to identify in which case are the parent-offspring haplotypes
        identifyCase=function(fam,sortedAlleles,childrenPair){
          group=fam$reImpAlls[,sortedAlleles,c(1,2,childrenPair+2),drop=FALSE]
          haps=matrix(NA,nrow=4,ncol=2)
          for (i in 1:4) for (j in 1:2) haps[i,j]=writeHap(group[j,,i],1)
          haps=matrix(haps,ncol=2)
          colnames(haps)=c("hap1","hap2")
          commChHap=intersect(haps[3,],haps[4,])         # Do the children have one common haplotype?
          ncch=length(commChHap)
          nChHaps=length(unique(as.vector(haps[3:4,])))  # How many unique haplotypes have the children?
          commParHap=intersect(haps[1,],haps[2,])        # Do the parents have one common haplotype?
          ncph=length(commParHap)
          nParHaps=length(unique(as.vector(haps[1:2,]))) # How many unique haplotypes have the parents?
          caso=0
          if (ncph==0){
            if (nChHaps==2&nParHaps==3){
              if (haps[1,1]==haps[1,2]) caso=3
              else if(haps[2,1]==haps[2,2]) caso=4
            } else if (nChHaps==3){
              if (nParHaps==4){
                if (commChHap%in%haps[1,]) caso=1
                else caso=2
              }
            }
          } else if (ncph==1){
            uch=unique(as.vector(haps[3:4,]))
            fcomm=length(intersect(uch,haps[1,]))
            mcomm=length(intersect(uch,haps[2,]))
            if (fcomm==2&mcomm==1&haps[2,1]!=haps[2,2]) caso=6
            else if (mcomm==2&fcomm==1&haps[1,1]!=haps[1,2]) caso=5
          }
          return(caso)
        }
        
        # Identification of such markers containing IDS=1 and missing parental values 
        sortedParentsMarkers=which(apply(fam$IDS[,,1:2],2,sum)==4)
        if (is.null(sortedParentsMarkers)) return(fam)
        
        # Markers with two heterozygous children (at least) and both parents with missing data
        fam$HMZ=t(apply(fam$reImpAlls,3,function(x) apply(x,2,function(a) as.integer(a[1]==a[2]))))
        twoHetChMkr=which(apply(fam$HMZ[-(1:2),],2,function(x) length(which(x==0)))>1)
        naParents=which(apply(fam$HMZ[1:2,],2,function(x) if (any(is.na(x))) TRUE else FALSE))
        mk2ch=intersect(twoHetChMkr,naParents) # Markers where the makeHapsFromTwoChildren can be applied
        
        # Search for those sibs pairs containing one common allele and two alleles from each other
        imputedMarkers=NULL
        for (mk in mk2ch){
          chPairs=combn(which(fam$HMZ[-(1:2),mk]==0),2)
          selec=which(apply(chPairs,2,function(cp){
            alleles=fam$reImpAlls[,mk,cp+2,drop=FALSE]
            (length(unique(as.vector(alleles)))==3)&(length(intersect(alleles[,,1],alleles[,,2]))==1)
          }))
          if (length(selec)>0){
            chp=chPairs[,selec,drop=FALSE] # Heterozygous children pairs in this marker with a common allele
            for (s in 1:length(selec)){
              sortedAlleles=sortedParentsMarkers[which(apply(fam$IDS[,sortedParentsMarkers,chp[,s]+2],2,sum)==4)]
              if (!is.null(sortedAlleles)){
                case=identifyCase(fam,sortedAlleles,chp[,s])
                if (case>0){
                  chAlle=fam$reImpAlls[,mk,chp[,s]+2,drop=FALSE]
                  commAllele=intersect(chAlle[,,1],chAlle[,,2])
                  if (case%in%c(1,4,5)) ref=2 else ref=1
                  for (i in 1:2) if (chAlle[ref,,i]==commAllele) chAlle[,,i]=chAlle[2:1,,i]
                  fam$reImpAlls[,mk,chp[,s]+2]=chAlle
                  fam$IDS[,mk,chp[,s]+2]=1
                  fam=findHaps(fam) 
                  imputedMarkers=c(imputedMarkers,mk)
                  break
                }
              }
            }
          }
        }
        return(fam)
      }
      whichCase5Mkrs=function(fam){ # Identifies those markers in case 5 from fam$markers matrix
        isCase5=function(mrkPos){
          mrk=fam$imputedAlls[,mrkPos]
          if (all(is.na(mrk))) return(TRUE)
          globAlNr=length(unique(as.vector(mrk[!is.na(mrk)])))  # Global number of alleles in marker
          indivAlNr=apply(mrk,1,function(x) if(any(is.na(x))) 2 else length(unique(x)))
          if ((globAlNr==2) & all(indivAlNr==2)) return(TRUE) else return(FALSE)
        }
        mrksPos=matrix(1:ncol(fam$imputedAlls),nrow=2)
        case5Mks=which(apply(mrksPos,2,isCase5))
        return(list(mrk=case5Mks,mrkCols=mrksPos[,case5Mks]))
      }
      writeHap=function(hap,ids){
        hap[ids==0]=NAsymbol
        return(paste(hap,collapse=alleSep))
      }
      writeFamilyHaps=function(fam){  
        haps=matrix(NA,nrow=fam$famSize,ncol=2)
        for (i in 1:fam$famSize)
          for (j in 1:2)
            haps[i,j]=writeHap(fam$reImpAlls[j,,i],fam$IDS[j,,i])
          haps=matrix(haps,ncol=2)
          colnames(haps)=c("hap1","hap2")
          return(haps)
      }
      
      #=== INITIALIZATION ===#
      fam <- list(infoFam=family[,1:6],imputedAlls=family[,-(1:6)],
                  nMarkers=(ncol(family)-6)/2, famSize=nrow(family))  
      fam$reImpAlls <- as.familyArray(fam$imputedAlls)
      
      #=== CASE 5 CHECKING ===#
      case5Mks <- whichCase5Mkrs(fam)
      if (length(case5Mks$mrk)>0){
        fam0 <- fam
        fam <- excludeCase5Mkrs(fam,case5Mks)
      }
      
      #=== HAPLOTYPING ===#
      if (fam$nMarkers==0){
        haps=paste(rep(NAsymbol,fam0$nMarkers),collapse=alleSep)
        haps=matrix(haps,ncol=2,nrow=fam0$famSize)
        colnames(haps)=c("hap1","hap2")
        reImputedAlls=as.matrix(fam0$imputedAlls)
        dimnames(reImputedAlls)=NULL
        outFam=list(reImputedAlls=reImputedAlls, haplotypes=cbind(fam$infoFam,haps), 
                       IDS=array(0,dim=c(fam0$famSize,fam0$nMarkers*2)),
                       haploLog=data.frame(famID=fam0$infoFam[1,1],hapIncidences="Not enough informative markers"), 
                                           missingParentsLog=NULL)
      } 
      else if (fam$nMarkers==1){
        fam=initializeIDS(fam)
        fam$haps=writeFamilyHaps(fam)
        if (length(case5Mks$mrk)>0) fam=includeCase5Mkrs(fam,fam0,case5Mks)
        outFam=list(reImputedAlls=as.familyMatrix(fam$reImpAlls), haplotypes=cbind(fam$infoFam,fam$haps), 
                       IDS=as.familyMatrix(fam$IDS), haploLog=NULL, missingParentsLog=NULL)
      }
      else {
        fam=initializeIDS(fam)
        fam=findHaps(fam)
        missingParentsMarkers=apply(fam$reImpAlls[,,1:2],2,function(x) any(is.na(x)))
        if (any(missingParentsMarkers)) if (fam$famSize>4) fam=makeHapsFromThreeChildren(fam)  # Tries to solve for three or more children
        missingParentsMarkers=apply(fam$reImpAlls[,,1:2],2,function(x) any(is.na(x)))
        if (any(missingParentsMarkers)) if (fam$famSize>3) fam=makeHapsFromTwoChildren(fam) # Tries to solve for two children
        # Two children function should be applied after the three children function, 
        # since it solves some cases that the other may have left behind
        if (length(case5Mks$mrk)>0) fam=includeCase5Mkrs(fam,fam0,case5Mks)
        haploLog <- if (exists("hapIncidences",fam)) data.frame(famID=fam$infoFam[1,1], hapIncidences=fam$hapIncidences) else NULL
        if (exists("alignedParents",fam)){
          if (!fam$alignedParents&!is.null(fam$unAlignedParentsHaps)){ 
            if (length(case5Mks$mrk)>0) {
              uaph=array(NA,dim=dim(fam$reImpAlls[,,1:2]))
              uaph[,-case5Mks$mrk,]=fam$unAlignedParentsHaps
              fam$unAlignedParentsHaps=uaph
            }
            missingParentsLog <- data.frame(famID=fam$infoFam[1:2,1],indID=c("Parent1","Parent2"),
                                            as.familyMatrix(fam$unAlignedParentsHaps))
          } else missingParentsLog <- NULL
        } else missingParentsLog <- NULL
        outFam=list(reImputedAlls=as.familyMatrix(fam$reImpAlls), haplotypes=cbind(fam$infoFam,fam$haps), 
                       IDS=as.familyMatrix(fam$IDS), haploLog=haploLog, missingParentsLog=missingParentsLog)
      }
      return(outFam)
    }
    # Haplotyping multiple families
    famsHaplotyper=function(fams){  
      #=== INITIALIZATION ===#
      mrksPos <- matrix(1:(ncol(fams$imputedMkrs)-6),nrow=2)
      nMrks <- ncol(mrksPos)
      idFams <- unique(fams$imputedMkrs$famID)
      fams$reImputedAlls <- NULL 
      fams$IDS <- NULL
      fams$haplotypes <- NULL                                   
      fams$haploLog <- NULL
      fams$missingParentsLog <- NULL                                 
      fams$Vdata <- NULL                                        
     
      #=== FAMILIES SCAN ===#
      for (f in idFams) {
        family <- subset(fams$imputedMkrs,fams$imputedMkrs$famID==f)
        fam <- famHaplotyper(family)
        fams$reImputedAlls <- rbind(fams$reImputedAlls,fam$reImputedAlls)                            
        fams$IDS <- rbind(fams$IDS,fam$IDS)
        fams$haplotypes <- rbind(fams$haplotypes,fam$haplotypes)                             
        if (!is.null(fam$haploLog)) fams$haploLog=rbind(fams$haploLog,fam$haploLog)
        if (!is.null(fam$missingParentsLog)) 
          fams$missingParentsLog=rbind(fams$missingParentsLog,fam$missingParentsLog)
        idsCount=c(rowSums(fam$IDS[,mrksPos[1,]]),rowSums(fam$IDS[,mrksPos[2,]]))
        nHaps=2*nrow(family)
        nFullHaps=length(which(idsCount==nMrks))
        nEmptyHaps=length(which(idsCount==0))
        nPartialHaps=nHaps-nFullHaps-nEmptyHaps
        initialMissingAlls=length(which(is.na(family[,-(1:6)])))
        nImputedAlls=initialMissingAlls-length(which(is.na(fam$reImputedAlls)))
        fams$Vdata <- rbind(fams$Vdata,
                            data.frame(famID=family$famID[1],nHaps=nHaps,nFullHaps=nFullHaps,
                                       nEmptyHaps=nEmptyHaps,nPartialHaps=nPartialHaps,
                                       initialMissingAlls,nImputedAlls=nImputedAlls))
      }
        
      #=== FAMILIES BOUND ===# 
      fams$reImputedAlls <- cbind(fams$imputedMkrs[,1:6],fams$reImputedAlls)          
      colnames(fams$IDS) <- names(fams$reImputedAlls)[-(1:6)] <- paste("Mk",gl(nMrks,2),"_",gl(2,1,nMrks*2),sep="")
      fams$haplotypes <- data.frame(fams$haplotypes)
      return(fams)
    }  
    # Data Summary 
    summarizeData=function(fams){
      nAlls <- nrow(fams$IDS)*ncol(fams$IDS)                          # Total number of alleles
      pPhasAlls <- sum(fams$IDS)/nAlls                                # Proportion of IDentified/Sorted alleles
      pNonPhasAlls <- 1-pPhasAlls                                     # Proportion of non-IDentified/Sorted alleles
      nHapsStats <- colSums(fams$Vdata[,-1])                          # Haplotyping Stats
      nHaps <- nrow(fams$haplotypes)*2                                # Number of haplotypes
      pFullHaps <- nHapsStats["nFullHaps"]/nHaps                      # Proportion of full reconstructed haplotypes
      pEmptyHaps <- nHapsStats["nEmptyHaps"]/nHaps                    # Proportion of empty haplotypes
      pPartialHaps <- nHapsStats["nPartialHaps"]/nHaps                # Proportion partially reconstructed haplotypes
      newImputedAlleles <- nHapsStats["nImputedAlls"]                 # New imputed alleles 
      totalImputedAlleles <- nHapsStats["nImputedAlls"]+fams$imputationSummary$nImputedAlls
      reImputationRate <- as.vector(totalImputedAlleles/fams$imputationSummary$nMissAlls)
      if (is.na(reImputationRate)) reImputationRate <- 0
      fams$haplotypingSummary <- data.frame(nAlls,pPhasAlls,pNonPhasAlls,nHaps,pFullHaps,pEmptyHaps,pPartialHaps,
                                            newImputedAlleles,reImputationRate,haplotypingTime) # Haplotyping summary
      fams[["Vdata"]] <- NULL; row.names(fams[["haplotypingSummary"]]) <- NULL
      return(fams)
    }
  }
  ############################################## II. RE-IMPUTATION ##############################################
  {
    # Preliminary imputation (marker by marker)
    fams <- alleImputer(data,invisibleOutput=invisibleOutput, dataSummary=dataSummary)
  }
  ############################################### III. HAPLOTYPING ##############################################
  {
    # Haplotyping and time computing
    haplotypingTime <- system.time(fams <- famsHaplotyper(fams))[[3]] 
  }
  ############################################### IV. DATA SUMMARY  #############################################
  {
    fams <- summarizeData(fams)
  }
  ###############################################  V. DATA STORING  #############################################
  {
    if (is.character(data)) {     
      baseName <- file_path_sans_ext(data)
      outName1 <- paste(baseName,"reImputed.ped",sep="_")
      outName2 <- paste(baseName,"haplotypes.txt",sep="_")
      write.table(fams$reImputedAlls,sep=" ",quote=FALSE,file=outName1)
      if (any(fams$hapIncidences[,-1]!="")) {
        outName3 <- paste(baseName,"hapIncidences.txt",sep="_")
        write.table(fams$hapIncidences,sep=" ",quote=FALSE,file=outName3)     
      }
    }    
  }
  ############################################## VI. FUNCTION OUTPUT  ###########################################
  {
    if (dataSummary==TRUE) {
      
      ### Data summary printing
      cat("\n========= HAPLOTYPING SUMMARY ==========")
      cat(paste("\nRe-imputation rate:",round(fams$haplotypingSummary$reImputationRate,4)))
      cat(paste("\nProportion of phased alleles:",round(fams$haplotypingSummary$pPhasAlls,4)))
      cat(paste("\nProportion of non-phased alleles:",round(fams$haplotypingSummary$pNonPhasAlls,4)))
      cat(paste("\nProportion of missing haplotypes:",round(fams$haplotypingSummary$pEmptyHaps,4)))
      cat(paste("\nProportion of partial haplotypes:",round(fams$haplotypingSummary$pPartialHaps,4)))
      cat(paste("\nProportion of full haplotypes:",round(fams$haplotypingSummary$pFullHaps,4)))
      cat(paste("\nHaplotyping time:",round(fams$haplotypingSummary$haplotypingTime,4)))
      cat("\n========================================\n")
      
      ### Haplotyping message printing
      if (!is.null(fams$haploLog)) {
        cat("Incidences were detected. Some haplotypes could not be generated properly.")
      }
      
      ### Storage path
      if (is.character(data)) {
        cat("\nGenerated data have been stored in:")
        cat(paste("\n",outName1,"\n",sep="")); cat(outName2)
        if (any(fams$hapIncidences[,2]!="")) cat(paste("\n",outName3,sep=""))
      } 
    }     
    
    ### Cleaning empty haplotype Incidences
    if (is.null(fams$haploLog)) fams$haploLog <- NULL

    ### Returning (or not) the output
    if (invisibleOutput) return(invisible(fams)) else return(fams) 
  }

}