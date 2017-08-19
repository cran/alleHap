#' @title Data loading of nuclear families (in .ped format)
#' @description The data to be loaded must be structured in .ped format and families must comprise by parent-offspring pedigrees.
#' @param data Data to be loaded.
#' @param invisibleOutput Data are not shown by default.
#' @param dataSummary A summary of the data is shown by default.
#' @param missingValues Specification of the character/numerical values which may be missing. 
#' @return Loaded dataset.
#' @import tools
#' @import utils
#' @export alleLoader
#' @references Medina-Rodriguez, N. Santana A. et al. (2014) alleHap: an efficient algorithm to reconstruct zero-recombinant haplotypes from parent-offspring pedigrees. BMC Bioinformatics, 15, A6 (S-3).
#' @examples
#' 
#' ## Loading of a dataset in .ped format with alphabetical alleles (A,C,G,T)
#' example1 <- file.path(find.package("alleHap"), "examples", "example1.ped")
#' example1Alls <- alleLoader(example1)
#' head(example1Alls)
#' 
#' ## Loading of a dataset in .ped format with numerical alleles
#' example2 <- file.path(find.package("alleHap"), "examples", "example2.ped")
#' example2Alls <- alleLoader(example2)
#' head(example2Alls)
#' 
alleLoader <- function(data, invisibleOutput=TRUE, dataSummary=TRUE, missingValues=c(-9,-99)){
  
  ######################################### I. INTERNAL FUNCTIONS #######################################
  {  
    # Recodification of pre-specified missing values as NA values
    recodeNA=function(data){             
      chMissing=c(as.character(missingValues),"NA","<NA>")
      dataClasses=lapply(data,class)
      cnumeric=which(dataClasses=="integer"|dataClasses=="numeric")
      ccharacter=which(dataClasses=="character")
      cfactor=which(dataClasses=="factor")
      if (length(cnumeric)>0) for (k in cnumeric) data[,k][data[,k]%in%missingValues] = NA
      if (length(ccharacter)>0) for (k in ccharacter) data[,k][data[,k]%in%chMissing] = NA
      if (length(cfactor)>0){
        for (k in cfactor) data[,k][data[,k]%in%chMissing]=NA
        data[,cfactor]=apply(data[,cfactor,drop=FALSE],2,function(v) factor(v))
      } 
      return(data)
    }
  }
  ################################ II. FILE EXTENSION CHECK & DATA READ #################################
  {  
    if (is.character(data)) {
      filetype <- file_ext(data)
      if (filetype=="ped") datasetAlls <- read.table(data)
      else {
        cat("===========================================")
        cat("\n===== alleHap package: version",
            sessionInfo("alleHap")$otherPkgs$alleHap$Version,"======")
        cat("\n===========================================")
        cat("\n\nYour data could not be loaded!!!") 
        return(stop("The file must have a .ped extension"))
      }
    } else if (is.data.frame(data)|is.matrix(data)) {
      datasetAlls <- data
    } else {
      cat("===========================================")
      cat("\n===== alleHap package: version",
          sessionInfo("alleHap")$otherPkgs$alleHap$Version,"======")
      cat("\n===========================================")
      cat("\n\nYour data could not be loaded!!!")
      return(stop("Your data must be in data.frame or matrix format"))
    }
  }
  ########################################## III. DATA CHECK ############################################
  {
    ### Variable identification and NA values assignment
    
    famIDs <- datasetAlls[,1]                                   # Family identifiers
    indIDs <- datasetAlls[,2]                                   # Individual identifiers
    patIDs <- datasetAlls[,3]                                   # Paternal identifiers
    matIDs <- datasetAlls[,4]                                   # Maternal identifiers
    sexIDs <- datasetAlls[,5]                                   # Sex identifiers
    phenot <- datasetAlls[,6]                                   # Phenotype values
    
    names(datasetAlls)[1:6]<- c("famID","indID","patID",
                                "matID","sex","phen")           # dataset names
    nMrk <- (ncol(datasetAlls)-6)/2
    names(datasetAlls)[-(1:6)] <- paste("Mk",gl(nMrk,2),"_",gl(2,1,2*nMrk),sep="")
    famAllelesCols <- (1:ncol(datasetAlls))[-(1:6)]             # Markers Alleles columns
    datasetAlls <- recodeNA(datasetAlls)                        # recode -9, -99, "NA" and "<NA>" as NA values
    
    ### Data Counting
    nFams <- length(unique(famIDs))                             # Number of families
    nIndivs <- length(indIDs)                                   # Number of individuals
    nPars <- length(which(indIDs<3))                            # Number of parents (founders)
    nOffs <- length(which(indIDs>2))                            # Number of children (offspring)
    nMals <- length(which(sexIDs%in%c(1,"male")))               # Number of males
    nFems <- length(which(sexIDs%in%c(2,"female")))             # Number of females
    nMrks <- length(famAllelesCols)/2                           # Number of markers
    
    ### Data Ranges
    if (length(unique(famIDs))>1) {                             # If there is more than one family
      minFamID <- min(famIDs); maxFamID <- max(famIDs)          # Range of Family IDs
    }
    minIndID <- min(indIDs); maxIndID <- max(indIDs)            # Range of Individual IDs
    levsPatID <- levels(as.factor(patIDs))                      # Paternal IDs
    levsMatID <- levels(as.factor(matIDs))                      # Maternal IDs
    levsSex <- levels(as.factor(sexIDs))                        # Sex values
    if (nlevels(as.factor(phenot))<=5) {                        # Phenotype values
      levsPhen <- levels(as.factor(phenot))                     # Phenotype levels
    } else { 
      minPhen <- min(phenot); maxPhen <- max(phenot)            # Min and max phenotype values 
    }
  }
  ####################################### IV. MISSING DATA COUNT ########################################
  { 
    noPars  <- 2*nFams-nPars                                        # Missing Founders
    noID    <- length(which(is.na(indIDs)))                         # Subjects without ID number
    noPatID <- length(which(is.na(patIDs)))                         # Unknown Paternal IDs
    noMatID <- length(which(is.na(matIDs)))                         # Unknown Maternal IDs
    noSex   <- length(which(is.na(sexIDs)))                         # Unknown Sex
    noPhen  <- length(which(is.na(phenot)))                         # Unknown Phenotypes
    missAlls  <- length(which(is.na(datasetAlls[,famAllelesCols]))) # Unknown Alleles
    NAmkrs  <- sum(sapply(1:nMrks, function(i) 
      any(is.na(datasetAlls[,2*i+c(5,6)]))))                        # Markers containing missing values
  }
  ######################################## V. FUNCTION OUTPUT ###########################################
  {
    if (dataSummary==TRUE) {
      cat("===========================================")
      cat("\n===== alleHap package: version",
          sessionInfo("alleHap")$otherPkgs$alleHap$Version,"======")
      cat("\n===========================================")
      cat("\n\nData have been successfully loaded from: \n")
      if (is.character(data)) cat(data) else cat(getwd())
      cat("\n\n===== DATA COUNTING ======")
      cat(paste("\nNumber of families:",nFams))
      cat(paste("\nNumber of individuals:",nIndivs))
      cat(paste("\nNumber of founders:",nPars))
      cat(paste("\nNumber of children:",nOffs))
      cat(paste("\nNumber of males:",nMals))
      cat(paste("\nNumber of females:",nFems))
      cat(paste("\nNumber of markers:",nMrks))
      cat("\n===========================")
      cat("\n\n======== DATA RANGES =========")
      if (length(unique(famIDs))==1) cat(paste("\nFamily ID: ",unique(famIDs),sep=""))
      else cat(paste("\nFamily IDs: [",minFamID,",...,",maxFamID,"]",sep=""))
      cat(paste("\nIndividual IDs: [",minIndID,",...,",maxIndID,"]",sep=""))
      cat(paste("\nPaternal IDs: [",paste(levsPatID,collapse=","),"]",sep=""))
      cat(paste("\nMaternal IDs: [",paste(levsMatID,collapse=","),"]",sep=""))
      cat(paste("\nSex values: [",paste(levsSex,collapse=","),"]",sep=""))
      if (nlevels(as.factor(phenot))<=5) 
        cat(paste("\nPhenotype values: [",paste(levsPhen,collapse=","),"]",sep=""))
      else cat(paste("\nRange of Phenotypes: [",minPhen,",...,",maxPhen,"]",sep=""))
      cat("\n==============================")
      cat("\n\n========= MISSING DATA =========")
      cat(paste("\nMissing founders:",noPars))
      cat(paste("\nMissing ID numbers:",noID))
      cat(paste("\nMissing paternal IDs:",noPatID))
      cat(paste("\nMissing maternal IDs:",noMatID))
      cat(paste("\nMissing sex:",noSex))
      cat(paste("\nMissing phenotypes:",noPhen))
      cat(paste("\nMissing alleles:",missAlls))
      cat(paste("\nMarkers with missing values:",NAmkrs))
      cat("\n================================\n")
    }
    if (invisibleOutput) return(invisible(datasetAlls)) else return(datasetAlls)
  }
}