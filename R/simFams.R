#' @title Generation of a group of simulated families
#' @description Database simulation can be performed taking into account many different factors such as number of families to generate, number of markers (allele pairs), number of different alleles per marker, type of alleles (numeric or character), number of different haplotypes in the population, probability of parent/offspring missing genotypes, probability of being affected by disease and recombination rate.
#' @param nfams Number of families to generate (integer: 1..10000+)
#' @param nMarkers Number of markers or allele pairs to generate (integer: 1..7+)
#' @param numAllperMrk Number of different alleles per marker (vector or NULL)
#' @param chrAlleles Should alleles be expressed as characters A,C,G,T ? (boolean: FALSE, TRUE)
#' @param nHaplos Number of different haplotypes in the population (numeric)
#' @param missParProb Probability of parents' missing genotype (numeric: 0..1)
#' @param missOffProb Probability of offspring missing genotype (numeric: 0..1)
#' @param pAff Probability of being affected by disease (numeric: 0..1)
#' @param recombRate Recombination rate (numeric: 0..1)
#' @return family alleles dataframe
#' @export simFams
#' @references Medina-Rodriguez, N. Santana A. et al. (2014) alleHap: an efficient algorithm to reconstruct zero-recombinant haplotypes from parent-offspring pedigrees. BMC Bioinformatics, 15, A6 (S-3).
#' @examples
#' ## Generation of 100 simulated families with 12 markers
#' dataset <- simFams(100,12)         # Creates a list with simulated alleles and haplotypes
#' datasetAlls <- dataset[[1]]        # Extracts familiar alleles
#' datasetHaps <- dataset[[2]]        # Extracts familiar haplotypes
#' 
#' ## Not run: 
#' # nMarkers>=15 and chrAlleles=TRUE   (May cause repetition in the allelic generation)
#' # nMarkers>=25 and chrAlleles=FALSE  (May cause repetition in the allelic generation)
simFams=function(nfams=2,nMarkers=6,numAllperMrk=NULL,chrAlleles=TRUE,nHaplos=1200,missParProb=0,missOffProb=0,pAff=0.2,recombRate=0){  
  if (is.null(numAllperMrk)){                                     # Simulate the number of alleles per marker if they are not supplied by user
    if (!chrAlleles)                                              # Checks if alleles are not character type
      numAllperMrk=sample(1:9,nMarkers,replace=TRUE)              # Assigns the alleles range per marker
    else numAllperMrk=rep(4,nMarkers)                             # Repeats 4 times if alles are character type
  }
  else nMarkers=length(numAllperMrk)                              # Number of markers
  popHaplos=simHapSelection(nHaplos,numAllperMrk)                 # Generate the haplotypes which exist in the population
  if (chrAlleles){                                                # Checks if alleles are character type
    chrMrkLabel=c('A','C','G','T')                                # Creates the 'A','C','G','T' character labels
    popHaplos=apply(popHaplos,2,function(k) chrMrkLabel[k])       # Assigns the 'A','C','G','T' characters labels to population haplotypes
  } 
  #else popHaplos=popHaplos+100                                   # Adds 100 if character type is numeric
  families=sapply(1:nfams,simOneFamily,popHaplos,
                  pAff,recombRate,simplify=FALSE)                 # Simulate one family
  families=do.call(rbind, families)                               # Concatenates the simulated families
  families=data.frame(families)                                   # Saves the generated families as data.frame
  i1=gl(nMarkers,2)                                               # Generates grade levels to name the first allele column
  i2=gl(2,1,2*nMarkers)                                           # Generates grade levels to name the second allele column
  markerCols=4+1:(2*nMarkers)                                     # Creates the markers' columns
  hapCols=ncol(families)-(1:0)                                    # Creates the haplotypes' columns 
  names(families)=c("famid","indid","sex","affected",             # markerrates the names of family columns 
                    paste("Mrk",i1,"_",i2,sep=""),"recombNr",
                    "hap_Parent1","hap_Parent2")   
  toInteger=function(f) as.numeric(levels(f))[f]                  # Converts to integer the factor levels
  toCharacter=function(f) as.character(f)                         # Converts to character the factor levels
  for (k in c(1:4,(ncol(families)-2)))                            # Browses the first 4 columns and recombNr column
    families[[k]]=toInteger(families[[k]])                        # Converts to integer famid, indid, sex, affection and recombNr
  if (chrAlleles) for (k in markerCols)                           # Checks if alleles are of character type
    families[[k]]=toCharacter(families[[k]])                      # Converts to character the markers' columns
  else for (k in markerCols)                                      # Checks if alleles are of numeric type
    families[[k]]=toInteger(families[[k]])                        # Converts to integer the markers' columns
  for (k in hapCols) families[[k]]=as.character(families[[k]])    # Browses the the haplotypes' columns and converts to character
  npLost=rbinom(1,2*nfams,missParProb)                            # Number of parents randomly selected to have missing alleles
  noffLost=rbinom(1,(nrow(families)-2*nfams),missOffProb)         # Number of offspring randomly selected to have missing alleles
  parentsRows=which(families$indid<3)                             # Counts the number of parents
  offspringRows=which(families$indid>2)                           # Counts the number of offspring
  lostParents=sample(parentsRows,npLost)                          # Samples the listed parents to create the missing ones
  lostoffspring=sample(offspringRows,noffLost)                    # Samples the listed offspring to create the missing ones
  lost=c(lostParents,lostoffspring)                               # Saves the lost parents and offspring
  if (length(lost)) families[lost,markerCols] <- NA               # Assigns to 0 values the lost parents and offspring alleles
  datasetAlls=families[-((length(families)-2):length(families))]  # Data.frame composed by: famid,indid,sex,affected,markers
  datasetHaps=families[-(5:(4+nMarkers*2))]                       # Data.frame composed by: famid,indid,sex,affected,recombNr,haplotypes
  return(list(datasetAlls=datasetAlls,datasetHaps=datasetHaps))   # Generates a families' dataset containing alleles and other with haplotypes (list: 2)
}