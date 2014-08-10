# Selects n different haplotypes between the total number of possible haplotypes
simHapSelection=function(n,numAllperMarker){        
  tnh=prod(numAllperMarker)                 # Total number of possible haplotypes
  if (n>tnh) n=tnh                          # nh is truncated if it is greater than the number of possible haplotypes
  if (tnh<1e7) hh=sample(0:(tnh-1),n)       # Randomly selects n haplotypes indexes
  else {                                    # This is really only needed for R version prior to 3.0.0. Previous versions could not sample for n>1e7
    pwr=(1+floor(log10(tnh)))               
    cvt=10^pwr
    mx=(tnh-1)/cvt
    u=round(runif(n,0,mx),pwr)
    hh=u*cvt
  } 
  asignAl=function(gNumber) {                # Indexes are converted to haplotypes
    (hh%/%prod(numAllperMarker[-(1:gNumber)]))%%numAllperMarker[gNumber]
  }
  1+sapply(1:length(numAllperMarker),asignAl) # Returns a matrix of indexes
}
# Generates n offspring by selecting randomly an haplotype from each parent
simOffspring=function(nOffs,father,mother,recombRate){  
  simChild=function(father,mother){                       # Assigns haplotypes randomly to a child from parents
    rcb=rbinom(2,1,recombRate)                            # Randomly decides if recombination occurs in parents
    if (rcb[1]) father=simRecombHap(father)               # Simulates the recombination of father's haplotypes
    if (rcb[2]) mother=simRecombHap(mother)               # Simulates the recombination of mother's haplotypes
    pos=sample(1:2,2,replace=TRUE)                        # Randomly decides which haplotye from each parent is inherited by the child
    offspring=rbind(father[pos[1],],mother[pos[2],])      # Generates the offspring haplotypes
    haps=apply(offspring,1,paste,collapse="")             # Converts the child's haplotypes in a character vector
    markers=as.vector(apply(offspring,2,sort))            # Sorts the alleles in alphanumeric order
    child=c(markers,sum(rcb),haps)                        # Generates a vector with alleles, haplotypes and number of recombinant haplotypes
    return(child)
  }
  offspring=t(replicate(nOffs,simChild(father,mother)))   # Concatenates the generated children
  return(offspring)
}
# Simulates one family from a population containing the haplotypes 'popHaplos'
simOneFamily=function(famid,popHaplos,pAff,recombRate){     
  dnof=c(0.1,0.45,0.25,0.15,0.044,0.005,0.001)              # Distribution of the number of offspring per family
  #dnof=c(0,0.45,0.3,0.2,0.044,0.005,0.001)                 # Excluding the trios
  nOffs=sample(1:length(dnof),1,prob=dnof)                  # Generates the number of offspring
  nHaplos=nrow(popHaplos)                                   # Number of haplotypes in population
  father=popHaplos[sample(1:nHaplos,2,replace=TRUE),]       # Simulates father's haplotypes
  mother=popHaplos[sample(1:nHaplos,2,replace=TRUE),]       # Simulates mother's haplotypes
  offspring=simOffspring(nOffs,father,mother,recombRate)    # Simulates offspring's haplotypes
  fhaps=sort(apply(father,1,paste,collapse=""))             # Father's haplotypes as character vector
  mhaps=sort(apply(mother,1,paste,collapse=""))             # Mother's haplotypes haplotypes by columns
  fmarkers=as.vector(apply(father,2,sort))                  # Sorts father's alleles alphanumerically
  mmarkers=as.vector(apply(mother,2,sort))                  # Sorts mother's alleles alphanumerically
  parents=rbind(c(fmarkers,0,fhaps),c(mmarkers,0,mhaps),
                deparse.level=0)                            # Array with parents markers and haplotypes
  markers=rbind(parents,offspring,deparse.level=0)          # Array with parents and offspring markers and haplotypes.
  nf=nOffs+2                                                # Number of family members
  indid=1:nf                                                # Each individual is assigned an identification index indid
  sex=c(1,2,sample(1:2,nOffs,replace=TRUE))                 # Random assignation of sex to offsprings
  affected=rbinom(nf,1,pAff)                                # Random assignation of affection status to each subject
  fam=cbind(famid,indid,sex,affected,markers,
            deparse.level=0)                                # Data.frame composed by: famid,indid,sex,affected,markers (ncol: 5..nmarkers*2)
  return(fam)
}
# Simulates the recombination of haplotypes
simRecombHap=function(haps){                           
  hl=ncol(haps)                                        # Number of alleles in haplotypes
  breakPoint=sample(2:hl,1)                            # Position in which recombination takes place
  newHaps=haps                                         # Initialize the new recombinated haplotype as the original haplotype
  newHaps[1,breakPoint:hl]=haps[2,breakPoint:hl]       # Recombination
  newHaps[2,breakPoint:hl]=haps[1,breakPoint:hl]       # Recombination
  return(newHaps)                                      # Returns the recombinated haplotypes
}