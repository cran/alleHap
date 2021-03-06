% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alleImputer.R
\name{alleImputer}
\alias{alleImputer}
\title{Imputation of missing alleles from a dataset composed by families.}
\usage{
alleImputer(data, invisibleOutput = TRUE, dataSummary = TRUE)
}
\arguments{
\item{data}{Data containing the families' identifiers and the corresponding genetic data (or the path of the PED file).}

\item{invisibleOutput}{Data are not shown by default.}

\item{dataSummary}{A summary of the data is shown by default.}
}
\value{
Imputed markers, Homozygosity (HMZ) matrix, marker messages and number of unique alleles per marker.
}
\description{
By analyzing all possible combinations of a parent-offspring pedigree in which parental and/or offspring genotypes may be missing; as long as one child was genotyped, in certain cases it is possible an unequivocal imputation of the missing genotypes both in parents and children.
}
\examples{

## Imputation of families containing parental missing data
simulatedFams <- alleSimulator(10,4,6,missParProb=0.2) 
famsAlls <- simulatedFams[[1]]       # Original data 
alleImputer(famsAlls)                # Imputed alleles (genotypes)

## Imputation of families containing offspring missing data
datasetAlls <- alleSimulator(10,4,6,missOffProb=0.2)
famsAlls <- simulatedFams[[1]]       # Original data 
alleImputer(famsAlls)                # Imputed alleles (genotypes)

## Imputation of a family marker containing missing values in one parent and one child
infoFam <- data.frame(famID="FAM03",indID=1:5,patID=c(0,0,1,1,1),
                      matID=c(0,0,2,2,2),sex=c(1,2,1,2,1),phenot=0)
mkr <- rbind(father=c(NA,NA),mother=c(1,3),child1=c(1,1),child2=c(2,3),child3=c(NA,NA))
colnames(mkr) <- c("Mkr1_1","Mkr1_2")
famMkr <- cbind(infoFam,mkr)         # Original data 
alleImputer(famMkr)                  # Imputed alleles (genotypes)

}
\references{
Medina-Rodriguez, N. Santana A. et al. (2014) alleHap: an efficient algorithm to reconstruct zero-recombinant haplotypes from parent-offspring pedigrees. BMC Bioinformatics, 15, A6 (S-3).
}
