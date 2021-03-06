% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alleSimulator.R
\name{alleSimulator}
\alias{alleSimulator}
\title{Simulation of genetic data (alleles) and non-genetic data (family identifiers)}
\usage{
alleSimulator(nFams = 2, nChildren = NULL, nMarkers = 3,
  numAllperMrk = NULL, chrAlleles = TRUE, nHaplos = 1200,
  missParProb = 0, missOffProb = 0, ungenotPars = 0, ungenotOffs = 0,
  phenProb = 0.2, recombRate = 0, invisibleOutput = TRUE)
}
\arguments{
\item{nFams}{Number of families to generate (integer: 1..1000+)}

\item{nChildren}{Number of children of each family (integer: 1..7 or NULL)}

\item{nMarkers}{Number of markers or allele pairs to generate (integer: 1..1000+)}

\item{numAllperMrk}{Number of different alleles per marker (vector or NULL)}

\item{chrAlleles}{Should alleles be expressed as characters A,C,G,T ? (boolean: FALSE, TRUE)}

\item{nHaplos}{Number of different haplotypes in the population (numeric)}

\item{missParProb}{Probability of parents' missing genotype (numeric: 0..1)}

\item{missOffProb}{Probability of offspring' missing genotype (numeric: 0..1)}

\item{ungenotPars}{Proportion of ungenotyped parents (numeric: 0..1)}

\item{ungenotOffs}{Proportion of ungenotyped offspring (numeric: 0..1)}

\item{phenProb}{Phenotype probability, e.g. being affected by disease (numeric: 0..1)}

\item{recombRate}{Recombination rate (numeric: 0..1)}

\item{invisibleOutput}{Data are not shown by default.}
}
\value{
Families' genotypes and haplotypes.
}
\description{
Data simulation can be performed taking into account many different factors such as number of families to generate, number of markers (allele pairs), number of different alleles per marker, type of alleles (numeric or character), number of different haplotypes in the population, probability of parent/offspring missing genotypes, proportion of missing genotypes per individual, probability of being affected by disease and recombination rate.
}
\examples{

## Generation of 5 simulated families with 2 children per family and 10 markers
simulatedFams <- alleSimulator(5,2,10)   # List with simulated alleles and haplotypes
simulatedFams[[1]]                       # Alleles (genotypes) of the simulated families
simulatedFams[[2]]                       # Haplotypes of the simulated families

}
\references{
Medina-Rodriguez, N. Santana A. et al. (2014) alleHap: an efficient algorithm to reconstruct zero-recombinant haplotypes from parent-offspring pedigrees. BMC Bioinformatics, 15, A6 (S-3).
}
