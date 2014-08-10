\name{alleHap-package}
\alias{alleHap-package}
\alias{alleHap}
\docType{package}
\title{
Allele simulation and Haplotype reconstruction from pedigree databases
}
\description{
This package contains tools for simulating alphanumeric alleles from parent-offspring pedigree databases to later reconstruct the corresponding non-recombinant haplotypes. Allelic simulations can be implemented taking into account many different factors (such as number of families, markers, alleles per marker, probability of missing genotypes, recombination rate, etc). Haplotype reconstruction can be performed even with missing data because of the algorithm is based on a preliminary analysis of all possible combinations that may exist in the genotyping of a family. All this considering that each member (due to meiosis) should unequivocally have two alleles, one from each parent.
}
\details{
\tabular{ll}{
Package: \tab alleHap\cr
Type: \tab Package\cr
Version: \tab 0.7.2\cr
Date: \tab 2014-08-10\cr
Depends: \tab gtools\cr
License: \tab GPL (>=2)\cr
}
}
\author{
Nathan Medina-Rodriguez and Angelo Santana

Maintainer: Nathan Medina-Rodriguez <natan.medina101@alu.ulpgc.es>
}
\examples{
## Generation of 100 simulated families with 12 markers
dataset <- simFams(100,12)         # Creates a list with simulated alleles and haplotypes
datasetAlls <- dataset[[1]]        # Extracts familiar alleles
datasetHaps <- dataset[[2]]        # Extracts familiar haplotypes

## Haplotype Reconstruction of a dataset composed by 100 families and 12 markers 
reconstHaps <- buildHaps(simFams(100,12)[[1]])

## Haplotype Reconstruction of a dataset composed by 100 families and 12 markers
# with missing parental information (Reconstruction Rate may lead a value lower than one, RR<1)
reconstHaps <- buildHaps(simFams(100,12,missParProb=0.1)[[1]])
 
## Haplotype Reconstruction of a dataset composed by 100 families and 12 markers 
# with missing offspring information (Reconstruction Rate may lead a value lower than one, RR<1)
reconstHaps <- buildHaps(simFams(100,12,missOffProb=0.05)[[1]])
}
\references{
Medina-Rodriguez, N. Santana A. et al. (2014) alleHap: an efficient algorithm to reconstruct zero-recombinant haplotypes from parent-offspring pedigrees. BMC Bioinformatics, 15, A6 (S-3).
}
