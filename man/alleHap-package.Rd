\name{alleHap-package}
\alias{alleHap-package}
\alias{alleHap}
\docType{package}
\title{
Allele Imputation and Haplotype Reconstruction from Pedigree Databases
}
\description{
Tools to simulate alphanumeric alleles, impute genetic missing data and reconstruct non-recombinant haplotypes from pedigree databases in a deterministic way. Allelic simulations can be implemented taking into account many factors (such as number of families, markers, alleles per marker, probability and proportion of missing genotypes, recombination rate, etc). Genotype imputation can be used with simulated datasets or real databases (previously loaded in .ped format). Haplotype reconstruction can be carried out even with missing data, since the program firstly imputes each family genotype (without a reference panel), to later reconstruct the corresponding haplotypes for each family member. All this considering that each individual (due to meiosis) should unequivocally have two alleles per marker (one inherited from each parent), and thus imputation and reconstruction results can be deterministically calculated.
}
\details{
\tabular{ll}{
Package: \tab alleHap\cr
Type: \tab Package\cr
Version: \tab 0.9.9\cr
Date: \tab 2017-08-19\cr
Depends: \tab abind, stats, tools, utils\cr
License: \tab GPL (>=2)\cr
}
}
\author{
Nathan Medina-Rodriguez and Angelo Santana

Maintainer: Nathan Medina-Rodriguez <nathan.medina@ulpgc.es>
}
\examples{
## Generation of 10 simulated families with 2 children per family and 20 markers
dataset <- alleSimulator(10,2,20)  # List with simulated alleles and haplotypes
datasetAlls <- dataset[[1]]        # Dataset containing alleles
datasetHaps <- dataset[[2]]        # Dataset containing haplotypes

## Loading of a dataset in .ped format with alphabetical alleles (A,C,G,T)
example1 <- file.path(find.package("alleHap"), "examples", "example1.ped")
datasetAlls1 <- alleLoader(example1)

## Loading of a dataset in .ped format with numerical alleles
example2 <- file.path(find.package("alleHap"), "examples", "example2.ped")
datasetAlls2 <- alleLoader(example2)

## Allele imputation of families with parental missing data
datasetAlls <- alleSimulator(10,4,6,missParProb=0.2)[[1]]
famsImputed <- alleImputer(datasetAlls)

## Allele imputation of families with offspring missing data
datasetAlls <- alleSimulator(10,4,6,missOffProb=0.2)[[1]]
famsImputed <- alleImputer(datasetAlls)

## Haplotype reconstruction for 3 families without missing data.
simulatedFams <- alleSimulator(3,3,6)  
(famsAlls <- simulatedFams[[1]])      # Original data 
famsList <- alleHaplotyper(famsAlls)  # List containing families' alleles and haplotypes
famsList$reImputedAlls                # Re-imputed alleles
famsList$haplotypes                   # Reconstructed haplotypes

## Haplotype reconstruction from a PED file
pedFamPath <- file.path(find.package("alleHap"), "examples", "example3.ped") # PED file path
pedFamAlls <- alleLoader(pedFamPath,dataSummary=FALSE) 
pedFamList <- alleHaplotyper(pedFamAlls)
pedFamAlls                # Original data 
pedFamList$reImputedAlls  # Re-imputed alleles 
pedFamList$haplotypes     # Reconstructed haplotypes
}
\references{
Medina-Rodriguez, N. Santana A. et al. (2014) alleHap: an efficient algorithm to reconstruct zero-recombinant haplotypes from parent-offspring pedigrees. BMC Bioinformatics, 15, A6 (S-3).
}
