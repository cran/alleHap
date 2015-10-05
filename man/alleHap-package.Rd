\name{alleHap-package}
\alias{alleHap-package}
\alias{alleHap}
\docType{package}
\title{
Allele Imputation and Haplotype Reconstruction from Pedigree Databases
}
\description{
Simulation of alphanumeric alleles, imputation of genetic missing data and reconstruction of non-recombinant haplotypes from pedigree databases in a deterministic way. Allelic simulations can be implemented taking into account many different factors (such as number of families, markers, alleles per marker, probability and proportion of missing genotypes, recombination rate, etc). Genotype imputation can be used with simulated datasets or with loaded databases (in .ped file biallelic format). Haplotype reconstruction can be performed even with missing data because of alleHap, firstly imputes the genotypes. All this considering that each member (due to meiosis) should unequivocally have two alleles, one inherited from each parent.
}
\details{
\tabular{ll}{
Package: \tab alleHap\cr
Type: \tab Package\cr
Version: \tab 0.9.2\cr
Date: \tab 2015-10-05\cr
Depends: \tab abind, tools\cr
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

## Haplotype reconstruction of a family containing missing values.
infoFam <- data.frame(famID="FAM001",indID=1:5,patID=c(0,0,1,1,1),
                     matID=c(0,0,2,2,2),sex=c(1,2,1,2,1),phenot=0)
Mkrs <- rbind(c(1,2, NA,NA, 1,2),c(3,4, 1,2, 3,4),
              c(2,4, 1,2, NA,NA),c(NA,NA, 1,1, 2,4),c(1,3, 2,1, 2,3))
colnames(Mkrs) <- c("Mk1_1","Mk1_2","Mk2_1","Mk2_2","Mk3_1","Mk3_2")
family <- cbind(infoFam,Mkrs)
hpList <- hapPhaser(family)

## Haplotype reconstruction of a family containing missing data in a parent. 
infoFam <- data.frame(famID="FAM002",indID=1:6,patID=c(0,0,1,1,1,1),
                     matID=c(0,0,2,2,2,2),sex=c(1,2,1,2,1,2),phenot=0)
Mkrs <- rbind(c(1,4,2,5,3,6),rep(NA,6),c(1,7,2,3,3,2),
              c(4,7,5,3,6,2),c(1,1,2,2,3,3),c(1,4,2,5,3,6))
colnames(Mkrs) <- c("Mk1_1","Mk1_2","Mk2_1","Mk2_2","Mk3_1","Mk3_2")
family <- cbind(infoFam,Mkrs)
hpList <- hapPhaser(family)
}
\references{
Medina-Rodriguez, N. Santana A. et al. (2014) alleHap: an efficient algorithm to reconstruct zero-recombinant haplotypes from parent-offspring pedigrees. BMC Bioinformatics, 15, A6 (S-3).
}
