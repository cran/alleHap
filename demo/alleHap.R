############################################################################
    ##################### alleHap DEMONSTRATION #######################
                 ######################################

# This demo shows the usual workflow of the alleHap in which examples  
# of the main functions usage will be secuentially introduced.

library(alleHap)

readline(prompt="Press [enter] to see the DATA LOADING stage")

############################# Ia. DATA LOADING #############################

## The data to be loaded must be structured in .ped format and families must 
## comprise by parent-offspring pedigrees.
## As default, alleLoader shows a summary about the loaded data.

readline(prompt="Press [enter] to see the loading of a dataset with alphabetical alleles")

# The first example is a dataset with alphabetical alleles (A,C,G,T)
example1 <- file.path(find.package("alleHap"), "examples", "example1.ped")
example1Alls <- alleLoader(example1)
head(example1Alls[,1:12])

readline(prompt="Press [enter] to see the loading of a dataset containing numerical alleles")

## As second example, a dataset with numerical alleles is shown
example2 <- file.path(find.package("alleHap"), "examples", "example2.ped")
example2Alls <- alleLoader(example2)
head(example2Alls)

readline(prompt="Press [enter] to see the DATA SIMULATION stage")

########################## Ib. DATA SIMULATION #############################

## Alternatively, you can use alleHap to simulate biallelic pedigree databases.
## For this purpose you only need the alleSimulator function.
## As default, alleSimulator shows a summary about the simulated data.

readline(prompt="Press [enter] to see the generation of simulated dataset")

## Generation of 5 simulated families with 2 children per family and 10 markers
simulatedFams <- alleSimulator(5,2,10)    # List with simulated alleles and haplotypes
head(simulatedFams[[1]])                  # Alleles (genotypes) of the simulated families
head(simulatedFams[[2]])                  # Haplotypes of the simulated families

readline(prompt="Press [enter] to see the DATA IMPUTATION stage")

########################## II. DATA IMPUTATION #############################

## Once data is properly loaded (or simulated), the imputation of missing alleles 
## from a dataset composed by nuclear families, can be deployed.
## As default, alleImputer does not show a summary about the imputed data.

readline(prompt="Press [enter] to see the imputation of families containing parental missing data")

## Imputation of families containing parental missing data
simulatedFams <- alleSimulator(10, 4, 6, missParProb=0.2) 
famsAlls <- simulatedFams[[1]]       # Original data 
alleImputer(famsAlls)                # Imputed alleles (genotypes)

readline(prompt="Press [enter] to see the imputation of families containing offspring missing data")

## Imputation of families containing offspring missing data
datasetAlls <- alleSimulator(10, 4, 6, missOffProb=0.2)
famsAlls <- simulatedFams[[1]]       # Original data 
alleImputer(famsAlls)                # Imputed alleles (genotypes)

readline(prompt="Press [enter] to see the Imputation of a family marker containing missing values in one parent and one child")

## Imputation of a family marker containing missing values in one parent and one child
infoFam <- data.frame(famID="FAM03",indID=1:5,patID=c(0,0,1,1,1), 
                      matID=c(0,0,2,2,2),sex=c(1,2,1,2,1),phenot=0)
mkr <- rbind(father=c(NA,NA),mother=c(1,3),child1=c(1,1),child2=c(2,3),child3=c(NA,NA))
colnames(mkr) <- c("Mkr1_1","Mkr1_2")
famMkr <- cbind(infoFam,mkr)             # Original data 
alleImputer(famMkr, dataSummary = TRUE)  # Imputed alleles (genotypes)

readline(prompt="Press [enter] to see the DATA HAPLOTYPING stage")

########################## III. DATA HAPLOTYPING #############################

## The final step is the haplotyping of a dataset composed by several families.
## As default, alleHaplotyper does not show a summary about the imputed data.

readline(prompt="Press [enter] to see the haplotype reconstruction for 3 families without missing data")

## Haplotype reconstruction for 3 families without missing data.
simulatedFams <- alleSimulator(3, 3, 6)  
famsList <- alleHaplotyper(simulatedFams[[1]])  # List containing families' alleles and haplotypes
famsList$haplotypes                             # Reconstructed haplotypes

readline(prompt="Press [enter] to see the haplotype reconstruction of a family containing missing data in a parent")

## Haplotype reconstruction of a family containing missing data in a parent. 
infoFam <- data.frame(famID="FAM002",indID=1:6,patID=c(0,0,1,1,1,1),matID=c(0,0,2,2,2,2),
                      sex=c(1,2,1,2,1,2),phenot=c(2,1,1,2,1,2))
Mkrs <- rbind(c(1,4,2,5,3,6),rep(NA,6),c(1,7,2,3,3,2),c(4,7,5,3,6,2),c(1,1,2,2,3,3),c(1,4,2,5,3,6))
colnames(Mkrs) <- c("Mk1_1","Mk1_2","Mk2_1","Mk2_2","Mk3_1","Mk3_2")
family <- cbind(infoFam,Mkrs)      
famList <- alleHaplotyper(family)  # List containing family's alleles and haplotypes
family                             # Original data 
famList$reImputedAlls              # Re-imputed alleles
famList$haplotypes                 # Reconstructed haplotypes

readline(prompt="Press [enter] to see the haplotype reconstruction from the previous loaded PED file")

## Haplotype reconstruction from a PED file
pedFamList <- alleHaplotyper(example1Alls, dataSummary = TRUE)
head(example1Alls[,1:20])              # Original data 
head(pedFamList$reImputedAlls[,1:20])  # Re-imputed alleles 
head(pedFamList$haplotypes)            # Reconstructed haplotypes
