# GAPIT - Genomic Association and Prediction Integrated Tool
#######################################################################################
library(multtest)
library('gplots')
library(LDheatmap)
library(genetics)
library(MASS)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R library("scatterplot3d")
library("scatterplot3d")
setwd("~/Desktop/del15k/MAF5%/")
source("emma.R")
source("gapit_function_new.R")
#############################################################################################
##############################################################################################
### start GAPIT ###
## Phenotype
pheno=read.csv('phenotype_for_GWAS.csv',sep=',',header=T)
sc.num  <-pheno[,c(1,3)] #number of sclerotia in one plate
sc.size  <-pheno[,c(1,4)] #average size of the sclerotia in one plate
## Genotype: SNP calling by bcftools; maf and missing rate were filtered by vcftools
myG <- read.delim("rs164.del15k.bialle.MAF5%.NA10%.taxa10%.maf5%.hmp.txt", header = F,sep = "\t")
## Calculate Kinship and PCA
myGAPIT <- GAPIT(
  Y=sc.size,
  G=myG,
  PCA.total = 10,
  #PCA.col = PCA.colors$color,
  SNP.MAF = 0.05,
  Model.selection=T,
  output.numerical=T)

#Step 2: Run GAPIT

myGAPIT <- GAPIT(
  Y=sc.size,
  G=myG,
  PCA.total = 2,
  SNP.MAF = 0.05,
  model=c('FarmCPU','BLINK'))
## EMMAxP3D estimates variance component ##
size.Vg=0.00202937037385262
size.Ve=0.00385870721317774
size.Vg/(size.Vg+size.Ve)

num.Vg= 48.7473 
num.Ve= 4.2589
num.Vg/(num.Vg+num.Ve)
