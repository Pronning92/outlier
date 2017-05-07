#!/usr/bin/env Rscript

# Simple_BoxCox_Outlier_Analysis: BoxCox transform a vector and calculate Z-scores
# author: Amila Weerasinghe (amila@wustl.edu)
# version: v0.0 - May 2017

args = commandArgs(trailingOnly=TRUE) # 1 - file name

fn = args[1] # file name

library(ggplot2) # for plotting
library(data.table)
library("EnvStats") # for calculating optimized lambda-parameter in BoxCox transformation

# read the input table (One column will be BoxCox transformed)
SMR = read.table(fn, header = FALSE, sep = "\t")

###########################################################
#####################  FUNCTIONS  #########################
###########################################################

# calculate the optimal lambda parameter for BoxCox transformation. PPCC is a measure for how Normal the transformed distribution will be.
calculate_lambda = function(vec) {
  vec = vec+0.00000001
  res = boxcox(vec, optimize = TRUE)
  lambda = res[1]
  PPCC = res[2]
  cat(paste("Lambda = ", lambda, "\t PPCC = ", PPCC, "\n"))
  return(res)
}

# perform BoxCox transformation
boxcox_transform = function(vec, lambda) {
  cat("\n##### BoxCox Transformation-Starting #####\n")
  transformedVec = vec
  if (lambda<0.1 & lambda>-0.1){
    transformedVec = log(vec + 0.0001)
  }
  else {
    transformedVec = (vec^lambda - 1)/lambda
  }
  cat("\n##### BoxCox Transformation-Done #####\n")
  
  return(transformedVec)
}

# calculate Z-Scores
calculate_outlier_zscores = function(vec){ 
  cat("\n##### OUTLIER ANALYSIS-Starting #####\n")
  
  outlierVec = ( vec - mean(vec, na.rm = T) )/sd(vec, na.rm = T)
    
  cat("\n Z-scores have being calculated...\n")
  return(outlierVec)
}

###########################################################
#######################  MAIN  ############################
###########################################################

# first get BoxCox lambda 
res = calculate_lambda(as.vector(SMR$V2))
lambda = as.numeric(res[1]) 

# do the BoxCox transformation
transformedVec = boxcox_transform(as.vector(SMR$V2), lambda)

# calculate Z-scores
zscoresVec = calculate_outlier_zscores(transformedVec)

# add to the table and print
SMR[, "BoxCoxFreq"] = transformedVec
SMR[, "SMRZScore"] = zscoresVec

# make gene name
#SMR[, "GeneName"] = sapply(strsplit(as.vector(SMR$V1),"[.]"), `[`, 1)
SMR[, "GeneName"] = as.vector(SMR$V1)

# determine top 10% cutoff
cutoff = as.numeric(quantile(as.vector(as.numeric(SMR[, "SMRZScore"])), 0.90))
numberOfSMR = length(unique(SMR[SMR$SMRZScore>=cutoff,]$V1))
totalRegions = length(unique(SMR$V1))
numberOfSMRgenes = length(unique(SMR[SMR$SMRZScore>=cutoff,]$GeneName))
cat( paste( "cutoff = ", cutoff , "Number of SMRs = ", numberOfSMR, "out of total regions = ", totalRegions, "and number of SMR genes =", numberOfSMRgenes, "\n" ) )

write.table(SMR, file = "Transformed_example_file.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(SMR[SMR$SMRZScore>=cutoff,], file = "SMR_exmaple_file.txt", quote = FALSE, sep = "\t", row.names = FALSE)


# Plot Zscore distribution
p = ggplot(SMR)+ theme_bw()
p = p + geom_freqpoly(aes(x = as.vector(SMR$SMRZScore)), bins=90, colour="blue" )
p = p + ggtitle("BoxCox Zscores of Mutation Frequency")
p = p + geom_vline(xintercept = cutoff)
p = p + ylab("Number of Regions") + xlab("Z-Score of the Frequency")
ggsave("./freqPlot_Zscore.pdf", width=10, useDingbat=F)

# Plot original distribution
p = ggplot(SMR)+ theme_bw()
p = p + geom_freqpoly(aes(x = as.vector(SMR$V2)), bins=90, colour="red" )
p = p + ggtitle("Mutation Frequency")
p = p + ylab("Number of Regions") + xlab("Mutation Frequency")
ggsave("./freqPlot_regular.pdf", width=10, useDingbat=F)

