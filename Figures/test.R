#!/usr/bin/env Rscript

# import detect chunks into R and make frequency distribution plots of BdrU probabilities
# for thymidiine positions

# usage Rscript figure1a.R <working_dir> <save_directory> <comma,separated,list,of,barcodes,to,be,used,together> <as,prev>

args = commandArgs(trailingOnly = TRUE)
#working_directory <- "/Users/rose/Data/Nanopore/2018_09_10_RW_ONT_RPE_24h_2h"

working_directory <- args[1]

working_directory

save_directory <- args[2]
#save_directory <- "/Users/rose/Data/Nanopore/2018_09_10_RW_ONT_RPE_24h_2h"

save_directory

vargs <- strsplit(args, ",")

set1 <- vargs[[3]]
#set2 <- vargs[[4]]

cat("set1 = ", set1, "\n", sep = " ")
