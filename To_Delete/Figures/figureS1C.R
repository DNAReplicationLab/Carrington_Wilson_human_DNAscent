#!/usr/bin/env Rscript

# import detect chunks into R and make frequency distribution plots of BdrU probabilities
# for thymidiine positions

# usage Rscript figureS1C.R <working_dir> <save_directory> <comma,separated,list,of,barcodes,to,be,used,together> <names> <set_name>

#rm(list = ls())

#Load libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)

args = commandArgs(trailingOnly = TRUE)

working_directory <- args[1]

save_directory <- args[2]

vargs <- strsplit(args, ",")

set1 <- vargs[[3]]
set1_names <- vargs[[4]]

set1 <- factor(set1, levels = set1)
set1_names <- factor(set1_names, levels = set1_names)

set_name <- args[5]

cat("set1 = ", set1, "\n", sep = " ")

barcode_dirs <-dir(path = working_directory, pattern = "^barcode[0-9]+.detect_chunks", all.files = FALSE, full.names = TRUE, recursive = FALSE) 

barcodes <- gsub(".detect_chunks", "", basename(barcode_dirs))

#### functions ############

# loading each barcode data
loadBarcode <- function(chunk_directory, barcode, name) {
  df <- tibble()
  
  file_list <- list.files(chunk_directory, pattern = "xx", full.names = TRUE)
  for (i in seq_along(file_list)) {
    read_df <- read.table(file_list[i], sep = "\t", skip = 1)
    header <- scan(file_list[i], what = "character", nlines = 1, sep = " ", quiet = TRUE)
    readID <- gsub(">", "", as.character(header[1]))
    colnames(read_df) <- c("position", "probability", "kmer")
    read_df <- read_df %>% mutate(readID = readID)
    sum_df <- read_df %>% group_by(readID) %>% summarise(n = n(), "GreaterThan0.5" = sum(probability >= 0.5), "FractionGreaterThan0.5" = GreaterThan0.5 / n)
    df <- bind_rows(df, sum_df)
  }
  df <- df %>% mutate(barcode = barcode, sample = name)
  return(df)
}

#########################
#set1

# for loop for barcodes in set1

set1_df <- tibble()

for (i in seq_along(set1)) {
  print(paste0("Set1[i] = ", set1[i]))
  barcode_df <- loadBarcode(paste0(working_directory, set1[i], ".detect_chunks"), set1[i], set1_names[i])
  set1_df <- bind_rows(set1_df, barcode_df)
  print(paste0("summary df made"))
}
  
write_csv(set1_df, paste0(save_directory, "/", set_name, "_FractionGreaterThan0.5.csv"))

print(paste0("Set1 .csv saved"))

########### 


