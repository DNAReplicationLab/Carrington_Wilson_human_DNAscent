#!/usr/bin/env Rscript

# import detect chunks into R and make frequency distribution plots of BdrU probabilities
# for thymidiine positions

# usage Rscript figure1a.R <working_dir> <save_directory> <comma,separated,list,of,barcodes,to,be,used,together> <names> <set_name>

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
  file_list <- list.files(chunk_directory, pattern = "xx", full.names = TRUE)
  reads <- lapply(file_list, read.table, sep = "\t", skip = 1)
  headers <- sapply(file_list, scan, what = "character", nlines = 1, sep = " ", quiet = TRUE)
  readIDs <- sapply(strsplit(as.character(headers[1,]), split = ">"), "[[", 2)
  names(reads) <- readIDs
  reads_df <- bind_rows(reads, .id = "readID")
  colnames(reads_df) <- c("readID", "position", "probability", "kmer")
  reads_df <- reads_df %>% mutate(barcode = barcode) %>% mutate(sample = name)
  return(reads_df)
}

# plot frequency plot of all probabilities
plotProb_freq <- function(df, barcode, name) {
  plot1 <- ggplot(df, aes(x = probability, y = after_stat(density))) +
    geom_histogram(binwidth = 0.005) +
    theme_bw() +
    xlab("D-NAscent BrdU probability") +
    ylab("Frequency") +
    ggtitle(paste0(set_name, "_", name, " BrdU probability frequency"))
  
  ggsave(paste0(save_directory, set_name, "_", name, "_prob_freq.pdf"), plot = plot1 + coord_cartesian(xlim=c(0,1)), width = 7, height = 4)
  ggsave(paste0(save_directory, set_name, "_", name, "_prob_freq_trim.pdf"), plot = plot1 + coord_cartesian(xlim=c(0.1,1)), width = 7, height = 4)
  
  plot_freq <- plot_freq_empty + geom_histogram(data = df, binwidth = 0.005, alpha = 0.1)
  plot_list <- ggplot_build(plot_freq)
  freq_df <- plot_list[[1]][[1]]
  write_csv(freq_df, paste0(save_directory, set_name, "_", name, "freq_data.csv"))
}

# single value for fraction of Ts >= 0.5
fractionGreaterThan0.5 <- function(df, name) {
  summary_df <- df %>% 
    group_by(barcode) %>%
    summarise(n = n(), "GreaterThan0.5" = sum(probability >= 0.5), "FractionGreaterThan0.5" = GreaterThan0.5 / n) %>%
    mutate(sample = name)
  fractBrdU <- bind_rows(fractBrdU, summary_df)
  return(fractBrdU)
}

# per read fraction Ts >= 0.5
readFractionGreaterThan0.5 <- function(df, barcode, name) {
  summary_df <- df %>% 
    group_by(readID) %>%
    summarise(Ts = n(), "GreaterThan0.5" = sum(probability >= 0.5), "FractionGreaterThan0.5" = GreaterThan0.5 / Ts) %>%
    mutate(barcode = as.character(barcode)) %>% mutate(sample = as.character(name))
  readFractBrdU <- bind_rows(readFractBrdU, summary_df)
  return(readFractBrdU)
}

#########################
#set1

# empty frequency and histogram plots for adding barcode data to

df <- data.frame(readID = character(0), position = integer(0), probability = numeric(0), kmer = character(0), barcode = character(0), sample = character(0))
plot_freq_empty <- ggplot(df, aes(x = probability, y = after_stat(density), colour = sample, fill = sample)) +
  theme_bw() +
  xlab("D-NAscent BrdU probability") +
  ylab("Frequency") +
  ggtitle(paste0(set_name, " BrdU probability frequency"))

# empty fraction greater than 0.5 dfs for adding to

fractBrdU <- tibble()
readFractBrdU <- data.frame(readID = character(0), Ts = integer(0), GreaterThan0.5 = integer(0), FractionGreaterThan0.5 = double(0), barcode = character(0), sample = character(0))

###############################

# for loop for barcodes in set1

for (i in seq_along(set1)) {
  print(paste0("Set1[i] = ", set1[i]), set1_names[i])
  barcode_df <- loadBarcode(paste0(working_directory, set1[i], ".detect_chunks"), set1[i], set1_names[i])
  print(paste0("reads imported"))
  plotProb_freq(barcode_df, set1[i], set1_names[i])
  fractBrdU <- fractionGreaterThan0.5(barcode_df, set1_names[i])
  readFractBrdU <- readFractionGreaterThan0.5(barcode_df, set1[i], set1_names[i])
  print(paste0("Individual plots made"))
}
  
rm(barcode_df)

###############################

# save fractBrdU and readFractBrdU

write_csv(fractBrdU, paste0(save_directory, set_name, "_FractionGreaterThan0.5.csv"))
write_csv(readFractBrdU, paste0(save_directory, set_name, "_readFractionGreaterThan0.5.csv"))

rm(readFractBrdU)

print(paste0("Set1 .csv saved"))

########### 


