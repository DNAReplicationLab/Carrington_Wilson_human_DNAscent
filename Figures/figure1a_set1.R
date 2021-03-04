#!/usr/bin/env Rscript

# import detect chunks into R and make frequency distribution plots of BdrU probabilities
# for thymidiine positions

# usage Rscript figure1a.R <working_dir> <save_directory> <comma,separated,list,of,barcodes,to,be,used,together> <as,prev>

#rm(list = ls())

#Load libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)

args = commandArgs(trailingOnly = TRUE)
#working_directory <- "/Users/rose/Data/Nanopore/2018_09_10_RW_ONT_RPE_24h_2h/"

working_directory <- args[1]

save_directory <- args[2]
#save_directory <- "/Users/rose/Data/Nanopore/2018_09_10_RW_ONT_RPE_24h_2h/test/"

vargs <- strsplit(args, ",")

set1 <- vargs[[3]]
# set2 <- vargs[[4]]

#set1 <- c("barcode06", "barcode07", "barcode05", "barcode09")
#set1 <- c("barcode06", "barcode05")
set1 <- factor(set1, levels = set1)
#set2 <- c("barcode06", "barcode10", "barcode11", "barcode12")
#set2 <- factor(set2, levels = set2)

cat("set1 = ", as.character(set1), "\n", sep = " ")
#cat("set2 = ", as.character(set2), "\n", sep = " ")

barcode_dirs <-dir(path = working_directory, pattern = "^barcode[0-9]+.detect_chunks", all.files = FALSE, full.names = TRUE, recursive = FALSE) 

barcodes <- gsub(".detect_chunks", "", basename(barcode_dirs))

#### functions ############

# loading each barcode data
loadBarcode <- function(chunk_directory, barcode) {
  file_list <- list.files(chunk_directory, pattern = "xx", full.names = TRUE)
  reads <- lapply(file_list, read.table, sep = "\t", skip = 1)
  headers <- sapply(file_list, scan, what = "character", nlines = 1, sep = " ", quiet = TRUE)
  readIDs <- sapply(strsplit(as.character(headers[1,]), split = ">"), "[[", 2)
  names(reads) <- readIDs
  reads_df <- bind_rows(reads, .id = "readID")
  colnames(reads_df) <- c("readID", "position", "probability", "kmer")
  reads_df <- reads_df %>% mutate(barcode = barcode)
  return(reads_df)
}

# plot frequency plot of all probabilities
plotProb_freq <- function(df, barcode) {
  plot1 <- ggplot(df, aes(x = probability, y = after_stat(density))) +
    geom_histogram(binwidth = 0.005) +
    theme_bw() +
    xlab("BrdU probability") +
    ylab("Frequency") +
    ggtitle(paste0(barcode, " BrdU probability frequency"))
  
#  ggsave(paste0(save_directory, barcode, "_prob_freq.pdf"), plot = plot1 + scale_x_continuous(limits = c(-0.01,1)), width = 7, height = 4)
#  ggsave(paste0(save_directory, barcode, "_prob_freq_trim.pdf"), plot = plot1 + scale_x_continuous(limits = c(0.1, 1)), width = 7, height = 4)
  
  plot_freq <- plot_freq + geom_histogram(data = df, binwidth = 0.005, alpha = 0.1)
  return(plot_freq)
}

# plot histogram plot of all probabilities
plotProb_hist <- function(df, barcode) {
  plot1 <- ggplot(df, aes(x = probability)) +
    geom_histogram(binwidth = 0.005) +
    theme_bw() +
    xlab("BrdU probability") +
    ylab("Count") +
    ggtitle(paste0(barcode, " BrdU probability histogram"))
  
#  ggsave(paste0(save_directory, barcode, "_prob_hist.pdf"), plot = plot1 + scale_x_continuous(limits = c(-0.01,1)), width = 7, height = 4)
#  ggsave(paste0(save_directory, barcode, "_prob_hist_trim.pdf"), plot = plot1 + scale_x_continuous(limits = c(0.1, 1)), width = 7, height = 4)
  
  plot_hist <- plot_hist + geom_histogram(data = df, binwidth = 0.005, alpha = 0.1)
  return(plot_hist)
}

# single value for fraction of Ts > 0.5
fractionGreaterThan0.5 <- function(df) {
  summary_df <- df %>% 
    group_by(barcode) %>%
    summarise(n = n(), "GreaterThan0.5" = sum(probability > 0.5), "FractionGreaterThan0.5" = GreaterThan0.5 / n)
  fractBrdU <- bind_rows(fractBrdU, summary_df)
  return(fractBrdU)
}

# per read fraction Ts > 0.5
readFractionGreaterThan0.5 <- function(df, barcode) {
  summary_df <- df %>% 
    group_by(readID) %>%
    summarise(Ts = n(), "GreaterThan0.5" = sum(probability > 0.5), "FractionGreaterThan0.5" = GreaterThan0.5 / Ts) %>%
    mutate(barcode = as.character(barcode))
  readFractBrdU <- bind_rows(readFractBrdU, summary_df)
  return(readFractBrdU)
}

#########################
#set1

# empty frequency and histogram plots for adding barcode data to

df <- data.frame(readID = character(0), position = integer(0), probability = numeric(0), kmer = character(0), barcode = character(0))
plot_freq <- ggplot(df, aes(x = probability, y = after_stat(density), colour = barcode, fill = barcode)) +
  theme_bw() +
  xlab("BrdU probability") +
  ylab("Frequency") +
  ggtitle("BrdU probability frequency")

plot_hist <- ggplot(df, aes(x = probability, colour = barcode, fill = barcode)) +
  theme_bw() +
  xlab("BrdU probability") +
  ylab("Count") +
  ggtitle("BrdU probability histogram")

# empty fraction greater than 0.5 dfs for adding to

fractBrdU <- tibble()
readFractBrdU <- data.frame(readID = character(0), Ts = integer(0), GreaterThan0.5 = integer(0), FractionGreaterThan0.5 = double(0), barcode = character(0))

###############################

# for loop for barcodes in set1

for (i in seq_along(set1)) {
  print(paste0("Set1[i] = ", set1[i]))
  barcode_df <- loadBarcode(paste0(working_directory, set1[i], ".detect_chunks"), set1[i])
  print(paste0("reads imported"))
  plot_freq <- plotProb_freq(barcode_df, set1[i])
  plot_hist <- plotProb_hist(barcode_df, set1[i])
  fractBrdU <- fractionGreaterThan0.5(barcode_df)
  readFractBrdU <- readFractionGreaterThan0.5(barcode_df, set1[i])
  print(paste0("Individual plots made"))
}
  

###############################
# save fractBrdU and readFractBrdU

write_csv(fractBrdU, paste0(save_directory, "set1_FractionGreaterThan0.5.csv"))
write_csv(readFractBrdU, paste0(save_directory, "set1_readFractionGreaterThan0.5.csv"))

print(paste0("Set1 .csv saved"))

# make set1 plots

ggsave(paste0(save_directory, "set1_prob_freq.pdf"), plot = plot_freq + scale_x_continuous(limits = c(-0.01,1)), width = 7, height = 4)
ggsave(paste0(save_directory, "set1_prob_freq_trim.pdf"), plot = plot_freq + scale_x_continuous(limits = c(0.1,1)), width = 7, height = 4)
ggsave(paste0(save_directory, "set1_prob_hist.pdf"), plot = plot_hist + scale_x_continuous(limits = c(-0.01,1)), width = 7, height = 4)
ggsave(paste0(save_directory, "set1_prob_hist_trim.pdf"), plot = plot_hist + scale_x_continuous(limits = c(0.1,1)), width = 7, height = 4)

ggsave(paste0(save_directory, "set1_prob_freq_facet.pdf"), plot = plot_freq + scale_x_continuous(limits = c(-0.01,1)) + facet_wrap(~ barcode), width = 7, height = 4)
ggsave(paste0(save_directory, "set1_prob_freq_facet_trim.pdf"), plot = plot_freq + scale_x_continuous(limits = c(0.1,1)) + facet_wrap(~ barcode), width = 7, height = 4)
ggsave(paste0(save_directory, "set1_prob_hist_facet.pdf"), plot = plot_hist + scale_x_continuous(limits = c(-0.01,1)) + facet_wrap(~ barcode), width = 7, height = 4)
ggsave(paste0(save_directory, "set1_prob_hist_facet_trim.pdf"), plot = plot_hist + scale_x_continuous(limits = c(0.1,1)) + facet_wrap(~ barcode), width = 7, height = 4)

print(paste0("Set1 plots made"))

# save fractBrdU and readFractBrdU

write_csv(fractBrdU, paste0(save_directory, "set1_FractionGreaterThan0.5.csv"))
write_csv(readFractBrdU, paste0(save_directory, "set1_readFractionGreaterThan0.5.csv"))

print(paste0("Set1 .csv saved"))

########### 




