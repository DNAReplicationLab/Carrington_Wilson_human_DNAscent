#!/usr/bin/env Rscript

# import detect chunks into R and make frequency distribution plots of BdrU probabilities
# for thymidiine positions

#Load libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)

args = commandArgs(trailingOnly = TRUE)

chunk_directory <- args[1]

#print(paste("chunk_directory = ", chunk_directory, sep = ""))

save_directory <- args[2]

#print(paste("save_directory = ", save_directory, sep = ""))

save_name <- args[3]

paste(save_directory, save_name, ".pdf", sep = "")

#window <- as.numeric(args[4])
#print(paste("window = ", window, sep = ""))

file_list <- list.files(chunk_directory, pattern = "xx", full.names = TRUE)
print("Imported filelist")

# pulls out header line and splits by sep into matrix
headers <- sapply(file_list, scan, what = "character", nlines = 1, sep = " ", quiet = TRUE)
readIDs <- sapply(strsplit(as.character(headers[1,]), split = ">"), "[[", 2)

#chromosomes <- as.character(headers[2,])
#starts <- as.numeric(headers[3,])
#ends <- as.numeric(headers[4,])
#strands <- as.character(headers[5,])

#make list of data frames, one for each read
reads <- lapply(file_list, read.table, sep = "\t", skip = 1)
print("Imported files")

# filter out any reads not long enough, from list and all vectors to iterate over
#long_enough <- sapply(reads, nrow) >= window
#reads <- reads[long_enough]
#file_list <- file_list[long_enough]
#readIDs <- readIDs[long_enough]
#chromosomes <- chromosomes[long_enough]
#starts <- starts[long_enough]
#ends <- ends[long_enough]
#strands <- strands[long_enough]

# make reads list into tibble
names(reads) <- readIDs

print("made dataframe")

reads_df <- bind_rows(reads, .id = "readID")
colnames(reads_df) <- c("readID", "position", "probability", "kmer")

plot1 <- ggplot(reads_df, aes(x = probability, y = after_stat(density))) +
  geom_histogram(binwidth = 0.005) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("BrdU probability") +
  ylab("Frequency")

ggsave(paste(save_directory, save_name, ".pdf", sep = ""), plot = plot1, width = 7, height = 4)

