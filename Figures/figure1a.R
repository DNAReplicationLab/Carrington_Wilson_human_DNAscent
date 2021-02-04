#!/usr/bin/env Rscript

# import detect chunks into R and make frequency distribution plots of BdrU probabilities
# for thymidiine positions

# usage Rscript figure1a.R <working_dir> <save_directory> <comma,separated,list,of,barcodes,to,be,used,together> <as,prev>

#Load libraries
library("dplyr", lib.loc="/home/nieduszynski/rosemary/R/x86_64-pc-linux-gnu-library/3.5")
library("ggplot2", lib.loc="/home/nieduszynski/rosemary/R/x86_64-pc-linux-gnu-library/3.5")
library("readr", lib.loc="/home/nieduszynski/rosemary/R/x86_64-pc-linux-gnu-library/3.5")
library("tibble", lib.loc="/home/nieduszynski/rosemary/R/x86_64-pc-linux-gnu-library/3.5")

args = commandArgs(trailingOnly = TRUE)
#working_directory <- "/Users/rose/Data/Nanopore/2018_09_10_RW_ONT_RPE_24h_2h"

working_directory <- args[1]

save_directory <- args[2]
#save_directory <- "/Users/rose/Data/Nanopore/2018_09_10_RW_ONT_RPE_24h_2h"

vargs <- strsplit(args, ",")

set1 <- vargs[[3]]
#set2 <- vargs[[4]]

cat("set1 = ", set1, "\n", sep = " ")
#cat("set2 = ", set2, "\n", sep = " ")

barcode_dirs <-dir(path = working_directory, pattern = "^barcode[0-9]+.detect_chunks", all.files = FALSE, full.names = TRUE, recursive = FALSE) 

barcodes <- gsub(".detect_chunks", "", basename(barcode_dirs))

file_list <- list()
i <- 1
for (i in seq_along(barcode_dirs)) {
  file_list[[i]] <- list.files(barcode_dirs[i], pattern = "xx", full.names = TRUE)
}

names(file_list) <- barcodes

print("Imported filelist")

# pulls out header line and splits by sep into matrix
headers <- list()
i <- 1
for (i in seq_along(barcodes)) {
  headers[[i]] <- sapply(file_list[[i]], scan, what = "character", nlines = 1, sep = " ", quiet = TRUE)
}

names(headers) <- barcodes

readIDs <- list()
i <- 1
for (i in seq_along(barcodes)) {
  readIDs[[i]] <- sapply(strsplit(as.character(headers[[i]][1,]), split = ">"), "[[", 2)
}

names(readIDs) <- barcodes

#make list of data frames, one for each read
reads <- list()
i <- 1
for (i in seq_along(barcodes)) {
  reads[[i]] <- lapply(file_list[[i]], read.table, sep = "\t", skip = 1)
}

names(reads) <- barcodes
for (i in seq_along(barcodes)) {
  names(reads[[i]]) <- readIDs[[i]]
}

print("Imported files")
str(reads)

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

barcode_list <- list()
i <- 1
for (i in seq_along(barcodes)) {
  barcode_list[[i]] <- bind_rows(reads[[i]], .id = "readID")
}
names(barcode_list) <- barcodes

str(barcode_list)

barcode_df <- bind_rows(barcode_list, .id = "barcode")
colnames(barcode_df) <- c("barcode", "readID", "position", "probability", "kmer")

print("made dataframe")

########### 

#make single probability plots
for (bar in barcodes) {
  plot1 <- ggplot(filter(barcode_df, barcode == bar), aes(x = probability, y = after_stat(density))) +
    geom_histogram(binwidth = 0.005) +
    theme_bw() +
    scale_x_continuous(limits = c(-0.01,1)) +
    xlab("BrdU probability") +
    ylab("Frequency") +
    ggtitle(bar)
  ggsave(paste(save_directory, "/", bar, "_probability_freq.pdf", sep = ""), plot = plot1, width = 7, height = 4)
}

# make single histogram plots
for (bar in barcodes) {
  plot2 <- ggplot(filter(barcode_df, barcode == bar), aes(x = probability)) +
    geom_histogram(binwidth = 0.005) +
    theme_bw() +
    scale_x_continuous(limits = c(-0.01,1)) +
    xlab("BrdU probability") +
    ylab("Count") +
    ggtitle(bar)
  ggsave(paste(save_directory, "/", bar, "_probability_hist.pdf", sep = ""), plot = plot2, width = 7, height = 4)
}

# Make facet frequency and histogram plots for set 1

plot3 <- ggplot(subset(barcode_df, barcode %in% set1), aes(x = probability, y = after_stat(density))) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

ggsave(paste(save_directory, "/set1_probability_freq.pdf", sep = ""), plot = plot3, width = 7, height = 4)

plot4 <- ggplot(subset(barcode_df, barcode %in% set1), aes(x = probability)) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")
  
ggsave(paste(save_directory, "/set1_probability_hist.pdf", sep = ""), plot = plot4, width = 7, height = 4)

# Make facet frequency and histogram plots for set 2

plot5 <- ggplot(subset(barcode_df, barcode %in% set2), aes(x = probability, y = after_stat(density))) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

ggsave(paste(save_directory, "/set2_probability_freq.pdf", sep = ""), plot = plot5, width = 7, height = 4)

plot6 <- ggplot(subset(barcode_df, barcode %in% set2), aes(x = probability)) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

ggsave(paste(save_directory, "/set2_probability_hist.pdf", sep = ""), plot = plot6, width = 7, height = 4)

# plot full and trim, freq and histogram overlay plots for set1

plot7 <-ggplot(subset(barcode_df, barcode %in% set1), aes(x = probability, y = after_stat(density), colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(-0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

ggsave(paste(save_directory, "/set1_overlay_probability_freq.pdf", sep = ""), plot = plot7, width = 7, height = 4)

plot8 <-ggplot(subset(barcode_df, barcode %in% set1), aes(x = probability, y = after_stat(density), colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

ggsave(paste(save_directory, "/set1_overlay_probability_freq_trim.pdf", sep = ""), plot = plot8, width = 7, height = 4)

plot9 <-ggplot(subset(barcode_df, barcode %in% set1), aes(x = probability, colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(-0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

ggsave(paste(save_directory, "/set1_overlay_probability_hist.pdf", sep = ""), plot = plot9, width = 7, height = 4)

plot10 <-ggplot(subset(barcode_df, barcode %in% set1), aes(x = probability, colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

ggsave(paste(save_directory, "/set1_overlay_probability_hist_trim.pdf", sep = ""), plot = plot10, width = 7, height = 4)

# plot full and trim, freq and histogram overlay plots for set2

plot11 <-ggplot(subset(barcode_df, barcode %in% set2), aes(x = probability, y = after_stat(density), colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(-0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

ggsave(paste(save_directory, "/set2_overlay_probability_freq.pdf", sep = ""), plot = plot11, width = 7, height = 4)

plot12 <-ggplot(subset(barcode_df, barcode %in% set2), aes(x = probability, y = after_stat(density), colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

ggsave(paste(save_directory, "/set2_overlay_probability_freq_trim.pdf", sep = ""), plot = plot12, width = 7, height = 4)

plot13 <-ggplot(subset(barcode_df, barcode %in% set2), aes(x = probability, colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(-0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

ggsave(paste(save_directory, "/set2_overlay_probability_hist.pdf", sep = ""), plot = plot13, width = 7, height = 4)

plot14 <-ggplot(subset(barcode_df, barcode %in% set2), aes(x = probability, colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

ggsave(paste(save_directory, "/set2_overlay_probability_hist_trim.pdf", sep = ""), plot = plot14, width = 7, height = 4)

####################
# single value for fraction of Ts > 0.5
SummaryTable <- barcode_df %>% 
                  group_by(barcode) %>%
                  summarise(n = n(), "GreaterThan0.5" = sum(probability > 0.5), "FractionGreaterThan0.5" = GreaterThan0.5 / n)
SummaryTable

write_csv(SummaryTable, paste0(save_directory, "/FractionGreaterThan0.5.csv"))
