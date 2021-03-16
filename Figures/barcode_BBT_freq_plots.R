#29.1.21 read BBT.tsv from Jamie's script and make frequency distribution plot 

#Clear R's brain
rm(list = ls())

#Load libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)

setwd("~/Data/Nanopore/2018_09_10_RW_ONT_RPE_24h_2h/")

#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly = TRUE)

#chunk_directory <- args[2]
file_directory <- "~/Data/Nanopore/2018_09_10_RW_ONT_RPE_24h_2h"

#save_directory <- args[3]
save_dir <- "~/Data/Nanopore/2018_09_10_RW_ONT_RPE_24h_2h"

file_list <- list.files(file_directory, pattern = ".tsv$", full.names = TRUE)

# read in *_BBT.tsv files into list
BBT_list <- list()
i <- 1
for (i in seq_along(file_list)) {
  BBT_list[[i]] <- read_tsv(file = file_list[i], col_names = TRUE, trim_ws = TRUE)
}

barcodes <- sapply(gsub("_BBT.tsv", "", file_list), basename)
names(BBT_list) <- barcodes

# make df
BBT_DF <- bind_rows(BBT_list, .id = "barcode")

BBT_DF$barcode <- factor(BBT_DF$barcode, levels = c("barcode06", "barcode07", "barcode05", "barcode09", "barcode10", "barcode11", "barcode12"))

plot_all <- ggplot(BBT_DF, aes(x = `B/(B+T)`, y = after_stat(density))) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_all
ggsave("barcodes_BBT_freq_all.pdf", plot = plot_all, width = 7, height = 4)

plot_hist_all <- ggplot(BBT_DF, aes(x = `B/(B+T)`)) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_hist_all
ggsave("barcodes_BBT_hist_all.pdf", plot = plot_hist_all, width = 7, height = 4)

plot_trim <- ggplot(BBT_DF, aes(x = `B/(B+T)`, y = after_stat(density))) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_trim
ggsave("barcodes_BBT_freq_trim.pdf", plot = plot_trim, width = 7, height = 4)

plot_hist_trim <- ggplot(BBT_DF, aes(x = `B/(B+T)`)) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_hist_trim
ggsave("barcodes_BBT_hist_trim.pdf", plot = plot_hist_trim, width = 7, height = 4)

# just plot 2hr
high2hr <- c("barcode09")

plot_50uM2hr <- ggplot(subset(BBT_DF, barcode %in% high2hr), aes(x = `B/(B+T)`, y = after_stat(density))) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_50uM2hr
ggsave("barcodes_BBT_50uM2hr_all.pdf", plot = plot_50uM2hr, width = 7, height = 4)

plot_50uM2hr_hist <- ggplot(subset(BBT_DF, barcode %in% high2hr), aes(x = `B/(B+T)`)) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_50uM2hr_hist
ggsave("barcodes_BBT_50uM2hr_all_hist.pdf", plot = plot_50uM2hr_hist, width = 7, height = 4)

plot_50uM2hrtrim_hist <- ggplot(subset(BBT_DF, barcode %in% high2hr), aes(x = `B/(B+T)`)) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_50uM2hrtrim_hist
ggsave("barcodes_BBT_50uM2hr_trim_hist.pdf", plot = plot_50uM2hrtrim_hist, width = 7, height = 4)

samples2hr <- c("barcode06", "barcode07", "barcode05", "barcode09")

plot_2hrtrim <-ggplot(subset(BBT_DF, barcode %in% samples2hr), aes(x = `B/(B+T)`, y = after_stat(density), colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_2hrtrim
ggsave("barcodes_BBT_2hr_trim.pdf", plot = plot_2hrtrim, width = 7, height = 4)

plot_2hrtrim_hist <-ggplot(subset(BBT_DF, barcode %in% samples2hr), aes(x = `B/(B+T)`, colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_2hrtrim_hist
ggsave("barcodes_BBT_2hr_trim_hist.pdf", plot = plot_2hrtrim_hist, width = 7, height = 4)

#just plot 24hr
high24hr <- c("barcode12")

plot_50uM24hr <- ggplot(subset(BBT_DF, barcode %in% high24hr), aes(x = `B/(B+T)`, y = after_stat(density))) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_50uM24hr
ggsave("barcodes_BBT_50uM24hr_all.pdf", plot = plot_50uM24hr, width = 7, height = 4)

plot_50uM24hr_hist <- ggplot(subset(BBT_DF, barcode %in% high24hr), aes(x = `B/(B+T)`)) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_50uM24hr_hist
ggsave("barcodes_BBT_50uM24hr_all_hist.pdf", plot = plot_50uM24hr_hist, width = 7, height = 4)

plot_50uM24hrtrim <- ggplot(subset(BBT_DF, barcode %in% high24hr), aes(x = `B/(B+T)`, y = after_stat(density))) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_50uM24hrtrim
ggsave("barcodes_BBT_50uM24hr_trim.pdf", plot = plot_50uM24hrtrim, width = 7, height = 4)

plot_50uM24hrtrim_hist <- ggplot(subset(BBT_DF, barcode %in% high24hr), aes(x = `B/(B+T)`)) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_50uM24hrtrim_hist
ggsave("barcodes_BBT_50uM24hr_trim_hist.pdf", plot = plot_50uM24hrtrim_hist, width = 7, height = 4)

samples24hr <- c("barcode06", "barcode10", "barcode11", "barcode12")

plot_24hrall <-ggplot(subset(BBT_DF, barcode %in% samples24hr), aes(x = `B/(B+T)`, y = after_stat(density), colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_24hrall
ggsave("barcodes_BBT_24hr_all.pdf", plot = plot_24hrall, width = 7, height = 4)

plot_24hrall_hist <-ggplot(subset(BBT_DF, barcode %in% samples24hr), aes(x = `B/(B+T)`, colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_24hrall_hist
ggsave("barcodes_BBT_24hr_all_hist.pdf", plot = plot_24hrall_hist, width = 7, height = 4)

plot_24hrtrim <-ggplot(subset(BBT_DF, barcode %in% samples24hr), aes(x = `B/(B+T)`, y = after_stat(density), colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_24hrtrim
ggsave("barcodes_BBT_24hr_trim.pdf", plot = plot_24hrtrim, width = 7, height = 4)

plot_24hrtrim_hist <-ggplot(subset(BBT_DF, barcode %in% samples24hr), aes(x = `B/(B+T)`, colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_24hrtrim_hist
ggsave("barcodes_BBT_24hr_trim_hist.pdf", plot = plot_24hrtrim_hist, width = 7, height = 4)
