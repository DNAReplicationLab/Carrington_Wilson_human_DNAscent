#29.1.21 read BBT.tsv from Jamie's script and make frequency distribution plot 

#Clear R's brain
rm(list = ls())

#Load libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)

setwd("~/Data/Nanopore/2019_07_17_JTC_ONT_HeRPBrdU/")

#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly = TRUE)

#chunk_directory <- args[2]
file_directory <- "~/Data/Nanopore/2019_07_17_JTC_ONT_HeRPBrdU/"

#save_directory <- args[3]
save_dir <- "~/Data/Nanopore/2019_07_17_JTC_ONT_HeRPBrdU/BBT_distribution_plots_2019/"

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

HeLa <- c("barcode01", "barcode02", "barcode03", "barcode04", "barcode05")
RPE <- c("barcode06", "barcode07", "barcode08", "barcode09", "barcode10")

# plot HeLa
plot_H_all <- ggplot(subset(BBT_DF, barcode %in% HeLa), aes(x = `B/(B+T)`, y = after_stat(density))) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_H_all
ggsave("barcodes_BBT_freq_H_all.pdf", plot = plot_H_all, width = 7, height = 4)

plot_hist_H_all <- ggplot(subset(BBT_DF, barcode %in% HeLa), aes(x = `B/(B+T)`)) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_hist_H_all
ggsave("barcodes_BBT_hist_H_all.pdf", plot = plot_hist_H_all, width = 7, height = 4)

plot_H_trim <- ggplot(subset(BBT_DF, barcode %in% HeLa), aes(x = `B/(B+T)`, y = after_stat(density))) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(0.05,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_H_trim
ggsave("barcodes_BBT_freq_H_trim.pdf", plot = plot_H_trim, width = 7, height = 4)

plot_hist_H_trim <- ggplot(subset(BBT_DF, barcode %in% HeLa), aes(x = `B/(B+T)`)) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(0.05,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_hist_H_trim
ggsave("barcodes_BBT_hist_H_trim.pdf", plot = plot_hist_H_trim, width = 7, height = 4)

HeLa_50uM <- c("barcode05")

plot_H50_all <- ggplot(subset(BBT_DF, barcode %in% HeLa_50uM), aes(x = `B/(B+T)`, y = after_stat(density))) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_H50_all
ggsave("barcodes_BBT_freq_H50_all.pdf", plot = plot_H50_all, width = 7, height = 4)

plot_H50_all_hist <- ggplot(subset(BBT_DF, barcode %in% HeLa_50uM), aes(x = `B/(B+T)`)) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_H50_all_hist
ggsave("barcodes_BBT_hist_H50_all.pdf", plot = plot_H50_all_hist, width = 7, height = 4)

#plot all HeLa
plot_HeLa <-ggplot(subset(BBT_DF, barcode %in% HeLa), aes(x = `B/(B+T)`, y = after_stat(density), colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(-0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_HeLa
ggsave("barcodes_BBT_freq_HeLa.pdf", plot = plot_HeLa, width = 7, height = 4)

plot_HeLa_trim <-ggplot(subset(BBT_DF, barcode %in% HeLa), aes(x = `B/(B+T)`, y = after_stat(density), colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_HeLa_trim
ggsave("barcodes_BBT_freq_HeLa_trim.pdf", plot = plot_HeLa_trim, width = 7, height = 4)

plot_HeLa_hist <-ggplot(subset(BBT_DF, barcode %in% HeLa), aes(x = `B/(B+T)`, colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(-0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_HeLa_hist
ggsave("barcodes_BBT_hist_HeLa.pdf", plot = plot_HeLa_hist, width = 7, height = 4)

plot_HeLa_hist_trim <-ggplot(subset(BBT_DF, barcode %in% HeLa), aes(x = `B/(B+T)`, colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_HeLa_hist_trim
ggsave("barcodes_BBT_hist_HeLa_trim.pdf", plot = plot_HeLa_hist_trim, width = 7, height = 4)

# plot RPE

plot_R_all <- ggplot(subset(BBT_DF, barcode %in% RPE), aes(x = `B/(B+T)`, y = after_stat(density))) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_R_all
ggsave("barcodes_BBT_freq_R_all.pdf", plot = plot_R_all, width = 7, height = 4)

plot_hist_R_all <- ggplot(subset(BBT_DF, barcode %in% RPE), aes(x = `B/(B+T)`)) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_hist_R_all
ggsave("barcodes_BBT_hist_R_all.pdf", plot = plot_hist_R_all, width = 7, height = 4)

plot_R_trim <- ggplot(subset(BBT_DF, barcode %in% RPE), aes(x = `B/(B+T)`, y = after_stat(density))) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(0.05,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_R_trim
ggsave("barcodes_BBT_freq_R_trim.pdf", plot = plot_R_trim, width = 7, height = 4)

plot_hist_R_trim <- ggplot(subset(BBT_DF, barcode %in% RPE), aes(x = `B/(B+T)`)) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(0.05,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_hist_R_trim
ggsave("barcodes_BBT_hist_R_trim.pdf", plot = plot_hist_R_trim, width = 7, height = 4)

RPE_50uM <- c("barcode10")

plot_R50_all <- ggplot(subset(BBT_DF, barcode %in% RPE_50uM), aes(x = `B/(B+T)`, y = after_stat(density))) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_R50_all
ggsave("barcodes_BBT_freq_R50_all.pdf", plot = plot_R50_all, width = 7, height = 4)

plot_R50_all_hist <- ggplot(subset(BBT_DF, barcode %in% RPE_50uM), aes(x = `B/(B+T)`)) +
  geom_histogram(binwidth = 0.005) +
  facet_wrap(~ barcode) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.01,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_R50_all_hist
ggsave("barcodes_BBT_hist_R50_all.pdf", plot = plot_R50_all_hist, width = 7, height = 4)

#plot all RPE
plot_RPE <-ggplot(subset(BBT_DF, barcode %in% RPE), aes(x = `B/(B+T)`, y = after_stat(density), colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(-0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_RPE
ggsave("barcodes_BBT_freq_RPE.pdf", plot = plot_RPE, width = 7, height = 4)

plot_RPE_trim <-ggplot(subset(BBT_DF, barcode %in% RPE), aes(x = `B/(B+T)`, y = after_stat(density), colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Frequency")

plot_RPE_trim
ggsave("barcodes_BBT_freq_RPE_trim.pdf", plot = plot_RPE_trim, width = 7, height = 4)

plot_RPE_hist <-ggplot(subset(BBT_DF, barcode %in% RPE), aes(x = `B/(B+T)`, colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(-0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_RPE_hist
ggsave("barcodes_BBT_hist_RPE.pdf", plot = plot_RPE_hist, width = 7, height = 4)

plot_RPE_hist_trim <-ggplot(subset(BBT_DF, barcode %in% RPE), aes(x = `B/(B+T)`, colour = barcode, fill = barcode)) +
  geom_histogram(binwidth = 0.01, alpha = 0.1, position = "identity") +
  theme_bw() +
  scale_x_continuous(limits = c(0.1,1)) +
  xlab("Fraction BrdU") +
  ylab("Count")

plot_RPE_hist_trim
ggsave("barcodes_BBT_hist_RPE_trim.pdf", plot = plot_RPE_hist_trim, width = 7, height = 4)
