#Using islands called by SICER2 as windows 4/3/21

#Clear R's brain
rm(list = ls())

#Load libraries
library(dplyr)
library(ggplot2)
library(readr)

#set working directory
setwd("/Users/pguest/Desktop/Rose_NewChanged/2019_11_04_sicer/SICER_peaks_beds/")

#create and populate list of files for Inputs and IPs
my_files <- list.files(pattern="\\peaks.bed$")
my_INPUTfiles <- grep("INPUT", my_files, value = TRUE)
my_IPfiles <- grep("IP", my_files, value = TRUE)

#vector for column names
bed.cols<-c("chrom","IslandStart","IslandEnd","IslandReadCount")

#import INPUT files to list of dfs
#Set the name of each list element to respective filename.
#convert to single df

myINPUTdata <- list()
for (i in seq_along(my_INPUTfiles)) {
  myINPUTdata[[i]] <- read_delim(file = my_INPUTfiles[i], "\t", escape_double = FALSE, col_names = bed.cols, trim_ws = TRUE)
}

names(myINPUTdata) <- my_INPUTfiles

dfINPUTall <- bind_rows(myINPUTdata, .id = "sample")

#add timepoint name for joining
dfINPUTall$sampleName = sub(".SICERpeaks.bed", "", dfINPUTall$sample)
dfINPUTall$sampleName = sub("INPUT_", "", dfINPUTall$sampleName)

#add column to join by
dfINPUTall <- dfINPUTall %>%
  mutate(ForJoin = paste(sampleName, chrom, IslandStart, sep = '_'))

#import IP files to list of dfs
#Set the name of each list element to respective filename.
#convert to single df

myIPdata <- list()
for (i in seq_along(my_IPfiles)) {
  myIPdata[[i]] <- read_delim(file = my_IPfiles[i], "\t", escape_double = FALSE, col_names = bed.cols, trim_ws = TRUE)
}

names(myIPdata) <- my_IPfiles

dfIPall <- bind_rows(myIPdata, .id = "sample")

#add timepoint name for joining
dfIPall$sampleName = sub(".SICERpeaks.bed", "", dfIPall$sample)
dfIPall$sampleName = sub("IP_", "", dfIPall$sampleName)

#add column to join by
dfIPall <- dfIPall %>%
  mutate(ForJoin = paste(sampleName, chrom, IslandStart, sep = '_'))

#Make df with INPUT and IP
dfall <- dfINPUTall %>%
          select(IslandReadCount, ForJoin) %>%
          inner_join(dfIPall, ., by = "ForJoin", suffix = c(".IP", ".INPUT")) %>%
          mutate(IPoverINPUT = IslandReadCount.IP/IslandReadCount.INPUT) %>%
          mutate(position = paste(chrom, IslandStart, sep = '_'))
dfall$sampleName <- factor(dfall$sampleName, 
                           levels = c("15min0uM",  "15min1uM", "15min5uM", "15min10uM", "30min_0uM", "30min_1uM", "30min_5uM", "30min_10uM", "60min_0uM", "60min_1uM", "60min_5uM", "60min_10uM", "120min_0uM", "120min_1uM", "120min_5uM", "120min_10uM"))
dfall <- dfall %>%
            mutate(Time = as.numeric(gsub("min[_]?[0-9]+uM", "", sampleName))) %>%
            mutate(BrdU = factor(gsub("[0-9]+min[_]?", "", sampleName), levels = c("0uM", "1uM", "5uM", "10uM")))

# Calculate average IP/INPUT for 0uM for each timepoint
noTreatAVG <- dfall %>%
                filter(BrdU == "0uM") %>%
                group_by(Time) %>%
                mutate(IPoverINPUT = na_if(IPoverINPUT, Inf)) %>%
                summarise(AverageIPoverINPUT = mean(IPoverINPUT, na.rm = TRUE))

# normallise IPoverINPUT to no treatment for each timepoint
dfall <- dfall %>%
          mutate(mean_0uMBrdU = case_when(Time == 15 ~ noTreatAVG$AverageIPoverINPUT[1], 
                                          Time == 30 ~ noTreatAVG$AverageIPoverINPUT[2], 
                                          Time == 60 ~ noTreatAVG$AverageIPoverINPUT[3], 
                                          Time == 120 ~ noTreatAVG$AverageIPoverINPUT[4])) %>%
          mutate(normalisedIPoverINPUT = IPoverINPUT / mean_0uMBrdU) %>%
          mutate(Timepoint = factor(paste0(Time, " min"), levels = c("15 min", "30 min", "60 min", "120 min")))

#boxplot of lengths
plot <- ggplot(dfall, aes(x = BrdU, y = normalisedIPoverINPUT)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2, size = 0.1) +
  geom_hline(yintercept = 1, colour = "lightseagreen") +
  facet_wrap(~ Timepoint) +
  ggtitle("Normalised IP/INPUT for peaks from 120min 10uM BrdU sample") +
  xlab("BrdU pulse") + ylab("Normalised IP/INPUT") +
  coord_cartesian(ylim = c(0.1, 50)) +
  scale_y_log10() +
  theme_bw()

plot

ggsave("fig2C_meta.pdf", plot = plot, width = 7, height = 7)

###############################

