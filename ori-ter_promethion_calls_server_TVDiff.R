#!/usr/bin/env Rscript

#--------------------------------------------------------------
# Copyright 2024-2025 Jamie T Carrington (University of Oxford)
# Written by Jamie T Carrington (University of Oxford)
# This software is licensed under GPL-3.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
#--------------------------------------------------------------

# This script is designed to identify DNA replication forks, and their orientation, from DNAscent 'detect' data.
# It also infers DNA replication initiation and termination zones from those fork calls.
# This script produces four BED3 files, "left_forks.bed", "right_forks.bed", "origins.bed", and "termini.bed".

# It has only been tested on DNAscent 2, which means it is restricted to BrdU detection only, and r9.4.1 ONT flow cells.
# The principles behind this script should be applicable to updated flow cells and software, but it is expected
# that the following parameters would have to optimised: window size, step size, probability threshold, and
# gradient threshold.

# It is assumed that human cells have been treated with a 'gradient' of BrdU: starting from a low concentration
# and ending on a high concentration (i.e. 0 uM to 12 uM over a period of an hour) and then chased with the highest BrdU concentration
# for an hour (i.e. 12 uM).
# In principle this script can work for non-human cells that can scavenge and then phosphorylate nucleotides from media
# (e.g. mouse, not S. cerevisiae), but this has not been tested and parameters may need optimising.
# It may also be possible to modify the 'gradient' of BrdU conditions (e.g. start gradient with a low background
# treatment of BrdU, or finishing the gradient on a much higher concentration of BrdU), but only the example above has been tested.
# We find that there are diminishing returns of BrdU detection somewhere between 10 uM BrdU and 50 uM BrdU treatments.

# DNAscent 2 produces a 'detect' file, a large table where each row describes the mapped coordinate of each thymidine contained in
# each read. Rows belonging to a single read are headed by an additional row starting with ">" followed by the readID, the
# coordinates of where the read mapped, and the strand of the genome to which it mapped. For the sake of memory management, a
# bash command was used to break the large detect table into individual tables, one for each read, each with its own header:
# csplit -s file.detect '/>/' '{*}'
# each output file pertains to a single read, and their names start with 'xx' following by a number. File xx00 contains the overall
# header information produced by the DNAscent detect, and is not imported by this script.

# The following arguments are required:
# Input --args in space separated list in this order:
# 1) path to directory that contains csplit detect files
# 2) path to directory where output files will be saved
# 3) B/(B+T) window length in basepairs
# 4) Sliding window step size in basepairs
# 5) Probability threshold at which a BrdU will be called
# 6) Gradient threshold for calling DNA replication forks
#example usage: R CMD ./ori-ter_promethion_calls_server_TVDiff.R --args /path/to/files/ /path/to/save/ 290 290 0.5 1

args = commandArgs(trailingOnly = TRUE)
#args = c("", "~/Desktop/gradients/fork_call_troubleshooting/", "~/Downloads/", 290, 290, 0.5, 1)

#load libraries
install.packages("data.table")
library(data.table)
install.packages("devtools")
devtools::install_github("natbprice/tvdiff")
library(tvdiff)

#set parameters
file_directory = args[2]
save_directory = args[3]
window = as.numeric(args[4])
step = as.numeric(args[5])
threshold = as.numeric(args[6])
gradient_threshold = as.numeric(args[7])

cat("\n---Command summary---\n")
cat(paste0("Input files location: ", args[2], "\n"))
cat(paste0("Output files location: ", args[3], "\n"))
cat(paste0("Window size (bp): ", args[4], "\n"))
cat(paste0("Step size (bp): ", args[5], "\n"))
cat(paste0("BrdU call probability threshold: ", args[6], "\n"))
cat(paste0("Gradient threshold for fork calling: ", args[7], "\n"))
cat("---------------------\n")

cat("\nImporting files\n")

#import files and metadata
files <- list.files(path = file_directory,
                    pattern = "xx",
                    full.names = TRUE)[-1]

headers <- matrix(nrow = 5, ncol = length(files))
invisible(sapply(X = 1:length(files),
                 FUN = function(file) {
                   headers[,file] <<- scan(file = files[file],
                                           what = "character",
                                           sep = " ",
                                           nlines = 1,
                                           quiet = TRUE)
                   cat(paste0("Processing metadata: ",
                              format(x = round(x = file / length(files) * 100,
                                               digits = 2),
                                     nsmall = 2),
                              "%\t\r"))
                 }))

cat("\n")

read_IDs <- gsub(x = as.character(headers[1,]),
                 pattern = ">",
                 replacement = "")
chromosomes <- as.character(headers[2,])
pre_strands <- as.character(headers[5,])
strands <- ifelse(test = pre_strands == "fwd",
                  yes = "+",
                  no = ifelse(test = pre_strands == "rev",
                              yes = "-",
                              no = "error"))

reads <- list()
invisible(lapply(X = 1:length(files),
                 FUN = function(file) {
                   reads[[file]] <<- fread(file = files[file],
                                           data.table = FALSE,
                                           sep = "\t",
                                           skip = 1,
                                           colClasses = c("numeric", "numeric", "NULL"))
                   cat(paste0("Importing detect data: ",
                              format(x = round(x = file / length(files) * 100,
                                               digits = 2),
                                     nsmall = 2),
                              "%\t\r"))
                 }))

cat("\nCalculating B / (B + T)\n")

#B/(B+T) function
BBT_fun <- function(read) {
  Ts <- reads[[read]][,1]
  Bs <- reads[[read]][,2] >= threshold
  centre_adjust <- floor((length(Ts) - (floor((length(Ts) - window) / step) * step + window)) / 2)
  
  slide_fun <- function(iteration) {
    sliding_window <- step * (iteration - 1) + centre_adjust
    c(Ts[window / 2 + sliding_window],
      sum(Bs[(1 + sliding_window):(window + sliding_window)]) / window,
      Ts[1 + sliding_window],
      Ts[window + sliding_window])
  }
  sapply(X = 1:((floor((nrow(reads[[read]]) - window) / step)) + 1),
         FUN = slide_fun)
}

read_BBTs <- list()
invisible(lapply(X = 1:length(reads),
                 FUN = function(read) {
                   read_BBTs[[read]] <<- as.data.frame(t(BBT_fun(read)))
                   cat(paste0("Calculating B/(B+T): ",
                              format(x = round(x = read / length(reads) * 100,
                                               digits = 2),
                                     nsmall = 2),
                              "%\t\r"))
                 }))
#reads that contain zero BrdU calls (rare) cause the TVDiff algorithm to produce an error, and so are removed here
filter_flat_reads <- sapply(X = 1:length(read_BBTs),
                            FUN = function(read) {
                              !all(read_BBTs[[read]][,2] == 0)
                              })

read_BBTs <- read_BBTs[filter_flat_reads]
read_IDs <- read_IDs[filter_flat_reads]
chromosomes <- chromosomes[filter_flat_reads]
strands <- strands[filter_flat_reads]

cat("Finding replication forks\n")

#The TVDiff algorithm is used to find gradients in the noisy B / (B+T) windowed data. Gradients greater than 1 are
#interpreted as consistent with a replication forks moving codirectional with the 5'->3' orientation of the top/forward strand of hg38
#(i.e. rightwards). Gradients less than -1 are interpreted as consistent with a replication forks moving codirectional with
# the 5'->3' orientation of the bottom/reverse strand of hg38 (i.e. leftwards).
gradient_calls <- list()
invisible(lapply(X = 1:length(read_BBTs),
                 FUN = function(read) {
                   est_derivative <- as.numeric(TVRegDiffR(data = read_BBTs[[read]][,2], iter = 10, alph = 0.01))
                   gradient_calls[[read]] <<- ifelse(test = est_derivative > gradient_threshold,
                                                     yes = 1,
                                                     no = ifelse(test = est_derivative < -gradient_threshold,
                                                                 yes = -1,
                                                                 no = 0))[-1]
                   cat(paste0("Calculating gradients: ",
                              format(x = round(x = read / length(read_BBTs) * 100,
                                               digits = 2),
                                     nsmall = 2),
                              "%\t\r"))
                 }))

fork_fun <- function(read, direction, coordinate = FALSE, index.mode = FALSE) {
  if(direction == "left") {direction = -1}
  if(direction == "right") {direction = 1}
  
  segment_end_index <- cumsum(rle(gradient_calls[[read]])$lengths)
  segment_start_index <- c(0, segment_end_index[-length(segment_end_index)]) + 1
  segment_gradient_call <- rle(gradient_calls[[read]])$values
  
  if(coordinate == "start") {coordinate = 3; side = segment_start_index}
  if(coordinate == "end") {coordinate = 4; side = segment_end_index}
  
  if(index.mode == TRUE) {
    sum(segment_gradient_call == direction)
  } else {
    if(any(segment_gradient_call == direction)) {
      lapply(X = which(segment_gradient_call == direction),
             FUN = function(x) {
               read_BBTs[[read]][side[x],coordinate]
             })
    } else {
      return(NA)
    }
  }
  
}

left_fork_start_coords <- unlist(sapply(X = 1:length(read_BBTs),
                                        FUN = fork_fun,
                                        direction = "left",
                                        coordinate = "start"))
left_fork_end_coords <- unlist(sapply(X = 1:length(read_BBTs),
                                      FUN = fork_fun,
                                      direction = "left",
                                      coordinate = "end"))
right_fork_start_coords <- unlist(sapply(X = 1:length(read_BBTs),
                                        FUN = fork_fun,
                                        direction = "right",
                                        coordinate = "start"))
right_fork_end_coords <- unlist(sapply(X = 1:length(read_BBTs),
                                      FUN = fork_fun,
                                      direction = "right",
                                      coordinate = "end"))
left_fork_index <- sapply(X = 1:length(read_BBTs),
                          FUN = fork_fun,
                          direction = "left",
                          index.mode = TRUE)
right_fork_index <- sapply(X = 1:length(read_BBTs),
                           FUN = fork_fun,
                           direction = "right",
                           index.mode = TRUE)
                          
cat("Finding replication origins and termini\n")

ori_ter_fun <- function(read, feature, coordinate = FALSE, index.mode = FALSE) {
  if(feature == "origin") {feature = 1}
  if(feature == "terminus") {feature = -1}
  
  segment_end_index <- cumsum(rle(gradient_calls[[read]])$lengths)
  segment_start_index <- c(0, segment_end_index[-length(segment_end_index)]) + 1
  segment_gradient_call <- rle(gradient_calls[[read]])$values
  
  vergence <- sapply(X = 1:(length(segment_gradient_call) - 1),
                     FUN = function(i) {
                       sum(segment_gradient_call[i] == -feature,
                           segment_gradient_call[i+1] == feature,
                           segment_gradient_call[i+2] == feature,
                           na.rm = TRUE)
                     })
  
  if(coordinate == "start") {coordinate = 3; side = segment_start_index}
  if(coordinate == "end") {coordinate = 4; side = segment_end_index}
  
  if(index.mode == TRUE) {
    sum(vergence == 2)
  } else {
    if(any(vergence == 2, na.rm = TRUE)) {
      lapply(X = which(vergence == 2) + 1,
             FUN = function(x) {
               if(segment_gradient_call[x] == 0) {
                 return(read_BBTs[[read]][side[x],coordinate])
               }
               if(segment_gradient_call[x] != 0 & coordinate == 3) {
                 return(read_BBTs[[read]][side[x]-1, coordinate])
               }
               if(segment_gradient_call[x] != 0 & coordinate == 4) {
                 return(read_BBTs[[read]][side[x-1]+1, coordinate])
               }
             })
    } else {
      return(NA)
    }
  }
  
}

origin_start_coords <- unlist(sapply(X = 1:length(read_BBTs),
                                     FUN = ori_ter_fun,
                                     feature = "origin",
                                     coordinate = "start"))
origin_end_coords <- unlist(sapply(X = 1:length(read_BBTs),
                                     FUN = ori_ter_fun,
                                     feature = "origin",
                                     coordinate = "end"))
terminus_start_coords <- unlist(sapply(X = 1:length(read_BBTs),
                                     FUN = ori_ter_fun,
                                     feature = "terminus",
                                     coordinate = "start"))
terminus_end_coords <- unlist(sapply(X = 1:length(read_BBTs),
                                   FUN = ori_ter_fun,
                                   feature = "terminus",
                                   coordinate = "end"))
origin_index <- sapply(X = 1:length(read_BBTs),
                       FUN = ori_ter_fun,
                       feature = "origin",
                       index.mode = TRUE)
terminus_index <- sapply(X = 1:length(read_BBTs),
                       FUN = ori_ter_fun,
                       feature = "terminus",
                       index.mode = TRUE)

cat("Writing files\n")
write.table(x = cbind(rep(x = chromosomes, times = left_fork_index),
                      na.omit(left_fork_start_coords),
                      na.omit(left_fork_end_coords),
                      rep(x = read_IDs, times = left_fork_index),
                      ".",
                      rep(x = strands, times = left_fork_index)),
            file = paste0(save_directory, "/left_forks.bed"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(x = cbind(rep(x = chromosomes, times = right_fork_index),
                      na.omit(right_fork_start_coords),
                      na.omit(right_fork_end_coords),
                      rep(x = read_IDs, times = right_fork_index),
                      ".",
                      rep(x = strands, times = right_fork_index)),
            file = paste0(save_directory, "/right_forks.bed"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(x = cbind(rep(x = chromosomes, times = origin_index),
                      na.omit(origin_start_coords),
                      na.omit(origin_end_coords),
                      rep(x = read_IDs, times = origin_index),
                      ".",
                      rep(x = strands, times = origin_index)),
            file = paste0(save_directory, "/origins.bed"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(x = cbind(rep(x = chromosomes, times = terminus_index),
                      na.omit(terminus_start_coords),
                      na.omit(terminus_end_coords),
                      rep(x = read_IDs, times = terminus_index),
                      ".",
                      rep(x = strands, times = terminus_index)),
            file = paste0(save_directory, "/termini.bed"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
cat("Done")
