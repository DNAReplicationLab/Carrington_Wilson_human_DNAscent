#!/usr/bin/env Rscript
#original code by Jamie Carrington (University of Oxford)
#Modifications by Conrad A. Nieduszynski (Earlham Instutute)
#Usage
args = commandArgs(trailingOnly = TRUE)
if (length(args)<5) {
  stop("At least five arguments must be supplied (input directory, output directory, window size, step size, threshold ).n", call.=FALSE)
} else if (length(args)==5) {
  # default export mode
  args[6] = "FALSE"
}
print(args)
chunk_directory <- args[1]
save_directory <- args[2]
window <- as.numeric(args[3])
step <- as.numeric(args[4])
threshold <- as.numeric(args[5])
export_mode <- args[6]
ROI_directory <- args[7]
file_list <- list.files(chunk_directory, pattern = "bedgraph", full.names = TRUE, recursive = TRUE)
headers <- sapply(file_list, scan, what = "character", nlines = 1, sep = " ", quiet = TRUE)
headers2 <- sapply(file_list, scan, what = "character", nlines = 1, skip=1, sep = " ", quiet = TRUE)
readInfo <- sapply(strsplit(as.character(headers[3,]), split = "="), "[[", 2)
readIDs <- as.character(substr(readInfo,1,nchar(readInfo)-3))
chromosomes <- as.character(headers2[1,])
strands <- as.character(substr(readInfo,nchar(readInfo)-2,nchar(readInfo)))
reads <- lapply(file_list, read.table, sep = " ", skip = 1)
long_enough <- sapply(reads, nrow) >= window
reads <- reads[long_enough]
file_list <- file_list[long_enough]
readIDs <- readIDs[long_enough]
chromosomes <- chromosomes[long_enough]
strands <- strands[long_enough]
BBT_fun <- function(read) {
	internal_fun <- function(iteration) {
		B_calls <- reads[[read]]$V4 >= threshold
		sum(B_calls[((iteration - 1) * step + 1):((iteration - 1) * step + window)]) / window}
	sapply(1:((floor((nrow(reads[[read]]) - window) / step)) + 1), internal_fun)}
calc_list <- sapply(1:length(reads), BBT_fun)
read_index_starts <- sapply(1:length(reads), function(read) {reads[[read]]$V2[step * (1:((floor((nrow(reads[[read]]) - window) / step)) + 1)) - step + 1]})
read_index_ends <- sapply(1:length(reads), function(read) {reads[[read]]$V2[step * (1:((floor((nrow(reads[[read]]) - window) / step)) + 1)) - step + window]})
read_index_mids <- sapply(1:length(reads), function(read) {reads[[read]]$V2[step * (1:((floor((nrow(reads[[read]]) - window) / step)) + 1)) - step + (window / 2)]})
if(export_mode == "FALSE" & is.na(ROI_directory)) {
export_list <- lapply(1:length(calc_list), function(read) cbind(chromosomes[read], read_index_starts[[read]], read_index_ends[[read]], round(calc_list[[read]],6), read_index_mids[[read]]))
invisible(sapply(1:length(export_list), function(read) {write.table(export_list[[read]], file = paste0(save_directory, "/", readIDs[read], "_", chromosomes[read], "_", strands[read], ".BBT.bedgraph"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)}))}
if(export_mode == "TRUE" & is.na(ROI_directory)) {
export_table <- cbind(rep(readIDs, sapply(calc_list, length)), Reduce(c, read_index_mids), Reduce(c, calc_list))
write.table(export_table, file = paste0(save_directory, "/BBT.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("read_IDs", "window_mid", "B/(B+T)"))}
