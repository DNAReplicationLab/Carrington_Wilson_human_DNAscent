library(dplyr)
library(IncDTW)

#useage: Rscript DNAscent_read_summary.R path_to_folder_with_bedgraphs read_summary_output_file barrier_summary_output_file

args = commandArgs(trailingOnly=TRUE)

my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)

options(width = 200)
rDNA_repeat_length = 9137
rDNA_barrier = 9000
rDNA_barrier_positions = seq(rDNA_barrier, 100000, rDNA_repeat_length)

print(args[1])

files <- list.files(path=args[1], pattern="*.detect.bedgraph", full.names=TRUE, recursive=FALSE)

threshold = 0.5
read_summary <- data.frame()
barrier_summary <-  data.frame()

for (i in 1:length(files)) {
#for (i in 1:10) {

	#bedgraph_file <- "BEDgraphs_with_strand/1/ec1cfba0-bafa-4fe4-83c1-38f336cf50f7.detect.bedgraph" 
	bedgraph_file <- files[i] 

	print(bedgraph_file)

	header_of_interest <- read.table(file = bedgraph_file, sep = " ", nrows = 1)
	read_info <- strsplit(header_of_interest[1,3], "=")[[1]][2]
	readID <- substr(read_info, 1, nchar(read_info) - 3)
	read_strand <- substr(read_info, nchar(read_info) - 2, nchar(read_info))

	read_of_interest <- read.table(file = bedgraph_file,
								   sep = " ",
								   skip = 1
	)

	read_of_interest <- subset(read_of_interest, select = -c(V1,V3))

	names(read_of_interest)[1] <- "position"
	names(read_of_interest)[2] <- "prob"

	read_of_interest$hard_call <- ifelse(read_of_interest$prob >= threshold, 1, 0)

	read_of_interest$hard_cumsum <- cumsum(read_of_interest$hard_call)

	read_of_interest$lag_cumsum_diff <- read_of_interest$hard_cumsum-lag(read_of_interest$hard_cumsum, n = 1000)
	read_of_interest$lead_cumsum_diff <- lead(read_of_interest$hard_cumsum, n = 1000)-read_of_interest$hard_cumsum

	read_of_interest$down_step <- (read_of_interest$lag_cumsum_diff - read_of_interest$lead_cumsum_diff)

	if(max(read_of_interest$position) > ((read_of_interest[1,1]%/%rDNA_repeat_length)*rDNA_repeat_length)+rDNA_barrier) {
		rDNA_barrier_positions <- seq(((read_of_interest[1,1]%/%rDNA_repeat_length)*rDNA_repeat_length)+rDNA_barrier,max(read_of_interest$position),rDNA_repeat_length)
		}
		
	barrier_stat <- data.frame(matrix(ncol = 6))
	colnames(barrier_stat) <- c("read_index", "file_name", "read_name", "position", "upstream", "downstream")
	for (j in 1:length(rDNA_barrier_positions)) {
		barrier_stat[j,1] <- i
		barrier_stat[j,2] <- bedgraph_file
		barrier_stat[j,3] <- readID
		barrier_stat[j,4] <- read_of_interest[which.min(abs(read_of_interest$position-rDNA_barrier_positions[j])),1]
		barrier_stat[j,5] <- read_of_interest[which.min(abs(read_of_interest$position-rDNA_barrier_positions[j])),5]
		barrier_stat[j,6] <- read_of_interest[which.min(abs(read_of_interest$position-rDNA_barrier_positions[j])),6]
	}
	
	
	step_pos <- filter(read_of_interest[find_peaks(read_of_interest$down_step, w = 500, get_min = FALSE, strict = FALSE),], down_step==my.max(read_of_interest$down_step))

	step_pos <- step_pos[!duplicated(step_pos[,c('down_step')]),]
	
	if(length(step_pos$position) > 0) {
	
		tmp_df <- data.frame(bedgraph_file,
								readID, 
								read_strand, 
								step_pos$position, 
								step_pos$down_step, 
								max(read_of_interest$hard_cumsum, na.rm = TRUE), 
								length(read_of_interest$position), 
								max(read_of_interest$hard_cumsum, na.rm = TRUE)/length(read_of_interest$position),
								max(read_of_interest$lag_cumsum_diff, na.rm = TRUE)/1000,
								min(read_of_interest$lag_cumsum_diff, na.rm = TRUE)/1000
		)

		barrier_summary <- rbind(barrier_summary, barrier_stat)
		read_summary <- rbind(read_summary, tmp_df)
		
	}
	
}

names(read_summary)[1] <- "file_name"
names(read_summary)[2] <- "read_name"
names(read_summary)[3] <- "read_strand"
names(read_summary)[4] <- "step_position"
names(read_summary)[5] <- "step_size"
names(read_summary)[6] <- "total_BrdU"
names(read_summary)[7] <- "total_positions"
names(read_summary)[8] <- "fraction_BrdU"
names(read_summary)[9] <- "most_pos_window"
names(read_summary)[10] <- "least_pos_window"

#rownames(read_summary) <- make.names(read_summary[,1], unique=TRUE)
#read_summary[,1] <- NULL

write.table(read_summary, args[2], sep = "\t")
write.table(barrier_summary, args[3], sep = "\t")

## clean up
rm(barrier_stat, barrier_summary, bedgraph_file, i, files, threshold, bedgraph_file, read_of_interest, read_summary)