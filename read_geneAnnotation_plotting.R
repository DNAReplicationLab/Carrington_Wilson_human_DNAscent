#R script intended to be used in Rstudio (for now).
#Inputs: 1) csplit .detect files where the first line of each chunk contains the read header starting with ">"
#        2) full gene annotation file (.gtf) e.g. https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz
#Output: pdf file containing B/(B+T) profile of a single read, with any overlapping gene annotations

#load data.table library for faster importing of read data
library(data.table)

#set parameters
file_directory <- "~/Desktop/fork direction metagene/old stuff/gradient2_sample_reads/hela1/"
gtf_file <- "~/Desktop/gradients/gencode.v38.primary_assembly.annotation.gtf.gz"
fork_calls_directory <- "/Users/jamie/Desktop/fork direction metagene/old stuff/fork_origin_calls"
OKseq_initiation_zone_file <- "~/Desktop/gradients/InitiationZones.Hela.hg38.bed"
ORM_initiation_zone_file <- "~/Desktop/gradients/ORM_initiation_zones_liftOver_sorted_hg38.bed"
save_directory <- "~/Downloads/test_figures/"
window = 290
step = 290
threshold = 0.5

#tests for existing save directory and creates one if it doesn't exist
ifelse(test = dir.exists(save_directory),
       yes = cat("WARNING: Your output directory already exists, figures created by this script will be written there\n"),
       no = dir.create(save_directory))

#load inputs
files <- list.files(path = file_directory,
                    pattern = "xx",
                    full.names = TRUE)

  #csplit puts the header lines of the whole .detect into the first chunk, typically called 'xx00000', this removes that file
if(as.numeric(strsplit(x = files[1], split = "xx")[[1]][2]) == 0) {
  files <- files[-1]
}

  #load in full gene annotation set
gtf <- fread(file = gtf_file,
             data.table = FALSE,
             skip = 5,
             select = c(1,3,4,5,7,9))

  #remove redundant features (that would be plotted on top of each other)
gtf_unique <- !duplicated(x = gtf[,1:5])
gtf <- gtf[gtf_unique,]
  #retain feature types that will be plotted (the rest of the features aren't informative)
gtf <- gtf[gtf[,2] == "transcript" | gtf[,2] == "exon" | gtf[,2] == "CDS",]

  #load in replication fork, origin, and termini calls
fork_call_files <- list.files(path = fork_calls_directory, pattern = "fork", full.names = TRUE)
left_forks <- rbind(fread(file = fork_call_files[1], data.table = FALSE),
                    fread(file = fork_call_files[2], data.table = FALSE))
right_forks <- rbind(fread(file = fork_call_files[3], data.table = FALSE),
                     fread(file = fork_call_files[4], data.table = FALSE))

origin_call_files <- list.files(path = fork_calls_directory, pattern = "origin", full.names = TRUE)
origin_calls <- rbind(fread(file = origin_call_files[1], data.table = FALSE),
                      fread(file = origin_call_files[2], data.table = FALSE))

termini_call_files <- list.files(path = fork_calls_directory, pattern = "termini", full.names = TRUE)
termini_calls <- rbind(fread(file = termini_call_files[1], data.table = FALSE),
                       fread(file = termini_call_files[2], data.table = FALSE))

OKseq_initiation_zones <- fread(file = OKseq_initiation_zone_file, data.table = FALSE)
ORM_initiation_zones <- fread(file = ORM_initiation_zone_file, data.table = FALSE)

  #load in .detect read header for each read
headers <- sapply(X = files,
                  FUN = scan,
                  what = "character",
                  sep = " ",
                  nlines = 1,
                  quiet = TRUE)

  #extract read ID, chromosome, read start, read end, and read strand from header lines
chromosomes <- as.character(headers[2,])
starts <- as.numeric(headers[3,])
ends <- as.numeric(headers[4,])
strands <- as.character(headers[5,])
strands <- ifelse(test = strands == "fwd",
                  yes = "+",
                  no = ifelse(test = strands == "rev",
                              yes = "-",
                              no = "ERROR"))
read_IDs <- sapply(X = strsplit(x = headers[1,], split = ">"),
                   FUN = '[[', 2)

  #load in read data. Each element of the list corresponds to a read
reads <- lapply(X = files,
                FUN = fread,
                data.table = FALSE,
                skip = 1,
                select = c(1,2))

#Calculate B/(B+T)
BBT_fun <- function(read) {
  #record location of all thymidines and BrdUs per read
  Ts <- reads[[read]][,1]
  Bs <- reads[[read]][,2] >= threshold
  
  #adjustment calculated so that trimming due to windowing occurs equally at each end of the read
  centre_adjust <- floor((length(Ts) - (floor((length(Ts) - window) / step) * step + window)) / 2)
  
  #function that calculates B/(B+T) by sliding over reads, per read
  slide_fun <- function(iteration) {
    #window that iterates (slides) along the read, within which the BrdUs are counted and dividied by the count of thymidines
    sliding_window <- step * (iteration - 1) + centre_adjust
    #output from sliding function, 1) coordinate of central thymidine in window, 2) B/(B+T) in window, 3) left-most thymidine in window, 4) right-might thymidine in window
    c(Ts[window / 2 + sliding_window],
      sum(Bs[(1 + sliding_window):(window + sliding_window)]) / window,
      Ts[1 + sliding_window],
      Ts[window + sliding_window])
  }
  #sliding function is applied per read
  sapply(X = 1:((floor((nrow(reads[[read]]) - window) / step)) + 1),
         FUN = slide_fun)
}

  #B/(B+T) function is applied to all reads
read_BBTs <- lapply(X = 1:length(reads),
                    FUN = function(read) {
                      as.data.frame(t(BBT_fun(read)))
                    })

#Define functions to plot graphics
  #plot a blank canvas upon which elements such as lines, rectangles (gene annotations), and axes can be drawn.
  #x limits are defined by the read start and end, taken from .detect per read header line
  #y limits are 0 -> 1, as this represents the proportion of B(B+T)
plot_blank <- function(read) {
  plot(NULL,
       xlim = c(starts[read], ends[read]),
       ylim = c(0,1),
       xlab = NA,
       ylab = NA,
       xaxt = 'n',
       yaxt = 'n',
       frame.plot = FALSE)
}

plot_probs <- function(read) {
  points(x = reads[[read]][,1],
         y = reads[[read]][,2],
         pch = 16,
         cex = 0.5,
         col = rgb(0,0,0,0.1))
}

  #plots simple curve, x values are the coordinates of the middle thymidine in each window, y values are the B(B+T) values calculated by the BBT function
plot_BBT <- function(read) {
  lines(x = read_BBTs[[read]][,1],
        y = read_BBTs[[read]][,2],
        col = 'blue')
}

plot_steps <- function(read) {
  length <- nrow(read_BBTs[[read]])
  segments(x0 = read_BBTs[[read]][,3],
           x1 = read_BBTs[[read]][,4],
           y0 = read_BBTs[[read]][,2],
           y1 = read_BBTs[[read]][,2],
           col = 'blue',
           lwd = 2)
  segments(x0 = read_BBTs[[read]][-length,4],
           x1 = read_BBTs[[read]][-1,3],
           y0 = read_BBTs[[read]][-length,2],
           y1 = read_BBTs[[read]][-1,2],
           col = 'blue',
           lwd = 2)
}

plot_fork_calls <- function(read) {
  chevron_width <- (ends[read] - starts[read]) / 100
  
  if(sum(left_forks[,4] == read_IDs[read]) > 0) {
    for(fork in which(left_forks[,4] == read_IDs[read])) {
      polygon(x = c(left_forks[fork, 2],
                    left_forks[fork, 2] + chevron_width,
                    left_forks[fork, 3],
                    left_forks[fork, 3] - chevron_width,
                    left_forks[fork, 3],
                    left_forks[fork, 2] + chevron_width),
                    y = c(0.5, 1, 1, 0.5, 0, 0))
    }
  }
  
  if(sum(right_forks[,4] == read_IDs[read]) > 0) {
    for(fork in which(right_forks[,4] == read_IDs[read])) {
      polygon(x = c(right_forks[fork, 2] + chevron_width,
                    right_forks[fork, 2],
                    right_forks[fork, 3] - chevron_width,
                    right_forks[fork, 3],
                    right_forks[fork, 3] - chevron_width,
                    right_forks[fork, 2]),
              y = c(0.5, 1, 1, 0.5, 0, 0))
    }
  }
}

plot_origins_termini <- function(read) {
  left_bound_origin <- origin_calls[origin_calls[,4] == read_IDs[read],2]
  right_bound_origin <- origin_calls[origin_calls[,4] == read_IDs[read],3]
  left_bound_termini <- termini_calls[termini_calls[,4] == read_IDs[read],2]
  right_bound_termini <- termini_calls[termini_calls[,4] == read_IDs[read],3]
  if(length(left_bound_origin != 0)) {
    rect(xleft = left_bound_origin,
         xright = right_bound_origin,
         ytop = 0.7,
         ybottom = 0.5,
         col = 'green')
  }
  if(length(left_bound_termini != 0)) {
    rect(xleft = left_bound_termini,
         xright = right_bound_termini,
         ytop = 0.7,
         ybottom = 0.5,
         col = 'red')
  }
}

plot_initiation_zones <- function(read) {
  OK_zones_left_bounds <- OKseq_initiation_zones[OKseq_initiation_zones[,1] == chromosomes[read],2]
  OK_zones_right_bounds <- OKseq_initiation_zones[OKseq_initiation_zones[,1] == chromosomes[read],3]
  ORM_zones_left_bounds <- ORM_initiation_zones[ORM_initiation_zones[,1] == chromosomes[read],2]
  ORM_zones_right_bounds <- ORM_initiation_zones[ORM_initiation_zones[,1] == chromosomes[read],3]
  if(length(OK_zones_left_bounds) != 0) {
    rect(xleft = OK_zones_left_bounds,
         xright = OK_zones_right_bounds,
         ytop = 1,
         ybottom = 0.5,
         col = "orange")
  }
  if(length(ORM_zones_left_bounds) != 0) {
    rect(xleft = ORM_zones_left_bounds,
         xright = ORM_zones_right_bounds,
         ytop = 0.5,
         ybottom = 0,
         col = "blue")
  }
}

  #adds x axis, with tick-marks every 10 kbp, and tick-mark labels every 100 kbp. x axis label shows chromosome
plot_x_axis <- function(read) {
  intervals <- seq(from = floor(x = starts[read] / 10000) * 10000,
                   to = ceiling(x = ends[read] / 10000) * 10000,
                   by = 10000)
  axis(side = 1,
       at = intervals,
       labels = ifelse(test = intervals %% 100000 == 0,
                       yes = formatC(format = 'f',
                                     x = intervals / 1000000,
                                     digits = 1),
                       no = NA))
  mtext(text = paste0("Genomic coordinates, ", chromosomes[read], " (Mbp). Read aligned to ", strands[read], "ve strand"),
        side = 1,
        line = 3)
}

  #adds y axis
plot_y_axis <- function() {
  axis(side = 2,
       at = seq(from = 0, to = 1, by = 0.2))
  mtext(text = "Proportion BrdU (B / (B + T)",
        side = 2,
        line = 3)
}

  #plots gene annotations
plot_genes <- function(read) {
  #gathers all annotations features that overlap with the read to be plotted
  overlapping_annotations <- gtf[chromosomes[read] == gtf[,1] &
                                 starts[read] < gtf[,4]       &
                                 ends[read] > gtf[,3],]
  #subsets overlapping features for 'transcript' features only
  overlapping_transcripts <- overlapping_annotations[overlapping_annotations[,2] == "transcript",]
  #gathers all unique gene IDs in overlapping features
  overlapping_gene_IDs <- unique(sapply(X = strsplit(x = overlapping_transcripts[,6],
                                                     split = "\""),
                                        FUN = "[[", 2))
  #gathers all common gene names in overlapping features
  overlapping_gene_names <- unique(sapply(strsplit(x = overlapping_transcripts[,6],
                                                   split = "\""),
                                          FUN = "[[", 8))
  #selects longest transcript (aka isoform) from each gene, output as transcript IDs.
  #selecting longest transcript per gene is up for debate, would be more valid to display most used isoform but this requires additional RNA-seq analysis per cell line.
  longest_transcripts <- sapply(X = overlapping_gene_IDs,
                                FUN = function(i) {
                                  #finds rows in overlapping transcripts that share a gene ID (i.e. all isoforms per gene)
                                  rows <- grep(pattern = i, x = overlapping_transcripts[,6])
                                  #calculates transcript length and selects largest value
                                  longest_transcript <- which.max((overlapping_transcripts[,4] - overlapping_transcripts[,3])[rows])
                                  #extracts transcript ID of longest transcript
                                  sapply(X = strsplit(x = overlapping_transcripts[,6][rows[longest_transcript]],
                                                      split = "\""),
                                         FUN = "[[", 4)
                                })
  #outputs list object per read, each element of list corresponds to feature details require to plot per gene, each row of each list element corresponds to transcript span, exons, and CDSs
  plotting_annotations <- lapply(X = longest_transcripts,
                                 FUN = function(transcript) {
                                   overlapping_annotations[grep(pattern = transcript,
                                                                x = overlapping_annotations[,6]),c(2,3,4,5)]
                                 })
  #if statement that asks if the read overlaps with any genes. If yes then graphical commands are executed, if not then they're skipped
  if(length(plotting_annotations) > 0) {
    #for loop iterates over number of genes that overlap read. Each loop plots a gene annotation
    for(gene in 1:length(plotting_annotations)) {
      #gene annotations frequently overlap, this object means that gene annotations are plotted on four different lines to reduce overlapping
      #four lines seem sufficient in nearly all cases, but if isn't then alter the rep command (e.g. rep(x = 1:5 / 6, times = 100) to distribute annotations over five lines)
      y_jitter <- rep(x = 1:4 / 5, times = 100)[gene]
      #this just stores the iteration of the loop, just to recall the gene name labels below
      i <- gene
      #this changes the iterations from a number series to list elements, just to improve readability of much of the below code
      gene <- plotting_annotations[[gene]]
      #this plots a simple rectangle to represent the span of the transcript/isoform
      rect(xleft = gene[gene[,1] == "transcript",2],
           xright = gene[gene[,1] == "transcript",3],
           ytop = y_jitter + 0.02,
           ybottom = y_jitter - 0.02,
           col = 'black',
           border = NA)
      #this adds the gene name label to each gene annotation
      #the median command is a neat trick to make sure the gene name label is plotted within the field of view. Some very long genes overlapped with reads, but extended well beyond the field of view, so their gene name labels would be as well.
      text(x = median(c((gene[gene[,1] == "transcript",2] + gene[gene[,1] == "transcript",3]) / 2,
                 starts[read],
                 ends[read])),
           y = y_jitter - 0.1,
           labels = overlapping_gene_names[i])
      #this if statement asks what strand the gene is on, and then plots the white chevrons that indicate the gene direction
      if(gene[gene[,1] == "transcript",4] == "+") {
        for(chevron in seq(from = gene[gene[,1] == "transcript",2],
                           to = gene[gene[,1] == "transcript",3],
                           by = (ends[read] - starts[read]) / 200)) {
          polygon(x = rep(x = c(chevron, chevron + ((ends[read] - starts[read]) / 1000)), times = 2),
                  y = c(y_jitter - 0.025, y_jitter, y_jitter + 0.025, y_jitter),
                  border = 'white',
                  lwd = 0.5)
        }
      }
      #this if statement does the same as above but for genes on the reverse strand
      if(gene[gene[,1] == "transcript",4] == "-") {
        for(chevron in rev(seq(from = gene[gene[,1] == "transcript",2],
                               to = gene[gene[,1] == "transcript",3],
                               by = (ends[read] - starts[read]) / 200))) {
          polygon(x = rep(x = c(chevron, chevron - ((ends[read] - starts[read]) / 1000)), times = 2),
                  y = c(y_jitter - 0.025, y_jitter, y_jitter + 0.025, y_jitter),
                  border = 'white',
                  lwd = 0.5)
        }
      }
      #this asks if the transcript has exonic features and then plots them as rectangles that are twice as thick as the transcript rectangle
      if(any(gene[,1] == "exon")) {
        rect(xleft = gene[gene[,1] == "exon",2],
             xright = gene[gene[,1] == "exon",3],
             ytop = y_jitter + 0.03,
             ybottom = y_jitter - 0.03,
             col = 'black',
             border = NA)
      }
      #this asks if the transcript has CDS features and then plots them as rectangles that are three times as thick as the transcript rectangle (exonic features without CDS features are UTRs)
      if(any(gene[,1] == "CDS")) {
        rect(xleft = gene[gene[,1] == "CDS",2],
             xright = gene[gene[,1] == "CDS",3],
             ytop = y_jitter + 0.04,
             ybottom = y_jitter - 0.04,
             col = 'black',
             border = NA)
      }
    }
  } 
}
#this command supplies the indexes of the reads which have a mean B(B+T) value of greater than 5%, and are 100 kbp or longer.
#I used this command to reduce the number of reads to plot for troubleshooting purposes, but have left it in in case its useful
#I recommend selecting reads for plotting before running this script on them
long_nascent_reads <- which(sapply(X = read_BBTs,
                                   FUN = function(read) {
                                     mean(read[,2])
                                   }) > 0.05 &
                            (ends - starts) > 100000)

#this for loop iterates, per read, over indexes that select elements of the read_BBT list object for plotting
#one file is written per read
for(read in long_nascent_reads) {
  #output is vector, adjusting width and height values here affectively changes line widths and text size
  pdf(file = paste0(save_directory, "read_", read, ".pdf"), width = 15, height = 9)
  #this establishes what proportion of the figure is taken up by gene annotations and read, currently set to 1:3. Increase 'times' argument to increase proportion of figure taken up by the read.
  layout(mat = matrix(c(rep(x = 1, times = 5), 2, 3, 4, rep(x = 5, times = 25))))
  #margin parameters set for gene annotation graphics. Important that the left and right values (the '5' and '2' in the vector) are the same as in the margin arguments for the read plot.
  par(mar = c(0,5,0,2) + 0.1)
  #plots blank canvas for gene annotation graphics
  plot_blank(read)
  #plots gene annotations that overlap with read of interest
  plot_genes(read)
  
  plot_blank(read)
  plot_initiation_zones(read)
  
  plot_blank(read)
  plot_origins_termini(read)
  
  plot_blank(read)
  plot_fork_calls(read)
  
  #margin parameters for read plot
  par(mar = c(5,5,0,2) + 0.1)
  #blank canvas for read plot
  plot_blank(read)
  #plots scatter plot where each point represents a thymidine and DNAscents probability score for BrdU
  plot_probs(read)
  #plots line graph coordinates of middle thymidine of windows against B/(B+T) values of that window
  #plot_BBT(read)
  #plots line graph in 'steps', using the 1st and last thymidine of each 290 thymidine window as the edges of each step
  plot_steps(read)
  #adds x axis graphics
  plot_x_axis(read)
  #adds y axis graphics
  plot_y_axis()
  #commits graphics to be written as .pdf file
  dev.off()
}

