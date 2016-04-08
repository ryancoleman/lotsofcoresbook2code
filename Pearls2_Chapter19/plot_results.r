# Copyright (c) 2013-2014 Matthias Noack (ma.noack.pr@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Usage
usage <- "Rscript <this_file.r> <data_file> <output.pdf> (<column_name=average>)"
DEFAULT_COLUMN = "average";

FONT_FAMILY <- "Helvetica" # "Times" # for cairo pdf
HEIGHT_DATA = 0.25 # plot height per bar
X_LIM_STEP <- 10 # x-axis length, 5 more than max rounded to next 5
ROUND_DIGITS <- 2
TITLE <- "Hexciton Runtime per Iteration [ms]"
COL_BAR <- "grey"
DATA_SCALE <- 1 / 1000000 # input data is in ns, output data is in ms

# command line argument handling
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
	print("Not enough command line arguments. Usage:")
	print(usage)
	quit()
}

file <- args[1]
output_file <- args[2]

if (length(args) < 3) {
	column <- DEFAULT_COLUMN
} else {
	column <- args[3]
}

if (!file.exists(file)) {
	print("Error: input file not found.")
	quit()
}

# read data
data <- read.table(file, header=TRUE, sep="\t")

# check if column exists
if (!(column %in% colnames(data))) {
	print("Error: Column name not found in input file.")
	quit()
}

# projection to two columns: name and the one to plot
data_plot <- data[, c("name", column)]

# preprocess data
num_data <- nrow(data_plot)
bar_sizes <- round(data_plot[dim(data_plot)[1]:1,2] * DATA_SCALE, ROUND_DIGITS)
X_LIM_STEP <- max(bar_sizes) * 0.25;
x_max <- round(max(bar_sizes) / X_LIM_STEP) * X_LIM_STEP + X_LIM_STEP

# setup pdf output
cairo_pdf(output_file, width=20, height=20, family=FONT_FAMILY) # size = large enough
par(pin=c(3.1, 1.5))
par(pin=c(4, HEIGHT_DATA * num_data))

# create a horizontal bar plot
bplt <- barplot(bar_sizes, col=COL_BAR, horiz=TRUE, names.arg=data_plot[dim(data_plot)[1]:1,1], xlim=c(0, x_max), las=1)
title(TITLE, line=-0.25)

# add value labels, rounded to two digits
text(x=x_max, y=bplt, labels=paste(sprintf("%.2f", bar_sizes), "ms"), xpd=TRUE, adj=c(1,0.5))

# write to file
useless_output <- dev.off()

# crop into same file using the pdfcrop command
system(paste("pdfcrop " , output_file, " ", output_file))

