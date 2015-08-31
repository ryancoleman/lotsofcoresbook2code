# Copyright (c) 2013-2014 Matthias Noack (ma.noack.pr@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Usage
usage <- "Rscript <this_file.r> <data_file_a> <data_file_b> <output.pdf> (<column_name=average>)"
DEFAULT_COLUMN = "average";

FONT_FAMILY <- "Helvetica" # "Times" # for cairo pdf
HEIGHT_DATA = 0.25 # plot height per bar
X_LIM_STEP <- 10 # x-axis length, 5 more than max rounded to next 5
ROUND_DIGITS <- 2
COL_BAR <- "grey"
DATA_SCALE <- 1 / 1000000 # input data is in ns, output data is in ms

# command line argument handling
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
	print("Not enough command line arguments. Usage:")
	print(usage)
	quit()
}

file_a <- args[1]
file_b <- args[2]
output_file <- args[3]

if (length(args) < 4) {
	column <- DEFAULT_COLUMN
} else {
	column <- args[4]
}

if (!file.exists(file_a) || !file.exists(file_b)) {
	print("Error: input file not found.")
	quit()
}

TITLE <- paste("Runtime diff: \"", file_a, "\" vs. \"", file_b, "\"", sep="")

# read data
data_a <- read.table(file_a, header=TRUE, sep="\t")
data_b <- read.table(file_b, header=TRUE, sep="\t")

# check if column exists
if (!(column %in% colnames(data_a)) || !(column %in% colnames(data_b))) {
	print("Error: Column name not found in input file.")
	quit()
}

# projection to two columns: name and the one to plot
data_plot_a <- data_a[, c("name", column)]
data_plot_b <- data_b[, c("name", column)]

# preprocess data
num_data_a <- nrow(data_plot_a)
num_data_b <- nrow(data_plot_b)
# check if same number of rows
if (num_data_a != num_data_b ) {
	print("Error: input files have differing row counts.")
	quit()
}

absolute_diff <- data_plot_b[dim(data_plot_b)[1]:1,2] - data_plot_a[dim(data_plot_a)[1]:1,2]
relative_diff <- 1 - data_plot_b[dim(data_plot_b)[1]:1,2] / data_plot_a[dim(data_plot_a)[1]:1,2]
bar_sizes <- round(absolute_diff * DATA_SCALE, ROUND_DIGITS)


X_LIM_STEP <- max(abs(bar_sizes)) * 0.5;
x_min <- round(min(bar_sizes) / X_LIM_STEP) * X_LIM_STEP - X_LIM_STEP
x_max <- round(max(bar_sizes) / X_LIM_STEP) * X_LIM_STEP + X_LIM_STEP

# setup pdf output
cairo_pdf(output_file, width=20, height=20, family=FONT_FAMILY) # size = large enough
par(pin=c(3.1, 1.5))
par(pin=c(4, HEIGHT_DATA * max(num_data_a, num_data_b)))

# create a horizontal bar plot
bplt <- barplot(bar_sizes, col=COL_BAR, horiz=TRUE, names.arg=data_plot_a[dim(data_plot_a)[1]:1,1], xlim=c(x_min, x_max), las=1)
title(TITLE, line=-0.25)

# add value labels
# ‘adj’ allows _adj_ustment of the text with respect to ‘(x, y)’.
# Values of 0, 0.5, and 1 specify left/bottom, middle and right/top
# alignment, respectively.  The default is for centered text, i.e.,
# ‘adj = c(0.5, NA)’
#text(x=x_max, y=bplt, labels=paste(sprintf("%.2f", bar_sizes), "ms", "=", sprintf("%.1f", relative_diff * 100), "%"), xpd=TRUE, adj=c(1,0.5))
text(x=x_max, y=bplt, labels=paste(sprintf("%.2f", bar_sizes), "ms"), xpd=TRUE, adj=c(1,0.5))
text(x=x_min, y=bplt, labels=paste(sprintf("%.1f", relative_diff * 100), "%"), xpd=TRUE, adj=c(0,0.5))

# write to file
useless_output <- dev.off()

# crop into same file using the pdfcrop command
system(paste("pdfcrop " , output_file, " ", output_file))

