args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
data_in <- read.table(input_file)
out <- p.adjust(unlist(data_in), method='BH')
write.table(out, file= output_file, row.names = FALSE, col.names = FALSE)
