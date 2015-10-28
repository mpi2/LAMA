args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
data_in <- read.csv(input_file, header=FALSE)[,1] # Convert to vector with [,1]
out <- p.adjust(unlist(data_in), method='BH')
write.table(out, file= output_file, row.names = FALSE, col.names = FALSE)
