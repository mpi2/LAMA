library(RcppCNPy)
args <- commandArgs(trailingOnly = TRUE)
input_file <- '/tmp/pvals.npy'
output_file <- args[2]
py <- npyLoad(input_file)
out <- p.adjust(data_in, method='BH')
npySave(output_file, out, mode='w')
