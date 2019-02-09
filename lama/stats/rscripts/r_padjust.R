args <- commandArgs(trailingOnly = TRUE);

pval_file <- args[1]
dim <- as.integer(args[2])
qvals_out <- args[3]

print('test=====')
print(pval_file)
print(dim)
print(qvals_out)

con <- file(pval_file, "rb")
mat <- readBin(con, "double", n=dim, size=4)
close(con)
out <- p.adjust(mat, method='BH')
writeBin(out, qvals_out)
