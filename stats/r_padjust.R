args <- commandArgs(trailingOnly = TRUE);

pval_file <- args[1]
dim <- as.integer(args[2])
qvals_out <- args[3]
pvals_figure_out <- args[4]

con <- file(pval_file, "rb")
mat <- readBin(con, "double", n=dim, size=4)
close(con)
out <- p.adjust(mat, method='BH')
writeBin(out, qvals_out)

png(filename=pvals_figure_out)
hist(c(mat), nclass=100, main = "p-value distribution", xlab = 'p-value')
dev.off()
