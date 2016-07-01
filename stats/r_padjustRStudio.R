args <- commandArgs(trailingOnly = TRUE);

args <- commandArgs(trailingOnly = TRUE);
pval_file <- '/home/neil/sig/LAMA_results/E14.5/300616_testing_cbx2_old_inputs/mut/output/stats/intensity/tempPvals.bin'
dim <- 4575567

con <- file(pval_file, "rb")
mat <- readBin(con, "double", n=dim, size=4)
close(con)
out <- p.adjust(mat, method='BH')
#writeBin(out, qvals_out)


png(filename='/home/neil/share/registration_projects/190315_CBX2/with_origins_set_to_0/stats/stats_300616/intensity/pvals.png')
hist(c(mat), nclass=100, main = "p-value distribution", xlab = 'p-value')
dev.off()