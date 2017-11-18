args <- commandArgs(trailingOnly = TRUE);

args <- commandArgs(trailingOnly = TRUE);
pval_file <- '/home/neil/sig/LAMA_results/E14.5/120716_E14.5_14um_test_set/mutant_runs/ATP2a1/output/stats/jacobians_brain_norm/tempPvals.bin'
dim <- 36605038

con <- file(pval_file, "rb")
mat <- readBin(con, "double", n=dim, size=4)
close(con)
out <- p.adjust(mat, method='BH')

min_p = min(mat)
min_q = min(out)
