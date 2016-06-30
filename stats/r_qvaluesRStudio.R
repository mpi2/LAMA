library("qvalue")

args <- commandArgs(trailingOnly = TRUE);

#pval_file <- '/home/neil/sig/LAMA_results/E14.5/28um_wt_test_set_240616_multiple_defs/mutant_runs/cbx2_240616/output/stats/intensity/tempPvals.bin'
pval_file <- '/home/neil/work/290616_cbx2_tcp_testing/intensity/tempPvals.bin'
dim <- 4575567


con <- file(pval_file, "rb")
mat <- readBin(con, "double", n=dim, size=4)
close(con)
qobj <- qvalue(p = mat )
qvals <- qobj$qvalues

#writeBin(qvals, qvals_out)


