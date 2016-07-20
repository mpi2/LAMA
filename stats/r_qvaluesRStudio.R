library("qvalue")

args <- commandArgs(trailingOnly = TRUE);

#pval_file <- '/home/neil/sig/LAMA_results/E14.5/28um_wt_test_set_240616_multiple_defs/mutant_runs/cbx2_240616/output/stats/intensity/tempPvals.bin'
pval_file <- '/home/neil/sig/LAMA_results/E14.5/120716_E14.5_14um_test_set/mutant_runs/ATP2a1/output/stats/jacobians_4_SPECIMENS/tempPvals.bin'

dim=36605038

con <- file(pval_file, "rb")
mat <- readBin(con, "double", n=dim, size=4)
close(con)
qobj <- qvalue(p = mat )
qvals <- qobj$qvalues
q_min = min(qvals)
p_min = min(mat)

#writeBin(qvals, qvals_out)


