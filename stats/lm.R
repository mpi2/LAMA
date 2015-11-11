library(RcppArmadillo)
args <- commandArgs(trailingOnly = TRUE);
pixels_file <- args[1];
groups_file <- args[2];
#pvals_out <- args[3];
tvals_out <- args[4];
#pixels_file <- '/tmp/raw_data_for_r.csv';
#groups_file <- '/tmp/groups_for_liear_model.csv';


#infile <- args[1];
#tvals_out <-'/tmp/tvals_lm.csv'

#infile <- '/tmp/tempdata.dat'
con <- file(infile, "rb")
dim <- readBin(con, "integer", 2)
mat <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])


close(con)

groups = factor(c(rep('wt', 10), rep('mut', 3)))
lm1 <- lm(mat ~ groups)

sum_lm1 <- summary(lm1)

num_pixels <- dim(mat)[2]

result_pvalues <- numeric(num_pixels)

# This works
result_pvalues <- vapply(sum_lm1, function(x) x$coefficients[,4][2], FUN.VALUE = 1)


result_pvalues <- lapply(sum_lm1, function(x) anova)

write.table(result_pvalues, tvals_out, sep=',');




