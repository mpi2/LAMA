library(RcppArmadillo)
args <- commandArgs(trailingOnly = TRUE);
#pixels_file <- args[1];
#groups_file <- args[2];
#pvals_out <- args[3];
#tvals_out <- args[4];
#pixels_file <- '/tmp/raw_data_for_r.csv';
#groups_file <- '/tmp/groups_for_liear_model.csv';

startTime <- Sys.time()

#infile <- args[1];
tvals_out <-'/tmp/tvals_lm.csv'

infile <- '/tmp/tempdata.dat'
con <- file(infile, "rb")
dim <- readBin(con, "integer", 2)
mat <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])

print('reading took')
print(Sys.time() - startTime)
close(con)

groups = factor(c(rep('wt', 10), rep('mut', 3)))
startTime <- Sys.time()
lm1 <- lm(mat ~ groups)
print('doing lm took')
print(Sys.time() - startTime)



startTime <- Sys.time()
sum_lm1 <- summary(lm1)
print('summary took')
print(Sys.time() - startTime)

num_pixels <- dim(mat)[2]

result_pvalues <- numeric(num_pixels)

count <- 0;


startTime <- Sys.time()
#for (idx in 1:length(sum_lm1)){
#    val <- sum_lm1[[idx]]$coefficients[,4][2]
#    result_pvalues[idx] <- val;
#};

# This works
#result_pvalues <- vapply(sum_lm1, function(x) x$coefficients[,4][2], FUN.VALUE = 1)


result_pvalues <- lapply(sum_lm1, function(x) anova)

print('t loop took')
print(Sys.time() - startTime)

write.table(result_pvalues, tvals_out, sep=',');




