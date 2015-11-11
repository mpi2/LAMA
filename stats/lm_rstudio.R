library(broom)
library(microbenchmark)
#args <- commandArgs(trailingOnly = TRUE);
pixels_file <- '/tmp/tmp_data_for_lm';
groups_file <- '/tmp/tmp_groups_for_lm';
#pvals_out <- args[3];
tvals_out <- '/tmp/tmp_tvals_out';
#pixels_file <- '/tmp/raw_data_for_r.csv';
#groups_file <- '/tmp/groups_for_liear_model.csv';


#infile <- args[1];
#tvals_out <-'/tmp/tvals_lm.csv'

#infile <- '/tmp/tempdata.dat'
con <- file(pixels_file, "rb")
dim <- readBin(con, "integer", 2)
m <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
mat <- m[, 0:200]

close(con)

g <- read.csv(groups_file, colClasses = "factor", header = FALSE)

groups <- apply(g, 2, factor) 

#groups = factor(c(rep('wt', 10), rep('mut', 3)))
fit <- lm(mat ~ groups)

#result_pvalues <- do.call("cbind", lapply(summary(fit), function(f) coef(f)[, 4]))

broomFit <- function(fit){
  result_pvalues <- subset(tidy(fit), term == "groupswt")[, c(1,6)]
  
}

doCall <- function(fit){
  result_pvalues <- do.call("cbind", lapply(summary(fit), function(f) coef(f)[, 4]))
}

pvalOnly <- function(fit){
    pt(abs(coef(fit) / sqrt(diag(vcov(fit)))), 
       df = fit$df.residual, 
       lower.tail = FALSE) * 2
  
}

sm = summary(fit)

myMethod <- function(fit){
  return(result_pvalues <- lapply(summary(fit), function(x) x$coefficients[2,3:4]))
}

#sum_lm1 <- summary(lm1)

num_pixels <- dim(mat)[2]



# This works
#result_pvalues <- vapply(sum_lm1, function(x) x$coefficients[,4][2], FUN.VALUE = 1)


#result_pvalues <- lapply(sum_lm1, function(x) anova)



#outCon <- file(tvals_out, "wb")
#writeBin(result_pvalues[2, ], outCon)
#close(outCon)




