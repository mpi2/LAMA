library(broom)
args <- commandArgs(trailingOnly = TRUE);
pixels_file <- args[1];
groups_file <- args[2];
#tvals_out <- args[3];
pvals_out <- args[3];

con <- file(pixels_file, "rb")
dim <- readBin(con, "integer", 2)
mat <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
close(con)

g <- read.csv(groups_file, colClasses = "factor", header = FALSE)

groups <- apply(g, 2, factor) 

fit <- lm(mat ~ groups)

sm = summary(fit)


result_tvalues <- lapply(sm, function(x) x$coefficients[2,3])
tvals <- unname(unlist(lapply(result_tvalues, function(x) x[[1]][[1]])))

result_pvalues <- lapply(sm, function(x) x$coefficients[2,4])
pvals <- unname(unlist(lapply(result_pvalues, function(x) x[[1]][[1]])))

num_pixels <- dim(mat)[2]

#result_pvalues <- numeric(num_pixels)

# This works
#result_pvalues <- vapply(sum_lm1, function(x) x$coefficients[,4][2], FUN.VALUE = 1)

#result_pvalues <- do.call("cbind", lapply(summary(fit), function(f) coef(f)[, 4]))

toutCon <- file(tvals_out, "wb")
writeBin(tvals, toutCon)
close(toutCon)

poutCon <- file(pvals_out, "wb")
writeBin(pvals, poutCon)
close(poutCon)



