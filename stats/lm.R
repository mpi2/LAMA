library(broom)
args <- commandArgs(trailingOnly = TRUE);
pixels_file <- args[1];
groups_file <- args[2];
tvals_out <- args[3];

con <- file(pixels_file, "rb")
dim <- readBin(con, "integer", 2)
mat <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
close(con)

g <- read.csv(groups_file, colClasses = "factor", header = FALSE)

groups <- apply(g, 2, factor) 

fit <- lm(mat ~ groups)

sum_lm1 <- summary(fit)

num_pixels <- dim(mat)[2]

result_pvalues <- numeric(num_pixels)

# This works
result_pvalues <- vapply(sum_lm1, function(x) x$coefficients[,4][2], FUN.VALUE = 1)

#result_pvalues <- do.call("cbind", lapply(summary(fit), function(f) coef(f)[, 4]))

outCon <- file(tvals_out, "wb")
writeBin(result_pvalues, outCon)
close(outCon)




