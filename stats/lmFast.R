args <- commandArgs(trailingOnly = TRUE);

pixels_file <- args[1];
groups_file <- args[2];
pvals_out <- args[3];
formula <- args[4]; # Just foing one formula at a time for now


# Create a data frame of the groups
g <- read.table(groups_file, header=TRUE, sep=',')
groups <- data.frame(g)

modelPvalOnly <- function(fit) {
  # This get the p-value for the whole model
  f <- t(fit$fitted.values)
  if (attr(fit$terms, "intercept"))  {
    mss <- rowSums((f - rowMeans(f)) ^ 2)
    numdf <- fit$rank - 1
  } else {
    mss <- rowSums(f ^ 2)
    numdf <- fit$rank
  }
  
  resvar <- colSums(fit$residuals^2) / fit$df.residual
  fstat <- (mss / numdf) / resvar
  pval <- pf(fstat, numdf, fit$df.residual, lower.tail = FALSE)
  pval
}

con <- file(pixels_file, "rb")
dim <- readBin(con, "integer", 2)
mat <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
close(con)

g <- read.csv(groups_file, colClasses = "factor", header = FALSE)

formula_elements <- strsplit(formula, split=',')
fit <- lm(mat ~., data=groups[, unlist(formula_elements)])

#fit <- lm(mat ~ groups)

p <- modelPvalOnly(fit)

#toutCon <- file(tvals_out, "wb")
#writeBin(tvals, toutCon)
#close(toutCon)

poutCon <- file(pvals_out, "wb")
writeBin(p, poutCon)
close(poutCon)



