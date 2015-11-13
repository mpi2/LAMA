
pixels_file <- '/tmp/tmp_data_for_lm';
groups_file <- '/tmp/tmp_groups_for_lm'
pvals_out <- '/tmp/pvals_out.dat'

modelPvalOnly <- function(fit) {
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

groups <- apply(g, 2, factor) 

fit <- lm(mat ~ groups)

p <- modelPvalOnly(fit)

#toutCon <- file(tvals_out, "wb")
#writeBin(tvals, toutCon)
#close(toutCon)

#poutCon <- file(pvals_out, "wb")
#writeBin(p, poutCon)
#close(poutCon)




