args <- commandArgs(trailingOnly = TRUE);

pixels_file <- args[1];
groups_file <- args[2];
pvals_out <- args[3];
tvals_out <- args[4];
formula <- args[5]; # Just foing one formula at a time for now


# Create a data frame of the groups
g <- read.table(groups_file, header=TRUE, sep=',')
groups <- data.frame(g)

pandt_vals <- function(fit) {
  # get estimates
  est <- fit$coefficients[fit$qr$pivot, ]
  
  # get R: see stats:::summary.lm to see how this is calculated
  p1 <- 1L:(fit$rank)
  R <- diag(chol2inv(fit$qr$qr[p1, p1, drop = FALSE]))
  
  # get residual sum of squares for each
  resvar <- colSums(fit$residuals^2) / fit$df.residual
  # R is same for each coefficient, resvar is same within each model 
  se <- sqrt(outer(R, resvar))
  
  tvals <- est / se
  pvals <- pt(abs(est / se), df = fit$df.residual, lower.tail = FALSE) * 2
  
  return(list(pvals=pvals, tvals=tvals))
}

con <- file(pixels_file, "rb")
dim <- readBin(con, "integer", 2)
mat <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
close(con)

#g <- read.csv(groups_file, colClasses = "factor", header = FALSE)

formula_elements <- strsplit(formula, split=',')
print('lm formula elements');
print(formula_elements)
fit <- lm(mat ~., data=groups[, unlist(formula_elements)])

#fit <- lm(mat ~ groups)

results <- pandt_vals(fit)

#toutCon <- file(tvals_out, "wb")
#writeBin(tvals, toutCon)
#close(toutCon)

poutCon <- file(pvals_out, "wb")
writeBin(results$pvals[2,], poutCon)
close(poutCon)

toutCon <- file(tvals_out, "wb")
writeBin(0 - results$tvals[2,], toutCon)
close(toutCon)


