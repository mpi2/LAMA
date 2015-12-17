

pvalOnly2 <- function(fit) {
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



#wt <- c(10,11,15,16,15,12,14,11,19,15,17)
#mut <- c(50, 55,48)

wt <- c(rep(37, 11))
mut <- c(rep(38, 3))

pixels <- c(wt, mut)

mat <- matrix(c(pixels, pixels), nrow = 14)


genotype <- c(rep('wt', 3), rep('mut',11))
#sex <- c('m', 'f', 'm', 'f', 'm', 'f', 'm', 'm', 'm', 'f', 'm', 'f', 'f', 'f')
#sex <- c(rep('m', 11), rep('f', 3))

df <- data.frame(genotype=genotype)

fit <- lm(mat ~ df$genotype)

results <- pvalOnly2(fit)

tvals <- results$tvals[2,]
pvals <- results$pvals[2,]