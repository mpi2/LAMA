

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


mut <- c(21306.957, 23608.543, 22996.174, 20319.906, 26505.662, 26277.314, 22667.098, 21594.791)
wt <- c(10446.006, 7997.1055, 7916.9165, 7543.1392, 8095.1958, 8084.9414, 8735.0625, 9735.5205)


#wt <- c(rep(37, 11))
#mut <- c(rep(38, 3))

pixels <- c(mut, wt)

mat <- matrix(c(pixels, pixels), nrow = 16)


genotype <- c(rep('mutant', 8), rep('wildtype',8))
#sex <- c('m', 'f', 'm', 'f', 'm', 'f', 'm', 'm', 'm', 'f', 'm', 'f', 'f', 'f')
#sex <- c(rep('m', 11), rep('f', 3))

df <- data.frame(genotype=genotype)

fit <- lm(mat ~ df$genotype)

results <- pandt_vals(fit)

tvals <- results$tvals[2,]
pvals <- results$pvals[2,]

t.test(wt, mut)