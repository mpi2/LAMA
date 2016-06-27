

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


mut <- c(95.99295, 82.011131, 110.04256, 95.991943)
wt <- c(10.910349, 14.313369, 48.914097, 21.896374, 20.247526, 5.489655, -0.1712565, 15.862488, 38.272747, 17.654718, 8.3408623, 20.005785, 20.468494, 91.264404, 10.043176, 32.852016, 15.724044, 14.639626, 8.5696468, 9.4951048)


#wt <- c(rep(37, 11))
#mut <- c(rep(38, 3))

pixels <- c(mut, wt)

mat <- matrix(c(pixels, pixels), nrow = 24)


genotype <- c(rep('mutant', 20), rep('wildtype',4))
#sex <- c('m', 'f', 'm', 'f', 'm', 'f', 'm', 'm', 'm', 'f', 'm', 'f', 'f', 'f')
#sex <- c(rep('m', 11), rep('f', 3))

df <- data.frame(genotype=genotype)

fit <- lm(mat ~ df$genotype)

#fit <- aov(mat ~ df$genotype)

results <- pandt_vals(fit)

tvals <- results$tvals[2,]
pvals <- results$pvals[2,]


#t.test(wt, mut)


summary(fit)