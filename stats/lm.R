args <- commandArgs(trailingOnly = TRUE)
pixels_file <- args[1]
groups_file <- args[2]

pixels <- read.csv(pixels_file, sep=',',header=FALSE) # Convert to vector with [,1]
groups <- read.csv(groups_file, header=TRUE)

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


df <- data.frame(pixels, groups)

pvals <- c()

for (i in 2:length(df) -1){
  lm1 = lm(df[,i] ~ df$Sex);
  pval = lmp(lm1);
  pvals[i] <- pval;
};

print(pvals);
