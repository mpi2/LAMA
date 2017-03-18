infile <- '/tmp/tmp7FsdFV'
con <- file(infile, "rb")
dim <- readBin(con, "integer", 2)
mat <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
close(con)


genotype <- c(rep('mutant', 14), rep('wildtype',3))

df <- data.frame(genotype=genotype)

fit <- lm(mat ~ df$genotype)