library(circular)
library(compiler)

args <- commandArgs(trailingOnly = TRUE);

pixels_file <- args[1];
groups_file <- args[2];
pvals_out <- args[3];
tvals_out <- args[4];
formula <- args[5]; # Just foing one formula at a time for now

# pixels_file <- '/home/neil/work/circular_stats/test_cir.dat';
# groups_file <- '/home/neil/work/circular_stats/groups.csv';
# pvals_out <- '/home/neil/work/circular_stats/pvals';
# tvals_out <- '/home/neil/work/circular_stats/tvals';
# formula <- 'genotype';


# Create a data frame of the groups
g <- read.table(groups_file, header=TRUE, sep=',')
groups <- data.frame(g)


con <- file(pixels_file, "rb")
dim <- readBin(con, "integer", 2)
mat <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
close(con)

#g <- read.csv(groups_file, colClasses = "factor", header = FALSE)

formula_elements <- strsplit(formula, split=',')
print('lm formula elements');
print(formula_elements)

circ = circular(mat, units = 'degrees')

f = factor(groups$genotype)

wwtest <- function(data){
  ww = watson.wheeler.test.default(data, group = f);
  c(ww$p.value, ww$statistic)
}

wwtest_comp <- cmpfun(wwtest)

result <- apply(mat, 2, wwtest_comp)
#results <- p,andt_vals(fit)



poutCon <- file(pvals_out, "wb")
writeBin(unlist(result[1,]), poutCon);
close(poutCon);

toutCon <- file(tvals_out, "wb");
writeBin(unlist(result[2,]), toutCon);
close(toutCon)


