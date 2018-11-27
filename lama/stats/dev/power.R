library('pwr')


print(pw)

# U will always be 2 (2 predictors)




#try altering the effect size
power = c()
spec_nums = c()

for (n_specimens in seq(6, 100)){
  v = n_specimens - 3
  pw = pwr.f2.test(u=2, v=v, f2=0.8, sig.level = 0.05, power=NULL)
  power = c(power, pw$power)
  spec_nums = c(spec_nums, n_specimens)
}

plot(x=spec_nums, y=power)
abline(h=0.8)
title('organ volume analysis power. sig=0.05, d=0.8')