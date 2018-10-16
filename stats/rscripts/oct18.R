# 081018
# Test whether crown-rump or whole embryo volume is the best variable to account for developmental differences at E14.5

df_crl <- read.csv('/home/neil/bit/LAMA_results/E14.5/paper_runs/output/staging_info.csv')

df_whole <- read.csv('/home/neil/bit/LAMA_results/E14.5/paper_runs/output/mask_volumes_08oct.csv')

df_organs <- read.csv('/home/neil/bit/LAMA_results/E14.5/paper_runs/output/organ_volumes_08oct.csv')

colnames(df_organs)[1] <- 'vol'
colnames(df_whole)[1] <- 'vol'
colnames(df_whole)[2] <- 'whole'
colnames(df_crl)[2] <- 'crl'


# Loop over the labels
df_test = merge(x=df_crl, y=df_whole, by.x='vol', by.y='vol')
df_test1 = merge(x=df_test, y=df_organs, by.x='vol', by.y='vol')

#fit_both <- lm(X100 ~ crl + whole, data=df_test1)
fit_crl <- lm(log(X100) ~ crl, data=df_test1)
fit_whole_log <- lm(log(X100) ~ log(whole), data=df_test1)
fit_whole_cubed <- lm(X100^(1/3) ~ I(whole^(1/3)), data=df_test1)

crl_bic <- BIC(fit_crl)
whole_bic_log <- BIC(fit_whole_log)
whole_biz_cub <- BIC(fit_whole_cubed)


