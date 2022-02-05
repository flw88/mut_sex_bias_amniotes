#!/usr/bin/env Rscript

source("alpha_regression.funcs.R")
load("lm_res/Mammals.Homo_sapiens.RData")

mt <- "mod"

alpha.samp <- c()
for(sp in names(lm.list[[mt]])){
  cur.lm <- lm.list[[mt]][[sp]]
  samp.coef <- samp.coef.list[[mt]][[sp]]
  
  # Compute 
  mean.gc <- cur.lm$data[region == "XZ", mean(gc)]
  a.sub <- samp.coef[,"(Intercept)"] + samp.coef[,"gc"]*mean.gc + samp.coef[,"I(gc^2)"]*(mean.gc^2)
  x.sub <- a.sub + samp.coef[,"regionXZ"]
  
  cur.alpha <- CalcAlpha(exp(x.sub) / exp(a.sub))

  alpha.samp <- cbind(alpha.samp, matrix(cur.alpha, ncol=1, dimnames=list(c(), sp)))
}

alpha.true <- alpha.dat[mut_type==mt, alpha]
names(alpha.true) <- alpha.dat[mut_type==mt, species]
alpha.true <- alpha.true[colnames(alpha.samp)]

rsq <- apply(alpha.samp, 1, function(x) summary(lm(alpha.true ~ x))$r.squared)

q <- quantile(rsq, c(0.025, 0.5, 0.975))

p <- ggplot(data=data.frame(rsq=rsq)) + geom_density(aes(x=rsq))
ggsave(p, filename="pdfs/alpha_resample.Mammals.rsq.pdf", width=4, height=4)

cat(str_interp("Mean r2 between estimated alpha and resamples: $[0.3f]{mean(rsq)}\n"))
cat(str_interp("Median + 95% CI: $[0.3f]{q[2]} ($[0.3f]{q[1]}, $[0.3f]{q[3]})\n"))
### OUTPUT ###
# Mean r2 between estimated alpha and resamples: 0.875
# Median + 95% CI: 0.885 (0.755, 0.949)


