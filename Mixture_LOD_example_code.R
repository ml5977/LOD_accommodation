###R packages for use
rpackage.list <- c("corrplot", "ggplot2","kableExtra","nlme","MASS","Matrix","mgcv",
                   "qgam","quantreg","splines2","survival","tableone","dplyr","gWQS","gglasso",
                   "qgcomp","bkmr","glmnet","caret","mice","tmvtnorm","tidyverse","devtools",
                   "janitor","Hmisc","Matrix","mvtnorm","magic","truncnorm","ranger","tmvtnorm")
lapply(rpackage.list, require, character.only = TRUE)

#set.seed
set.seed(020123)


###True & simulation setting
n <- 500 #sample size
p <- 10 #number of exposures
lod.p <- 0.30 #Percentage of missing due to LOD
true.beta <- c(rep(1,3), 1, 0.8, 0, 0.6, 0.4, 0, 0.2, 0.1, 0, 0)
corr1 <- 0.25 #correlation between exposures within a group; group1
corr2 <- 0.75 #group2
corr3 <- 0.50 #group3
sig.aft <- (1/2)^2 #variance of AFT error term
sig.lr <- 2 #variance of outcome model error term
true.alpha <- matrix(c(-0.20, -0.35, -0.30, -0.25, -0.35, -0.25, -0.25, -0.40, -0.25, -0.30,
                        -0.25, -0.50, -0.25, -0.05, -0.03, -0.10, -0.25, -0.25, -0.50, -0.25,
                        -0.05, -0.02, -0.00, -0.25, -0.25, -0.25, -0.25, -0.50, -0.25, -0.25), ncol=10, byrow=T)


###Data Generation
#covariates X & error terms generation
mat1 <- matrix(corr1, ncol=3, nrow=3)
mat2 <- matrix(corr2, ncol=3, nrow=3)
mat3 <- matrix(corr3, ncol=4, nrow=4)
diag(mat1) <- diag(mat2) <- diag(mat3) <- 1
mat.t <- adiag(mat1, mat2, mat3);mat.t
ee1 <- rmvnorm(n, sigma=sig.aft*mat1)
ee2 <- rmvnorm(n, sigma=sig.aft*mat2)
ee3 <- rmvnorm(n, sigma=sig.aft*mat3)
ee.t <- cbind(ee1, ee2, ee3) #error 
#head(ee.t)

x1 <- rnorm(n, mean=1, sd=1)
x2 <-  rbinom(n, 1, prob=0.5)
x.t <- cbind(rep(1,n), x1, x2) #covariates X including intercept
#head(x.t)

#covariates Z generation
gen.t <- x.t %*% true.alpha + ee.t
true.z <- exp(-gen.t)
cor(true.z)
#cor(cbind(true.z, x1, x2))

#LOD value
lod.cut <- c(1.639771, 1.294595, 1.265859, 1.379562, 1.652190)
log.lod.cut <- log(lod.cut)

#Indicator for values < LOD
gen.z <- log(true.z)
gen.z[,2] <- ifelse(gen.z[,2]<=log.lod.cut[1],log.lod.cut[1],gen.z[,2])
gen.z[,3] <- ifelse(gen.z[,3]<=log.lod.cut[2],log.lod.cut[2],gen.z[,3])
gen.z[,5] <- ifelse(gen.z[,5]<=log.lod.cut[3],log.lod.cut[3],gen.z[,5])
gen.z[,7] <- ifelse(gen.z[,7]<=log.lod.cut[4],log.lod.cut[4],gen.z[,7])
gen.z[,9] <- ifelse(gen.z[,9]<=log.lod.cut[5],log.lod.cut[5],gen.z[,9])

ind.z <- matrix(1, ncol=p, nrow=n)
ind.z[,2] <- ifelse(gen.z[,2]<=log.lod.cut[1],0,ind.z[,2]) #0 if missing due to LOD, 1 if observed
ind.z[,3] <- ifelse(gen.z[,3]<=log.lod.cut[2],0,ind.z[,3])
ind.z[,5] <- ifelse(gen.z[,5]<=log.lod.cut[3],0,ind.z[,5])
ind.z[,7] <- ifelse(gen.z[,7]<=log.lod.cut[4],0,ind.z[,7])
ind.z[,9] <- ifelse(gen.z[,9]<=log.lod.cut[5],0,ind.z[,9])

#outcome Y generation
gen.y <- cbind(x.t, log(true.z)) %*% as.matrix(true.beta) + rnorm(n, 0, sqrt(sig.lr))

#generated dataset: gold.data and gen.data
gold.original.data <- data.frame(y=gen.y, x=x.t[,-1], z=true.z, ind=ind.z)
colnames(gold.original.data) <- c("y","x1","x2","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                                  "ind1","ind2","ind3","ind4","ind5","ind6","ind7","ind8","ind9","ind10")
gen.data <- data.frame(y=gen.y, x=x.t[,-1], z=gen.z, ind=ind.z)
colnames(gen.data) <- c("y","x1","x2","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                        "ind1","ind2","ind3","ind4","ind5","ind6","ind7","ind8","ind9","ind10")
covariate.names <- c("x1","x2","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10")


###LOD accommodation dataset
##gold analysis data
gold.data <- data.frame(y=gen.y, x=x.t[,-1], z=log(true.z), ind=ind.z)
colnames(gold.data) <- c("y","x1","x2","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                         "ind1","ind2","ind3","ind4","ind5","ind6","ind7","ind8","ind9","ind10")

##complete case data
complete.data <- gen.data
complete.data$ind.t <- complete.data$ind1*complete.data$ind2*complete.data$ind3*complete.data$ind4*
  complete.data$ind5*complete.data$ind6*complete.data$ind7*complete.data$ind8*complete.data$ind9*complete.data$ind10
complete.data <- subset(complete.data, complete.data$ind.t==1)

##LOD/sqrt(2) data
sqrt2.data <- gen.data
sqrt2.data$z2 <- ifelse(sqrt2.data$ind2==0, log(lod.cut[1]/sqrt(2)), sqrt2.data$z2)
sqrt2.data$z3 <- ifelse(sqrt2.data$ind3==0, log(lod.cut[2]/sqrt(2)), sqrt2.data$z3)
sqrt2.data$z5 <- ifelse(sqrt2.data$ind5==0, log(lod.cut[3]/sqrt(2)), sqrt2.data$z5)
sqrt2.data$z7 <- ifelse(sqrt2.data$ind7==0, log(lod.cut[4]/sqrt(2)), sqrt2.data$z7)
sqrt2.data$z9 <- ifelse(sqrt2.data$ind9==0, log(lod.cut[5]/sqrt(2)), sqrt2.data$z9)

##Conventional MI with outcome
imput.data <- gen.data
imput.data$z2 <- ifelse(imput.data$ind2==0, NA, imput.data$z2)
imput.data$z3 <- ifelse(imput.data$ind3==0, NA, imput.data$z3)
imput.data$z5 <- ifelse(imput.data$ind5==0, NA, imput.data$z5)
imput.data$z7 <- ifelse(imput.data$ind7==0, NA, imput.data$z7)
imput.data$z9 <- ifelse(imput.data$ind9==0, NA, imput.data$z9)
temp.imput.data <- as.matrix(imput.data[,c(covariate.names, "y")])
#norm.boot
imput.data <- mice(temp.imput.data, m=5, maxit=100, metho="norm.boot", seed=100, printFlag = FALSE)

###truncated MI
ind.data <- gen.data[,c("ind1","ind2","ind3","ind4","ind5","ind6","ind7","ind8","ind9","ind10")]
imput.trunc.data <- list()
kk <- 1
loop.num <- c(2, 3, 5, 7, 9)
for(kk in 1:5){
  imput.trunc.data[[kk]] <- complete(imput.data, kk)
  imput.trunc.data[[kk]] <- cbind(imput.trunc.data[[kk]], ind.data)
  ll <- 1
  for(ll in 1:100){
    k <- 2
    for(k in loop.num){
      fit.temp <- lm(imput.data$formulas[[k+2]], data=imput.trunc.data[[kk]])
      imput.trunc.data[[kk]][,(k+2)][imput.trunc.data[[kk]][,(k+13)]==0] <- 
        rtruncnorm(1, a=-Inf, b=log.lod.cut[which(loop.num==k)], 
                   mean=fit.temp$fitted.values[imput.trunc.data[[kk]][,(k+13)]==0], 
                   sd=sd(fit.temp$residuals))
    }
    #print(paste0(ll," th iteration at ",kk))
  }
}
print("Truncated MI finished")

#AFT without outcome
aft.data <- gen.data
aft.data[,covariate.names[c(3:12)]] <-  exp(aft.data[,covariate.names[c(3:12)]])
lod.var.names <- c("z2","z3","z5","z7","z9")
lod.n <- length(lod.var.names)
lod.values <- lod.cut
#AFT model fit
n <- nrow(aft.data)
aft.coef.mat <- matrix(NA, ncol=3, nrow=lod.n)
aft.resi.mat <- matrix(NA, ncol=lod.n, nrow=n)
aft.mean.mat <- matrix(NA, ncol=lod.n, nrow=n)
lod.var.num <- c(2, 3, 5, 7, 9)
for(k in lod.var.num){
  TT <- -log(aft.data[,paste("z",k,sep="")])
  delta <- aft.data[,paste("ind",k,sep="")]
  aft.reg <- survreg(Surv(exp(TT), delta) ~ x1 + x2, dist="lognormal", data=aft.data)
  aft.coef.mat[which(k==lod.var.num),] <- aft.reg$coefficients
  aft.resi.mat[,which(k==lod.var.num)] <- TT-as.matrix(cbind(1, aft.data[,c("x1","x2")]))%*%matrix(aft.reg$coefficients)
  aft.mean.mat[,which(k==lod.var.num)] <- as.matrix(cbind(1, aft.data[,c("x1","x2")]))%*%matrix(aft.reg$coefficients)
}

#estimated covariance matrix of residuals: pairwise.complete.obs setting
aft.resi.mat2 <- aft.resi.mat
aft.resi.mat2 <- ifelse(aft.data[,c("ind2","ind3","ind5","ind7","ind9")]==1,aft.resi.mat2, NA)
sigma <- cov(aft.resi.mat2, use="pairwise.complete.obs")
#imputation
imputed.res.mat <- aft.resi.mat2
delta.mat <- aft.data[,c("ind2","ind3","ind5","ind7","ind9")]
for(k in 1:n){
  ind.delta <- delta.mat[k,]
  ind.res <- aft.resi.mat[k,]
  obs.d <- which(ind.delta==1)
  mis.d <- which(ind.delta==0)
  if(sum(ind.delta)==0){
    mu.cond <- rep(0, lod.n)
    cov.cond <- sigma
    tmu.cond <- mtmvnorm(mean=c(mu.cond), sigma=cov.cond, lower=ind.res[mis.d])$tmean
    imputed.res.mat[k,mis.d] <- tmu.cond
  }else if(sum(ind.delta)!=lod.n){
    mu.cond <- sigma[mis.d, obs.d]%*%solve(sigma[obs.d, obs.d])%*%(ind.res[obs.d])
    cov.cond <- sigma[mis.d, mis.d] - sigma[mis.d, obs.d]%*%solve(sigma[obs.d, obs.d])%*%sigma[obs.d, mis.d]
    tmu.cond <- mtmvnorm(mean=c(mu.cond), sigma=cov.cond, lower=ind.res[mis.d])$tmean
    tmu.cond <- ifelse(is.na(tmu.cond)==TRUE, ind.res[mis.d], tmu.cond)
    imputed.res.mat[k,mis.d] <- tmu.cond
  }
}
imputed.Z <- exp(-(aft.mean.mat + imputed.res.mat))
aft.data[,lod.var.names] <- imputed.Z
aft.data[,covariate.names[c(3:12)]] <- log(aft.data[,covariate.names[c(3:12)]])
print("AFT finished")


###Mixture analysis code
##Note: here I only provide example codes for mixture analysis (WQS, Elastic, BKMR) for gold standard dataset.

###ANALYSIS: WQS regression
##gold data
exposures <- covariate.names[c(3:length(covariate.names))]
wqs.gold <- gwqs(y ~ wqs + x1 + x2, mix_name=exposures, data=gold.data, q=4, 
                 validation=0.6, b=200, b1_pos=T, b1_constr=F, family="gaussian", seed=100)

###ANALYSIS: Elastic net
##gold data
x <- as.matrix(gold.data[,exposures])
z <- as.matrix(gold.data[,covariate.names[c(1:2)]]) 
xz <- cbind(x,z)
y <- c(gold.data$y)
# For penalty factors
is_penalized <- c(rep(1, ncol(xz[,1:10])), rep(0, ncol(xz[,11:12])))
# Use CV to select best alpha value
# Create a tuning grid of alpha and lambda values
cv.elnet.lambda <- cv.glmnet(xz, y, 
                             penalty.factor = is_penalized,
                             type.measure = "mse", alpha = 0.5)
egrid <- expand.grid(alpha = (1:10) * 0.1, 
                     lambda = cv.elnet.lambda$lambda)
# Create a tuning control for cv
control <- trainControl(method = "repeatedcv", 
                        repeats = 3, 
                        verboseIter = TRUE)
# Use the tuning grid and control to find best alpha
garbage <- capture.output( elnetFit <- train(x = xz, 
                                             y = y,
                                             method = "glmnet",
                                             penalty.factor = is_penalized,
                                             tuneGrid = egrid,
                                             trControl = control))
# For the optimal 'alpha', identify best lamdba and explore the results
# Use CV to find best lambda value
cv.elnet <- cv.glmnet(xz, y, 
                      penalty.factor = is_penalized,
                      type.measure = "mse", alpha = elnetFit$bestTune$alpha)
elnet.mod <- glmnet(xz, y, 
                    penalty.factor = is_penalized,
                    alpha = elnetFit$bestTune$alpha, lambda = cv.elnet$lambda.min)


###ANALYSIS: BKMR
##gold data
covar.gold <- data.frame(gold.data[,covariate.names[c(1:2)]])
expos.gold <- data.frame(gold.data[,covariate.names[c(3:length(covariate.names))]])
Y <- gold.data[,"y"]
colnames(covar.gold) <- covariate.names[c(1:2)]
colnames(expos.gold) <- covariate.names[c(3:length(covariate.names))]
bkmr_m.gold <- kmbayes(Y, Z=expos.gold, X=covar.gold, iter=500, verbose=FALSE, varsel=TRUE)

