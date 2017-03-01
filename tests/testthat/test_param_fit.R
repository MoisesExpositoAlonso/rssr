context("Fitting parameters")

data(betahat)

data(se)

data(alpha_test)
data(mu_test)
data(R_shrink)

betahat <- c(betahat)
se <- c(se)
mu_test <- c(mu_test)
alpha_test <- c(alpha_test)
t_SiRiS <- SiRSi(R_shrink,1/se)
SiRiS_f <- as.matrix(SiRSi(R_shrink,1/se))
p <- length(betahat)
SiRiSr=c(SiRiS_f%*%(alpha_test*mu_test))


sigb <- 1
log10oddsvec <- seq(-6,-1,0.5)
logodds <- log10oddsvec*log(10)
# input_h5files <- dir("/media/nwknoblauch/Data/GTEx/1kg_LD",full.names = T)[1]
# output_h5files <- gsub("1kg_LD","1kg_IBD",input_h5files)
# file.remove(output_h5files)

SiRiS  <-as(SiRiS_f,"dgCMatrix")
toutfile <- tempfile()
options <- list(plan=list(engine="MC",resources=list(nodes=2)),
                datafile=toutfile,
                betahat=betahat,
                se=se,
                logodds=logodds,
                toFile=toutfile,
                alpha=alpha_test,
                mu=mu_test,
                SiRiS=SiRiS,
                sigb=sigb,
                datafile=tempfile(),
                tolerance=1e-4)
datafiles <- tempfile()
mmat <- rss_varbvsr_parallel_future(options=options)



lnzm <- gen_lnzmat(mmat[[1]],logoddsvec = logodds,sigbvec = sigb)
test_that("We can estimate a known value of pi by grid search",
          expect_equal(c(marg_pi(log10oddsvec,c(lnzm))),(10/p),tolerance=0.01)
)
lnzm
#file.remove(toutfile)
sigb <- seq(0.5,1.5,0.1)
log10oddsvec <- log((10/p)/(1-(10/p)),base = 10)
logodds <- log10oddsvec*log(10)
options <- list(plan=list(engine="MC",resources=list(nodes=2)),
                betahat=betahat,
                se=se,
                logodds=logodds,
                toFile=toutfile,
                alpha=alpha_test,
                mu=mu_test,
                SiRiS=SiRiS,
                sigb=sigb,
                datafile=tempfile(),
                tolerance=1e-4)
datafiles <- tempfile()
mmat <- rss_varbvsr_parallel_future(options=options)
lnzm <- gen_lnzmat(mmat[[1]],logoddsvec = logodds,sigbvec=sigb)
gsigb <- marg_param(c(lnzm),param = sigb)
test_that("We can estimate a value of sigb by grid search",
          expect_equal(gsigb,1,tolerance=0.1))
tlnzmat <- numeric(length(logodds))
ttSiRiS <- SiRiS
t_R_shrink <- R_shrink

SiRiS = SiRSi(R_shrink,1/se)

SiRiSr=c((SiRiS%*%(alpha_test*mu_test))@x)

update_alpha_mu_alt <- function(SiRiS,
                                sigma_beta,logodds,betahat,se,alpha0,mu0,SiRiSr,reverse){
  sigma_beta_square <- sigma_beta[1]*sigma_beta[1]
  se_square <- se*se
  sigma_square <- (se_square * sigma_beta_square) / (se_square + sigma_beta_square)
  # tSiRiS <- SiRiS
  # diag(tSiRiS) <- 0
  sigmoid <- function(x) {
    return(1/(1 + exp(-x)))
  }
  alpha = sigmoid(logodds[1] + 0.5 * (log(sigma_square/(sigma_beta_square)) +  mu0 * mu0 / sigma_square))
  r=alpha0*mu0
  mu <- sigma_square*(betahat/se_square-(SiRiS%*%r-r/se_square))@x
  SiRiSr=(SiRiS%*%(alpha*mu))@x
  return(list(alpha1=alpha,mu1=mu,SiRiSr=SiRiSr))
}
i <- 1
# fam_alt <- update_alpha_mu_alt(alpha_test,mu_test,se,sigb,logodds,betahat,SiRiS)
fstra <- update_alpha_mu_alt(SiRiS = SiRiS,sigma_beta = 1,
                             logodds = logodds[i],betahat = betahat,se = se,
                             alpha = alpha_test,mu = mu_test,SiRiSr = SiRiSr,reverse = F)


fstra <- wrap_rss_varbvsr_iter(SiRiS = SiRiS,sigma_beta = 1,
                      logodds = logodds[i],betahat = betahat,se = se,
                      alpha = alpha_test,mu = mu_test,SiRiSr = SiRiSr,reverse = F)

stra <- wrap_rss_varbvsr_iter_alt(SiRiS = SiRiS,sigma_beta = 1,
                                  logodds = logodds[i],betahat = betahat,se = se,
                                  alpha = alpha_test,mu = mu_test,SiRiSr = SiRiSr,reverse = F)

expect_equal(fstra$mu1,stra$mu1)
expect_equal(fstra$SiRiSr,stra$SiRiSr)
expect_equal(fstra$alpha1,stra$alpha1)



for(i in 1:length(tlnzmat)){

  tres  <- rss_varbvsr_squarem_alt(SiRiS = SiRiS,
                          sigma_beta = 1,
                          logodds = logodds[i],
                          betahat = betahat,
                          se = se,
                          talpha0 = alpha_test,
                          tmu0 = mu_test,
                          tSiRiSr0 = SiRiSr,
                          tolerance = 1e-4,
                          itermax = 100,
                          verbose = T,
                          lnz_tol = T)
  
  
  stres <-  rss_varbvsr_squarem(SiRiS = SiRiS,
                                    sigma_beta = 1,
                                    logodds = logodds[i],
                                    betahat = betahat,
                                    se = se,
                                    talpha0 = alpha_test,
                                    tmu0 = mu_test,
                                    tSiRiSr0 = SiRiSr,
                                    tolerance = 1e-4,
                                    itermax = 200,
                                    verbose = T,
                                    lnz_tol = T)

}



