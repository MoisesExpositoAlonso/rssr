#include <RcppEigen.h>
#include "rssvarbvsr.hpp"
#include "sigmoid.hpp"
#include <math.h>



// [[Rcpp::export]]
void rss_varbvsr_update(const double betahat,
                         const double se_square,
                         const double sigma_beta_square,
                         const double sigma_square,
                         const Eigen::ArrayXd &SiRiS_snp,
                         Eigen::ArrayXd &SiRiSr,
                         const double SiRiSr_snp,
                         const double logodds,
                         double &alpha,
                         double &mu) {
  
  
  // Update the variational estimate of the posterior mean.
  double r = alpha * mu;
  mu = sigma_square * (betahat / se_square + r / se_square - SiRiSr_snp);
  
  // Update the variational estimate of the posterior inclusion probability.
  double SSR = mu * mu / sigma_square;
  alpha = sigmoid(logodds + 0.5 * (log(sigma_square/(sigma_beta_square)) + SSR));
  
  // Update SiRiSr = inv(S)*R*inv(S)*r
  double r_new = alpha * mu;
  SiRiSr+=(SiRiS_snp*(r_new-r));
}


//In this method, each SNP can be updated in parallel, as they do not depend is update 

void rss_varbvsr_iter_alt(const Eigen::SparseMatrix<double> SiRiS,
                      const double sigma_beta,
                      const double logodds,
                      const Eigen::ArrayXd &betahat,
                      const Eigen::ArrayXd &se_square,
                      const Eigen::ArrayXd &sigma_square,
                      Eigen::ArrayXd &alpha,
                      Eigen::ArrayXd &mu,
                      Eigen::ArrayXd &SiRiSr,
                      bool reverse){
  
  double sigma_beta_square = sigma_beta*sigma_beta;
  
  size_t p=betahat.size();
  
  Eigen::ArrayXd  SiRiS_snp(p);
  Eigen::VectorXd  SiRiS_snp_v(p);
  

  Eigen::ArrayXd r(p);
  
  
//  double r = alpha * mu;
//  mu = sigma_square * (betahat / se_square + r / se_square - SiRiSr_snp);
  
  // Update the variational estimate of the posterior inclusion probability.
//  double SSR = mu * mu / sigma_square;
//  alpha = sigmoid(logodds + 0.5 * (log(sigma_square/(sigma_beta_square)) + SSR));
  
  // Update SiRiSr = inv(S)*R*inv(S)*r
//  double r_new = alpha * mu;
//  SiRiSr+=(SiRiS_snp*(r_new-r));
  Eigen::ArrayXd alpha0=alpha;
  Eigen::ArrayXd  mu0=mu;
  Eigen::ArrayXd SiRiSr0=SiRiSr;
  if(!reverse){
    for(size_t i=0; i<p; i++){
      r[i]=alpha0[i]*mu0[i];
      mu[i]=sigma_square[i]*(betahat[i]/se_square[i]+r[i]/se_square[i]-SiRiSr0[i]);
      //    mu = sigma_square * (betahat / se_square + r / se_square - SiRiSr_snp);
      // Update the variational estimate of the posterior inclusion probability.
      double SSR = mu0[i] * mu0[i] / sigma_square[i];
      
      
      alpha[i] = sigmoid(logodds + 0.5 * (log(sigma_square[i]/(sigma_beta_square)) + SSR));
      
      // Update SiRiSr = inv(S)*R*inv(S)*r
      double r_new = alpha[i] * mu[i];
      SiRiS_snp_v=(SiRiS.col(i));
//      SiRiS_snp=SiRiS_snp_v.array();
      SiRiSr+=(SiRiS_snp_v.array()*(r_new-r[i]));
    }
  }else{
    for(size_t i= p-1; i>=0; i--){
      r[i]=alpha0[i]*mu0[i];
      mu[i]=sigma_square[i]*(betahat[i]/se_square[i]+r[i]/se_square[i]-SiRiSr0[i]);
      //    mu = sigma_square * (betahat / se_square + r / se_square - SiRiSr_snp);
      // Update the variational estimate of the posterior inclusion probability.
      double SSR = mu0[i] * mu0[i] / sigma_square[i];
      
      
      alpha[i] = sigmoid(logodds + 0.5 * (log(sigma_square[i]/(sigma_beta_square)) + SSR));
      
      // Update SiRiSr = inv(S)*R*inv(S)*r
      double r_new = alpha[i] * mu[i];
      SiRiS_snp_v=(SiRiS.col(i));
      //      SiRiS_snp=SiRiS_snp_v.array();
      SiRiSr+=(SiRiS_snp_v.array()*(r_new-r[i]));
    }
    
    
    
    
    
  }
    
    
 
// //  Eigen::ArrayXd alpha0=alpha;
//   Eigen::ArrayXd r=alpha*mu;
//   
//   alpha = sigmoid(logodds+0.5*(log(sigma_square/(sigma_beta_square))+(mu*mu)/sigma_square));
//   
//   
//  // sigma_square * (betahat / se_square + r / se_square - SiRiSr_snp)
// 
//   // mu <- sigma_square*(betahat/se_square-(SiRiS%*%r-r/se_square))@x
//   mu=sigma_square * (betahat / se_square -((SiRiS*r.matrix()).array()-r/se_square));
//   SiRiSr=(SiRiS*(alpha*mu).matrix()).array();

}



void rss_varbvsr_iter(const Eigen::SparseMatrix<double> SiRiS,
                      const double sigma_beta,
                      const double logodds,
                      const Eigen::ArrayXd &betahat,
                      const Eigen::ArrayXd &se_square,
                      const Eigen::ArrayXd &sigma_square,
                      Eigen::ArrayXd &alpha,
                      Eigen::ArrayXd &mu,
                      Eigen::ArrayXd &SiRiSr,
                      bool reverse){
  
  
  double sigma_beta_square = sigma_beta*sigma_beta;
  // Get the number of SNPs (p) and coordinate ascent updates (m).
  const size_t p = betahat.size();
  
  // Initialize outputs.
  
  // Store a single column of matrix inv(S)*R*inv(S).
  
  Eigen::ArrayXd  SiRiS_snp(p);
  Eigen::VectorXd  SiRiS_snp_v(p);
  // Run coordinate ascent updates.
  // Repeat for each coordinate ascent update.
  size_t i=0;
  for (size_t j = 0; j < p; j++) {
    if(reverse){
      i=p-1-j;
    }else{
      i=j;
    }
    SiRiS_snp_v=(SiRiS.col(i));
    SiRiS_snp=SiRiS_snp_v.array();
    // Copy the kth column of matrix inv(S)*R*inv(S).
    // copySparseColumn(SiRiS_snp, k, SiRiS.elems, Ir, Jc, p);
    
    // Copy the kth element of vector inv(S)*R*inv(S)*r.
    double SiRiSr_snp = SiRiSr(i);
    
    // Perform the mean-field variational update.
    rss_varbvsr_update(betahat(i),
                       se_square(i),
                       sigma_beta_square, 
                       sigma_square(i),
                       SiRiS_snp,
                       SiRiSr,
                       SiRiSr_snp,
                       logodds,
                       alpha(i),
                       mu(i));
  }
}