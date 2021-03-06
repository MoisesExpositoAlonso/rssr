#ifndef INCLUDE_BVSR
#define INCLUDE_BVSR
#include <RcppEigen.h>
#include "rssr_types.h"
using namespace Eigen;





void rss_varbvsr_update (const double betahat, const double se, const double sigma_beta, const c_arrayxd_internal SiRiS_snp, 
                         arrayxd_internal SiRiSr, const double SiRiSr_snp, const double logodds, double &alpha, double &mu);
void rss_varbvsr_iter(const c_sparseMatrix_internal SiRiS,
                      const double sigma_beta,
                      const double logodds,
                      const c_arrayxd_internal betahat,
                      const c_arrayxd_internal se,
                      arrayxd_internal alpha,
                      arrayxd_internal mu,
                      arrayxd_internal SiRiSr,
                      bool reverse);


Rcpp::List rss_varbvsr_naive(const c_sparseMatrix_internal SiRiS,
                             const double sigma_beta,
                             const double logodds,
                             const c_arrayxd_internal betahat,
                             const c_arrayxd_internal se,
                             const c_arrayxd_internal talpha0,
                             const c_arrayxd_internal tmu0,
                             const c_arrayxd_internal tSiRiSr0,
                             double tolerance,
                             int itermax,
                             Rcpp::LogicalVector verbose,
                             Rcpp::LogicalVector lnz_tol);

#endif
