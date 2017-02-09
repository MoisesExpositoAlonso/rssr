#ifndef KL_HPP
#define KL_HPP
#include <RcppEigen.h>
#include <limits>

typedef std::numeric_limits<double> double_lim;

Eigen::ArrayXd betavar(const Eigen::ArrayXd &p,const Eigen::ArrayXd &mu,Eigen::ArrayXd &s);

double intklbeta_rssbvsr(const Eigen::ArrayXd &alpha,const Eigen::ArrayXd &mu,const Eigen::ArrayXd &sigma_square,const double &sigma_beta_square);

double intgamma(const Eigen::ArrayXd &logodds, const Eigen::ArrayXd &alpha);

double find_maxerr(const Eigen::ArrayXd &alpha,
                   const Eigen::ArrayXd &alpha0,
                   const Eigen::ArrayXd &r,
                   const Eigen::ArrayXd &r0);

double calculate_lnZ(const Eigen::VectorXd &q,
                     const Eigen::VectorXd &r,
                     const Eigen::VectorXd &SiRiSr,
                     const Eigen::ArrayXd &logodds,
                     const Eigen::VectorXd &sesquare,
                     const Eigen::VectorXd &alpha,
                     const Eigen::VectorXd &mu,
                     const Eigen::VectorXd &s,
                     const double sigb);

#endif