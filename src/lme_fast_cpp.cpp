/*
BrainStat-style mass-univariate random-intercept fitting, RcppArmadillo.
Independent specialized implementation of the algorithm in:
https://github.com/MICA-MNI/BrainStat/blob/master/brainstat/stats/_linear_model.py
and coefficient standard errors in stats/_t_test.py (source retrieved 2026-09-12).

Model: Y[,j] = X beta[,j] + subject random intercept + independent noise.
Covariance bases are [K, I], K[i,k]=1 if observations i,k share a subject.
This is not a full BrainStat port: no other covariance bases, random slopes,
multivariate outcomes, surface processing, multiple-testing correction, or
contrast-specific degrees of freedom. Full-column-rank X is required.

R USAGE (install Rcpp and RcppArmadillo first):
  Rcpp::sourceCpp("brainstat_cpp.cpp")
  fit <- fast_rint_reg_multi_y_cpp(Y, X, group_offsets, group_sizes)
Y is observations x outcomes; X is observations x fixed predictors. Group rows
must be contiguous and group_offsets ZERO-based, as in lme_fast_cpp.cpp.
For unsorted data:
  stopifnot(!anyNA(dat$PTID))
  ord <- order(dat$PTID)
  Y2 <- as.matrix(Y[ord,,drop=FALSE])
  dat2 <- dat[ord,,drop=FALSE]
  X2 <- model.matrix(~ Age + sex, dat2)
  stopifnot(nrow(X2)==nrow(Y2))
  sizes <- rle(as.character(dat2$PTID))$lengths
  offsets <- c(0, head(cumsum(sizes),-1))
  fit <- fast_rint_reg_multi_y_cpp(Y2,X2,offsets,sizes)

The exported name, argument types/defaults, and seven output fields match
lme_fast_cpp.cpp. Replace that source when using this in a package; do not
compile two definitions of the exported function. Regenerate Rcpp attributes.

CONTROLS:
 gamma_max: upper bound on the random/residual ratio (compatibility extension).
 tol, max_iter, grid_size: validated but otherwise unused legacy REML controls.
 BrainStat controls are internal constants in fit(): niter=1 (extra update),
 thetalim=0.01, drlim=0.1. Bins are adaptive, NOT a fixed number of bins.

METHOD:
 For residual precision R = V^-1 - V^-1 X (X'V^-1 X)^-1 X'V^-1,
 E[a,j]=y_j' R V_a R y_j; M[a,b]=trace(R V_a R V_b).
 theta=pinv(M) E, then apply BrainStat's variance lower limits and normalize.
 An initial estimate uses V=I. The reference's initial dr expression is kept,
 including its zero-initialized slm_r term. One further Fisher-scoring update
 uses mean variance proportions in rounded r/dr bins, then final binned GLS.
 Ties in bin rounding use nearest-even, matching NumPy.

OUTPUT:
 coefficients, std_error, t_stat: P x M numeric matrices.
 var_random, var_residual, gamma: length M numeric vectors. These describe
 the covariance USED for final GLS, not the unbinned scoring theta estimates:
 total = final whitened residual SSE/(N-P), var_random=total*r_bin,
 var_residual=total*(1-r_bin), gamma=var_random/var_residual.
 std_error^2=total*diag((X'Vshape^-1 X)^-1); t_stat=coefficient/std_error.
 status: integer M-vector, 0 fitted, 1 effective gamma ceiling reached,
 2 failed/degenerate. Legacy status 3 (optimizer limit) is never returned.
 No p-values are returned. BrainStat's contrast-specific dfs are not returned;
 do not assume N-P reproduces BrainStat's p-values.

Differences from Python: specialized block covariance algebra; double precision;
 outcome/design rescaling; ill-conditioned designs rejected instead of pinv;
 degenerate outcomes excluded before global bin-width estimation; finite gamma
 ceiling; explicit NA failures. Removing outcomes can change global bin widths.
 No full N x N x number-of-variance-components array is formed. The scoring
 information calculation still uses O(N^2) memory; it is not a huge-N solver.

Define BRAINSTAT_STANDALONE only to compile the core outside R with Armadillo.
*/
#ifdef BRAINSTAT_STANDALONE
#include <armadillo>
#else
#include <RcppArmadillo.h>
#endif
#include <cmath>
#include <limits>
#include <map>
#include <vector>
#include <stdexcept>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
namespace brainstat_detail {
using arma::uword;
inline void interrupt(){
#ifndef BRAINSTAT_STANDALONE
 Rcpp::checkUserInterrupt();
#endif
}
struct Result {
 arma::mat coefficients,std_error,t_stat;
 arma::vec var_random,var_residual,gamma;
 arma::ivec status;
 // Internal diagnostics, not exported by the compatible R interface.
 arma::vec rho,rho_used; double dr;
};
// Vshape=(1-r)I+rK. Centering avoids subtracting nearly equal quantities.
arma::mat weight(const arma::mat& A,double r,const arma::uvec& off,const arma::uvec& sz){
 arma::mat B=A;
 for(uword g=0;g<sz.n_elem;++g){
  uword lo=off[g],hi=lo+sz[g]-1;
  arma::rowvec mean=arma::mean(A.rows(lo,hi),0);
  B.rows(lo,hi).each_row()-=mean;
  B.rows(lo,hi)/=(1-r);
  B.rows(lo,hi).each_row()+=mean/(1-r+r*double(sz[g]));
 }
 return B;
}
struct Context {arma::mat cov,map,mi;};
Context context(const arma::mat& X,double r,const arma::uvec& off,const arma::uvec& sz,bool scoring){
 Context c;
 arma::mat W=weight(X,r,off,sz),A=X.t()*W;
 A=0.5*(A+A.t());
 if(!arma::inv_sympd(c.cov,A))throw std::runtime_error("GLS design is numerically singular.");
 c.map=c.cov*W.t();
 if(!scoring)return c;
 const uword N=X.n_rows,G=sz.n_elem;
 arma::mat R=weight(arma::eye<arma::mat>(N,N),r,off,sz)-W*c.map;
 R=0.5*(R+R.t());
 arma::mat RZ(N,G);
 for(uword g=0;g<G;++g)RZ.col(g)=arma::sum(R.cols(off[g],off[g]+sz[g]-1),1);
 double m00=0;
 for(uword g=0;g<G;++g){
  arma::rowvec row=arma::sum(RZ.rows(off[g],off[g]+sz[g]-1),0);
  m00+=arma::dot(row,row);
 }
 arma::mat M(2,2);
 M(0,0)=m00; M(0,1)=M(1,0)=arma::accu(arma::square(RZ));
 M(1,1)=arma::accu(arma::square(R));
 // Match NumPy pinv's default relative singular-value threshold.
 if(!arma::pinv(c.mi,M,1e-15*arma::norm(M,2)))
  throw std::runtime_error("Variance information pseudoinverse failed.");
 return c;
}
arma::rowvec score(const arma::mat& Y,const arma::mat& X,const Context& c,double r,
                  const arma::uvec& off,const arma::uvec& sz,double ceiling){
 arma::mat e=Y-X*(c.map*Y),u=weight(e,r,off,sz);
 arma::mat E(2,Y.n_cols,arma::fill::zeros);
 E.row(1)=arma::sum(arma::square(u),0);
 for(uword g=0;g<sz.n_elem;++g){
  arma::rowvec sums=arma::sum(u.rows(off[g],off[g]+sz[g]-1),0);
  E.row(0)+=arma::square(sums);
 }
 arma::mat theta=c.mi*E;
 arma::rowvec ans(Y.n_cols);
 for(uword j=0;j<Y.n_cols;++j){
  double total=theta(0,j)+theta(1,j);
  if(!std::isfinite(total)||total<=0){ans[j]=arma::datum::nan;continue;}
  for(uword a=0;a<2;++a){
   double floor=std::sqrt(std::max(0.0,2*c.mi(a,a)))*0.01*total;
   theta(a,j)=std::max(theta(a,j),floor);
  }
  ans[j]=std::min(ceiling,theta(0,j)/arma::accu(theta.col(j)));
 }
 return ans;
}
// Explicit nearest-even, independent of the process floating-point round mode.
double round_even(double x){double f=std::floor(x),d=x-f;
 if(d<0.5)return f;if(d>0.5)return f+1;return std::fmod(f,2)==0?f:f+1;}
Result fit(const arma::mat& Y,const arma::mat& X,const arma::uvec& off,
           const arma::uvec& sz,double gamma_max,double tol,int max_iter,int grid_size){
 uword N=Y.n_rows,M=Y.n_cols,P=X.n_cols;
 if(!M||!P||X.n_rows!=N||N<=P)throw std::invalid_argument("Invalid dimensions; require N>P.");
 if(!X.is_finite()||!Y.is_finite())throw std::invalid_argument("X and Y must be finite.");
 if(sz.n_elem<2||off.n_elem!=sz.n_elem)throw std::invalid_argument("Need at least two groups and matching offsets/sizes.");
 if(!std::isfinite(gamma_max)||gamma_max<=0||!std::isfinite(tol)||tol<=0||max_iter<1||grid_size<5)
  throw std::invalid_argument("Invalid optimization controls.");
 uword expected=0;bool repeated=false;
 for(uword g=0;g<sz.n_elem;++g){
  if(!sz[g]||off[g]!=expected||sz[g]>N-expected)throw std::invalid_argument("Invalid contiguous grouping; offsets are zero-based.");
  expected+=sz[g]; repeated=repeated||sz[g]>1;
 }
 if(expected!=N||!repeated)throw std::invalid_argument("Groups must cover rows and include repeated observations.");
 double ceiling=gamma_max<=1?gamma_max/(1+gamma_max):1-1/(1+gamma_max);
 ceiling=std::min(ceiling,std::nextafter(1.0,0.0));
 arma::rowvec xs=arma::max(arma::abs(X),0);
 if(arma::any(xs==0))throw std::invalid_argument("X has an all-zero column.");
 arma::mat D=X;D.each_row()/=xs;
 arma::mat Q,R;
 if(!arma::qr_econ(Q,R,D)||arma::rcond(R)<1e-12)throw std::invalid_argument("X is rank deficient or ill-conditioned.");
 Q.reset();R.reset();
 const double nan=arma::datum::nan,df=double(N-P);
 const uword block=256;
 const int extra_updates=1; // BrainStat SLM default niter
 Result z;
 z.coefficients.set_size(P,M);z.coefficients.fill(nan);
 z.std_error=z.coefficients;z.t_stat=z.coefficients;
 z.rho.set_size(M);z.rho.fill(nan);z.rho_used=z.rho;
 z.var_random=z.rho;z.var_residual=z.rho;z.gamma=z.rho;
 z.status=arma::ivec(M,arma::fill::value(2));z.dr=nan;
 arma::vec scale(M,arma::fill::ones);
 auto scaled=[&](const std::vector<uword>& ids,uword start,uword count){
  arma::mat yy(N,count);
  for(uword k=0;k<count;++k)yy.col(k)=Y.col(ids[start+k])/scale[ids[start+k]];
  return yy;
 };
 Context init=context(D,0,off,sz,true);
 std::vector<uword> active;
 for(uword lo=0;lo<M;lo+=block){
  interrupt();uword count=std::min(block,M-lo);
  std::vector<uword> ids(count);
  for(uword k=0;k<count;++k){ids[k]=lo+k;double s=arma::abs(Y.col(lo+k)).max();scale[lo+k]=s>0?s:1;}
  arma::mat yy=scaled(ids,0,count),e=yy-D*(init.map*yy);
  arma::rowvec rr=score(yy,D,init,0,off,sz,ceiling);
  for(uword k=0;k<count;++k){
   double energy=arma::dot(yy.col(k),yy.col(k));
   if(energy==0||arma::dot(e.col(k),e.col(k))<=1e-24*energy||!std::isfinite(rr[k]))continue;
   z.rho[lo+k]=rr[k];active.push_back(lo+k);
  }
 }
 if(active.empty())return z;
 double mean_r2=0;for(uword j:active)mean_r2+=z.rho[j]*z.rho[j]/double(active.size());
 arma::mat Vt=2*init.mi;
 // Reference slm_r is zero-initialized, so its middle term vanishes.
 z.dr=std::sqrt(Vt(0,0)+arma::accu(Vt)*mean_r2)*0.1;
 if(!std::isfinite(z.dr)||z.dr<=0)throw std::runtime_error("Variance ratio is numerically unidentified.");
 init=Context();
 for(int pass=0;pass<=extra_updates;++pass){
  std::map<double,std::vector<uword>> bins;
  for(uword j:active)if(std::isfinite(z.rho[j]))bins[round_even(z.rho[j]/z.dr)].push_back(j);
  for(const auto& bin:bins){
   interrupt();const auto& ids=bin.second;
   double r=0;for(uword j:ids)r+=z.rho[j]/double(ids.size());r=std::min(r,ceiling);
   Context c;
   try{c=context(D,r,off,sz,pass<extra_updates);}catch(const std::runtime_error&){
    for(uword j:ids)z.rho[j]=nan;continue;
   }
   for(uword start=0;start<ids.size();start+=block){
    interrupt();uword count=std::min(block,uword(ids.size())-start);
    arma::mat yy=scaled(ids,start,count);
    if(pass<extra_updates){
     arma::rowvec rr=score(yy,D,c,r,off,sz,ceiling);
     for(uword k=0;k<count;++k)z.rho[ids[start+k]]=rr[k];
    }else{
     arma::mat b=c.map*yy,e=yy-D*b,we=weight(e,r,off,sz);
     for(uword k=0;k<count;++k){
      uword j=ids[start+k];double s=scale[j],q=arma::dot(e.col(k),we.col(k));
      if(!std::isfinite(q)||q<=1e-24*arma::dot(yy.col(k),yy.col(k)))continue;
      double total=q/df;
      arma::vec se=arma::sqrt(total*c.cov.diag()),coef=b.col(k)/xs.t()*s,serr=se/xs.t()*s;
      if(!coef.is_finite()||!serr.is_finite()||arma::any(serr<=0)||!std::isfinite(total*s*s)||total*s*s<=0)continue;
      z.coefficients.col(j)=coef;z.std_error.col(j)=serr;z.t_stat.col(j)=b.col(k)/se;
      z.rho_used[j]=r;z.gamma[j]=std::min(gamma_max,r/(1-r));
      z.var_residual[j]=total*s*s*(1-r);z.var_random[j]=z.gamma[j]*z.var_residual[j];
      z.status[j]=(r>=ceiling*(1-1e-12)?1:0);
     }
    }
   }
  }
 }
 return z;
}
} // namespace brainstat_detail
#ifndef BRAINSTAT_STANDALONE
// [[Rcpp::export]]
Rcpp::List fast_rint_reg_multi_y_cpp(
    const arma::mat& Y,
    const arma::mat& X,
    const arma::uvec& group_offsets,
    const arma::uvec& group_sizes,
    double gamma_max = 1e8,
    double tol = 1e-8,
    int max_iter = 200,
    int grid_size = 41) {
 auto z=brainstat_detail::fit(Y,X,group_offsets,group_sizes,gamma_max,tol,max_iter,grid_size);
 for(arma::uword j=0;j<Y.n_cols;++j)if(z.status[j]==2){
  z.coefficients.col(j).fill(NA_REAL);z.std_error.col(j).fill(NA_REAL);z.t_stat.col(j).fill(NA_REAL);
  z.var_random[j]=z.var_residual[j]=z.gamma[j]=NA_REAL;
 }
 return Rcpp::List::create(
  Rcpp::Named("coefficients")=z.coefficients,
  Rcpp::Named("std_error")=z.std_error,
  Rcpp::Named("t_stat")=z.t_stat,
  Rcpp::Named("var_random")=Rcpp::NumericVector(z.var_random.begin(),z.var_random.end()),
  Rcpp::Named("var_residual")=Rcpp::NumericVector(z.var_residual.begin(),z.var_residual.end()),
  Rcpp::Named("gamma")=Rcpp::NumericVector(z.gamma.begin(),z.gamma.end()),
  Rcpp::Named("status")=Rcpp::IntegerVector(z.status.begin(),z.status.end())
 );
}
#endif
