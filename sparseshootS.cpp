#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma; 

// [[Rcpp::export]]
double rhobi(double z, double k){
  // rho function Tukey's biweight

  double absz=std::abs(z);
  double value=0;

  if (absz <= k) {
    value = (1-pow(  (1 - pow(z/k,2.0))   ,3.0))*pow(k,2.0)/6;
  } else {
    value=pow(k,2.0)/6;
  }
  return value;
}


// [[Rcpp::export]]
double scaleC(arma::vec u, double kp,double cc, double initialsc, int maxit, double tol){
  // u: n-dimensional vector of residuals
  // kp: constant for consistency of M-scale : E(z)=kp
  // cc: tuning parameter for rho function
  // initialsc: initial scale esitmate
  // maxit: maximum number of iterations
  // tol: numerical tolerance for convergence

  // Return: M-scale estimate of the residuals u
  double sc=initialsc;
  int n=u.size();
  double total=0;
  double sc2=1;
  double err=10*tol;
  int j=0;
  while((j<maxit) & (err>tol)){
    for( int i=0; i<n; ++i){
      total+= rhobi(u[i]/sc,cc)/n;
    }
    sc2=sqrt(pow(sc,2.0)*total/kp);
    err=std::abs(sc2/sc-1);
    sc=sc2;
    total=0;
    j+=1;
  }
  return sc;
}

// [[Rcpp::export]]
double wbi(double z, double k){
// weights Tukey's biweight
  double absz=std::abs(z);
  double value=0;

  if (absz <= k) {
    value = pow(1-pow(z/k,2.0),2.0);
  } else {
    value=0;
  }
  return value;
}

// [[Rcpp::export}]]
double sign(double x){
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}

// [[Rcpp::export]]
arma::mat betawls( arma::mat X,  arma::mat Y,  arma::mat W){
  // estimated regression coefficient for weighted least squares 
  
  // X: nxp matrix (in shooting loops, always nx2: one predictor and constant)
  // Y: nx1 matrix
  // W: nxn matrix
  arma::mat betahat = inv(X.t()*W*X)*X.t()*W*Y;

  return(betahat);
}

// [[Rcpp::export]]
List univest(arma::mat x, arma::mat y, double k, arma::vec wt, double delta, int maxiter, double tol,
             arma::vec prevresiduals, double tolscale, int maxitscale, double sinit){
  // Simple regression
  
  // x: nx1 matrix
  // y: nx1 matrix
  // k: tuning parameter of rho-function
  // wt: n-dimensional vector of regression weights
  // delta: constant of consistency of M-scale (E[Z]=delta)
  // maxiter: maximum number of iterations in IRWLS algorithm
  // tol: numerical tolerance for convergence
  // prevresiduals: n-dimensional vector of residuals
  // tolscale: numerical tolerance in fixed point iterations
  // maxitscale: maximum number of iterations in fixed point iterations
  // sinit: initial scale estimator

  // Initialization
  int n=x.n_rows;
  arma::mat mW=diagmat(wt);
  arma::mat bwls;
  arma::vec residuals;
  double err=10*tol;
  double s=1;
  int j=0;

  // Add column for the intercept
  arma::mat I;
  I.ones(n,1); 
  arma::mat xI = join_rows(I, x); 


  while((j<maxiter) & (err>tol)){

    // weighted least squares
    bwls=betawls(xI,y,mW);
    // residuals WLS
    residuals=y-xI*bwls;

    // compute residual scale
    if(j==0){
      s=scaleC(residuals,delta,k,sinit, maxitscale, tolscale);
    }else{
      s=scaleC(residuals,delta,k,s, maxitscale, tolscale);
    }

    // update weights
    for(int iw=0; iw<n;++iw){
	     wt[iw]=wbi(residuals[iw]/s,k);
    }
    mW= diagmat(wt);
    err=max(abs(prevresiduals-residuals));
    prevresiduals=residuals;
    j+=1;
  }

  List results=List::create(Named("s")=s,Named("coef")=bwls,Named("err")=err,Named("iter")=j,
                            Named("wt")=wt,Named("residuals")=residuals);
  return(results);
}


// [[Rcpp::export]]
List univest_sparse(arma::mat x, arma::mat y, double k, arma::vec wt, double delta, int maxiter, double tol,
             arma::vec prevresiduals,double tolscale, int maxitscale, double lambda, int post, double sinit){
  // Sparse simple regression
  
  // x: nx1 matrix
  // y: nx1 matrix
  // k: tuning parameter of rho-function
  // wt: n-dimensional vector of regression weights
  // delta: constant of consistency of M-scale (E[Z]=delta)
  // maxiter: maximum number of iterations in IRWLS algorithm
  // tol: numerical tolerance for convergence
  // prevresiduals: n-dimensional vector of residuals
  // tolscale: numerical tolerance in fixed point iterations
  // maxitscale: maximum number of iterations in fixed point iterations
  // lambda is penalty term
  // sinit: initial scale estimator

  // Initialization
  int n=x.n_rows;
  arma::mat mW=diagmat(wt);
  arma::mat bwls;
  arma::vec residuals;
  double err=10*tol;
  double s=1;
  int j=0;

  // Add column for the intercept
  arma::mat I;
  I.ones(n,1); 
  arma::mat xI = join_rows(I, x); 


  while((j<maxiter) & (err>tol)){

    // weighted least squares
    bwls= betawls(xI,y,mW);

	  // Sparsify 
	  if(bwls[1] > (lambda/2) || bwls[1]< -(lambda/2) ){
		if(post == 1){
		     bwls[1] = bwls[1];
		}else{
			if(bwls[1] > (lambda/2) ){
	    		bwls[1] = bwls[1]-(lambda/2) ;
			} else{
			bwls[1] = bwls[1]+ (lambda/2);
			}
		}

    }else{
	    bwls[1] = 0;
    }
    
    // residuals WLS
    residuals=y-xI*bwls;

    // compute residual scale
    if(j==0){
      s=scaleC(residuals,delta,k,sinit, maxitscale, tolscale);
    }else{
      s=scaleC(residuals,delta,k,s, maxitscale, tolscale);
    }

    // update weights
    for(int iw=0; iw<n;++iw){
	    wt[iw]=wbi(residuals[iw]/s,k);
    }
    mW= diagmat(wt);
    err=max(abs(prevresiduals-residuals));
    prevresiduals=residuals;
    j+=1;
  }

  List results=List::create(Named("s")=s,Named("coef")=bwls,Named("err")=err,Named("iter")=j,
                            Named("wt")=wt,Named("residuals")=residuals);
  return(results);
}

// [[Rcpp::export]]
List shootloop(arma::mat ytilde, arma::mat xtilde, arma::mat betaEst, arma::mat x, arma::mat wt,
               arma::vec scaleVar,double k, double delta, int maxIteration,double tol,
               arma::mat intercept, arma::mat Xexp, arma::vec scaleVarOLD,double tolout,
               double tolscale, int maxitscale, arma::mat y, int maxituniv , double wcut){
  // Shooting loop
  
  // ytilde: nx1 matrix
  // xtilde: nxp matrix
  // betaEst: px1 matrix
  // x: nxp matrix
  // wt: nxp matrix
  // scaleVar:p-dimensional vector
  // k: tuning parameter for rho-function
  // delta: constant of consistency of M-scale (E[Z]=delta)
  // maxIteration: maximum number of iterations shooting loop
  // tol: numerical tolerance for convergence IRWLS
  // intercept:px1 matrix
  // Xexp: nxp matrix xhat matrix
  // scaleVarOLD: p-dimensional vector
  // tolout: numerical tolerance for convergence shooting loop
  // tolscale: numerical tolerance for congervence fixed point iterations
  // maxitscale: maximum number of iterations fixed point iterations
  // y: nx1 matrix
  // maxituniv : maximum number of iterations in reweighted least squares
  // wcut: cut-off value for outlier flagging

  // Initializations
  int n=wt.n_rows; 
  arma::vec auxvec=ytilde;
  double medaux=1;
  arma::mat auxmat=ytilde;
  arma::vec resid=ytilde;
  arma::vec wjt=ytilde;
  arma::vec betaEstold = betaEst;
  List outputuniv;
  arma::vec residuniv=ytilde;
  arma::vec absresid=ytilde;
  double shat=1;
  int p=xtilde.n_cols+1;
  double err=10*tolout;
  vec lastresidscale = zeros(maxIteration,1);
  int flagstuck = 0;
  
  int j=0;
 
  while((j<maxIteration) & (err>tolout)& (flagstuck == 0)){

    for(int indexcol=1; indexcol<p; ++indexcol){
      
      // Regression step
      ytilde = y - xtilde*betaEst + xtilde.col(indexcol-1)*betaEst.row(indexcol-1);
 
      // Residuals
      auxvec=ytilde - x.col(indexcol-1)*betaEst.row(indexcol-1);
      medaux=median(auxvec);
      auxmat=ytilde;
      auxmat.fill(medaux);
      resid=ytilde - x.col(indexcol-1)*betaEst.row(indexcol-1) - auxmat;

      // Regression weights
      for(int iw=0; iw<n;++iw){
	      wt(iw,indexcol-1)=wbi(resid[iw]/scaleVar[indexcol-1],k);
      }
      wjt=wt.col(indexcol-1);



      // Simple regression
      outputuniv=univest(x.col(indexcol-1), ytilde, k, wjt, delta, maxituniv, tol, resid,tolscale, maxitscale, scaleVar(indexcol-1));

      // Assign output of simple regression
      betaEst(indexcol-1,0)=as<arma::mat>(outputuniv["coef"])(1,0);
      intercept(indexcol-1,0)=as<arma::mat>(outputuniv["coef"])(0,0);
      scaleVar(indexcol-1)=outputuniv["s"];
      wt.col(indexcol-1)=as<arma::vec>(outputuniv["wt"]);
      residuniv=as<arma::vec>(outputuniv["residuals"]);
      absresid=abs(residuniv);
      shat=outputuniv["s"];

      // Regression weights
      for(int il=0; il<n; ++il){
        if(absresid(il)/shat<=wcut){
          wt(il,indexcol-1)=1;
          xtilde(il,indexcol-1)=x(il,indexcol-1);
        }else{
          wt(il,indexcol-1)=0;
          xtilde(il,indexcol-1)=Xexp(il,indexcol-1);
        }
      }

      // Update Response
      ytilde=ytilde - xtilde.col(indexcol-1)*betaEst(indexcol-1);
    }
	
    arma::vec error_vec = abs(betaEst-betaEstold);	
    err=*std::max_element(error_vec.begin(), error_vec.end());
    lastresidscale(j) = err;
    for(int l = 0;l < j; ++l){
      double value = lastresidscale(l);
      if((value == lastresidscale(j))){
         flagstuck = 1;
         }
     } 
    scaleVarOLD=scaleVar;
    betaEstold = betaEst;
    j+=1;
  }

  arma::vec ytildevec = y-xtilde*betaEst;
  double alpha=median(ytildevec);

  List results=List::create(Named("resid")=resid,Named("betaEst")=betaEst,Named("alpha")=alpha,Named("scaleVar")=scaleVar,Named("weights")=wt,
                                  Named("iter")=j, Named("xtilde")= xtilde, Named("flagstuck")= flagstuck, Named("lastresidscale") = lastresidscale(j-1));
  return(results);

}


// [[Rcpp::export]]
List shootloop_sparse(arma::mat ytilde, arma::mat xtilde, arma::mat betaEst, arma::mat x, arma::mat wt,
               arma::vec scaleVar,double k, double delta, int maxIteration,double tol,
               arma::mat intercept, arma::mat Xexp, arma::vec scaleVarOLD,double tolout,
               double tolscale, int maxitscale, arma::mat y, int maxituniv, double lambda, int post, double wcut){
  // Sparse Shooting loop
  
  // ytilde: nx1 matrix
  // xtilde: nxp matrix
  // betaEst: px1 matrix
  // x: nxp matrix
  // wt: nxp matrix
  // scaleVar:p-dimensional vector
  // k: tuning parameter for rho-function
  // delta: constant of consistency of M-scale (E[Z]=delta)
  // maxIteration: maximum number of iterations shooting loop
  // tol: numerical tolerance for convergence IRWLS
  // intercept:px1 matrix
  // Xexp: nxp matrix of xhat
  // scaleVarOLD: p-dimensional vector
  // tolout: numerical tolerance for convergence shooting loop
  // tolscale: numerical tolerance for congervence fixed point iterations
  // maxitscale: maximum number of iterations fixed point iterations
  // y: nx1 matrix
  // maxituniv : maximum number of iterations in reweighted least squares
  // lambda: penalty term
  // post: 1 if post-lasso approach should be used, 0 else
  // wcut: cut-off value for outlier flagging

  // Initializations
  int n=wt.n_rows; 
  arma::vec auxvec=ytilde;
  double medaux=1;
  arma::mat auxmat=ytilde;
  arma::vec resid=ytilde;
  arma::vec wjt=ytilde;
  arma::vec betaEstold = betaEst;
  List outputuniv;
  arma::vec residuniv=ytilde;
  arma::vec absresid=ytilde;
  double shat=1;
  int p=xtilde.n_cols+1;
  double err=10*tolout;
  vec lastresidscale = zeros(maxIteration,1);
  int flagstuck = 0;
 
  int j=0;
  
  while((j<maxIteration) & (err>tolout)& (flagstuck == 0)){
      for(int indexcol=1; indexcol<p; ++indexcol){

      // Regression step
      ytilde = y - xtilde*betaEst + xtilde.col(indexcol-1)*betaEst.row(indexcol-1);

      // Residuals
      auxvec=ytilde - x.col(indexcol-1)*betaEst.row(indexcol-1);
      medaux=median(auxvec);
      auxmat=ytilde;
      auxmat.fill(medaux);
      resid=ytilde - x.col(indexcol-1)*betaEst.row(indexcol-1) - auxmat;

      for(int iw=0; iw<n;++iw){
        wt(iw,indexcol-1)=wbi(resid[iw]/scaleVar[indexcol-1],k);
      }
      wjt=wt.col(indexcol-1);



      // Simple regression
      outputuniv=univest_sparse(x.col(indexcol-1), ytilde, k, wjt, delta, maxituniv, tol, resid,tolscale, maxitscale, lambda, post, scaleVar(indexcol-1));

      // Assign output of simple regression
      betaEst(indexcol-1,0)=as<arma::mat>(outputuniv["coef"])(1,0);
      intercept(indexcol-1,0)=as<arma::mat>(outputuniv["coef"])(0,0);
      scaleVar(indexcol-1)=outputuniv["s"];
      wt.col(indexcol-1)=as<arma::vec>(outputuniv["wt"]);
      residuniv=as<arma::vec>(outputuniv["residuals"]);
      absresid=abs(residuniv);
      shat=outputuniv["s"];
	
      // Regression weights
      for(int il=0; il<n; ++il){
        if(absresid(il)/shat<=wcut){
          wt(il,indexcol-1)=1;
          xtilde(il,indexcol-1)=x(il,indexcol-1);
        }else{
          wt(il,indexcol-1)=0;
          xtilde(il,indexcol-1)=Xexp(il,indexcol-1);
        }
      }
      
      // Update Response
      ytilde=ytilde - xtilde.col(indexcol-1)*betaEst(indexcol-1);
    }
     
    arma::vec error_vec = abs(betaEst-betaEstold);	
    err=*std::max_element(error_vec.begin(), error_vec.end());
    lastresidscale(j) = err;
    for(int l = 0;l < j; ++l){
      double value = lastresidscale(l);
      if((value == lastresidscale(j))){
         flagstuck = 1;
         }
     } 
    scaleVarOLD=scaleVar;
    betaEstold = betaEst;
    j+=1;
  }

  arma::vec ytildevec = y-xtilde*betaEst;
  double alpha=median(ytildevec);
  arma::vec residall = ytildevec-alpha;
  double sguess = median(abs(residall))/.6745;
  double sigmahat = scaleC(residall, delta, k,sguess,maxitscale,tolscale);

  List results=List::create(Named("resid")=resid,Named("betaEst")=betaEst,Named("alpha")=alpha,Named("scaleVar")=scaleVar,Named("weights")=wt,
                                  Named("iter")=j, Named("sigmahat")=sigmahat, Named("flagstuck")= flagstuck, Named("lastresidscale") = lastresidscale(j-1));
  return(results);

}


// [[Rcpp::export]]
List shootloop_sparse_lambdas(arma::mat ytilde, arma::mat xtilde, arma::mat betaEst, arma::mat x, arma::mat wt,
                                arma::vec scaleVar,double k, double delta, int maxIteration,double tol,
                                  arma::mat intercept, arma::mat Xexp, arma::vec scaleVarOLD,double tolout,
                                    double tolscale, int maxitscale, arma::mat y, int maxituniv, arma::vec lambda, int post, double wcut){
  // Select sparsity parameter
  // ytilde: nx1 matrix
  // xtilde: nxp matrix
  // betaEst: px1 matrix
  // x: nxp matrix
  // wt: nxp matrix
  // scaleVar:p-dimensional vector
  // k: tuning parameter for rho-function
  // delta: constant of consistency of M-scale (E[Z]=delta)
  // maxIteration: maximum number of iterations shooting loop
  // tol: numerical tolerance for convergence IRWLS
  // intercept:px1 matrix
  // Xexp: nxp matrix of xhat
  // scaleVarOLD: p-dimensional vector
  // tolout: numerical tolerance for convergence shooting loop
  // tolscale: numerical tolerance for congervence fixed point iterations
  // maxitscale: maximum number of iterations fixed point iterations
  // y: nx1 matrix
  // maxituniv : maximum number of iterations in reweighted least squares
  // lambda: vector of sparsity parameters
  // post: 1 if post-lasso approach should be used, 0 else
  // wcut: cut-off value for outlier flagging
  
  int llength = lambda.size();
  List fit;
  Rcpp::List resid(llength);
  Rcpp::List beta(llength);
  Rcpp::List alpha(llength);
  Rcpp::List scales(llength);
  Rcpp::List weights(llength);
  Rcpp::List iter(llength);
  Rcpp::List sigmahat(llength);
  Rcpp::List flagstuck(llength);
  Rcpp::List lastresidscale(llength);

  for(int il=0; il < llength; ++ il){ 
    fit = shootloop_sparse(ytilde, xtilde, betaEst,x, wt, scaleVar, k,  delta,  maxIteration, tol,
                            intercept, Xexp, scaleVarOLD, tolout, tolscale,  maxitscale,  y,  maxituniv,  lambda[il],  post,  wcut);
    resid[il] = fit["resid"];
    beta[il] = fit["betaEst"];
    alpha[il] = fit["alpha"];
    scales[il] = fit["scaleVar"];
    weights[il] = fit["weights"];
    iter[il] = fit["iter"];
    sigmahat[il] = fit["sigmahat"];
    flagstuck[il] = fit["flagstuck"];
    lastresidscale[il] = fit["lastresidscale"];
  }
  
  List results=List::create(Named("resid") = resid, Named("betaEst") = beta, Named("alpha") = alpha, Named("scaleVar") = scales, Named("weights") = weights,
                                  Named("iter") = iter, Named("sigmahat") = sigmahat, Named("flagstuck") = flagstuck, Named("lastresidscale") = lastresidscale);
  return(results);
  
}
