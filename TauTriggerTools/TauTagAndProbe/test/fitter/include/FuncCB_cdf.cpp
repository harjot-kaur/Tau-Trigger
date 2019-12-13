1;95;0c /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * The skeleton of this code is taken from FuncCB class                      * 
  *****************************************************************************/ 

 // Your description goes here... 

// #include "Riostream.h" 

#include "FuncCB_cdf.h"

 
ClassImp(FuncCB_cdf) 
  
  FuncCB_cdf::FuncCB_cdf(const char *name, const char *title, 
		 RooAbsReal& _m,
		 RooAbsReal& _m0,
		 RooAbsReal& _sigma,
		 RooAbsReal& _alpha,
		 RooAbsReal& _n,
		 RooAbsReal& _norm,
		 RooAbsReal& _yrise) :
    RooAbsReal(name,title), 
    m("m","m",this,_m),
    m0("m0","m0",this,_m0),
    sigma("sigma","sigma",this,_sigma),
    alpha("alpha","alpha",this,_alpha),
    n("n","n",this,_n),
    norm("norm","norm",this,_norm),
    yrise("yrise","yrise",this,_yrise)
{ 
} 


FuncCB_cdf::FuncCB_cdf(const FuncCB_cdf& other, const char* name) :  
  RooAbsReal(other,name), 
  m("m",this,other.m),
  m0("m0",this,other.m0),
  sigma("sigma",this,other.sigma),
  alpha("alpha",this,other.alpha),
  n("n",this,other.n),
  norm("norm",this,other.norm),
  yrise("yrise",this,other.yrise)
{ 
} 


Double_t FuncCB_cdf::valeur(Double_t et)
{
  m = et;
  return evaluate();
}

Double_t FuncCB_cdf::evaluate() const
{ 
	if (n <= 1.){
	 	MATH_ERROR_MSG("crystalball_cdf","CrystalBall cdf not defined for n <=1");
  return std::numeric_limits<double>::quiet_NaN();
	}
   	   static const double kSqrt2 = 1.41421356237309515; // sqrt(2.)
       double abs_alpha = std::abs(alpha);
       double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);
       double D = std::sqrt(M_PI/2.)*(1.+ RooMath::erf(abs_alpha/std::sqrt(2.)));
       double totIntegral = sigma*(C+D);
 
       double integral = crystalball_integral(-m, alpha, n, sigma, m0);
       return (alpha > 0) ? yrise -(1. - integral/totIntegral)*(norm) : yrise - (integral/totIntegral)*(norm);
  }


//_____________________________________________________________________________
Double_t FuncCB_cdf::crystalball_integral(Double_t m, Double_t alpha, Double_t n, Double_t sigma, Double_t m0) const
    {
       // compute the integral of the crystal ball function (ROOT::Math::crystalball_function)
       // If alpha > 0 the integral is the right tail integral.
       // If alpha < 0 is the left tail integrals which are always finite for finite x.     
       // parameters:
       // alpha : is non equal to zero, define the # of sigma from which it becomes a power-law function (from mean-alpha*sigma)
       // n > 1 : is integrer, is the power of the low  tail
       // add a value xmin for cases when n <=1 the integral diverges 
       if (sigma == 0)   return 0;
       if (alpha==0)
       {
          MATH_ERROR_MSG("crystalball_integral","CrystalBall function not defined at alpha=0");
          return 0.;
       }
       bool useLog = (n == 1.0); 
       if (n<=0)   MATH_WARN_MSG("crystalball_integral","No physical meaning when n<=0");
 
       double z = (m-m0)/sigma;
       if (alpha < 0 ) z = -z;
       
       double abs_alpha = std::abs(alpha);
       
       //double D = *(1.+ROOT::Math::erf(abs_alpha/std::sqrt(2.)));
       //double N = 1./(sigma*(C+D));
       double intgaus = 0.;
       double intpow  = 0.;
  
       const double sqrtpiover2 = std::sqrt(M_PI/2.);
       const double sqrt2pi = std::sqrt( 2.*M_PI); 
       const double oneoversqrt2 = 1./sqrt(2.);
       if (z <= -abs_alpha)
       {
          double A = std::pow(n/abs_alpha,n) * std::exp(-0.5 * alpha*alpha);
          double B = n/abs_alpha - abs_alpha;
 
          if (!useLog) {
             double C = (n/abs_alpha) * (1./(n-1)) * std::exp(-alpha*alpha/2.);
             intpow  = C - A /(n-1.) * std::pow(B-z,-n+1) ;
          }
          else {
             // for n=1 the primitive of 1/x is log(x)
             intpow = -A * std::log( n / abs_alpha ) + A * std::log( B -z );
          }
          intgaus =  sqrtpiover2*(1.+RooMath::erf(abs_alpha*oneoversqrt2));
       }
       else
       {
          intgaus = gaussian_cdf_c(z, 1, 0);
          intgaus *= sqrt2pi;
          intpow  =  0;  
       }
       return sigma * (intgaus + intpow);
    }
 

//_____________________________________________________________________________
Double_t FuncCB_cdf::gaussian_cdf_c(double m, double sigma, double m0) const 
    {
	   static const double kSqrt2 = 1.41421356237309515; // sqrt(2.)    
       double z = (m-m0)/(sigma*kSqrt2);
       if (z > 1.)  return 0.5*RooMath::erfc(z);
       else         return 0.5*(1.-RooMath::erf(z));
    }
 

//_____________________________________________________________________________
Double_t FuncCB_cdf::ApproxErf(Double_t arg) const 
{
  static const double erflim = 5.0;
  if( arg > erflim )
    return 1.0;
  if( arg < -erflim )
    return -1.0;
  
  return RooMath::erf(arg);
}
