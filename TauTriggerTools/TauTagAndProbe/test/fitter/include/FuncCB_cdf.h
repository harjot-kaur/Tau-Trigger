#ifndef FUNCCB_CDF
#define FUNCCB_CDF

#include "TMath.h" 
 #include "Math/Error.h"
#include <math.h> 

#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooMath.h"

class FuncCB_cdf : public RooAbsReal {
public:
  FuncCB_cdf() {} ; 
  FuncCB_cdf(const char *name, const char *title,
	      RooAbsReal& _m,
	      RooAbsReal& _m0,
	      RooAbsReal& _sigma,
	      RooAbsReal& _alpha,
	      RooAbsReal& _n,
	      RooAbsReal& _norm,
	      RooAbsReal& _yrise);
  FuncCB_cdf(const FuncCB_cdf& other, const char* name=0) ;
 
  Double_t evaluate() const;
  Double_t valeur(Double_t et);

  virtual TObject* clone(const char* newname) const { return new FuncCB_cdf(*this,newname); }
  inline virtual ~FuncCB_cdf() { }

protected:
  Double_t ApproxErf(Double_t arg) const ;
  Double_t crystalball_integral(Double_t m, Double_t alpha, Double_t n, Double_t sigma, Double_t m0) const;
  Double_t gaussian_cdf_c(double m, double sigma, double m0) const;
  RooRealProxy m ;
  RooRealProxy m0 ;
  RooRealProxy sigma ;
  RooRealProxy alpha ;
  RooRealProxy n ;
  RooRealProxy norm ; 
  RooRealProxy yrise;  

private:

  ClassDef(FuncCB_cdf,1) // Your description goes here...
};
 
#endif
