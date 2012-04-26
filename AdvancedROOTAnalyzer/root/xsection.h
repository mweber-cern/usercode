#ifndef _xsection_h
#define _xsection_h

void xsratio(Double_t xs[4], // xs, parabolic error, -error, +error
	     Int_t order = 1, Bool_t draw = kFALSE);

void correctfactor(const char * pname, Double_t xlow, Double_t xup, bool setWeight = true);

void BinnedPoissonianLikelihood(Int_t    & npar,  // # of internal parameters
				Double_t * grad,  // array of first derivatives
				Double_t & fval,  // the function value
				Double_t * xvar,  // the input variable values
				Int_t      iflag);// status flag

#endif

