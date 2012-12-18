#ifndef __fakerate_h
#define __fakerate_h

void tightlooseplots(int start = 1, int end = 5);
TH1D * fake_estimate_1d(const char * sel, const char * hname);
TH2D * fake_estimate_2d(const char * sel, const char * hname);

#endif
