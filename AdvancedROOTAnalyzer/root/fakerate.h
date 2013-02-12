#ifndef __fakerate_h
#define __fakerate_h

void tightlooseplots(int start = 1, int end = 5);
void tight_loose_ratioplot();
TH1D * get_subtracted_tight_loose_ratio(bool save=true, bool draw=true);
TH1D * fake_estimate_1d(const char * sel, const char * hname);
TH2D * fake_estimate_2d(const char * sel, const char * hname);

#endif
