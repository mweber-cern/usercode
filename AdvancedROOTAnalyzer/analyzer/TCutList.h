#ifndef _TCutList_h_
#define _TCutList_h_

#include <utility>
#include <map>
#include <string>

class TH1D;

typedef std::pair<bool, double> TCutB;

class TCutList
{
protected:
  // cuts
  std::map <std::string, TCutB *> fCuts;
  std::map <std::string, TH1D *> & fHisto;
  double fWeight;

public:
  TCutList(std::map <std::string, TH1D *> & histo, double weight); 
  ~TCutList();
  void Set(std::string name, bool cut, double value);
  bool PassesAll();
  bool PassesAllBut(std::string name);
  void FillHistograms();
};

#endif
