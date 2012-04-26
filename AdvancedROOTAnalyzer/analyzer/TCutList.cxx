#include "TCutList.h"
#include "Utilities.h"

#include "TH1D.h"

TCutList::TCutList(std::map <std::string, TH1D *> & histo, double weight)
  : fHisto(histo), fWeight(weight)
{
  // nothing to do
}

TCutList::~TCutList()
{
  // free head from created cuts
  std::map<std::string, TCutB *>::const_iterator p;
  for (p = fCuts.begin(); p != fCuts.end(); p++) {  
    delete p->second;
  }
}


void TCutList::Set(std::string name, bool cut, double value)
{
  if (fCuts[name] == 0) {
    fCuts[name] = new TCutB(cut, value);
  }
  else {
    fCuts[name]->first = cut;
    fCuts[name]->second = value;
  }
}

bool TCutList::PassesAll()
{
  std::map<std::string, TCutB *>::const_iterator p;
  bool passes = true;
  for (p = fCuts.begin(); p != fCuts.end(); p++) {
    passes = passes && p->second->first;
  }
  return passes;
}

bool TCutList::PassesAllBut(std::string name)
{
  std::map<std::string, TCutB *>::const_iterator p;
  bool passes = true;
  bool found = false;
  for (p = fCuts.begin(); p != fCuts.end(); p++) {
    if (p->first != name)
      passes = passes && p->second->first;
    else 
      found = true;
  }
  if (!found)
    THROW("cut with name " + name + " not found!");
  return passes;
}

void TCutList::FillHistograms()
{
  std::map<std::string, TCutB *>::const_iterator p;
  for (p = fCuts.begin(); p != fCuts.end(); p++) {
    if (!PassesAllBut(p->first)) {
      continue;
    }
    TH1D * h = fHisto[p->first];
    if (h != 0)
      h->Fill(p->second->second, fWeight);
    else {
      THROW(std::string("Histogram \"") + p->first + std::string("\" not existing. Did you misspell or forgot to create?"));
    }
  }
}
