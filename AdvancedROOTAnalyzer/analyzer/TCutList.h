#ifndef _TCutList_h_
#define _TCutList_h_

#include <map>
#include <string>

class TCutList
{
protected:
  // cuts
  std::map <std::string, bool> fCuts;

public:
  void Set(std::string name, bool value);
  bool PassesAll();
  bool PassesAllBut(std::string name);
};

#endif
