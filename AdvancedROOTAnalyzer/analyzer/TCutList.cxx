#include "TCutList.h"

void TCutList::Set(std::string name, bool value)
{
  fCuts[name] = value;
}

bool TCutList::PassesAll()
{
  std::map<std::string, bool>::const_iterator p;
  bool passes = true;
  for (p = fCuts.begin(); p != fCuts.end(); p++) {
    passes = passes && p->second;
  }
  return passes;
}

bool TCutList::PassesAllBut(std::string name)
{
  std::map<std::string, bool>::const_iterator p;
  bool passes = true;
  for (p = fCuts.begin(); p != fCuts.end(); p++) {
    if (p->first != name)
      passes = passes && p->second;
  }
  return passes;
}
