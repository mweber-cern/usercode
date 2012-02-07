#include "TCutList.h"
#include "Utilities.h"

void TCutList::Set(std::string name, bool value)
{
  // check if already set
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
  bool found = false;
  for (p = fCuts.begin(); p != fCuts.end(); p++) {
    if (p->first != name)
      passes = passes && p->second;
    else 
      found = true;
  }
  if (!found)
    THROW("cut with name " + name + " not found!");
  return passes;
}
