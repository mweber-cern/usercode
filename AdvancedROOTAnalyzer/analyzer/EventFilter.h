#ifndef EventFilter_h
#define EventFilter_h

// system include files
#include <string>
#include <vector>

class EventFilter {
  
public :
  explicit EventFilter(std::string eventFileName, bool verbose = false);
  ~EventFilter();

  bool filter(int run, unsigned int lumisection, unsigned int event);

private:
  void readEventListFile(std::string eventFileName);
  void addEventString(const std::string & eventString);
        
  typedef std::vector< std::string >::iterator strVecI;

  // vector of strings representing bad events, with each string in "run:LS:event" format
  std::vector< std::string > EventList_;

  bool verbose_;

  // Set run range of events in the BAD LASER LIST.  
  // The purpose of these values is to shorten the length of the EventList_
  // vector when running on only a subset of data
  int minrun_;
  int maxrun_;

  // if specified (i.e., values > -1), then only events in the given range will be filtered
  int minRunInFile, maxRunInFile;
};
#endif 
