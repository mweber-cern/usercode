#include "Utilities.h"
#include "EventFilter.h"

#include <stdlib.h>
#include <zlib.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

EventFilter::EventFilter(string eventFileName, bool verbose)
{
  minrun_  = -1;
  maxrun_  = -1;
  minRunInFile=999999;
  maxRunInFile=1;
  readEventListFile(eventFileName);
  std::sort(EventList_.begin(), EventList_.end());

  if (minrun_ == -1 || minrun_ < minRunInFile)
    minrun_ = minRunInFile;
  if (maxrun_ == -1 || maxrun_ > maxRunInFile)
    maxrun_ = maxRunInFile;
  verbose_ = verbose;
}

EventFilter::~EventFilter()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

#define LENGTH 0x2000
#define CHUNK 16384
void EventFilter::readEventListFile(string eventFileName)
{
#ifdef NO_GZCOMPRESS
  ERROR("NO_GZCOMPRESS -- gz* functions cannot compress");
#else
  
  // code does nothing when string is empty
  if ( eventFileName.empty() ) {
    return;
  }


  gzFile file = gzopen(eventFileName.c_str(), "r");
  if (!file) {
    WARNING("Unable to open event list file " << eventFileName);
    return;
  }
  string b2;
  int err;
  int bytes_read;
  char buffer[LENGTH];
  unsigned int i;
  char * pch;

  while (1) {
    bytes_read = gzread(file, buffer, LENGTH - 1);
    buffer[bytes_read] = '\0';
    i = 0;
    pch = strtok(buffer,"\n");
    if (buffer[0] == '\n' ) {
      addEventString(b2);
      ++i;
    } 
    else 
      b2 += pch;

    while (pch != NULL)
    {
      i += strlen(pch)+1;
      if (i > b2.size()) 
	b2 = pch;
      if (i == (LENGTH-1)) {
        if ((buffer[LENGTH-2] == '\n' )|| (buffer[LENGTH-2] == '\0' )){
          addEventString(b2);
          b2="";
        }
      } 
      else if (i < LENGTH) {
        addEventString(b2);
      }
      pch = strtok (NULL, "\n");
    }
    if (bytes_read < LENGTH - 1) {
      if (gzeof(file)) 
	break;
      else {
	const char * error_string;
	error_string = gzerror(file, & err);
	if (err) {
	  ERROR("Error while reading gzipped file:  " << error_string);
	  return;
	}
      }
    }
  }
  gzclose (file);
#endif
  return;
}

void EventFilter::addEventString(const string & eventString)
{
  // Loop through list of bad events, and if run is in allowed range, add bad event to EventList
  int run = 0;
  unsigned int ls = 0;
  unsigned int event = 0;
  // Check that event list object is in correct form
  size_t found = eventString.find(":");  // find first colon
  if (found != std::string::npos) {
    run = atoi((eventString.substr(0,found)).c_str());  // convert to run
  }
  else
  {
    WARNING("Unable to parse event list input '" << eventString <<"' for run number!");
    return;
  }
  size_t found2 = eventString.find(":", found+1);  // find second colon
  if (found2 != std::string::npos) {
    // Some event numbers are less than 0?  \JetHT\Run2012C-v1\RAW:201278:2145:-2130281065
    // -- due to events being dumped out as ints, not uints!
    ls = atoi((eventString.substr(found+1,(found2-found-1))).c_str());  // convert to ls
    event = atoi((eventString.substr(found2+1)).c_str()); // convert to event
    /// Some event numbers are less than 0?  \JetHT\Run2012C-v1\RAW:201278:2145:-2130281065
    if (ls==0 || event==0 || event<0)
      WARNING("Strange lumi, event numbers for input '" << eventString << "'");
  }
  else {
    WARNING("Unable to parse Event list input '" << eventString << "' for run number!");
    return;
  }
  // If necessary, check that run is within allowed range
  if (minrun_ > -1 && run < minrun_) {
    DEBUG("Skipping Event list input '" << eventString << "' because it is less than minimum run # "<< minrun_);
    return;
  }
  if (maxrun_> -1 && run > maxrun_) {
    DEBUG("Skipping Event list input '" << eventString << "' because it is greater than maximum run # "<< maxrun_);
    return;
  }
  if (minRunInFile > run) 
    minRunInFile = run;
  if (maxRunInFile < run) 
    maxRunInFile = run;
  // Now add event to Event List
  EventList_.push_back(eventString);
}


bool EventFilter::filter(int run, unsigned int lumisection, unsigned int event)
{
  // if run is outside filter range, then always return true
  if (minrun_ > -1 && run < minrun_) 
    return true;
  if (maxrun_ > -1 && run > maxrun_) 
    return true;

  // Okay, now create a string object for this run:ls:event
  std::stringstream thisevent;
  thisevent << run << ":" << lumisection << ":" << event;

  strVecI it = std::lower_bound(EventList_.begin(), EventList_.end(), thisevent.str());
  // Event not found in bad list; it is a good event
  if (it == EventList_.end() || thisevent.str() < *it) 
    return true;

  // Otherwise, this is a bad event
  DEBUG("Filtering bad event " << thisevent.str());

  return false;
}
