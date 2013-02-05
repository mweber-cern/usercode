######################################################################
# quite generic Makefile for ROOT-based programs with automatic
# dependency scanning
#
# (C) Martin Weber 2011

# Name of the executable
PROGRAM        = analyzer

# external includes / libraries
#CLHEP_VERSION = 2.0.2.1
#EXTINCDIR     = -I/usr/local/include/ # external header files should be here, leading slash is mandatory!
#EXTLIBS       = -L/usr/local/lib -lCLHEP-Matrix-$(CLHEP_VERSION) -lCLHEP-Vector-$(CLHEP_VERSION)

######################################################################
# You should not modify the lines below here
#

# Setup for ROOT libraries
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTGLIBS   += -lX3d

# Linux with gcc / g++ / g77
FC            = g77
FCFLAGS       = -I/cern/pro/include -fPIC -fno-automatic -fdollar-ok -fno-backslash \
		-finit-local-zero -fno-second-underscore -fugly-logint -ftypeless-boz
# -pg for gprof

CXX           = g++
CXXFLAGS      = -O2 -Wall -fPIC -fsignaling-nans -g # -DNDEBUG # -pg for gprof
CXXFLAGS     += $(EXTINCDIR)   # external header files should be here
CXXFLAGS     += $(ROOTCFLAGS)

LD            = g++
LDFLAGS       = -O2 -fsignaling-nans -g -lz # -pg for gprof
LDFLAGS      += $(ROOTGLIBS)
LDFLAGS      += $(EXTLIBS)

SOURCES       = $(wildcard *.cxx)
HEADERS       = $(SOURCES:.cxx=.h)
OBJECTS       = $(SOURCES:.cxx=.o)
DEPENDS       = $(SOURCES:.cxx=.d)

DICT	      = Dict

CUTFLOW	      = CutFlow

DATE	     := $(shell date "+%F")

######################################################################
# default targets

all:            $(PROGRAM)

clean:
		@rm -f $(OBJECTS) $(PROGRAM) $(DICT).* $(CUTFLOW).* $(DEPENDS)

distclean:      clean
		@rm -f *~ core *.def *.exp *.root *.ps .def gmon.out

tarball:
		tar cvjf $(PROGRAM)-$(DATE).tar.bz2 Makefile *.h *.cxx

######################################################################
# rules

.SUFFIXES: .cxx .$(ExeSrcSuf) .C .f .o .so .d

.cxx.o:
	$(CXX) $(CXXFLAGS) -c $<

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

.cc.o:
	$(CXX) $(CXXFLAGS) -c $<

.C.o:
	$(CXX) $(CXXFLAGS) -c $<

.f.o:
	$(FC) $(FCFLAGS) -c $< -o $@

.cpp.d:
	$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< > $@'

.cxx.d:
	$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< > $@'

.cc.d:
	$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< > $@'

include $(DEPENDS)

######################################################################
# ROOT dictionary

$(DICT).C: $(HEADERS) LinkDef.h
	@echo "Generating dictionary $(DICT)..."
	rootcint -f $@ -c $(EXTINCDIR) $+

$(CUTFLOW).C: $(SOURCES)
	@echo "Generating cutflow..."
	gencutflow.py -f $@ $+

######################################################################
# targets

$(PROGRAM):  $(OBJECTS) $(DICT).o $(CUTFLOW).o
	$(LD) $(LDFLAGS) $+ -o $@
		@echo "$@ done"
