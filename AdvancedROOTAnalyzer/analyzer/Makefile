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
CXXFLAGS      = -O2 -Wall -fPIC -fsignaling-nans -g -lz # -DNDEBUG # -pg for gprof
CXXFLAGS     += $(EXTINCDIR)   # external header files should be here
CXXFLAGS     += $(ROOTCFLAGS)

# -DVERSION=88

LD            = g++
LDFLAGS       = -O2 -fsignaling-nans -g # -pg for gprof
LDFLAGS      += $(ROOTGLIBS)
LDFLAGS      += $(EXTLIBS)

SOURCES_COMMON	= $(shell ls *.cxx | grep -v _v[0-9][0-9].cxx )
SOURCES_78	= $(shell ls *_v78.cxx )
SOURCES_88	= $(shell ls *_v88.cxx )
SOURCES_92	= $(shell ls *_v92.cxx )
SOURCES_96	= $(shell ls *_v96.cxx )
HEADERS_COMMON  = $(SOURCES_COMMON:.cxx=.h)
HEADERS_78     	= $(SOURCES_78:.cxx=.h)
HEADERS_88     	= $(SOURCES_88:.cxx=.h)
HEADERS_92     	= $(SOURCES_92:.cxx=.h)
HEADERS_96     	= $(SOURCES_96:.cxx=.h)
OBJECTS_78	= $(SOURCES_78:.cxx=.o) $(SOURCES_COMMON:.cxx=_v78.o)
OBJECTS_88	= $(SOURCES_88:.cxx=.o) $(SOURCES_COMMON:.cxx=_v88.o)
OBJECTS_92	= $(SOURCES_92:.cxx=.o) $(SOURCES_COMMON:.cxx=_v92.o)
OBJECTS_96	= $(SOURCES_96:.cxx=.o) $(SOURCES_COMMON:.cxx=_v96.o)
DEPENDS_78	= $(SOURCES_78:.cxx=.d) $(SOURCES_COMMON:.cxx=_v78.d)
DEPENDS_88	= $(SOURCES_88:.cxx=.d) $(SOURCES_COMMON:.cxx=_v88.d)
DEPENDS_92	= $(SOURCES_92:.cxx=.d) $(SOURCES_COMMON:.cxx=_v92.d)
DEPENDS_96	= $(SOURCES_96:.cxx=.d) $(SOURCES_COMMON:.cxx=_v96.d)

DICT_78	      = Dict_v78
DICT_88	      = Dict_v88
DICT_92	      = Dict_v92
DICT_96	      = Dict_v96

CUTFLOW	      = CutFlow

DATE	     := $(shell date "+%F")

######################################################################
# default targets

.PHONY:		all clean distclean tarball

all:            $(PROGRAM)_v78 $(PROGRAM)_v88 $(PROGRAM)_v92 $(PROGRAM)_v96

clean:
		rm -f $(OBJECTS_COMMON) $(OBJECTS_78) $(OBJECTS_88) $(OBJECTS_92) $(OBJECTS_96) \
		 $(PROGRAM)_v78 $(PROGRAM)_v88 $(PROGRAM)_v92 $(PROGRAM)_v96 \
		 $(DEPENDS_78) $(DEPENDS_88) $(DEPENDS_92) $(DEPENDS_96) \
		 $(DICT_78).* $(DICT_88).* $(DICT_92).* $(DICT_96).* $(CUTFLOW).*

distclean:      clean
		@rm -f *.o *.d *~ core *.def *.exp *.root *.ps .def gmon.out

tarball:
		tar cvjf $(PROGRAM)-$(DATE).tar.bz2 Makefile *.h *.cxx

######################################################################
# rules

.SUFFIXES: .cxx .$(ExeSrcSuf) .C .f .o .so .d

%_v78.o: %.cxx 
	$(CXX) $(CXXFLAGS) -DVERSION=78 -c $< -o $@

%_v78.o: %_v78.cxx 
	$(CXX) $(CXXFLAGS) -DVERSION=78 -c $< -o $@

%_v78.o: %_v78.C
	$(CXX) $(CXXFLAGS) -DVERSION=78 -c $<

%_v88.o: %.cxx 
	$(CXX) $(CXXFLAGS) -DVERSION=88 -c $< -o $@

%_v88.o: %_v88.cxx 
	$(CXX) $(CXXFLAGS) -DVERSION=88 -c $< -o $@

%_v88.o: %_v88.C
	$(CXX) $(CXXFLAGS) -DVERSION=88 -c $<

%_v92.o: %.cxx 
	$(CXX) $(CXXFLAGS) -DVERSION=92 -c $< -o $@

%_v92.o: %_v92.cxx 
	$(CXX) $(CXXFLAGS) -DVERSION=92 -c $< -o $@

%_v92.o: %_v92.C
	$(CXX) $(CXXFLAGS) -DVERSION=92 -c $<

%_v96.o: %.cxx 
	$(CXX) $(CXXFLAGS) -DVERSION=96 -c $< -o $@

%_v96.o: %_v96.cxx 
	$(CXX) $(CXXFLAGS) -DVERSION=96 -c $< -o $@

%_v96.o: %_v96.C
	$(CXX) $(CXXFLAGS) -DVERSION=96 -c $<


%.o: %.f
	$(FC) $(FCFLAGS) -c $< -o $@

%_v78.d: %.cxx
	$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) -DVERSION=78 $< > $@'

%_v78.d: %_v78.cxx
	$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) -DVERSION=78 $< > $@'

%_v88.d: %_v88.cxx
	$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) -DVERSION=88 $< > $@'

%_v88.d: %.cxx
	$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) -DVERSION=88 $< > $@'

%_v92.d: %_v92.cxx
	$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) -DVERSION=92 $< > $@'

%_v92.d: %.cxx
	$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) -DVERSION=92 $< > $@'

%_v96.d: %_v96.cxx
	$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) -DVERSION=96 $< > $@'

%_v96.d: %.cxx
	$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) -DVERSION=96 $< > $@'

include $(DEPENDS_78) $(DEPENDS_88) $(DEPENDS_92) $(DEPENDS_96)

######################################################################
# special procedures

######################################################################
# ROOT dictionary

$(DICT_78).C: $(HEADERS_COMMON) $(HEADERS_78) LinkDef.h
	@echo "Generating dictionary $(DICT)..."
	rootcint -f $@ -c $(EXTINCDIR) -DVERSION=78 $+

$(DICT_88).C: $(HEADERS_COMMON) $(HEADERS_88) LinkDef.h
	@echo "Generating dictionary $(DICT)..."
	rootcint -f $@ -c $(EXTINCDIR) -DVERSION=88 $+

$(DICT_92).C: $(HEADERS_COMMON) $(HEADERS_92) LinkDef.h
	@echo "Generating dictionary $(DICT)..."
	rootcint -f $@ -c $(EXTINCDIR) -DVERSION=92 $+

$(DICT_96).C: $(HEADERS_COMMON) $(HEADERS_96) LinkDef.h
	@echo "Generating dictionary $(DICT)..."
	rootcint -f $@ -c $(EXTINCDIR) -DVERSION=96 $+

$(CUTFLOW).C: Analysis.cxx
	@echo "Generating cutflow..."
	../bin/gencutflow.py -f $@ $+

######################################################################
# targets

$(PROGRAM)_v78:  $(OBJECTS_COMMON) $(OBJECTS_78) $(DICT_78).o $(CUTFLOW).o
	$(LD) $(LDFLAGS) $+ -o $@
		@echo "$@ done"

$(PROGRAM)_v88:  $(OBJECTS_COMMON) $(OBJECTS_88) $(DICT_88).o $(CUTFLOW).o
	$(LD) $(LDFLAGS) $+ -o $@
		@echo "$@ done"

$(PROGRAM)_v92:  $(OBJECTS_COMMON) $(OBJECTS_92) $(DICT_92).o $(CUTFLOW).o
	$(LD) $(LDFLAGS) $+ -o $@
		@echo "$@ done"

$(PROGRAM)_v96:  $(OBJECTS_COMMON) $(OBJECTS_96) $(DICT_96).o $(CUTFLOW).o
	$(LD) $(LDFLAGS) $+ -o $@
		@echo "$@ done"
