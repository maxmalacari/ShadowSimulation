#
# General Makefile for ADST analysis examples
#
USER_SRCS = $(wildcard *.cc)
#
# Executable names come from the .cc sources in this directory.
# Replace the wildcard expression with .cc file list if you do
# not want to compile all .cc files in this directory
#
EXE = $(patsubst %.cc,%, $(wildcard *.cc))
#
#############################################################

## You should not need to change anything below this line ###

.PHONY: all depend clean

ifdef AUGEROFFLINEROOT
ADSTROOT = $(AUGEROFFLINEROOT)
endif

CXXFLAGS    = $(shell root-config --cflags)
CXXFLAGS   += -I$(ADSTROOT)/include/adst
LDFLAGS     = $(shell root-config --ldflags)
LDFLAGS    += $(shell root-config --libs)
all: $(EXE)

%: %.cc
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) -lMinuit 

#############################################################
# gcc can generate the dependency list

depend: Make-depend

Make-depend: $(USER_SRCS)
	$(CPP) $(CXXFLAGS) -MM $^ > $@

clean:
	- rm -f *.o  *.so *.ps core Make-depend $(EXE)

-include Make-depend
