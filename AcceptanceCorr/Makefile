SHELL = /bin/bash

.DELETE_ON_ERROR:

.PHONY: all clean


ROOTCONFIG  := root-config
ROOTCFLAGS  := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS    := $(shell $(ROOTCONFIG) --libs)
ROOTINCDIR  := $(shell $(ROOTCONFIG) --incdir)

CXX       := g++
CXXFLAGS  += -std=c++11 -O2 -Wall -fPIC $(ROOTCFLAGS)
LD        = g++
LDFLAGS   = -O2 $(ROOTLDFLAGS) -Wl,--no-as-needed

INCLUDES  := -I/$(ROOTINCDIR) -I$(MYROOTLIB) -I$(CLASTOOL)/include -I$(ANALYSER)/include
LIBS      := $(ROOTLIBS) -L$(CLASTOOL)/slib/Linux -lClasTool -L$(ANALYSER)/slib -lTIdentificator -L$(MYROOTLIB) -lmyROOTLib

##############################################################################
all: AcceptanceCorr
	
AcceptanceCorr:AcceptanceCorr.o
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

AcceptanceCorr.o:AcceptanceCorr.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	@rm -rf AcceptanceCorr AcceptanceCorr.o
