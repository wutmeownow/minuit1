# This Makefile builds two stand alone programs that have access to 
# the ROOT libraries.  To use this with your own code just substitute
# the name of your program below.

# here we access the root configuration, include files, and libraries
ROOTCFLAGS=$(shell root-config --cflags)
ROOTINC=$(shell root-config --incdir)
ROOTLIBDIR=$(shell root-config --libdir)
ROOTLIBS=$(shell root-config --libs) -lMinuit
ROOTLDFLAGS=$(shell root-config --ldflags)

ROOTC=$(ROOTCFLAGS) 
#-I$(ROOTINC)
ROOTLINK=-L$(ROOTLIBDIR) $(ROOTLIBS) $(ROOTLDFLAGS)

CPP=g++

default: expFit rootExample dist1Fit simulExpsFit templBkgSubtr

templBkgSubtr: templBkgSubtr.cpp
	$(CPP) -O -Wall $(ROOTC) -o templBkgSubtr templBkgSubtr.cpp $(ROOTLINK)

simulExpsFit: simulExpsFit.cpp
	$(CPP) -O -Wall $(ROOTC) -o simulExpsFit simulExpsFit.cpp $(ROOTLINK) 

dist1Fit: dist1Fit.cpp
	$(CPP) -O -Wall $(ROOTC) -o dist1Fit dist1Fit.cpp $(ROOTLINK) 

expFit: expFit.cpp
	$(CPP) -O -Wall $(ROOTC) -o expFit expFit.cpp $(ROOTLINK) 

rootExample: rootExample.cpp
	$(CPP) -O -Wall $(ROOTC) -o rootExample rootExample.cpp $(ROOTLINK) 
# note: just replace the -O flag with -g to build a debug version

clean: 
	rm -f expFit rootExample *~ *.d *.so
