#

G77 = gfortran
CXX = g++
CC = gcc
LD = g++ 

FASTJETHOME	  = /wsu/home/groups/rhi/build/fastjet-install
FASTJETLIBS   = -Wl,-rpath,$(FASTJETHOME)/lib -lm \
				-L$(FASTJETHOME)/lib -lfastjettools -lfastjet -lfastjetplugins -lsiscone_spherical -lsiscone
THERMUSHOME   = /wsu/home/groups/rhi/build/thermus2020
THERMUSLIBS   = -L${THERMUSHOME}/build/lib -lFunctions -lTHERMUS
GSLLIBS      := $(shell gsl-config --libs) 

INCLUDES      = -I. -I$(shell root-config --incdir) 
ARCHFLAG      = -m64
ROOTCFLAGS   := $(shell root-config --cflags) 
ROOTLIBS     := $(shell root-config --libs --evelibs) -lMinuit
CERNLIBS      = -L$(CERN)/pro/lib -lkernlib -lpacklib -lmathlib
LIBS          = $(GSLLIBS) $(ROOTLIBS) -lGui -lHtml -lGuiHtml -lHistPainter -lMathCore -lMathMore 
# ROOTGLIBS    := $(shell root-config --new --glibs)
# ROOTGLIBS    += -lNew -lGui
# GLIBS         = $(ROOTGLIBS) $(SYSLIBS) 

#---- enable profiling? (if yes, use GPROF in following lines; if no, use NOGPROF)
GPROF = -pg -O3
NOGPROF = -ggdb
LDFLAGS = ${NOGPROF}
CFLAGS = $(ARCHFLAG) $(ROOTCFLAGS)
CXXFLAGS = $(CFLAGS) $(NOGPROF) -Wuninitialized

#---- enable fastjet?
# CXXFLAGS += -DUSEFASTJET -Wno-deprecated-declarations
# INCLUDES += -I$(FASTJETHOME)/include
# LIBS     += $(FASTJETLIBS)
#---- enable thermus?
# CXXFLAGS += -DUSETHERMUS
# INCLUDES += -I$(THERMUSHOME)/main
# LIBS     += $(THERMUSLIBS)

all: m

#----------------------------------------------------------------------------------

m: m.o MultCum.o
	${CXX} m.o MultCum.o \
		$(LDFLAGS) \
		$(ARCHFLAG) \
		$(INCLUDES) \
		$(LIBS) \
		-o m

clean:
	rm -fv *.o

m.o: m.cxx
	${CXX} -c $(CXXFLAGS) $< $(INCLUDES)

MultCum.o: MultCum.cxx
	${CXX} -c $(CXXFLAGS) $< $(INCLUDES)

m.o: MultCum.h

