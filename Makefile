#------------------------------------------------------------------------------
# Core library
SRC  =  particle.cxx planeHit.cxx plane.cxx target.cxx sft.cxx tof1.cxx ac.cxx mwpc.cxx tof2.cxx matrix.cxx fieldMap.cxx fieldStepper.cxx trackState.cxx trackSite.cxx trackSystem.cxx track.cxx trackFinder.cxx

EXTRAHDR = utility.h particle.h planeHit.h target.h sft.h tof1.h ac.h mwpc.h tof2.h matrix.h fieldMap.h fieldStepper.h trackState.h trackSite.h trackSystem.h track.h

CORE = trek
CORELIB  = lib$(CORE).so
COREDICT = $(CORE)Dict

LINKDEF = $(CORE)_LinkDef.h


#------------------------------------------------------------------------------
# Compile extra code for printing verbose messages (enabled with fDebug)
export VERBOSE = 1

# Architecture to compile for
ARCH          = linux

#------------------------------------------------------------------------------
ROOTCFLAGS   := $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS     := $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS    := $(shell $(ROOTSYS)/bin/root-config --glibs)
ROOTBIN      := $(shell $(ROOTSYS)/bin/root-config --bindir)
CXX          := $(shell $(ROOTSYS)/bin/root-config --cxx)
LD           := $(shell $(ROOTSYS)/bin/root-config --ld)

PKGINCLUDES  = -I$(shell pwd)
INCLUDES     = -I$(shell root-config --incdir) $(PKGINCLUDES)

LIBS          = 
GLIBS         = 

ifeq ($(ARCH),linux)
# Linux with gcc (RedHat)
ifdef DEBUG
  CXXFLAGS    = -g -O0
  LDFLAGS     = -g -O0
  DEFINES     =
else
  CXXFLAGS    = -O2 -g #-march=pentium4
  LDFLAGS     = -O -g
#  DEFINES     = -DNDEBUG
endif
DEFINES      += -DLINUXVERS -DHAS_SSTREAM
CXXFLAGS     += -Wall -Woverloaded-virtual -fPIC
DICTCXXFLG   :=
ifdef EXTRAWARN
#FIXME: should be configure'd:
CXXVER       := $(shell g++ --version | head -1 | sed 's/.* \([0-9]\)\..*/\1/')
ifeq ($(CXXVER),4)
CXXFLAGS     += -Wextra -Wno-missing-field-initializers
DICTCXXFLG   := -Wno-strict-aliasing 
endif
endif
SOFLAGS       = -shared
ifdef I387MATH
CXXFLAGS     += -mfpmath=387
else
CXXFLAGS     += -march=core2 -mfpmath=sse
endif
endif


ifdef VERBOSE
DEFINES      += -DVERBOSE
endif


CXXFLAGS     += $(DEFINES) $(ROOTCFLAGS) $(PKGINCLUDES)
LIBS         += $(ROOTLIBS) $(SYSLIBS)
GLIBS        += $(ROOTGLIBS) $(SYSLIBS)


#------------------------------------------------------------------------------
OBJSRC		= $(SRC:.cxx=.o)
OBJ           = $(SRC:.cxx=.o) $(COREDICT).o
HDR           = $(SRC:.cxx=.h) $(EXTRAHDR)
DEP           = $(SRC:.cxx=.d)

EXEC = trekTracker

all:		$(EXEC).o $(CORELIB) 
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(EXEC).o -L$(shell pwd) -l$(CORE) $(LIBS)

$(EXEC).o:	$(HDR) ${EXEC}.cxx 
	$(CXX) $(CXXFLAGS) -o $(EXEC).o -c ${EXEC}.cxx

$(CORELIB):	$(OBJ) 
		$(LD) $(LDFLAGS) $(SOFLAGS) -o $@ $^ $(LIBS)
		@echo "$@ done"


#dbconvert:	dbconvert.o
#		$(LD) $(LDFLAGS) $(LIBS) -o $@ $^

ifeq ($(ARCH),linux)
$(COREDICT).o:	$(COREDICT).cxx
	$(CXX) $(CXXFLAGS) $(DICTCXXFLG) -o $@ -c $^
endif

$(COREDICT).cxx: $(HDR) $(LINKDEF)
	@echo "Generating dictionary $(COREDICT)..."
	$(ROOTBIN)/rootcint -f $@ -c $(INCLUDES) $(DEFINES) $^

$(OBJSRC):	 %.o: %.cxx $(EXTRAHDR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
		rm -f *.o *~ $(CORELIB) $(COREDICT).*
