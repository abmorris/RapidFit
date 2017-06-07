# $Id: Makefile,v 1.30 2009/11/11 17:18:14 gcowan Exp $
SHELL=/bin/bash
UNAME=$(shell uname -s )
CC=g++ -fdiagnostics-color=always

# Location of compiled "common" libraries from ssh://git@gitlab.cern.ch:7999/admorris/common.git
COMMONDIR   = common
COMCXXFLAGS = $(shell make -sC $(COMMONDIR) cflags)
COMLIBS     = $(shell make -sC $(COMMONDIR) libs)
COMLIBDIR   = $(shell make -sC $(COMMONDIR) libdir)
COMLIBFLAGS = -L$(COMLIBDIR) -Wl,--as-needed $(COMLIBS) -Wl,-rpath,$(COMLIBDIR)

#		Include Root files as system headers as they're NOT standards complient and we do not want to waste time fixing them!
#		ROOT has some broken backwards compatability for OSX so won't claim to be a set of system headers
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTLIBDIR = $(shell root-config --libdir)
ROOTGLIBS  = $(shell root-config --glibs)

#               On some Systems with Mathmore compiled, sometimes things need to be resolved against it... I don't know why
EXTRA_ROOTLIBS=-lTreePlayer -lThread -lMinuit -lMinuit2 -lRooFit -lRooStats -lRooFitCore -lFoam -lMathMore

#		Command Line Tools
CXX          = $(CC) $(shell if [ "$(shell root-config --arch | grep 32)" = "" ]; then echo ""; else echo "--arch=i386"; fi)
RM           = rm -f

CXXFLAGS_BASE_MINIMAL = -D_GNU_SOURCE -D__USE_GNU -fPIC
CXXFLAGS_BASE_WARNINGS = -Werror -Wall -Wextra -Wno-non-virtual-dtor -Wno-reorder -Wshadow -Wmissing-noreturn -Wcast-align
#		Compiler Flags
CXXFLAGS_BASE_COMMON  = $(CXXFLAGS_BASE_MINIMAL) -D__ROOFIT_NOBANNER  $(CXXFLAGS_BASE_WARNINGS)
CXXFLAGS_BASE = $(CXXFLAGS_BASE_COMMON) -O3 -msse2 -msse3 -m3dnow -ftree-vectorize -finline-limit=2000 -fprefetch-loop-arrays -fmerge-all-constants

#		Some Useful global variables, makes this file MUCH easier to maintain
SRCEXT    = cpp
HDREXT    = h
SRCDIR    = framework/src
SRCPDFDIR = pdfs/src
INCDIR    = framework/include
INCPDFDIR = pdfs/include
INCGSL    = $(shell if command -v gsl-config >/dev/null 2>&1; then echo "$(shell gsl-config --cflags)"; else echo ""; fi )
LINKGSL   = $(shell if command -v gsl-config >/dev/null 2>&1; then echo "$(shell gsl-config --libs)"; else echo ""; fi )
USE_GSL   = $(shell if command -v gsl-config >/dev/null 2>&1; then echo "-D__RAPIDFIT_USE_GSL"; else echo ""; fi )
OBJDIR    = framework/build
OBJPDFDIR = pdfs/build
EXEDIR    = bin
LIBDIR    = lib
SRCDALITZEXT = cc
HDRDALITZEXT = hh
SRCDALITZDIR = pdfs/dalitz/src
INCDALITZDIR = pdfs/dalitz/include
OBJDALITZDIR = pdfs/dalitz/build

#	Source Files to be Built	ignoring all files in 'unused' and the RapidRun source for ROOT linking
SRCS    := $(shell find $(SRCDIR) -name '*.$(SRCEXT)' | grep -v 'unused' )
PDFSRCS := $(shell find $(SRCPDFDIR) -name '*.$(SRCEXT)' | grep -v 'unused' )
DALITZSRCS := $(shell find $(SRCDALITZDIR) -name '*.$(SRCDALITZEXT)' | grep -v 'unused' )

#	Absolute Paths of headers	ignoring the LinkDef written for ROOT and ignoring unused code
HEADERS := $(shell find $(INCDIR) -name '*.$(HDREXT)' | grep -v 'unused' | grep -v 'LinkDef' )
PDFHEAD := $(shell find $(INCPDFDIR) -name '*.$(HDREXT)' )
DALITZHEAD := $(shell find $(INCDALITZDIR) -name '*.$(HDRDALITZEXT)' )

#	Binary Objects to make in the build
OBJS    := $(patsubst $(SRCDIR)/%.$(SRCEXT),$(OBJDIR)/%.o,$(SRCS))
PDFOBJS := $(patsubst $(SRCPDFDIR)/%.$(SRCEXT),$(OBJPDFDIR)/%.o,$(PDFSRCS))
DALITZOBJS := $(patsubst $(SRCDALITZDIR)/%.$(SRCDALITZEXT),$(OBJDALITZDIR)/%.o,$(DALITZSRCS))

#################
##Dependencies

LINKFLAGS = -Wl,--no-undefined -Wl,--no-allow-shlib-undefined -lpthread

LIBS=-lstdc++

CXXFLAGS     = $(CXXFLAGS_BASE_COMMON) -I$(INCDIR) -I$(INCPDFDIR) -I$(INCDALITZDIR) -I$(INCGSL) $(ROOTCFLAGS) $(COMCXXFLAGS) -std=c++1y

#CHECKGCCABI := $(shell expr `gcc -dumpversion | cut -f1 -d.` \>= 5)
#ifeq ("$(CHECKGCCABI)","1")
#	CXXFLAGS+= -D_GLIBCXX_USE_CXX11_ABI=0
#	LINKFLAGS+= -D_GLIBCXX_USE_CXX11_ABI=0
#endif

LINKFLAGS += $(USE_GSL) $(LINKGSL) $(ROOTLIBS) $(EXTRA_ROOTLIBS) $(COMLIBFLAGS) -Wl,-rpath,$(COMLIBDIR):$(ROOTLIBDIR)

.PHONY : all clean extra lib $(COMMONDIR)

#	Default build command when someone asks for 'make'
all : $(EXEDIR)/fitting $(COMMONDIR)

$(COMMONDIR) :
	make -C $@

$(OBJDALITZDIR)/%.o : $(SRCDALITZDIR)/%.$(SRCDALITZEXT) $(INCDALITZDIR)/%.$(HDRDALITZEXT)
	@echo "Compiling $@"
	@$(CXX) $(CXXFLAGS) $(USE_GSL) $(INCGSL) -c $< -o $@

$(OBJPDFDIR)/%.o : $(SRCPDFDIR)/%.$(SRCEXT) $(INCPDFDIR)/%.$(HDREXT) $(DALITZOBJS)
	@echo "Compiling $@"
	@$(CXX) $(CXXFLAGS) $(USE_GSL) $(INCGSL) -c $< -o $@

$(OBJDIR)/%.o : $(SRCDIR)/%.$(SRCEXT) $(INCDIR)/%.$(HDREXT)
	@echo "Compiling $@"
	@$(CXX) $(CXXFLAGS) $(USE_GSL) $(INCGSL) -c $< -o $@

#	Main Build of RapidFit Binary
$(EXEDIR)/fitting : $(OBJS) $(PDFOBJS) $(DALITZOBJS) $(OBJDIR)/rapidfit_dict.o | $(COMMONDIR)
	@echo "Linking $@"
	@$(CXX) $(OBJDIR)/*.o $(OBJPDFDIR)/*.o $(OBJDALITZDIR)/*.o -o $@ $(LINKFLAGS)
	chmod +t $(EXEDIR)/fitting

#	Cleanup
clean :
	$(RM) $(EXEDIR)/* $(OBJDIR)/* $(OBJPDFDIR)/* $(OBJDALITZDIR)/* $(LIBDIR)/*
	make -C $(COMMONDIR) clean

extra:	$(EXEDIR)/Per-Event $(EXEDIR)/lifetime_tool $(EXEDIR)/weighted $(EXEDIR)/ApplyWeights $(EXEDIR)/Compare $(EXEDIR)/tupleDiff $(EXEDIR)/AngularDist $(EXEDIR)/plotDists

#	For Compiling RapidFit as a library to use within CINT which makes life easier on the grid... (supposedly)
#	make lib

lib: $(LIBDIR)/libRapidRun.so

#	This command will generate a C++ file which interfaces the rest of humanity with root...
#	It requires the explicit paths of all files, or that you remain in the same working directory at all times during the build process
#	We want to place the output dictionary in the Build directory as this is CODE that is NOT to be editted by the $USER!
$(OBJDIR)/rapidfit_dict.cpp: framework/include/RapidRun.h framework/include/LinkDef.h
	@echo "Generating $@"
	rootcint -f $(OBJDIR)/rapidfit_dict.cpp -c $^

#	Compile the class that root has generated for us which is the linker interface to root	(i.e. dictionaries & such)
$(OBJDIR)/rapidfit_dict.o: $(OBJDIR)/rapidfit_dict.cpp
	@echo "Compiling $@"
	@$(CXX) $(CXXFLAGS) -o $@ -I. -c $<

#	Class which has a dictionary generated for it, think of this as the equivalent to int main() in a CINT-y Universe
$(OBJDIR)/RapidRun.o: $(SRCDIR)/RapidRun.cpp
	@echo "Compiling $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

#	Finally, Compile RapidFit as a library making use of the existing binaries for other classes
$(LIBDIR)/libRapidRun.so: $(OBJDIR)/RapidRun.o $(OBJS) $(PDFOBJS) $(DALITZOBJS) $(OBJDIR)/rapidfit_dict.o
	@echo "Linking $@"
	@$(CXX) -shared $(OBJDIR)/*.o $(OBJPDFDIR)/*.o $(OBJDALITZDIR)/*.o -o $@ $(LINKFLAGS)

