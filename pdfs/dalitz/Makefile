# Compiler and shell
CC         = g++
SHELL      = /bin/bash
RM         = rm -f

# ROOT
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTLIBDIR = $(shell root-config --libdir)
#EXTRA_ROOTLIBS = -lRooFit -lRooStats -lRooFitCore

# Extensions
SRCEXT     = cc
HDREXT     = hh
OBJEXT     = o
LIBEXT     = so

# Directories
TSRCDIR    = test/src
BINDIR     = test/bin
TOBJDIR    = test/build
DSRCDIR    = src
OBJDIR     = build
HDRDIR     = include
LIBDIR     = test/lib

# Get files and make list of objects and libraries
TSRCS     := $(shell find $(TSRCDIR) -name '*.$(SRCEXT)')
DSRCS     := $(shell find $(DSRCDIR) -name '*.$(SRCEXT)')
HDRS      := $(shell find $(HDRDIR) -name '*.$(HDREXT)')
OBJS      := $(patsubst $(DSRCDIR)/%.$(SRCEXT),$(OBJDIR)/%.$(OBJEXT),$(DSRCS))
BINS      := $(patsubst $(TSRCDIR)/%.$(SRCEXT),$(BINDIR)/%,$(TSRCS))

# Where the output is
OUTPUT     = $(OBJDIR)/*.$(OBJEXT) $(LIBDIR)/*.$(LIBEXT) $(BINDIR)/*

# Compiler flags
CXXFLAGS   = -Wall -Werror -Wno-reorder -fPIC -I$(HDRDIR) $(ROOTCFLAGS)
LIBFLAGS   = -Wl,--as-needed $(ROOTLIBS) $(EXTRA_ROOTLIBS) -Wl,-rpath,$(ROOTLIBDIR)

# Make the libraries
all: $(BINS)
# Build binaries
$(BINDIR)/%: $(TOBJDIR)/%.$(OBJEXT) $(LIBDIR)/libDP.$(LIBEXT) | $(BINDIR)
	@echo "Linking $@"
	@$(CC) $< -o $@ $(LIBFLAGS) -Wl,--as-needed -L$(LIBDIR) -lDP -Wl,-rpath,$(shell pwd)/$(LIBDIR)
$(TOBJDIR)/%.$(OBJEXT): $(TSRCDIR)/%.$(SRCEXT) $(HDRS) | $(TOBJDIR)
	@echo "Compiling $@"
	@$(CC) $(CXXFLAGS) -c $< -o $@
$(LIBDIR)/libDP.$(LIBEXT): $(OBJS) | $(LIBDIR)
	@$(CC) -shared $^ -o $@ $(LIBFLAGS)
# Build objects
$(OBJDIR)/%.$(OBJEXT): $(DSRCDIR)/%.$(SRCEXT) $(HDRS)
	@echo "Compiling $@"
	@$(CC) $(CXXFLAGS) -c $< -o $@
$(BINDIR) $(LIBDIR) $(TOBJDIR):
	mkdir -p $@
# Remove all the output
clean :
	$(RM) $(OUTPUT)
.PHONY : all clean objects cflags libs libdir
