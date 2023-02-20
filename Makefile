# Set compiler
CC := g++

# global include
INCLDIR := include

# header only libraries
HOLIB := holib

# resource on flags: https://stackoverflow.com/questions/2855121/what-is-the-purpose-of-using-pedantic-in-the-gcc-g-compiler
# Set compiler flags
CFLAGS := -std=c++11 -pedantic -Wall -O3 -I$(INCLDIR) -I$(HOLIB) -I/usr/include/hdf5/serial #-I/user/include/boost
# Set linker flags
LDFLAGS := -std=c++11 -lncurses -lm -lboost_filesystem -lhdf5 -L/usr/lib/x86_64-linux-gnu/hdf5/serial #-L/usr/include/boost #-lboost_mpi -lboost_serialization 


# Define extensions of files
SRCEXT := cpp
OBJEXT := o
HEADEXT := h
DEPEXT := d

# Define directories to look for certain types of files
SRCDIR := src
OBJDIR := obj
DEPDIR := dep
OUTDIR := out

# Autogenerate lists for all files in the project 
SOURCES := $(wildcard $(SRCDIR)/*.$(SRCEXT))
DEPENDS := $(SOURCES:$(SRCDIR)/%.$(SRCEXT)=$(DEPDIR)/%.$(DEPEXT))
OBJECTS := $(SOURCES:$(SRCDIR)/%.$(SRCEXT)=$(OBJDIR)/%.$(OBJEXT))

# get kernel name to be able to run sed correctly on Darwin (MacOS) or Linux kernels
KERNEL := $(shell uname -s)
ifeq ($(KERNEL), Darwin) 
	SED := sed -i "~"
else
	SED := sed -i
endif

# Define executable directory and executable file
EXEDIR := bin
EXE := exe

# MAIN TARGET
all: $(OBJECTS) $(EXEDIR)/$(EXE)

# build all with debug flags
debug: CFLAGS := -std=c++11 -pedantic -Wall -O3 -I$(INCLDIR) -I$(HOLIB) -I/usr/include/hdf5/serial
# show linker invocation when building debug target
debug: LDFLAGS += -v
debug: all

# define phony target for cleaning project
# use this e.g. to rebuild for debugging
clean:
	@rm $(EXEDIR)/$(EXE) $(OBJECTS) $(DEPENDS)

cleaner:
	@rm $(EXEDIR)/$(EXE) $(OBJECTS) $(DEPENDS)
	@rm -r $(OUTDIR)
	@mkdir $(OUTDIR)


# use this to delete all autogenerated files and directories
superclean: 
	@rm -r $(EXEDIR) $(OBJDIR) $(DEPDIR) $(OUTDIR) $(HOLIB)

# autogenerate directories if they don't exist
$(EXEDIR):
	@mkdir -p $@

$(DEPDIR):
	@mkdir -p $@

$(OBJDIR):
	@mkdir -p $@

$(OUTDIR):
	@mkdir -p $@

$(HOLIB):
	@mkdir -p $@

prepare-run:
	@mkdir -p $(OUTDIR)

# build module as shared library and install it on your system
install: prepare-run | $(HOLIB)
	@echo "Downloading cxxopts.hpp into holib/ ..."
	curl https://raw.githubusercontent.com/jarro2783/cxxopts/master/include/cxxopts.hpp -o $(HOLIB)/cxxopts.hpp
	@echo "... done. Downloading HighFive as zip-file ..."
	@mkdir tmp
	curl https://github.com/BlueBrain/HighFive/archive/refs/heads/master.zip -L -o tmp/HighFive.zip
	@echo "... extracting and copying include/highfive to holib/ ..."
	@unzip -a tmp/HighFive.zip -d tmp/
	@cp -r tmp/HighFive-master/include/highfive $(HOLIB)/.
	@echo "... done. Deleting tmp/ ..."
	@rm -rf tmp
	@echo "... done. Ready to build via 'make all|debug'."

# generating dependency files for all sources
# sed changes '%.o: ...' to '%.o %.d: ...' in dependency file
$(DEPDIR)/%.$(DEPEXT): $(SRCDIR)/%.$(SRCEXT) | $(DEPDIR) $(OUTDIR)
	@echo "Generating dependency file '$@' ..."
	@$(CC) -MM $(CFLAGS) $< -MF $@
	@$(SED) 's,$(*F)\.$(OBJEXT),$*\.$(OBJEXT) $@,' $@
	@rm -f $@~
	@echo "... done."

# include targets from generated dependency files
-include $(DEPENDS)

# build main target
# check if target directory 'bin' already exist via prerequisite
$(EXEDIR)/$(EXE): $(OBJECTS) | $(EXEDIR)
	@echo "Linking binary '$@'..."
	@$(CC) $(CFLAGS) $(OBJECTS) -o $(EXEDIR)/$(EXE) $(LDFLAGS)
	@echo "... done."

# pattern rule to build object from source file
$(OBJDIR)/%.$(OBJEXT): $(DEPDIR)/%.$(DEPEXT) | $(OBJDIR)
	@echo "Compiling '$@' ..."
	@$(CC) -c $(CFLAGS) $(@:$(OBJDIR)/%.$(OBJEXT)=$(SRCDIR)/%.$(SRCEXT)) -o $@
	@echo "... done."

.PHONY: all debug clean cleaner superclean install 
