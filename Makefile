# All targets with # symbol are self-documenting, i.e. make help or simply make will
# show the targets among available options
#
# User targets are at the bottom
#
ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else

CC = g++
CMSROOT = ./
ROOFITINCLUDE = 
RM  = /bin/rm
ifdef CMSSW_VERSION
	ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
endif
INCLUDE = -I$(CMSROOT) $(ROOFITINCLUDE) 
CFLAGS = -Wall -Wno-unused-function -g -O2 -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS)

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags)

ifeq ($(shell root-config --platform),macosx)
ifdef CMSSW_RELEASE_BASE
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH)/ -L$(CMSSW_RELEASE_BASE)/external/$(SCRAM_ARCH)/lib $(shell root-config --libs) -lEG -lGenVector
else
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
endif
endif


SOURCES = $(wildcard *.cc)  $(wildcard foam/*.cc)  
OBJECTS = $(SOURCES:.cc=.o) LinkDef_out.o
LIB = libME.so


LIBS = $(LIB)

.PHONY: all help compile clean cms2env

libs:	$(LIBS)


$(LIB):	$(OBJECTS) 
	$(QUIET) echo "Linking $(LIB)"; \
	$(LINKER) $(LINKERFLAGS) -shared $(OBJECTS) -o $@ 2>&1|perl -ne 'print if(!/skipping incompatible/)'

LinkDef_out.cxx: LinkDef.h
	$(QUIET) echo "Making CINT dictionaries"; \
        rootcint -f LinkDef_out.cc -c -p $(INCLUDE) LinkDef.h; \
        cat LinkDef.h LinkDef_out.cc > LinkDef_out.cxx; rm LinkDef_out.cc

# General rule for making object files
%.d:	%.cc
	$(QUIET) echo "Checking dependencies for $<"; \
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@
%.o: 	%.cc 
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.cxx 
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

clean:   
	$(QUIET) rm -v -f \
	processed_data.* *.list \
	*.o *.d libME.so LinkDef_out* *_C.so; echo "Done"


endif

