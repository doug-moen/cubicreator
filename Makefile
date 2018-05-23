#
# Makefile for CuraEngine
#

# simplest working invocation to compile it
#g++ main.cpp modelFile/modelFile.cpp clipper/clipper.cpp -o CuraEngine

CXX ?= g++
CFLAGS += -I. -c -Wall -Wextra -O3 -fomit-frame-pointer
# also include debug symbols
#CFLAGS+=-ggdb
#CFLAGS+=$(LIBS)

LDFLAGS += -static
ODIR = obj
MINGWDIR = C:/MinGW
EXCUTE_DIR = ./../../cubicreator/Cubicreator_newGen/Application
LIBS    =    -L$(MINGWDIR)/lib -lgcc_ss -lstdc++

SOURCES  = bridge.cpp \
		comb.cpp gcodeExport.cpp \
		infill.cpp \
		inset.cpp \
		layerPart.cpp \
		main.cpp \
		optimizedModel.cpp \
		pathOrderOptimizer.cpp \
		polygonOptimizer.cpp \
		raft.cpp \
		settings.cpp \
		skin.cpp \
		skirt.cpp \
		slicer.cpp \
		support.cpp \
		timeEstimate.cpp \
		gcodeControl.cpp \
#		mesh.cpp \
#		meshGroup.cpp \		

SOURCES += clipper/clipper.cpp \
		modelFile/modelFile.cpp \
		utils/gettime.cpp \
		utils/logoutput.cpp \
		utils/util.cpp \
		utils/meshTopology.cpp \
#		utils/socket.cpp \

OBJECTS = $(SOURCES:.cpp=.o)
#EXECUTABLE = ./$(EXCUTE_DIR)/CubiEngine
EXECUTABLE = ./CubiEngine
UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
	OPEN_HTML=firefox
	LDFLAGS += --static
endif
ifeq ($(UNAME), Darwin)
	OPEN_HTML=open
	#For MacOS force to build
	CFLAGS += -force_cpusubtype_ALL -mmacosx-version-min=10.6 -arch x86_64 -arch i386
	LDFLAGS += -force_cpusubtype_ALL -mmacosx-version-min=10.6 -arch x86_64 -arch i386
endif
ifeq ($(UNAME), MINGW32_NT-6.1)
	#For windows make it large address aware, which allows the process to use more then 2GB of memory.
	EXECUTABLE := $(EXECUTABLE).exe
	CFLAGS += -march=pentium4
	LDFLAGS += -Wl,--large-address-aware -lm -lwsock32
endif

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CXX) $(CFLAGS) $< -o $@

test: $(EXECUTABLE)
	python _tests/runtest.py $(abspath $(EXECUTABLE))

## clean stuff
clean:
	rm -f $(EXECUTABLE) $(OBJECTS)

help:
	@cat Makefile |grep \#\#| grep \: |cut -d\# -f3
