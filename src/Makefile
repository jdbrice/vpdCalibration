
# project name and executable name
project = vpd
# source suffix
source = .cpp 
# object files to make
objects = vpd.o histoBook.o calib.o chainLoader.o TOFrPicoDst.o xmlConfig.o splineMaker.o utils.o reporter.o

# ROOT libs and includes
ROOTCFLAGS    	= $(shell root-config --cflags)
ROOTLIBS      	= $(shell root-config --libs)
ROOTGLIBS     	= $(shell root-config --glibs)

# Custom Project includes and libs
includes		= -I../include -I. -I$(ROOTSYS)/include -I$(ROOTDEV)/include 
libraries 			= 	$(ROOTLIBS) -lGui -lHtml -lGuiHtml -lMathCore -lMathMore $(CERNLIBS) $(SYSLIBS) \
#					$(ROOTDEV)/lib/libSplineFit.so \
#					$(ROOTDEV)/lib/libTwoPad.so

cxx 		= g++ -std=c++0x
flags 		= -Wall -g $(ROOTCFLAGS) $(includes) -Wno-write-strings -m32 -fno-inline
compile 	= $(cxx) $(flags) -c 
ldFlags  	= -m32 -g

all: $(project)

$(project): $(objects)
		$(cxx) $(ldFlags) -o  $(project) $(objects) $(ROOTGLIBS) $(libraries)
		@cp $(project) ../bin

%.o: %$(source)
		$(compile)  $<

clean:
		@rm -f $(objects) $(project)
		@rm -f ../bin/$(project)
		@echo "clean done"

full:
	clean
	$(project)
