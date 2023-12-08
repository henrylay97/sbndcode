# Automatically finds all .cxx files in a directory and compiles it with correct flags
CXX = g++
CXXFLAGS = -g -Wall -fPIC
ROOTFLAGS = `root-config --cflags --glibs`

LIBS = SelectionLibrary/*.a /exp/sbnd/app/users/hlay/plotting_utils/*.a

TARGETS   = $(wildcard *.o)
BINARIES  = $(TARGETS:.o=)
all : $(patsubst %.cxx, %.o, $(wildcard *.cxx))

#--- if making individual targets ---#
%.o : %.cxx $(DEPS)
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $(LDLIBS) -c $< #compile
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) $(LDLIBS) -o $* $@  $(LIBS)       #link 
.PHONY: clean

clean:
	rm -f $(BINARIES) $(TARGETS)
