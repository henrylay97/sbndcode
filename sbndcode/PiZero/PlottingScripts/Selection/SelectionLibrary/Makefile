# Automatically finds all .cxx files in a directory and compiles it with correct flags
CXX = g++
CXXFLAGS = -g -Wall -fPIC
ROOTFLAGS = `root-config --cflags --glibs`

TARGETS   = $(wildcard *.o)
BINARIES  = $(TARGETS:.o=)
all : $(patsubst %.cxx, %.o, $(wildcard *.cxx))

#--- if making individual targets ---#
%.o : %.cxx %.h $(DEPS)
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $(LDLIBS) -c $< #compile
	ar rvs $*.a $@
.PHONY: clean

clean:
	rm -f $(BINARIES) $(TARGETS)
