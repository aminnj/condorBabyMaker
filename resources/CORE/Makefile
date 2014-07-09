include Makefile.arch

#
# stuff to make
#

SOURCES=$(wildcard *.cc) $(wildcard jetsmear/*.cc) $(wildcard MT2/*.cc) $(wildcard lostlepton/*.cc)


OBJECTS=$(SOURCES:.cc=.o)
LIB=libCMS2NtupleMacrosCORE.so

#
# how to make it
#

$(LIB): $(OBJECTS) 
	$(LD) $(LDFLAGS) $(SOFLAGS) $(OBJECTS) $(ROOTLIBS) -o $@

%.o:	%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

#
# target to build
# likelihood id library
#

all: $(LIB) 
clean:
	rm -f *.o \ 
	rm -f jetsmear/*.o \ 
	rm -f MT2/*.o \
	rm -f lostlepton/*.o \
	rm -f *.d \
	rm -f *.so
