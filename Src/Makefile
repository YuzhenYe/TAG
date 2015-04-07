# Makefile
#
#objects := $(patsubst %.cpp, %.o, $(wildcard *.cpp))
objects0 :=  smallapp.o seq.o contigapp.o
objects1 :=  reads.o soap.o path.o seqgraph.o taag.o

CXX := g++

CXXFLAGs :=

BUILD_TYPE :=

ifeq ($(BUILD_TYPE), debug)
	CXXFLAGS += -g
else
	CXXFLAGS += -O2 -w
endif

ifeq ($(pg), 1)
	CXXFLAGS += -pg
else
	CXXFLAGS +=
endif

INC := -I ./

tag: taag_main.cpp $(objects0) $(objects1)
	$(CXX) $(CXXFLAGS) -o tag taag_main.cpp $(objects0) $(objects1) $(INC) 

$(objects0) : %.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(INC)
$(objects1) : %.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(INC)


clean:
	$(RM) *.o tag
