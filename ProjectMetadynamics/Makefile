SHELL=/bin/bash

##########################
## USER OPTION SETTINGS ##
##########################
# Set to yes if OpenMP or Boost library should be used
USE_OMP    = no
USE_GSL    = no
USE_BOOST  = no

# set system compiler and flags
CXX	 = clang++
CXXFLAGS = -c -w -std=c++11 -I$(INCLUDE) -O3

MAINFLAGS =

#########################
## OPTIONAL COMPONENTS ##
#########################

## Include OpenMP headers and flags
ifeq ($(USE_OMP), yes)
        CXXFLAGS  := $(CXXFLAGS) -DOMP -fopenmp
        MAINFLAGS := $(MAINFLAGS) -DOMP -fopenmp
endif
##


## Determine if boost is to be used
ifeq ($(USE_BOOST), yes)
        CXXFLAGS := $(CXXFLAGS) -DBOOST
        #
	# add directories where the boost headers and libraries were installed
        #
	CXXFLAGS  := $(CXXFLAGS) -I/opt/local/include
        MAINFLAGS := $(MAINFLAGS) -I/opt/local/include
endif
##





.PHONY: all
all: MD.o Init.o main.o
	$(CXX) -o run MD.o Init.o main.o $(MAINFLAGS)

Init.o: Init.hpp Init.cpp
	gcc $(CXXFLAGS) -c Init.cpp

MD.o: Init.hpp MD.hpp MD.cpp
	gcc $(CXXFLAGS) -c MD.cpp

main.o: Init.hpp MD.hpp main.cpp
	gcc $(CXXFLAGS) -c main.cpp

clean:
	rm *.o
	rm run
