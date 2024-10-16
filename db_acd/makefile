
OS_NAME := $(shell uname -o)

CPPFLAGS = -c -std=c++2b -O3 -Wall -pedantic -I"../acdR/src/"

src=acd.cpp 
acd.obj=../acdR/src/tree.o \
        ../acdR/src/parameters.o \
		../acdR/src/stand.o

OBJECTS=$(src:.cpp=.o)

EXECUTABLE=db_acd

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(acd.obj)
ifeq ($(OS_NAME),Linux)
	g++ $(OBJECTS) $(acd.obj) -lodbc -o $@ 
else
	g++ -static $(OBJECTS) $(acd.obj) -lodbc32 -o $@ 
endif

.cpp.o:
	g++ $(CPPFLAGS) $< -o $@

clean:
	rm *.o 
	rm ../acdR/src/*.o

