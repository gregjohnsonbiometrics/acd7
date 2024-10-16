
OS_NAME = $(shell uname -o)

CPPFLAGS = -c -std=c++2b -O3 -Wall -pedantic -I"../acdR/src/"

src=acd.cpp 
acd.obj=../acdR/src/tree.o \
        ../acdR/src/parameters.o \
		../acdR/src/stand.o

OBJECTS=$(src:.cpp=.o)

EXECUTABLE=acd

cl: $(EXECUTABLE) $(OBJECTS) $(acd.obj)

$(EXECUTABLE): $(OBJECTS) $(acd.obj) 
	g++ -static $(OBJECTS) $(acd.obj) -o $@

.cpp.o:
	g++ $(CPPFLAGS) $< -o $@

clean:
	rm *.o 
	rm ../acdR/src/*.o

