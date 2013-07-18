CC          = g++

CFLAGS      = -c -O3 -I.. -I../include

LDFLAGS     = -lCamorra

SOURCES     = example1.cpp SMseesaw.cpp example2.cpp

OBJECTS     = $(SOURCES:.cpp=.o)  

EXECUTABLES = example1 example2

all: $(SOURCES) $(EXECUTABLES)
	
example1: example1.o
	$(CC) -g example1.o -o $@ $(LDFLAGS)


example2: example2.o SMseesaw.o
	$(CC) -g example2.o SMseesaw.o -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -rf *.o
