OBJECTS = ttcaibifurcationmain.cpp
EXECUTABLE = caibif.out
CC = g++
#FLAGS = -gencode arch=compute_30,code=sm_30

numcai = 1
numvar = 1
BEATS = 400
REMOVEBEATS = 200

BIFMACROS = -D numcai=${numcai} -D numvar=${numvar} -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS}

default: all

all: TTbif

TTbif:
	$(CC) $(FLAGS) $(OBJECTS) -o TT$(EXECUTABLE) $(BIFMACROS) -D TT -D clampnai -D clampki

clean:
	-rm -f *$(BIFEXECUTABLE) *$(RESTEXECUTABLE)

