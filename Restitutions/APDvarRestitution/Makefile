OBJECTS = apdvarrestitutionmain.cpp
EXECUTABLE = rest.out
CC = g++
#FLAGS = -gencode arch=compute_30,code=sm_30

NUMDIS = 1
NUMVARS = 10001
yshift=-8.0
ikifac=1.0
ikfac=1.0
ibarcafac=1.0

VARIABLE = xr

MACROS = -D NUMDIS=${NUMDIS} -D NUMVARS=${NUMVARS} -D VARIABLE=${VARIABLE} -D yshift=${yshift} -D ikifac=${ikifac} -D ikfac=${ikfac} -D ibarcafac=${ibarcafac}

default: all

all: LR1 TT
	
LR1:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1$(EXECUTABLE) $(MACROS) -D LR1

TT:
	$(CC) $(FLAGS) $(OBJECTS) -o TT$(EXECUTABLE) $(MACROS) -D TT

clean:
	-rm -f *$(EXECUTABLE) *$(EXECUTABLE)

