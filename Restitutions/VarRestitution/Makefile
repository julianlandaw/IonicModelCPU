OBJECTS = varrestitutionmain.cpp
EXECUTABLE = rest.out
CC = g++
#FLAGS = -gencode arch=compute_30,code=sm_30

BEATS = 200

NCELLS = 10001

VARIABLE = xr
#yshift=0.0
#ibarcafac=1.0
#ikfac=1.0
#ikifac=1.0

#yshift=8.0
#ibarcafac=1.15
#ikfac=1.0
#ikifac=2.2

#yshift=-8.0
#ibarcafac=1.0
#ikfac=1.0
#ikifac=1.0

MACROS = -D BEATS=${BEATS} -D NCELLS=${NCELLS} -D VARIABLE=${VARIABLE}

TTMACROS = -D BEATS=${BEATS} -D NCELLS=${NCELLS} -D VARIABLE=${VARIABLE} -D clampki 

default: all

all: LR1 TT
	
LR1:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1$(EXECUTABLE) $(MACROS) -D LR1

TT:
	$(CC) $(FLAGS) $(OBJECTS) -o TT$(EXECUTABLE) $(TTMACROS) -D TT

clean:
	-rm -f *$(EXECUTABLE) *$(EXECUTABLE)

