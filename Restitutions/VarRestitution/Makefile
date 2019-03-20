OBJECTS = varrestitutionmain.cpp
EXECUTABLE = rest.out
CC = g++
#FLAGS = -gencode arch=compute_30,code=sm_30

BEATS = 200

NCELLS = 1001

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
MACROS2 = -D BEATS=${BEATS} -D NCELLS=${NCELLS}
TTMACROS = -D BEATS=${BEATS} -D NCELLS=${NCELLS} -D VARIABLE=${VARIABLE} -D clampki 
TP06MACROS = -D BEATS=${BEATS} -D NCELLS=${NCELLS} -D VARIABLE=${VARIABLE}

default: all

all: LR1 TT
	
LR1:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1$(EXECUTABLE) $(MACROS) -D LR1

TT:
	$(CC) $(FLAGS) $(OBJECTS) -o TT$(EXECUTABLE) $(TTMACROS) -D TT

TP06:
	$(CC) $(FLAGS) $(OBJECTS) -o TP06$(EXECUTABLE) $(TP06MACROS) -D TP06 -D EPI

OHara:
	$(CC) $(FLAGS) $(OBJECTS) -o OHara$(EXECUTABLE) $(MACROS2) -D OHara -D VARIABLE=${VARIABLE} -D stimulus=-34.0 -D stimduration=2.0 -D noslowito

LR2:
	$(CC) $(FLAGS) $(OBJECTS) -o LR2$(EXECUTABLE) $(MACROS2) -D LR2 -D Miyoshi -D VARIABLE=${VARIABLE} -D stimulus=-38.0 -D stimduration=2.0

clean:
	-rm -f *$(EXECUTABLE) *$(EXECUTABLE)

