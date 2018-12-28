OBJECTS = varrestitutionmain.cpp
EXECUTABLE = rest.out
CC = g++
#FLAGS = -gencode arch=compute_30,code=sm_30

BEATS = 2048

NUMVAR1 = 201
NUMVAR2 = 201

VARIABLE1 = casr
VARIABLE2 = nai
curfac = itofac

MACROS = -D curfac=${curfac} -D BEATS=${BEATS} -D NUMVAR1=${NUMVAR1} -D NUMVAR2=${NUMVAR2} -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2}

default: all

all: LR1 TT UCLA OHara LR2
	
LR1:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1$(EXECUTABLE) $(MACROS) -D LR1

TT:
	$(CC) $(FLAGS) $(OBJECTS) -o TT$(EXECUTABLE) $(MACROS) -D TT

UCLA:
	$(CC) $(FLAGS) $(OBJECTS) -o UCLA$(EXECUTABLE) $(MACROS) -D UCLA

OHara:
	$(CC) $(FLAGS) $(OBJECTS) -o OHara$(EXECUTABLE) $(MACROS) -D OHara -D noslowito

LR2:
	$(CC) $(FLAGS) $(OBJECTS) -o LR2$(EXECUTABLE) $(MACROS) -D LR2


clean:
	-rm -f *$(EXECUTABLE) *$(EXECUTABLE)

