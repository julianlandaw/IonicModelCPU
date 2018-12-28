OBJECTS = apdrestitutionmain.cpp
EXECUTABLE = rest.out
CC = g++
#FLAGS = -gencode arch=compute_30,code=sm_30

#NCELLS = 2000
#yshift= 8.0
#ibarcafac= 1.15
#ikfac= 1.0
#ikifac= 2.2

NCELLS = 1000
#yshift=-8.0
yshift=0.0
ibarcafac=1.0
ikfac=1.0
ikifac=1.0

BEATS=2000
curfac=itofac

#MACROS = -D NCELLS=${NCELLS} -D yshift=${yshift} -D ibarcafac=${ibarcafac} -D ikfac=${ikfac} -D ikifac=${ikifac}
MACROS = -D NCELLS=${NCELLS} -D curfac=${curfac}

MACROS2 = -D NCELLS=${NCELLS} -D curfac=${curfac}

default: all

all: LR1 LR2 TT LR1Modified TTMod
	
LR1:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1$(EXECUTABLE) $(MACROS) -D LR1 -D BEATS=${BEATS}

LR2:
	$(CC) $(FLAGS) $(OBJECTS) -o LR2$(EXECUTABLE) $(MACROS2) -D LR2 -D UCLAito -D BEATS=${BEATS}

TT:
	$(CC) $(FLAGS) $(OBJECTS) -o TT$(EXECUTABLE) $(MACROS) -D TT -D BEATS=${BEATS}

LR1Modified:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1Mod$(EXECUTABLE) $(MACROS) -D LR1CellModified

TTMod:
	$(CC) $(FLAGS) $(OBJECTS) -o TTMod$(EXECUTABLE) $(MACROS) -D TTMod

clean:
	-rm -f *$(EXECUTABLE) *$(EXECUTABLE)

