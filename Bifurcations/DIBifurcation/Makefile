OBJECTS = dibifurcationmain.cpp
EXECUTABLE = dibif.out
CC = g++
#FLAGS = -gencode arch=compute_30,code=sm_30

numdi = 1
numvar = 1
BEATS = 400
REMOVEBEATS = 200

#yshift= 0.0
#ibarcafac= 1.00
#ikfac= 1.0
#ikifac= 1.0

yshift=-8.0
ikifac=1.0
ibarcafac=1.0
ikfac=1.0

#yshift=8.0
#ibarcafac=1.15
#ikfac=1.0
#ikifac=2.2

BIFMACROS = -D numdi=${numdi} -D numvar=${numvar} -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D yshift=${yshift} -D ibarcafac=${ibarcafac} -D ikfac=${ikfac} -D ikifac=${ikifac}

BIFMACROS2 = -D numdi=${numdi} -D numvar=${numvar} -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS}

default: all

all: LR1bif LR1_taudbif TTbif TTModbif

LR1bif:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1$(EXECUTABLE) $(BIFMACROS) -D LR1

LR1_taudbif:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1_taud$(EXECUTABLE) $(BIFMACROS) -D LR1_taud -D VARIABLE=tauXd

TTbif:
	$(CC) $(FLAGS) $(OBJECTS) -o TT$(EXECUTABLE) $(BIFMACROS2) -D TT
TTModbif:
	$(CC) $(FLAGS) $(OBJECTS) -o TTMod$(EXECUTABLE) $(BIFMACROS2) -D TTMod
clean:
	-rm -f *$(BIFEXECUTABLE) *$(RESTEXECUTABLE)

