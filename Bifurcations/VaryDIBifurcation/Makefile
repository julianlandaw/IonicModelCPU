OBJECTS = varydibifurcationmain.cpp
#OBJECTS = varydibifurcationmain.cpp
EXECUTABLE = varydibif.out
CC = g++
#FLAGS = -gencode arch=compute_30,code=sm_30

numdi = 1
numvar = 1
numalpha = 1
BEATS = 100
REMOVEBEATS = 75

#yshift= 0.0
#ibarcafac= 1.00
#ikfac= 1.0
#ikifac= 1.0

#yshift=-8.0
#ikifac=1.0
#icafac=1.0
#ikfac=1.0
#taufac=1.0

#yshift=8.0
#icafac=1.15
#ikfac=1.0
#ikifac=2.2
#taufac=5.0

yshift=0.0
icafac=1.0
ikfac=1.0
ikifac=1.0
taufac=10.0

BIFMACROS = -D numdi=${numdi} -D numvar=${numvar} -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS}

BIFMACROS2 = -D numdi=${numdi} -D numvar=${numvar} -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS}

default: all

all: LR1bif LR1_taudbif TTbif TTModbif

LR1bif:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1$(EXECUTABLE) $(BIFMACROS) -D LR1 

LR1_taudbif:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1_taud$(EXECUTABLE) $(BIFMACROS) -D LR1_taud -D VARIABLE=tauXd

TTbif:
	$(CC) $(FLAGS) $(OBJECTS) -o TT$(EXECUTABLE) $(BIFMACROS2) -D TT -D clampnai -D clampki
TTModbif:
	$(CC) $(FLAGS) $(OBJECTS) -o TTMod$(EXECUTABLE) $(BIFMACROS2) -D TTMod
clean:
	-rm -f *$(BIFEXECUTABLE) *$(RESTEXECUTABLE)

