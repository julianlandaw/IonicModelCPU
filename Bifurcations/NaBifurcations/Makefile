OBJECTS = nabifurcationmain.cpp
EXECUTABLE = bif.out
CC = g++
#FLAGS = -gencode arch=compute_30,code=sm_30
FLAGS = -std=c++11

VARIABLE1 = iskfac
VARIABLE2 = skh

SKVARIABLE1 = iskfac
SKVARIABLE2 = skh

gnaped = 0

numpcl = 1
numvar1 = 1
numvar2 = 1
BEATS = 2000
REMOVEBEATS = 1500
PRECTYPE = double

yshift= 0.0
#ibarcafac= 1.00
ikfac= 1.0
ikifac= 1.0

#yshift=-8.0
#ikifac=1.0
#ibarcafac=1.0
#ikfac=1.0

#yshift=8.0
#ibarcafac=1.15
#ikfac=1.0
#ikifac=2.2

BIFMACROS = -D numpcl=${numpcl} -D numvar=${numvar} -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D YSHIFT=${yshift} -D IKFAC=${ikfac} -D IKIFAC=${ikifac}

BIFMACROS2 = -D numpcl=${numpcl} -D numvar1=${numvar1} -D numvar2=${numvar2} -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D PRECTYPE=${PRECTYPE}

DV_MAX = 0.1

default: all

all: LR1bif LR1_taudbif LR2bif LR2SKbif TTbif TTModbif UCLAbif OHarabif

LR1bif:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1$(EXECUTABLE) $(BIFMACROS2) -D LR1 -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2}

LR1_taudbif:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1_taud$(EXECUTABLE) $(BIFMACROS) -D LR1_taud -D VARIABLE=tauXd

LR2bif:
	$(CC) $(FLAGS) $(OBJECTS) -o LR2$(EXECUTABLE) $(BIFMACROS2) -D LR2 -D UCLAito -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2}

LR2SKbif:
	$(CC) $(FLAGS) $(OBJECTS) -o LR2SK$(EXECUTABLE) $(BIFMACROS2) -D LR2 -D VARIABLE1=${SKVARIABLE1} -D VARIABLE2=${SKVARIABLE2} -D SK -D UCLAito 

TTbif:
	$(CC) $(FLAGS) $(OBJECTS) -o TT$(EXECUTABLE) $(BIFMACROS2) -D TT -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D noslowito -D clampki -D gnaped=${gnaped} 
TTModbif:
	$(CC) $(FLAGS) $(OBJECTS) -o TTMod$(EXECUTABLE) $(BIFMACROS2) -D TTMod

#UCLAbif:
#	$(CC) $(FLAGS) $(OBJECTS) -o UCLA$(EXECUTABLE) $(BIFMACROS2) -D UCLA -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D xiaodong

UCLAbif:
	$(CC) $(FLAGS) $(OBJECTS) -o UCLA$(EXECUTABLE) $(BIFMACROS2) -D UCLA -D VARIABLE1=$(VARIABLE1) -D VARIABLE2=$(VARIABLE2) -D noslowito 

OHarabif:
	$(CC) $(FLAGS) $(OBJECTS) -o OHara$(EXECUTABLE) $(BIFMACROS2) -D OHara -D noslowito -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2}
clean:
	-rm -f *$(BIFEXECUTABLE) *$(RESTEXECUTABLE)

