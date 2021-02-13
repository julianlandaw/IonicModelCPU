OBJECTS = apdbifurcationmain4.cpp
EXECUTABLE = bif.out
CC = g++
#FLAGS = -gencode arch=compute_30,code=sm_30
FLAGS = -std=c++11

VARIABLE1 = icalfac
VARIABLE2 = ikrfac
VARIABLE3 = iksfac
VARIABLE4 = ik1fac

SKVARIABLE1 = iskfac
SKVARIABLE2 = skh

gnaped = 0
alphafcass = 0

numpcl = 1
numvar1 = 1
numvar2 = 1
BEATS = 2000
REMOVEBEATS = 1500
PRECTYPE = double

#yshift= 0.0
#ibarcafac= 1.00
#ikfac= 1.0
#ikifac= 1.0

#yshift=-8.0
#ikifac=1.0
#ibarcafac=1.0
#ikfac=1.0

yshift=8.0
ibarcafac=1.15
ikfac=1.0
ikifac=2.2

mshift=0
hshift=0

BIFMACROS = -D numpcl=${numpcl} -D numvar=${numvar} -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D YSHIFT=${yshift} -D IKFAC=${ikfac} -D IKIFAC=${ikifac}

BIFMACROS2 = -D numpcl=${numpcl} -D numvar1=${numvar1} -D numvar2=${numvar2} -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D PRECTYPE=${PRECTYPE}

DV_MAX = 0.1

default: all

all: LR1bif LR1_taudbif LR2bif LR2SKbif TTbif TTModbif UCLAbif OHarabif TP06bif

LR1bif:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1$(EXECUTABLE) $(BIFMACROS2) -D LR1 -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D VARIABLE3=${VARIABLE3} -D VARIABLE4=${VARIABLE4} -D stimulus=-30.0 -D stimduration=2.0 -D TAUDFAC=0.3 -D TAUFFAC=0.3

LR1_taudbif:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1_taud$(EXECUTABLE) $(BIFMACROS) -D LR1_taud -D VARIABLE=tauXd -D stimulus=-30.0 -D stimduration=2.0

LR2bif:
	$(CC) $(FLAGS) $(OBJECTS) -o LR2$(EXECUTABLE) $(BIFMACROS2) -D LR2 -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D VARIABLE3=${VARIABLE3} -D VARIABLE4=${VARIABLE4} -D stimulus=-38.0 -D stimduration=2.0

LR2SKbif:
	$(CC) $(FLAGS) $(OBJECTS) -o LR2SK$(EXECUTABLE) $(BIFMACROS2) -D LR2 -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D VARIABLE3=${VARIABLE3} -D VARIABLE4=${VARIABLE4} -D SK -D UCLAito -D stimulus=-38.0 -D stimduration=2.0

TTbif:
	$(CC) $(FLAGS) $(OBJECTS) -o TT$(EXECUTABLE) $(BIFMACROS2) -D TT -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D VARIABLE3=${VARIABLE3} -D VARIABLE4=${VARIABLE4} -D stimulus=-26.0 -D stimduration=2.0 -D EPI -D xiaodong

TP06bif:
	$(CC) $(FLAGS) $(OBJECTS) -o TP06$(EXECUTABLE) $(BIFMACROS2) -D TP06 -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D VARIABLE3=${VARIABLE3} -D VARIABLE4=${VARIABLE4} -D stimulus=-52.0 -D stimduration=2.0 -D EPI

TTModbif:
	$(CC) $(FLAGS) $(OBJECTS) -o TTMod$(EXECUTABLE) $(BIFMACROS2) -D TTMod -D stimulus=-26.0 -D stimduration=2.0

UCLAbif:
	$(CC) $(FLAGS) $(OBJECTS) -o UCLA$(EXECUTABLE) $(BIFMACROS2) -D UCLA -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D VARIABLE3=${VARIABLE3} -D VARIABLE4=${VARIABLE4} -D stimulus=-30.0 -D stimduration=2.0

#UCLAbif:
#	$(CC) $(FLAGS) $(OBJECTS) -o UCLA$(EXECUTABLE) $(BIFMACROS2) -D UCLA -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D VARIABLE3=${VARIABLE3} -D VARIABLE4=${VARIABLE4} -D noslowito 

OHarabif:
	$(CC) $(FLAGS) $(OBJECTS) -o OHara$(EXECUTABLE) $(BIFMACROS2) -D OHara -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D VARIABLE3=${VARIABLE3} -D VARIABLE4=${VARIABLE4} -D stimulus=-34.0 -D stimduration=2.0 -D mshift=${mshift} -D hshift=${hshift}

Hundbif:
	$(CC) $(FLAGS) $(OBJECTS) -o Hund$(EXECUTABLE) $(BIFMACROS2) -D Hund -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D VARIABLE3=${VARIABLE3} -D VARIABLE4=${VARIABLE4} -D stimulus=-80.0 -D stimduration=1.0

clean:
	-rm -f *$(BIFEXECUTABLE) *$(RESTEXECUTABLE)

