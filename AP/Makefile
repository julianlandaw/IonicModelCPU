OBJECTS = main.cpp
#OBJECTS = main_varypcl.cpp
EXECUTABLE = ap.out
CC = g++
#FLAGS = -gencode arch=compute_30,code=sm_30
yshift=-8.0
ikifac=1.0
#ibarcafac=1.0
ikfac=1.0
PRECTYPE=double
BEATS=2000
REMOVEBEATS=1990
VARIABLE1=itofac
VARIABLE2=iskfac

gnaped = 0.0

default: all

all: LR1ap LR2ap TTap TTHopen LR1SKap TTSKap LR1Mod TTMod UCLAap OHaraap TP06ap

LR1ap:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1$(EXECUTABLE) -D LR1 -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D PRECTYPE=${PRECTYPE} -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D stimulus=-30.0 -D stimduration=2.0

LR2ap:
	$(CC) $(FLAGS) $(OBJECTS) -o LR2$(EXECUTABLE) -D LR2 -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D PRECTYPE=${PRECTYPE} -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D stimulus=-38.0 -D stimduration=2.0

TTap:
	$(CC) $(FLAGS) $(OBJECTS) -o TT$(EXECUTABLE) -D TT -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D clampki -D xiaodong -D noslowito -D stimulus=-26.0 -D stimduration=2.0 -D gnaped=${gnaped}

TP06ap:
	$(CC) $(FLAGS) $(OBJECTS) -o TP06$(EXECUTABLE) -D TP06 -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D stimulus=-52.0 -D stimduration=2.0 -D EPI

TTHopenap:
	$(CC) $(FLAGS) $(OBJECTS) -o TTHopen$(EXECUTABLE) -D TTHopen -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D clampki -D stimulus=-26.0 -D stimduration=2.0 -D gnaped=0 -D EPI

LR1SKap:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1SK$(EXECUTABLE) -D LR1SK -D curfac=iskfac -D stimulus=-30.0 -D stimduration=2.0

LR2SKap:
	$(CC) $(FLAGS) $(OBJECTS) -o LR2SK$(EXECUTABLE) -D LR2 -D SK -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D PRECTYPE=${PRECTYPE} -D curfac=iskfac -D stimulus=-38.0 -D stimduration=2.0

TTSKap:
	$(CC) $(FLAGS) $(OBJECTS) -o TTSK$(EXECUTABLE) -D TTSK -D curfac=iskfac -D stimulus=-26.0 -D stimduration=2.0

#Fentonap:
#	$(CC) $(FLAGS) $(OBJECTS) -o Fenton$(EXECUTABLE) -D Fenton
LR1Mod:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1Mod$(EXECUTABLE) -D LR1Modified -D stimulus=-30.0 -D stimduration=2.0

TTMod:
	$(CC) $(FLAGS) $(OBJECTS) -o TTMod$(EXECUTABLE) -D TTMod -D stimulus=-26.0 -D stimduration=2.0

UCLAap:
	$(CC) $(FLAGS) $(OBJECTS) -o UCLA$(EXECUTABLE) -D UCLA -D BEATS=$(BEATS) -D REMOVEBEATS=$(REMOVEBEATS) -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D noslowito -D stimulus=-30.0 -D stimduration=2.0
#UCLAap:	
#	$(CC) $(FLAGS) $(OBJECTS) -o UCLA$(EXECUTABLE) -D UCLA -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D xiaodong
    
OHaraap:	
	$(CC) $(FLAGS) $(OBJECTS) -o OHara$(EXECUTABLE) -D OHara -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D noslowito -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D stimulus=-34.0 -D stimduration=2.0

clean:
	-rm -f *$(BIFEXECUTABLE) *$(RESTEXECUTABLE)

