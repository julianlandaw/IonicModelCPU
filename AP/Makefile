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
REMOVEBEATS=1900
curfac=itofac

default: all

all: LR1ap LR2ap TTap LR1SKap TTSKap LR1Mod TTMod UCLAap OHaraap

LR1ap:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1$(EXECUTABLE) -D LR1 -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D PRECTYPE=${PRECTYPE} -D curfac=${curfac} -D stimulus=-30.0 -D stimduration=2.0

LR2ap:
	$(CC) $(FLAGS) $(OBJECTS) -o LR2$(EXECUTABLE) -D LR2 -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D PRECTYPE=${PRECTYPE} -D curfac=${curfac} -D UCLAito -D stimulus=-38.0 -D stimduration=2.0

TTap:
	$(CC) $(FLAGS) $(OBJECTS) -o TT$(EXECUTABLE) -D TT -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D curfac=${curfac} -D clampki -D noslowito -D stimulus=-26.0 -D stimduration=2.0

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
	$(CC) $(FLAGS) $(OBJECTS) -o UCLA$(EXECUTABLE) -D UCLA -D BEATS=$(BEATS) -D REMOVEBEATS=$(REMOVEBEATS) -D curfac=$(curfac) -D noslowito -D xiaodong -D stimulus=-30.0 -D stimduration=2.0
#UCLAap:	
#	$(CC) $(FLAGS) $(OBJECTS) -o UCLA$(EXECUTABLE) -D UCLA -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D curfac=${curfac} -D xiaodong
    
OHaraap:	
	$(CC) $(FLAGS) $(OBJECTS) -o OHara$(EXECUTABLE) -D OHara -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D noslowito -D curfac=iskfac -D stimulus=-34.0 -D stimduration=2.0

clean:
	-rm -f *$(BIFEXECUTABLE) *$(RESTEXECUTABLE)

