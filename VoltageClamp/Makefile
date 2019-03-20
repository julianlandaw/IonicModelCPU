OBJECTS = voltageclampmain.cpp
EXECUTABLE = voltclamp.out
CC = g++
#FLAGS = -gencode arch=compute_30,code=sm_30

BEATS=2000
REMOVEBEATS=0

MACROS= -D BEATS=$(BEATS) -D REMOVEBEATS=$(REMOVEBEATS)

default: all

all: LR1ap TTap UCLAap

LR1ap:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1$(EXECUTABLE) $(MACROS) -D LR1

TTap:
	$(CC) $(FLAGS) $(OBJECTS) -o TT$(EXECUTABLE) $(MACROS) -D TT -D clampki -D xiaodong

UCLAap:
	$(CC) $(FLAGS) $(OBJECTS) -o UCLA$(EXECUTABLE) $(MACROS) -D UCLA -D noslowito

clean:
	-rm -f *$(BIFEXECUTABLE) *$(RESTEXECUTABLE)

