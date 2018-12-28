OBJECTS = cable1dmain.cpp
EXECUTABLE = cable.out
CC = g++

default: all

all: LR1 TT LR1nsca

LR1:
	$(CC) $(OBJECTS) -o LR1$(EXECUTABLE) -D LR1 -D DV_MAX=-1

TT:
	$(CC) $(OBJECTS) -o TT$(EXECUTABLE) $(MACROS) -D TT -D DV_MAX=-1
    
LR2:
	$(CC) $(OBJECTS) -o LR2$(EXECUTABLE) $(MACROS) -D LR2 -D Miyoshi

LR1nsca:
	$(CC) $(OBJECTS) -o LR1nsca$(EXECUTABLE) -D LR1nsca -D DV_MAX=-1

clean:
	-rm -f *$(EXECUTABLE)

