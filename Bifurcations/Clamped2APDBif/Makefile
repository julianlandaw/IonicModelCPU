OBJECTS = clampedapdbifmain.cpp
EXECUTABLE = bif.out
CC = g++
#FLAGS = -gencode arch=compute_30,code=sm_30

numpcl = 1
numvar1 = 1
numvar2 = 1
numitofac = 1
BEATS = 2200
REMOVEBEATS = 200
DV_MAX = 0.1

VARIABLE1 = cai
VARIABLE2 = casr
VARIABLE3 = nai
VARIABLE4 = ki

curfac = iskfac

MACROS = -D DV_MAX=${DV_MAX} -D numpcl=${numpcl} -D numvar1=${numvar1} -D numvar2=${numvar2} -D numitofac=${numitofac} -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D VARIABLE3=${VARIABLE3} -D VARIABLE4=${VARIABLE4}

default: all

all: LR1 TT

LR1:
	$(CC) $(OBJECTS) -o LR1$(EXECUTABLE) $(MACROS) -D LR1

TT: 	
	$(CC) $(OBJECTS) -o TT$(EXECUTABLE) $(MACROS) -D TT

TTSK:
	$(CC) $(OBJECTS) -o TTSK$(EXECUTABLE) $(MACROS) -D TTSK -D curfac=$(curfac)

clean:
	-rm -f *$(EXECUTABLE)
