OBJECTS = stability.cpp
EXECUTABLE = stab.out
CC = g++
#FLAGS = -gencode arch=compute_30,code=sm_30

numpcl=1
numitofac=301
yshift= 8.0
ibarcafac= 1.15
ikfac= 1.0
ikifac= 2.2

#yshift=-8.0
#ibarcafac=1.0
#ikfac=1.0
#ikifac=1.0

MACROS = -D numpcl=${numpcl} -D numitofac=${numitofac} -D yshift=${yshift} -D ibarcafac=${ibarcafac} -D ikfac=${ikfac} -D ikifac=${ikifac}

default: all

all: LR1 LR1_taud TT LR1SK TTSK
	
LR1:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1$(EXECUTABLE) $(MACROS) -D LR1
LR1_taud:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1_taud$(EXECUTABLE) $(MACROS) -D LR1_taud -D curfac=tauXd

TT:
	$(CC) $(FLAGS) $(OBJECTS) -o TT$(EXECUTABLE) $(MACROS) -D TT
LR1SK:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1SK$(EXECUTABLE) $(MACROS) -D LR1SK -D curfac=iskfac
TTSK:
	$(CC) $(FLAGS) $(OBJECTS) -o TTSK$(EXECUTABLE) $(MACROS) -D TTSK -D curfac=iskfac

clean:
	-rm -f *$(EXECUTABLE) *$(EXECUTABLE)

