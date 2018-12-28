OBJECTS = apdfixedmain.cpp
EXECUTABLE = apdfix.out
CC = g++
#FLAGS = -gencode arch=compute_30,code=sm_30

numpcl=1
numvar1=1

curfac1=nacafac
TYPEVAR1='"INaCa"'
curfac2=ikrfac
TYPEVAR2='"IKr"'

#yshift= 8.0
#ibarcafac= 1.15
#ikfac= 1.0
#ikifac= 2.2

#yshift=-8.0
#ibarcafac=1.0
#ikfac=1.0
#ikifac=1.0

#MACROS = -D numpcl=${numpcl} -D numitofac=${numitofac} -D yshift=${yshift} -D ibarcafac=${ibarcafac} -D ikfac=${ikfac} -D ikifac=${ikifac}
MACROS = -D numvar1=${numvar1} -D curfac1=${curfac1} -D TYPEVAR1=${TYPEVAR1} -D curfac2=${curfac2} -D TYPEVAR2=${TYPEVAR2}

default: all

#all: LR1 LR1_taud TT LR1SK TTSK
all: LR1 TT
	
LR1:
	$(CC) $(FLAGS) $(OBJECTS) -o LR1$(EXECUTABLE) $(MACROS) -D LR1
#LR1_taud:
#	$(CC) $(FLAGS) $(OBJECTS) -o LR1_taud$(EXECUTABLE) $(MACROS) -D LR1_taud -D curfac=tauXd

TT:
	$(CC) $(FLAGS) $(OBJECTS) -o TT$(EXECUTABLE) $(MACROS) -D TT
#LR1SK:
#	$(CC) $(FLAGS) $(OBJECTS) -o LR1SK$(EXECUTABLE) $(MACROS) -D LR1SK -D curfac=iskfac
#TTSK:
#	$(CC) $(FLAGS) $(OBJECTS) -o TTSK$(EXECUTABLE) $(MACROS) -D TTSK -D curfac=iskfac

clean:
	-rm -f *$(EXECUTABLE) *$(EXECUTABLE)

