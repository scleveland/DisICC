CC=g++ -g
CFLAGS=#-W -Wall
OPT=-O3 -funroll-loops
MATH =-lm
DEBUG=-g
BIOPHYL=-lbpp-phyl
BIONUC=-lbpp-numcalc
BIOUT=-lbpp-utils
BIOSEQ=-lbpp-seq
RM=rm -f
ARCH=#-arch i386
GSL=-lgsl -lgslcblas #-I/sw/include -L/sw/lib

all:	BCFasta.o file_manip.o Blosum.o create.o caps 

caps: caps.cpp file_manip.o create.o BCFasta.o
	$(CC) $(CFLAGS) $(OPT) $(GSL) $(ARCH) $(BIOPHYL) $(BIOUT) $(BIONUC) $(BIOSEQ) -o caps caps.cpp file_manip.o BCFasta.o create.o Blosum.o
Blosum.o: Blosum.cpp
	$(CC) $(CFLAGS) $(OPT) -c Blosum.cpp
create.o: create.cpp
	$(CC) $(OPT) $(CFLAGS) -c create.cpp
file_manip.o: file_manip.cpp
	$(CC) $(OPT) $(CFLAGS) -c file_manip.cpp
BCFasta.o: BCFasta.cpp
	$(CC) $(OPT) $(CFLAGS) -c BCFasta.cpp
clean:
	$(RM) *.o *~ caps
cleanout:
	$(RM) *.o *~ caps *.csv *.out 
