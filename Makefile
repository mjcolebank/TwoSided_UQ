#*************************************************************************#
#                                                                         #
#  Program: Makefile                                                      #
#  Version: 2.0                                                           #
#  By: Mette Olufsen                                                      #
#  Date: 14. Jan. 1997                                                    #
#                                                                         # 
#  A makefile that ensures that all modules are linked together in the    #
#  right order.                                                           #
#*************************************************************************#

# $Id: Makefile,v 1.9 2014/10/17 12:15:01 heine Exp $
# Last modified on October 23, 2014 by M. Umar Qureshi

# FOR TYPICAL MAC
CXX=g++
CXXFLAGS=-O2 -g -Wall -D_REENTRANT -fPIC

# FOR TYPICAL MAC (one or the other)
FC=gfortran
FFLAGS=-O2 -g -Wall
FLIBS=-lgfortran -lquadmath

LIBS=$(FLIBS) -lm

LDFLAGS=-O2

OBJS1=tools.o sor06.o arteries.o
OBJS2=impedance_sub.o new_match.o f90_tools.o

MAIN=sor06

all: $(MAIN)

$(MAIN): $(OBJS1) $(OBJS2) 
	$(CXX) -o $(MAIN) $(LDFLAGS) $(OBJS1) $(OBJS2) $(LIBS)
	
sor06.o: sor06.c sor06.h
	$(CXX) -c $(CXXFLAGS) sor06.c
	
arteries.o: arteries.c arteries.h tools.h sor06.h
	$(CXX) -c $(CXXFLAGS) arteries.c
	
tools.o: tools.C tools.h
	$(CXX) -c $(CXXFLAGS) tools.c
		
new_match.o: new_match.f90 f90_tools.o
	$(FC) -c $(FFLAGS) new_match.f90
	
f90_tools.o: f90_tools.f90
	$(FC) -c $(FFLAGS) f90_tools.f90
	
impedance_sub.o: impedance_sub.f90 f90_tools.o new_match.o
	$(FC) -c $(FFLAGS) impedance_sub.f90
		
clean:
	-rm -f *.o *.mod *.2d Admit* sor06
	
veryclean: clean
	-rm $(MAIN) a.out *~
