FC=/mingw64/bin/gfortran


#FFLAGS=-g -fcheck=all -Wall
FFLAGS=-O2 -g -Wall
FLIBS=-lgfortran -lquadmath

LIBS=$(FLIBS) -lm

MAIN=treebranch

all: $(MAIN) 

treebranch.o: treebranch.f90 f90_tools.mod tree_pres.mod
	$(FC) -c $(FFLAGS) treebranch.f90

f90_tools.o: f90_tools.f90
	$(FC) -c $(FFLAGS) f90_tools.f90

f90_tools.mod: f90_tools.f90
	$(FC) -c $(FFLAGS) f90_tools.f90

tree_pres.o: tree_pres.f90 f90_tools.mod
	$(FC) -c $(FFLAGS) tree_pres.f90

tree_pres.mod: tree_pres.f90 f90_tools.mod
	$(FC) -c $(FFLAGS) tree_pres.f90

$(MAIN): treebranch.o f90_tools.o tree_pres.o
	$(FC) -o treebranch treebranch.o f90_tools.o tree_pres.o $(LIBS)

clean: 
	rm -f *.o *.mod *.2d treebranch

veryclean: clean
	-rm -r alpha_beta/*.2d
