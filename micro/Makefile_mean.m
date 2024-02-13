FC=/mingw64/bin/gfortran


#FFLAGS=-g -fcheck=all -Wall
FFLAGS=-O2 -g -Wall
FLIBS=-lgfortran -lquadmath

LIBS=$(FLIBS) -lm

MAIN=treebranch

all: $(MAIN) 

treebranch_mean.o: treebranch_mean.f90 f90_tools.mod tree_pres_mean.mod
	$(FC) -c $(FFLAGS) treebranch_mean.f90

f90_tools.o: f90_tools.f90
	$(FC) -c $(FFLAGS) f90_tools.f90

f90_tools.mod: f90_tools.f90
	$(FC) -c $(FFLAGS) f90_tools.f90

tree_pres_mean.o: tree_pres_mean.f90 f90_tools.mod
	$(FC) -c $(FFLAGS) tree_pres_mean.f90

tree_pres_mean.mod: tree_pres_mean.f90 f90_tools.mod
	$(FC) -c $(FFLAGS) tree_pres_mean.f90

$(MAIN): treebranch_mean.o f90_tools.o tree_pres_mean.o
	$(FC) -o treebranch_mean treebranch_mean.o f90_tools.o tree_pres_mean.o $(LIBS)

clean: 
	rm -f *.o *.mod *.2d treebranch_mean

veryclean: clean
	-rm -r alpha_beta/*.2d
