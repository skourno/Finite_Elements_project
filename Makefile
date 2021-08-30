FC = ifort
LD = ifort
FFLAGS = -O3 #-traceback
LDFLAGS = -O3 #-traceback

.DEFAULT:
	-touch $@

all: $1
systemVars_2D.o: ./systemVars_2D.f90 numerics.o input.o
	$(FC) $(FFLAGS) -c	./src/systemVars_2D.f90
systemVars_1D.o: ./systemVars_1D.f90 input.o
	$(FC) $(FFLAGS) -c	./src/systemVars_1D.f90
numerics.o: ./numerics.f90
	$(FC) $(FFLAGS) -c	./src/numerics.f90
compMethodsProblems.o: ./compMethodsProblems.f90
	$(FC) $(FFLAGS) -c	./src/compMethodsProblems.f90
input.o: ./input.f90
	$(FC) $(FFLAGS) -c	./src/input.f90
1Dproblem.o: ./1Dproblem.f90 input.o systemVars_2D.o finElements.o numerics.o
	$(FC) $(FFLAGS) -c	./src/1Dproblem.f90
2Dproblem.o: ./2Dproblem.f90 input.o systemVars_2D.o finElements.o numerics.o
	$(FC) $(FFLAGS) -c	./src/2Dproblem.f90
2DproblemNL.o: ./2DproblemNL.f90 input.o systemVars_2D.o finElements.o numerics.o
	$(FC) $(FFLAGS) -c	./src/2DproblemNL.f90
finElements.o: ./finElements.f90 compMethodsProblems.o systemVars_1D.o input.o systemVars_2D.o
	$(FC) $(FFLAGS) -c	./src/finElements.f90

SRC  = ./src/systemVars_2D.f90 ./src/systemVars_1D.f90 ./src/numerics.f90 ./src/compMethodsProblems.f90 ./src/input.f90 ./src/2Dproblem.f90 ./src/finElements.f90
OBJ  = systemVars_2D.o systemVars_1D.o numerics.o compMethodsProblems.o input.o finElements.o
CMD  = $1
M1   = 1Dproblem.o
M2   = 2Dproblem.o
M3   = 2DproblemNL.o
OBJ1 = $(OBJ) $(M1)
OBJ2 = $(OBJ) $(M2)
OBJ3 = $(OBJ) $(M3)

clean: neat
	-rm -f .cppdefs $(OBJ) $(M1) $(M2) $(M3) *.mod exc*
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
exc1: $(OBJ1)
	$(LD) $(OBJ1) -o exc1 $(LDFLAGS)
exc2: $(OBJ2)
	$(LD) $(OBJ2) -o exc2 $(LDFLAGS)
exc3: $(OBJ3)
	$(LD) $(OBJ3) -o exc3 $(LDFLAGS)
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... exc1 - make 1D Linear FEM solver"
	@echo "... exc2 - make 2D Linear FEM solver"
	@echo "... exc3 - make 2D non-Linear FEM solver"
	@echo "... clean"
	@echo "... neat"
	@echo "... tags"


