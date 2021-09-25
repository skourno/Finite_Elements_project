# Finite_Elements_project
Finite elements project for the computational methods course (CompMech MSc).

Instructions:
*************************************************************
Just type make with one of the following arguments:

        exc1 - for 1D linear problems
        exc2 - for 2D linear problems
        exc3 - for 2D non-linear problems 
*************************************************************
 There exist 4 types of inputs used:
 
        input_fem.txt              : mandatory - contains essential problem info
        input_fem_bcs.txt          : mandatory - There exist two versions, pick the one corresponding to the specific problem dimensions
        input_fem_NL_Sol_Guess.txt : optional  - contains a prediction for the solution vector of 2D non-Linear problems. Usage is flagged in "input_fem.txt"
*************************************************************
This version has been tested with ifort compiler:
ifort (IFORT) 10.1 20080801
