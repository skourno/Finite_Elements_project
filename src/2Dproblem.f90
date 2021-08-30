!***********************************************************
!***  THIS IS THE MAIN PROGRAM FOR PROBLEM 2
!***  2-dimensional linear problem
!***  biquadratic finite element basis functions
!
!***  AUTHORS: S Kournopoulos, V Gkoutzioupa
!***  (Based on previous work with V Pappa)
!***********************************************************
Program problem2
    Use input,          Only: readInput, IBasis, nDim, IPr, AB, AD, theta
    Use SystemVars_2D
    Use FiniteElements, Only: SetupStiffnessMatrix
    Use Numerics,       Only: GaussElimination
    Implicit none
    Integer             :: INode , i , j

    call readInput

    call systemSetup ( IBasis , AB , AD , theta )

    call SetupStiffnessMatrix ( nDim , IBasis , IPr )  ! input: ( 2D problem , biquadratic basis , IPr)

    call GaussElimination ( nn , SMatrix , Sol , RHS )

    open(50,file='res.dat')
    write(50,*) nnx, nny, nn
    do INode = 1, nn
		write(50,*) Node(INode)%r(1), Node(INode)%r(2), Sol(INode)  
    enddo
    close(50)

End Program problem2
