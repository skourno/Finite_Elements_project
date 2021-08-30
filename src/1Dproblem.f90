!***********************************************************
!***  THIS IS THE MAIN PROGRAM FOR PROBLEM 1
!***  one dimensiom finite element problem with linear and
!***  quadratic basis functions
!
!***  AUTHORS: S Kournopoulos, V Gkoutzioupa
!***  (Based on previous work with V Pappa)
!***********************************************************
Program Problem_1
      Use SystemVars_1D
      Use input,          Only: readInput, IBasis, xEast, xWest, nDim, IPr
      Use FiniteElements, Only: SetupStiffnessMatrix
      Use Numerics,       Only: GaussElimination
      Implicit None
      Integer    :: i, j

      call readInput

      call SetupVars ( IBasis , xEast , xWest )      ! *finish variable setup*

      call SetupStiffnessMatrix ( nDim , IBasis , IPr )

      call GaussElimination ( nn , SMatrix , Sol , RHS )

      do i = 1, nn
        print*, Sol(i), 'at node #', i
      enddo
End Program Problem_1
