!**************************************************************************
!                     compMethodsProblems.f90  -  description
!                          ---------------------------
!       authors              : skou, vgoutzi
!       email                : skournopoulos@msn.com
!**************************************************************************
! This module includes all subroutines that are used to build a specific
! problem for the Computational Methods Class of 2014-2015. Specifically
! it is used to impose Boundary Conditions (BCs) and to define the problem
! solved by inserting the appropriate differential equation.
!
! IPr : 1:2 - Exc1, Computational Methods Class 2016
!       3   - Exc2, Computational Methods Class 2016
!       4   - a = 1, Exc3, Computational Methods Class 2016
!       5   - a = 0, Exc3, Computational Methods Class 2016
!**************************************************************************
Module compMethodsProblems

    ! Set this up to calibrate module based in the problem solved
    Integer :: IPr

    Real*8  :: lamda       ! lamda value

contains

    !********************************************************************
    ! Define differential equation. IEq is used so that many different
    ! equations can be called. Here:
    ! IEq = 1  ->  10.*exp(-200.*(x-0.5)**2)
    !********************************************************************
    Real*8 function Equation ( x , LDerivative )
        Implicit None
        Real*8               ::  x
        Logical              ::  LDerivative
        Optional             ::  LDerivative
        Intent(IN)           ::  x

        Select Case (IPr)
            Case(1:2)
                Equation  =  10.*exp(-200.*(x-0.5)**2)
            Case(3)
                Equation  =   0.
            Case(4:5)
                if (.NOT. LDerivative) then ; Equation = lamda*EXP(x/(1. + alfa()*x))
                                       else ; Equation = lamda*EXP(x/(1. + alfa()*x))*(1./(1. + alfa()*x)**2.)
                endif
            Case Default
                STOP '*** ERROR: tried to solve for an undefined equation ***'
        End Select
    End function Equation

    !********************************************************************
    ! Parameter that modifies the stiffness matrix depending on the
    ! defined problem
    ! IPr = 1  ->  alfa = 10. (parameter from problem 1)
    !********************************************************************
    Real*8 function alfa()
        Implicit None

        Select Case(IPr)
            Case(1:2)
                alfa  =  10.d0
            Case(4)
                alfa  =   1.d0
            Case(5)
                alfa  =   0.d0
            Case Default
                STOP '*** ERROR: undefined alfa function in compMethodsProblems ***'
        End Select
    End function alfa

    !********************************************************************
    ! Lamda is aParameter that modifies the stiffness matrix depending on the
    ! defined problem and running non-Linear Mode. Change lamda modifies it
    ! by using a step picked by the user in the input
    !********************************************************************
    Subroutine Set_lamda
        Use Input,   Only : lamdaStep, lamdaInitial
        Implicit None
        Logical          :: LFirstEntry = .TRUE.
        Save             :: LFirstEntry

        if (LFirstEntry) then ; lamda = lamdaInitial      ; LFirstEntry = .FALSE.
                         else ; lamda = lamda + lamdaStep
        endif

    End Subroutine Set_lamda

End Module compMethodsProblems
