!***********************************************************
!***  THIS IS THE MAIN PROGRAM FOR PROBLEM 3
!***  2-dimensional non-linear problem
!***  biquadratic finite element basis functions
!
!***  AUTHORS: S Kournopoulos, V Gkoutzioupa
!***  (Based on previous work with V Pappa)
!***********************************************************
Program problem3
    Use input,               Only: readInput, IBasis, nDim, IPr, AB, AD, theta, IModeNL, lamdaMax, lamdaStep
    Use SystemVars_2D
    Use FiniteElements,      Only: SetupStiffnessMatrix
    Use Numerics,            Only: GaussElimination
    Use CompMethodsProblems, Only: Set_lamda, lamda
    Implicit none
    Integer             :: INode , i , j, iter
    Real*8              :: residual

    call readInput

    call systemSetup ( IBasis , AB , AD , theta )

    do ! infinite loop
        call Set_lamda  ! update lamda value

        iter     = 0
        residual = 100.d0
        do while (residual > 1.D-6) ! solve with six decimal accuracy

            iter = iter + 1
            call SetupStiffnessMatrix ( nDim , IBasis , IPr )  ! input: ( 2D problem , biquadratic basis , IPr)

            call GaussElimination ( nn , SMatrix , delta , RHS )

            ! calculate the residual
            residual = 0.d0
            do i = 1, nn
                !print*, delta(i)
                residual = residual + delta(i)*delta(i)
            enddo
            !residual = residual / nn
            residual = SQRT(residual)
!            print*, "#"
!            print*, "#   Iter  -  residual"
!            print*,     iter, " ", residual
            if (residual > 1.D+10) STOP "ERROR - large residual"

            !update Sol before next iteration
            Sol(:)       = Sol(:) + delta(:)

            call Clean_System(1)
        enddo

!        Print*, "#"
!        Print*, "#####################"
!        Print*, "# lamda value     :", lamda
!        Print*, "# Sol at mid node :", Sol( nn / 2 )

        write(*,*) lamda, Sol(nn/2), iter, residual
        !20 FORMAT( 2x, 6F.3, 2x, 14F.6, 2x, I4, 2x, E3.1E3  )

        if (IModeNL == 1) EXIT ! perform analysis only for lamda initial

        Select Case (INT(SIGN(1.d0,lamdaStep)))
            Case(+1) ! step > 0
                If (lamda > lamdaMax) EXIT
            Case(-1) ! step < 0
                If (lamda < lamdaMax) EXIT
            Case Default
                STOP "Main: ERROR - Step Value can not be ZERO"
        End Select
    enddo


!    if (IModeNL == 1) then
!        open(50,file='res.dat')
!        write(50,*) nnx, nny, nn
!        do INode = 1, nn
!            write(50,*) Node(INode)%r(1), Node(INode)%r(2), Sol(INode)
!        enddo
!        close(50)
!    endif

End Program problem3
