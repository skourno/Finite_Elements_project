!**************************************************************************
!                         Numerics.f90  -  description
!                          ---------------------------
!       authors              : skou, vgoutzi
!       email                : skournopoulos@msn.com
!**************************************************************************
! This module includes various numerical analysis subroutines. For further
! information read the dedicated descriptions of each module procedure.
!**************************************************************************
!  CONTAINS :
! !---------------------------------!
! ! Gauss Elimination    - gElim    !
! ! Pi function          - PI       !
! !---------------------------------!
!**************************************************************************
Module Numerics
    private

    Interface GaussElimination
        module procedure gElim
    End Interface

    public GaussElimination
    public PI
contains

    !**********************************************************************
    ! gElim performs Gauss elimination to solve Au = r. The dimensions of
    ! each input variable are A(NxN), u(N), r(N).
    !**********************************************************************
    Subroutine gElim ( N , A , u , r )
        Implicit None
        Integer       :: N, i, j, k
        Real*8        :: A(N,N) , u(N) , r(N)
        Real*8        :: p
        Intent(INOUT) :: N, A, u, r

        do i = 1, N-1
            do j = i+1, N
                p    =  A(j,i) / A(i,i)
                r(j) =  r(j) - p * r(i)

                do k = 1, N
                    A(j,k) = A(j,k) - p * A(i,k)
                enddo
            enddo
        enddo

        do i = N, 1, -1
            u(i) = r(i) / A(i,i)
            do j = i-1, 1, -1
                p      = A(j,i) / A(i,i)
                A(j,i) = 0.
                r(j)   = r(j) - p * r(i)
            enddo
        enddo
    End Subroutine gElim

    !**********************************************************************
    ! PI returns a value of PI with high decimal accuracy
    !**********************************************************************
    Real*8 function PI
        PI = 3.14159265358979323846264338327950288419716939937510
    End function PI

End Module Numerics
