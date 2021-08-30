!**************************************************************************
!                       SystemVars_1D.f90  -  description
!                          ---------------------------
!       author               : skourno, vgkoutzi
!       email                : skournopoulos@msn.com
!**************************************************************************
! This module contains all the variables that need to be defined for 1D
! problems. Use by defining the necessary input (# of Elements, type of
! basis functions and boundary coordinates) and calling "SetupVars"
! subroutine to complete allocations and calulations of the all system
! parameters.
!**************************************************************************
Module SystemVars_1D
    Implicit none

    Type Point
        Real*8   ::  r      ! Coordinate
    End Type Point

    Type(Point), Allocatable ::  Node(:)       ! Node Array
    Real*8     , Allocatable ::  SMatrix(:,:)  ! Stiffness Matrix
    Real*8     , Allocatable ::  RHS(:)        ! Right Hand Side
    Real*8     , Allocatable ::  Sol(:)        ! Solution Vector
    Integer    , Allocatable ::  NOP(:,:)      ! Correlates node # and element #

    Integer                  ::  nn            ! # of Nodes

contains
!*********************************************************************
! Used to setup system variables based on IMode. If IMode = 1 uses
! linear basis, else if IMode = 2 a quadratic base is used
!*********************************************************************
    Subroutine SetupVars ( IMode , xEast , xWest )
        Use input,          Only: ne
        Implicit None
        Integer             ::  INode , IElement , IMode
        Integer             ::  nLocalNodes
        Real*8              ::  xEast , xWest    , dr
        Intent(IN)          ::  IMode , xEast    , xWest

        if     ( IMode == 1 ) then
            nn  =   ne + 1
        elseif ( IMode == 2 ) then
            nn  = 2*ne + 1
        endif

        NLocalNodes  =  IMode + 1

        Allocate ( Node(nn) , SMatrix(nn,nn) , RHS(nn) , Sol(nn) )
        Allocate ( NOP(ne,NLocalNodes) )

        ! Setup Node Coordinates
        dr         =  (xEast - xWest ) / (ne*IMode)
        Node(1)%r  =  xWest
        do INode   =  2, nn
            Node(INode)%r  =  Node(INode - 1)%r + dr
        enddo

        ! Setup NOP
        do IElement = 1, ne
            do INode = 1, NLocalNodes
                NOP( IElement , INode ) = IMode*(IElement - 1) + INode
            enddo
        enddo

!        do IElement = 1, ne
!            write(*,*)  ( NOP(IElement , INode ) , INode = 1, NLocalNodes )
!        enddo

        ! Initialize everything else
        SMatrix(:,:)  =  0.d0
        RHS(:)        =  0.d0
        Sol(:)        =  0.d0

    End Subroutine SetupVars

End Module
