!**************************************************************************
!                       SystemVars_2D.f90  -  description
!                          ---------------------------
!       authors              : skou, vgoutzi
!       email                : skournopoulos@msn.com
!**************************************************************************
! This module contains all the variables that need to be defined for 2D
! problems. Use by defining the necessary input (# of Elements, type of
! basis functions and boundary coordinates) and calling "SetupVars" subroutine
! to complete allocations and calulations of the all system parameters.
!**************************************************************************
Module SystemVars_2D
    Implicit none

    Type Point
        Real*8,  dimension(2)  ::  r                          ! Coordinates (x,y)
        Real*8                 ::  bc1    = 0., bc2    = 0.
        Integer                ::  iBC    = 0                 ! Boundary Condition
    End Type Point

    Type Chunk
        Real*8                 ::  bc1    = 0., bc2    = 0.
        Integer                ::  iBC    = 0                 ! Boundary Condition
    End Type Chunk

    Type(Point), Allocatable ::  Node(:)              ! Node Array
    Type(Chunk), Allocatable ::  Element(:)           ! Array saving boundary data for each Element
    Real*8     , Allocatable ::  SMatrix(:,:)         ! Stiffness Matrix
    Real*8     , Allocatable ::  RHS(:)               ! Right Hand Side
    Real*8     , Allocatable ::  Sol(:)               ! Solution Vector
    Real*8     , Allocatable ::  delta(:)             ! deviation from solution. used in non-Linear
    Integer    , Allocatable ::  NOP(:,:)             ! Correlates node # and element #

    Integer                  ::  ne , nn              ! total # of elements and nodes
    Integer                  ::  nnx, nny             ! # of Nodes in each direction

contains
!*********************************************************************
! Used to setup system variables based on IMode. If IMode = 1 uses
! linear basis, else if IMode = 2 a quadratic base is used
!*********************************************************************
    Subroutine systemSetup ( IMode , AB , AD , theta )
        Use Numerics,       Only:  PI
        Use input,          Only:  nex, ney, iBC_west, iBC_east, iBC_north, iBC_south
        Use input,          Only:  aw , bw , ae , be , an ,  bn , as , bs
        Use input,          Only:  Linearity, LGuessedSol, SolGuess
        Implicit None
        Integer           ::  INode , IElement , i , j , k , IMode
        Integer           ::  nLocalNodes
        Real*8            ::  AB , AD , theta   , thetaRad
        Real*8            ::  xorigin , yorigin , xlast , ylast , deltax , deltay
        Intent(IN)        ::  AB , AD , theta , IMode

        ne   =    nex*ney
        nnx  =  2*nex + 1
        nny  =  2*ney + 1
        nn   =    nnx*nny

        NLocalNodes  =  9

        Allocate ( Node(nn) , SMatrix(nn,nn) , RHS(nn) , Sol(nn) )
        Allocate ( NOP(ne,NLocalNodes) )
        Allocate ( Element(ne) )

        if (.NOT. Linearity) Allocate( delta(nn) )

        thetaRad   =  ( theta*PI() ) / 180.  ! deg to rad


        ! Setup Node Coordinates using (0,0) as origin
        xorigin      =  0.d0
        yorigin      =  0.d0
        xlast        =  AD
        ylast        =  AB
        Node(1)%r(1) =  xorigin
        Node(1)%r(2) =  yorigin
        deltay       =  (ylast - yorigin) / ney

        do i = 1 , nny
            Node(i)%r(2) =  Node(1)%r(2) + (i-1)*deltay/2.
            Node(i)%r(1) =  Node(1)%r(1)

          !  print*, i, Node(i)%r(1), Node(i)%r(2)

            if ( theta == 90 ) then
                deltax = (xlast - xorigin) / nex
            else
                deltax = ( (xlast - (i-1)*(deltay/2.)/DTAN(thetaRad) ) ) / nex
            endif

            do j = 2 , nnx
                INode = (j-1)*nny+1
                Node(INode+i-1)%r(1) =Node(i)%r(1) + (j-1)*(deltax/2.)
                Node(INode+i-1)%r(2) =Node(i)%r(2)

           !     print*, INode+i-1, Node(INode+i-1)%r(1), Node(INode+i-1)%r(2)
            enddo
        enddo

        ! Nodal numbering
        IElement = 0
        do i = 1, nex
            do j = 1, ney
                IElement = IElement + 1

                do k = 1, 3
                    INode                 =  3*k - 2
                    NOP(IElement,INode)   =  nny*(2*i+k-3)+2*j-1
                    NOP(IElement,INode+1) =  NOP(IElement,INode) + 1
                    NOP(IElement,INode+2) =  NOP(IElement,INode) + 2
                enddo
            enddo
        enddo

        ! prepare for essential boundary conditions

        ! AB segment
        do i = 1, nny
            Node(i)%iBC     = iBC_west
!           print*, "iBC_west=", iBC_west
            if (Node(i)%iBC > 1) CYCLE

            Node(i)%bc1     = aw
            Node(i)%bc2     = bw
!            print*, aw, bw
        enddo

        ! AD segment
        do i = 1, (nn-nny+1), nny
            Node(i)%iBC     = iBC_south
!           print*, "iBC_south=", iBC_south
            if (Node(i)%iBC > 1) CYCLE

            Node(i)%bc1     = as
            Node(i)%bc2     = bs
!            print*, as, bs
        enddo

        ! BC segment
        do i = nny, nn, nny
            Node(i)%iBC     = iBC_north
 !           print*, "iBC_north=", iBC_north
            if (Node(i)%iBC > 1) CYCLE

            Node(i)%bc1     = an
            Node(i)%bc2     = bn
!            print*, an, bn
        enddo

        ! CD segment
        do i = (nn-nny+1), nn
            Node(i)%iBC     = iBC_east
!          print*, "iBC_east=", iBC_east
            if (Node(i)%iBC > 1) CYCLE

            Node(i)%bc1     = ae
            Node(i)%bc2     = be
!            print*, ae, be
        enddo

        !*** prepare for natural boundary conditions
        do i = (ne-ney+1), ne
            Element(i)%iBC     = iBC_east

            if (Element(i)%iBC <= 1) CYCLE

            Element(i)%bc1     = ae
            Element(i)%bc2     = be
        enddo

        do i = 1, ney
            Element(i)%iBC     = iBC_west

            if (Element(i)%iBC <= 1) CYCLE

            Element(i)%bc1     = aw
            Element(i)%bc2     = bw
        enddo

        do i = ney, ne, ney
            Element(i)%iBC      = iBC_north

            if (Element(i)%iBC <= 1) CYCLE

            Element(i)%bc1      = an
            Element(i)%bc2      = bn
        enddo

        do i = 1, (ne-ney+1), ney
            Element(i)%iBC     = iBC_south

            if (Element(i)%iBC <= 1) CYCLE

            Element(i)%bc1     = as
            Element(i)%bc2     = bs
        enddo

        ! Initialize everything else
        SMatrix(:,:)  =  0.d0
        RHS(:)        =  0.d0
        delta(:)      =  0.d0

        if (LGuessedSol) then ; Sol(:) = SolGuess(:) ; Deallocate(SolGuess)
                         else ; Sol(:) = 0.d0        ;
        endif


    End Subroutine systemSetup

!*********************************************************************
! Use this routine to clean up your system variables.
!
! ICleanType: 1 - Zero everything
!             2 - Deallocate and nullify everything. Then allocate again
!*********************************************************************
    Subroutine Clean_System (ICleanType)
        Implicit None
        Integer    :: ICleanType

        if ((ICleanType .NE. 1).AND.(ICleanType .NE. 2)) STOP "Clean_System: ERROR"
        Select Case(ICleanType)
        Case(1)
            SMatrix(:,:)  = 0.d0  ; RHS(:) = 0.d0 ; delta(:) = 0.d0
        Case(2)
            Deallocate( SMatrix, RHS, delta     )
        End Select

    End Subroutine Clean_System

End Module
