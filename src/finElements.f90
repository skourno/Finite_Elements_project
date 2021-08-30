!**************************************************************************
!                         FinElements.f90  -  description
!                          ---------------------------
!       authors              : skou, vgoutzi
!       email                : skournopoulos@msn.com
!**************************************************************************
! This module includes all subroutines that are used for implementing finite
! element analysis. All components are meant to be private, except the
! interface. To use this module you have to define problem dimensions (1D or
! 2D problems are supported) and the type of basis functions that will be used
! (linear or quadratic - IMode = 1 or 2).
!**************************************************************************

Module FiniteElements
    PRIVATE

    ! 1D variables
    Real*8, allocatable :: ph(:), phc(:), phx(:)

    ! 2D variables
    Real*8, allocatable :: phi(:), phie(:), phic(:)

    ! # of Gauss Points & # of Nodes per Element
    Integer             :: NGaussPoints , NLocalNodes

    ! Gauss Points and Gauss Weights
    Real*8, allocatable :: gp(:) , gw(:)

    Interface SetupStiffnessMatrix
        module procedure ModDriver
    End Interface

    public SetupStiffnessMatrix

contains

    !********************************************************************
    ! ModDriver is a driver routine designed to coordinate the setup of
    ! the stiffeness Matrix based on DIMENSION (dd), basis function
    ! type (IMode), and problem specification (IPr).
    !********************************************************************
    Subroutine ModDriver ( dd , IMode , IProblem )
        Use input,               Only: Linearity
        Use compMethodsProblems, Only: alfa, IPr
        Implicit None
        Real*8     :: success
        Integer    :: dd , IMode , IProblem
        Intent(IN) :: dd , IMode , IProblem

        ! Clean module vars
        if (Allocated(ph))  Deallocate (ph, phc, phx)
        if (Allocated(phi)) Deallocate (phi, phie, phic, gp, gw)

        if ( (dd    <= 0).OR. (dd    >= 3) ) STOP ' *ERROR*: dimension should be 1 or 2'
        if ( (IMode /= 1).AND.(IMode /= 2) ) STOP ' *ERROR*: basis functions 1-lin , 2-quad'

        ! Set compMethodsProblems IPr variable
        IPr = IProblem

        Select case(dd)
            ! 1D Problem
            case(1)
                call SetupStiffMatrix_1D ( IMode )
                call ImposeBoundaryConditions_1D
            ! 2D Problem
            case(2)
                if (Linearity) then ; call SetupStiffMatrix_2D ( IMode )
                               else ; call SetupStiffMatrix_2D_NonLinear( IMode )
                endif
                call ImposeBoundaryConditions_2D
        End Select
    End Subroutine ModDriver

    !********************************************************************
    !********************************************************************
    Subroutine SetupStiffMatrix_1D (IMode)
        Use SystemVars_1D,       Only: SMatrix, RHS, Sol, Node, NOP
        Use input,               Only: ne
        Use compMethodsProblems, Only: Equation, alfa
        Implicit None
        Integer               :: IMode, ig
        Integer               :: IElement               ! Element     #
        Integer               :: i, j                   ! Node Local  #
        Integer               :: IGlobal(3), iGl, jGl   ! Node Global #
        Real*8                :: x, xc, p1, p2
        Intent(IN)            :: IMode

        call GenerateGaussPoints(1,IMode)

        NLocalNodes = IMode + 1

        do IElement = 1, ne

            do i = 1, NLocalNodes
                IGlobal(i) = NOP(IElement,i)
            enddo

            do ig = 1,NGaussPoints
                call TsFun ( 1 , IMode , gp(ig) )

                x      = 0.d0
                xc     = 0.d0
                do i = 1, NLocalNodes
                    x  = x  + Node(IGlobal(i))%r * ph(i)
                    xc = xc + Node(IGlobal(i))%r * phc(i)
                enddo

                do i = 1, NLocalNodes
                    phx(i) = phc(i) / xc
                enddo

                do i = 1, NLocalNodes
                    iGl      = IGlobal(i)
                    RHS(iGl) = RHS(iGl) + Equation(x) * gw(ig)*ph(i)*xc

                    do j = 1, NLocalNodes
                        jGl              = IGlobal(j)
                        SMatrix(iGl,jGl) = SMatrix(iGl,jGl) - gw(ig)*xc*phx(i)*phx(j) &
                                           - gw(ig)*ph(i)*xc*phx(j)*alfa()
                    enddo
                enddo
            enddo
        enddo
    End Subroutine SetupStiffMatrix_1D

    !********************************************************************
    !********************************************************************
    Subroutine SetupStiffMatrix_2D ( IMode )
        Use SystemVars_2D
        Use compMethodsProblems, Only: Equation, alfa
        Implicit None
        Integer             ::  IElement, IGlobal(9), ILocal, JLocal
        Integer             ::  i, j, k, iGl, jGl
        Real*8              ::  xc, xe, yc, ye, dett, x, y
        Real*8              ::  tphx(9) , tphy(9)
        Integer, Intent(IN) ::  IMode

        call GenerateGaussPoints(2,IMode)

        NLocalNodes  =  9

        do IElement  =  1, ne
            do ILocal = 1, 9
                IGlobal(ILocal)  =  NOP(IElement,ILocal)
            enddo

            do i = 1, 3
                do j = 1, 3
                    call tsfun( 2 , IMode , gp(i) , gp(j) )

                    ! isoparametric transformation
                    x   =  0.
		    y   =  0.
		    xc  =  0.
                    xe  =  0.
                    yc  =  0.
                    ye  =  0.
                    do k = 1, 9
			x  = x  + Node(IGlobal(k))%r(1)*phi(k)
                        y  = y  + Node(IGlobal(k))%r(2)*phi(k)
                        xc = xc + Node(IGlobal(k))%r(1)*phic(k)
                        xe = xe + Node(IGlobal(k))%r(1)*phie(k)
                        yc = yc + Node(IGlobal(k))%r(2)*phic(k)
                        ye = ye + Node(IGlobal(k))%r(2)*phie(k)
                    enddo
                    dett = xc*ye - xe*yc

                    do ILocal = 1, 9
                        tphx(ILocal) = (ye*phic(ILocal)-yc*phie(ILocal))/dett
                        tphy(ILocal) = (xc*phie(ILocal)-xe*phic(ILocal))/dett
                    enddo

                    !***** residuals
                    do ILocal = 1, 9
			iGl   = IGlobal(ILocal)
                        do JLocal = 1, 9
			    jGl              = IGlobal(JLocal)
                            SMatrix(iGl,jGl) = SMatrix(iGl,jGl) - gw(i)*gw(j)*dett*(tphx(ILocal)*tphx(JLocal) + tphy(ILocal)*tphy(JLocal))
                        enddo
                    enddo
                enddo
            enddo
        enddo
    End Subroutine SetupStiffMatrix_2D

    !********************************************************************
    !********************************************************************
    Subroutine SetupStiffMatrix_2D_NonLinear ( IMode )
        Use SystemVars_2D
        Use compMethodsProblems, Only: Equation, alfa
        Implicit None
        Integer             ::  IElement, IGlobal(9), ILocal, JLocal
        Integer             ::  i, j, k, iGl, jGl
        Real*8              ::  xc, xe, yc, ye, dett, solu, dsolux, dsoluy, x, y
        Real*8              ::  tphx(9) , tphy(9)
        Integer, Intent(IN) ::  IMode

        call GenerateGaussPoints(2,IMode)

        NLocalNodes  =  9

        do IElement  =  1, ne
            do ILocal = 1, 9
                IGlobal(ILocal)  =  NOP(IElement,ILocal)
            enddo

            do i = 1, 3
                do j = 1, 3
                    call tsfun( 2 , IMode , gp(i) , gp(j) )

                    ! isoparametric transformation
                    x   =  0.
		    y   =  0.
		    xc  =  0.
                    xe  =  0.
                    yc  =  0.
                    ye  =  0.
                    do k = 1, 9
			x  = x  + Node(IGlobal(k))%r(1)*phi(k)
                        y  = y  + Node(IGlobal(k))%r(2)*phi(k)
                        xc = xc + Node(IGlobal(k))%r(1)*phic(k)
                        xe = xe + Node(IGlobal(k))%r(1)*phie(k)
                        yc = yc + Node(IGlobal(k))%r(2)*phic(k)
                        ye = ye + Node(IGlobal(k))%r(2)*phie(k)
                    enddo
                    dett = xc*ye - xe*yc

                    do ILocal = 1, 9
                        tphx(ILocal) = (ye*phic(ILocal)-yc*phie(ILocal))/dett
                        tphy(ILocal) = (xc*phie(ILocal)-xe*phic(ILocal))/dett
                    enddo

                    solu   = 0.d0
                    dsolux = 0.d0
                    dsoluy = 0.d0
                    do ILocal = 1, 9
                        solu   = solu   + Sol(IGlobal(ILocal))*phi(ILocal)
                        dsolux = dsolux + Sol(IGlobal(ILocal))*tphx(ILocal)
                        dsoluy = dsoluy + Sol(IGlobal(ILocal))*tphy(ILocal)
                    enddo

                    !***** residuals
                    do ILocal = 1, 9
			iGl       =   IGlobal(ILocal)
                        RHS(iGl)  =   RHS(IGlobal(ILocal)) + gw(i)*gw(j)*dett*(tphx(ILocal)*dsolux + tphy(ILocal)*dsoluy) &
                                    - gw(i)*gw(j)*dett*phi(ILocal)*Equation(solu, .FALSE.) ! lamda*exp(sol/(1 + alpha*sol))
                        do JLocal = 1, 9
			    iGl               =   IGlobal(JLocal)
                            SMatrix(iGl,jGl)  =   SMatrix(iGl,jGl) - gw(i)*gw(j)*dett*(tphx(ILocal)*tphx(JLocal) + tphy(ILocal)*tphy(JLocal)) &
                                                + gw(i)*gw(j)*dett*phi(ILocal)*phi(JLocal)*Equation(solu, .TRUE.) ! flag true for derivative of equation
                        enddo
                    enddo
                enddo
            enddo
        enddo

    End Subroutine SetupStiffMatrix_2D_NonLinear

    !********************************************************************
    ! Basis function routine
    !********************************************************************
    Subroutine TsFun (dd , IMode , x , y )
        Implicit None
        Integer                       :: N
        Integer, dimension(2)         :: vcase
        Real*8                        :: c
        Real*8 , Intent(IN)           :: x
        Real*8 , Intent(IN), Optional :: y
        Integer, Intent(IN)           :: dd , IMode

        vcase = (/ dd , IMode /)

        ! case 1D Linear
        If ( ALL(vcase==(/1,1/)) ) then
            N      = NLocalNodes
            if ( .NOT. allocated(ph) ) Allocate ( ph(N), phc(N) , phx(N) )

            c      =  x
            ph(1)  =  1.- c
            ph(2)  =  c
            phc(1) = -1.
            phc(2) =  1.
        ! case 1D quadratic
        elseif( ALL(vcase==(/1,2/)) ) then
            N      = NLocalNodes
            if ( .NOT. allocated(ph) ) Allocate ( ph(N), phc(N) , phx(N) )

            c      =  x
            ph(1)  =  2.*c*c - 3.*c + 1.
            ph(2)  = -4.*c*c + 4.*c
            ph(3)  =  2.*c*c -    c
            phc(1) =           4.*c - 3.
            phc(2) =          -8.*c + 4.
            phc(3) =           4.*c - 1.
        ! case 2D biquadratic
        elseif ( ALL(vcase==(/2,2/)) ) then
            N  =  NLocalNodes
            if ( .NOT. allocated(phi) ) Allocate ( phi(N), phic(N) , phie(N) )
            if ( .NOT. PRESENT(y) )     STOP '*** ERROR: missing coordinate y in TsFun routine ***'

            phi(1)   =  l1(x)*l1(y)
            phi(2)   =  l1(x)*l2(y)
            phi(3)   =  l1(x)*l3(y)
            phi(4)   =  l2(x)*l1(y)
            phi(5)   =  l2(x)*l2(y)
            phi(6)   =  l2(x)*l3(y)
            phi(7)   =  l3(x)*l1(y)
            phi(8)   =  l3(x)*l2(y)
            phi(9)   =  l3(x)*l3(y)
            phic(1)  =  l1(y)*dl1(x)
            phic(2)  =  l2(y)*dl1(x)
            phic(3)  =  l3(y)*dl1(x)
            phic(4)  =  l1(y)*dl2(x)
            phic(5)  =  l2(y)*dl2(x)
            phic(6)  =  l3(y)*dl2(x)
            phic(7)  =  l1(y)*dl3(x)
            phic(8)  =  l2(y)*dl3(x)
            phic(9)  =  l3(y)*dl3(x)
            phie(1)  =  l1(x)*dl1(y)
            phie(2)  =  l1(x)*dl2(y)
            phie(3)  =  l1(x)*dl3(y)
            phie(4)  =  l2(x)*dl1(y)
            phie(5)  =  l2(x)*dl2(y)
            phie(6)  =  l2(x)*dl3(y)
            phie(7)  =  l3(x)*dl1(y)
            phie(8)  =  l3(x)*dl2(y)
            phie(9)  =  l3(x)*dl3(y)
        else
            STOP '***ERROR: TsFun Mode requested is not supported'
        endif
    End Subroutine TsFun

    !********************************************************************
    !
    !********************************************************************
    Subroutine ImposeBoundaryConditions_1D
        Use input,         Only : iBC_east, iBC_west, aw, bw, ae, be
        Use SystemVars_1D, Only : SMatrix, RHS, nn
        Implicit None

        ! west
        Select Case ( iBC_west )
        Case (1) ! Dirichlet
            SMatrix(1,:) = 0.d0
            SMatrix(1,1) = 1.d0
            RHS(1)       = aw
        Case (2) ! Neumann
            RHS(1)       = RHS(1) - aw
        Case (3) ! Robin
            SMatrix(1,1) = SMatrix(1,1) + aw
            RHS(1)       = RHS(1)       - bw
        End Select

        ! east
        Select Case ( iBC_east )
        Case (1) ! Dirichlet
            SMatrix(nn,:)  = 0.d0
            SMatrix(nn,nn) = 1.d0
            RHS(nn)        = ae
        Case (2) ! Neumann
            RHS(nn)        = RHS(nn) - ae
        Case (3) ! Robin
            SMatrix(nn,nn) = SMatrix(nn,nn) + ae
            RHS(nn)        = RHS(nn)        - be
        End Select
    End Subroutine ImposeBoundaryConditions_1D


    !********************************************************************
    ! IPr = 3 -> Problem 2 - 2D temperature distribution
    !     B _____north______C
    !      |                \
    !      |                 \ east
    ! west |                  \
    !      |___________________\
    !     A       south         D
    !
    !  west  BCs -> natural dT/dx = 0
    !  east  BCs -> dirichlet Sol = T1 = 200 C
    !  south BCs -> dirichlet Sol = T1 = 200 C
    !  north BCs -> Robin  dT/dy = -h*(T-T0)
    !
    !********************************************************************
    Subroutine ImposeBoundaryConditions_2D
        Use SystemVars_2D
        Use input, Only : Linearity
        Implicit None
        Integer    :: INode, JNode, iEl, i, j, k, ig, iGlobal(9) = 0
        Real*8     :: x1, x2, y1, y2, dett
        Real*8     :: tphx(9) , tphy(9)

        ! Robin BCs
        do iEl = 1, ne
            if ( Element(iEl)%iBC <= 1) CYCLE ! Only Robin and Neumann Boundary Elements, should pass this point

            ! loop over all element nodes
            do i = 1, 9
                iGlobal(i)  =  NOP(iEl,i)
            enddo

            do ig = 1,3
                call tsFun ( 2 , 2 , gp(ig) , 1.d0 )

                ! isoparametric transformation
                x1  =  0.
                x2  =  0.
                y1  =  0.
                y2  =  0.
                do k = 1, 9
                    x1 = x1 + Node(iGlobal(k))%r(1)*phic(k)
                    x2 = x2 + Node(iGlobal(k))%r(1)*phie(k)
                    y1 = y1 + Node(iGlobal(k))%r(2)*phic(k)
                    y2 = y2 + Node(iGlobal(k))%r(2)*phie(k)
                enddo
                dett = x1*y2 - x2*y1

                do i = 1, 9
                    tphx(i) = (y2*phic(i)-y1*phie(i))/dett
                    tphy(i) = (x1*phie(i)-x2*phic(i))/dett
                enddo

                ! loop over only north gauss points
                do i = 3, 9, 3
                    INode                =  iGlobal(i)
                    RHS(INode)           =  RHS(INode) - Element(iEl)%bc1*phi(i)*x1*gw(ig)

                    do j = 3, 9, 3
                        SMatrix(INode,INode) =  SMatrix(INode,INode) + Element(iEl)%bc2*phi(i)*phi(j)*x1*gw(ig)
                    enddo
                enddo
            enddo
        enddo

        Select Case (Linearity)
        Case(.TRUE.)
            ! Dirichlet BCs
            do INode = 1, nn
                if ( Node(INode)%iBC .NE. 1 ) CYCLE

                SMatrix(INode,    :) = 0.d0
                SMatrix(INode,INode) = 1.d0
                RHS(INode)           = Node(INode)%bc1
            enddo
        Case(.FALSE.)
            do INode = 1, nn
                if ( Node(INode)%iBC .NE. 1 ) CYCLE
                SMatrix(INode,    :) = 0.d0
                SMatrix(INode,INode) = 1.d0
                RHS(INode)           = Sol(INode) - Node(INode)%bc1
            enddo
        End Select

    End Subroutine ImposeBoundaryConditions_2D


    !********************************************************************
    ! Generates Gauss Points based on dimension and basis functions
    !********************************************************************
    Subroutine GenerateGaussPoints ( dd , IMode )
        Implicit None
        Integer               :: N
        Integer, Intent(IN)   :: dd , IMode
        Integer, dimension(2) :: vcase

        vcase = (/ dd , IMode /)

        ! case 1D Linear
        If ( ALL(vcase==(/1,1/)) )  then
            NGaussPoints  =  1
            N             =  NGaussPoints

            if ( .NOT. allocated(gp) ) Allocate ( gp(N) , gw(N) )

            gp(1) = 0.5   ;  gw(1) = 1.

        ! case 1D Quadratic
        elseif ( ALL(vcase==(/1,2/)) ) then
            NGaussPoints  =  3
            N             =  NGaussPoints

            if ( .NOT. allocated(gp) ) Allocate ( gp(N) , gw(N) )

            gp(:) = (/0.1127016654    , 0.5           , 0.8872983346    /)
            gw(:) = (/0.27777777777778, 0.444444444444, 0.27777777777778/)

        ! case 2D biquadratic
        elseif ( ALL(vcase==(/2,2/)) ) then
            NGaussPoints  =  9
            N             =  NGaussPoints / 3

            if ( .NOT. allocated(gp) ) Allocate ( gp(N) , gw(N) )

            gp(:) = (/0.1127016654    , 0.5           , 0.8872983346    /)
            gw(:) = (/0.27777777777778, 0.444444444444, 0.27777777777778/)
        endif
    End Subroutine GenerateGaussPoints

    !************************************************************************
    ! Functions supporting TsFun
    !************************************************************************
    Real*8 Function l1(c)
        Real*8 :: c
        l1     = 2.*c**2-3.*c+1.
    End function l1
    !************************************************************************
    Real*8 Function l2(c)
        Real*8 :: c
        l2     = -4.*c**2+4.*c
    End function l2
    !************************************************************************
    Real*8 Function l3(c)
        Real*8 :: c
        l3     = 2.*c**2-c
    End function l3
    !************************************************************************
    Real*8 Function dl1(c)
        Real*8 :: c
        dl1    = 4.*c-3.
    End function dl1
    !************************************************************************
    Real*8 Function dl2(c)
        Real*8 :: c
        dl2    = -8.*c+4.
    End function dl2
    !************************************************************************
    Real*8 Function dl3(c)
        Real*8 :: c
        dl3    = 4.*c-1.
    End function dl3

End Module FiniteElements
