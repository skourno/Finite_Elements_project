!**************************************************************************
!                          input.f90  -  description
!                          ---------------------------
!       authors              : skou, vgoutzi
!       email                : skournopoulos@msn.com
!**************************************************************************
! This module includes variables loaded through the input files and the
! necessary routine to read them
!**************************************************************************
Module Input
    Implicit None

    !Global input variables
    Integer    :: nDim                 ! problem dimensions
    Integer    :: IBasis               ! type of basis functions
    Integer    :: IPr                  ! computational methods class problem #

    !1-D input variables
    Integer    :: ne                   ! # of elements
    Real*8     :: xEast , xWest        ! boundary coordinates
    Real*8     :: aw, bw, ae, be       ! BC parameters
    Real*8     :: an, bn, as, bs       ! BC parameters
    Integer    :: iBC_east , iBC_west  ! type of BCs
    Integer    :: iBC_north, iBC_south ! type of BCs

    !2-D input variables
    Integer    :: nex                  ! # of elements in x axis
    Integer    :: ney                  ! # of elements in y axis
    Real*8     :: AB                   ! length of west  boundary
    Real*8     :: AD                   ! length of south boundary
    Real*8     :: theta                ! degrees of south-eastern angle
    Logical    :: Linearity            ! defines if the equation is linear or non-linear. TRUE if linear

    ! non-Linear 2D problem variables
    Integer    :: IModeNL              ! non-Linear problem mode (1 -> solve for lamdaInitial, 2 -> solve for many lamdas)
    Real*8     :: lamdaInitial         ! value of lamda
    Real*8     :: lamdaStep            ! IMode = 2 - change lamda with this step
    Real*8     :: lamdaMax             ! lamda threshold that stops the run. lamdaMax should be less than lamdaInitial for negative step
    Logical    :: LGuessedSol          ! TRUE if an input file has been used to guess the solution. If FALSE the initial guess is set ZERO everywhere
    Real*8, Allocatable  :: SolGuess(:)! temp var to save the initial solution guess
contains

    !******************************************************************************
    ! This routine is used to load the input files provided with the code.
    ! There exist 4 types of inputs used:
    ! input_fem.txt              : mandatory - contains essential problem info
    ! input_fem_bcs.txt          : mandatory - There exist two versions, pick the one
    !                              corresponding to the specific problem dimensions
    ! input_fem_NL_Sol_Guess.txt : optional  - contains a prediction for the solution
    !                              vector of 2D non-Linear problems. Usage is flagged
    !                              in "input_fem.txt"
    !******************************************************************************
    Subroutine readInput
        Implicit None
        Integer    :: INode, NNodes
        Real*8     :: x, y
        Character  :: CDump*16

        open(99, File="input_fem.txt", Form='Formatted', Status='OLD',ERR=66)

        call ReadCommentLines(99,4)

        read(99,*) nDim, IBasis

        If (nDim == 1) then
            !read 1D data
            call ReadCommentLines(99,5)
            read(99,*) ne, xWest, xEast, IPr

        else if (nDim == 2) then
            !read 2D data
            call ReadCommentLines(99,22)
            read(99,*) nex, ney, AB, AD, theta, IPr, Linearity

            if (.NOT. Linearity) then
                call ReadCommentLines(99,9)
                read(99,*) IModeNL, lamdaInitial, lamdaStep, lamdaMax, LGuessedSol
            endif
        endif

        close(99)

        open(98, File="input_fem_bcs.txt", Form='Formatted', Status='OLD',ERR=67)

        If (nDim == 1) then
            !read 1D data
            call ReadCommentLines(98,11)
            read(98,*) CDump, CDump, CDump, iBC_west
            read(98,*) CDump, CDump, aw
            read(98,*) CDump, CDump, bw
            call ReadCommentLines(98,4)
            read(98,*) CDump, CDump, CDump, iBC_east
            read(98,*) CDump, CDump, ae
            read(98,*) CDump, CDump, be
        else if (nDim == 2) then
            !read 2D data
            call ReadCommentLines(98,19)
            read(98,*) CDump, CDump, CDump, iBC_west
            read(98,*) CDump, CDump, aw
            read(98,*) CDump, CDump, bw
            call ReadCommentLines(98,4)
            read(98,*) CDump, CDump, CDump, iBC_north
            read(98,*) CDump, CDump, an
            read(98,*) CDump, CDump, bn
            call ReadCommentLines(98,4)
            read(98,*) CDump, CDump, CDump, iBC_east
            read(98,*) CDump, CDump, ae
            read(98,*) CDump, CDump, be
            call ReadCommentLines(98,4)
            read(98,*) CDump, CDump, CDump, iBC_south
            read(98,*) CDump, CDump, as
            read(98,*) CDump, CDump, bs
!            print*, iBC_west, iBC_north, iBC_east, iBC_south
!            print*, aw, bw, an, bn
!            print*,  ae, be, as, bs
        endif

        close(98)

        if (      Linearity  )  return
        if (.NOT. LGuessedSol)  return

        open(97, File="input_fem_NL_Sol_Guess.txt", Form='Formatted', Status='OLD',ERR=68)

        NNodes = (2*ney + 1)*(2*nex + 1)

        Allocate (SolGuess(NNodes))

        do INode = 1, NNodes
            read(97,*) x, y, SolGuess(INode)
        enddo

        close(97)

        return
        66 STOP "readInput : ERROR - input_fem.txt does not exist"
        67 STOP "readInput : ERROR - input_fem_bcs.txt does not exist"
        68 STOP "readInput : ERROR - input_fem_NL_Sol_Guess.txt does not exist"
    End Subroutine readInput

    !-----------------------------------------------------------
    !  Read a number of comment lines from a given file
    !-----------------------------------------------------------
    Subroutine ReadCommentLines(IFileUnit, NLines)
        Implicit None
        Integer      :: IFileUnit, NLines, i ; Character :: CDummy*1
        Intent(IN)   :: IFileUnit, NLines

        do i=1, NLines  ;  read(IFileUnit,*) CDummy  ;  enddo

    End Subroutine ReadCommentLines

End Module Input
