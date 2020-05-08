SUBROUTINE READKEY
! this subroutine reads in keywords from a parameter file
! it sets the various global variables defined in KEYS module
! name of the parameter file is param.* where * is a keyword argument
! if no keyword argument is supplied, the default is just a file called param
! The EXTRAPARAMFILES keyword will allow extra parameter files to be 
! read in as well

  USE KEYS
  USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
  USE GENUTIL

  IMPLICIT NONE

  ! ---- stuff for inputing the parameter file in free format --------
  CHARACTER*100 :: ARG ! command line argument
  INTEGER :: NUMARG ! number of command line arguments
  INTEGER :: NITEMS ! number of items on the line in the parameter file
  INTEGER :: PF ! input file unit
  LOGICAL :: FILEEND=.FALSE. ! done reading file?
  CHARACTER*100 :: WORD ! keyword
  ! -------------- for reading multiple parameter files --------  
  INTEGER, PARAMETER :: MAXNFILES = 10
  CHARACTER*100 :: PARAMFILES(MAXNFILES)
  INTEGER :: NPARAMFILES, NPARAMREAD
  ! ------ for initializing random number generator
  INTEGER :: TIMEVAL(8), SEED
  ! ---------------- temporary variables ---------------
  INTEGER :: I,J
  CHARACTER*100 :: DUMSTR
  LOGICAL :: LDUM

  ! ------------------------
  ! set variable defaults
  ! ------------------------

  ! general program control
  ACTION = 'NONE'
  RNGSEED = 0
  VERBOSE = .FALSE.

  ! simulation info
  NTRIALS = 1 ! number of trials

  ! particle properties
  NPART = 1 ! number of particles
  NTYPE = 1 ! number of particle types
  VEL = 1D0 ! velocity of particles
  PRAD = 1D-2 ! particle radius (different radii can be assigned to different particles)
  BUFLEN = 0D0 ! buffer length to determine contact

  KREV = 1D0 ! rate of direction reversal
  KSTOP = 0D0 ! rate of stopping
  D = 0D0 ! diffusion coefficient

  ! protein dynamics
  NPROT = 1
  PROTCONV = 0

  ! fusion details
  KFUSE = 0D0 ! matrix of fusion rates
  PFUSE = 0D0 ! matrix of fusion probabilities
  DOFUSION = .FALSE.
  USEPFUSE = .FALSE.

  ! particle production
  KPROD = 0D0
  STARTPOS = 0D0 ! starting position
  STARTPROT = -1D0 ! amount of protein to start with

  ! domain info
  DOMLEN = 1D0 ! domain length
  NREG = 0 ! number of regions
  REGPOS = 0D0 ! region positions


  ! input/output  
  SNAPFILE = '*.snap.out' ! snapshot file
  OUTFILE = '*.out' ! output file
  DUMPSNAPSHOTS = .FALSE. ! periodically dump snapshots
  APPENDSNAPSHOTS = .FALSE. ! append snapshots to file rather than replacing
  SNAPEVERY = 1 ! how often to dump snapshots
  OUTEVERY = 1 ! how often to write to output file
  
  ! dynamics
  DELT = 1D-4 ! time step
  NSTEPS = INT(1D6) ! number of brownian steps to run
  PRINTEVERY = 1 ! how often to print output
    
  ! -------------------------
  ! Read in all parameter files, starting with the ones specified on command line
  ! --------------------------

  PF = 55 ! i/o unit number to be used for parameter files

  ! get input parameter files from command line
  NPARAMFILES = 0
  NUMARG = COMMAND_ARGUMENT_COUNT()  
  IF (NUMARG==0) THEN
    NPARAMFILES = 1
    PARAMFILES(1) = 'param'
    ARG = ''
  ELSE
    DO I = 1,NUMARG
      CALL GETARG(I, ARG)
      NPARAMFILES = NPARAMFILES + 1
      WRITE(DUMSTR,'(A)') 'param.' //TRIM(ADJUSTL(ARG))
      PARAMFILES(NPARAMFILES) = DUMSTR
    END DO
    ! reset arg to its original value
    IF (NUMARG.GT.1) CALL GETARG(1,ARG)
  ENDIF

  NPARAMREAD = 0 ! keep track of how many files have been read
  DO WHILE (NPARAMREAD.LT.NPARAMFILES)
    NPARAMREAD = NPARAMREAD + 1

    PRINT*, 'Reading parameter file: ', PARAMFILES(NPARAMREAD)
    INQUIRE(FILE=PARAMFILES(NPARAMREAD),EXIST=LDUM)
    IF (.NOT.LDUM) THEN
      PRINT*, 'ERROR in READKEY: Parameter file ', TRIM(ADJUSTL(PARAMFILES(NPARAMREAD))), ' does not exist.'
      STOP 1
    ENDIF
    OPEN(UNIT=PF, FILE=PARAMFILES(NPARAMREAD), STATUS='OLD')

    ! read in the keywords one line at a time
    DO 
      CALL READLINE(PF,FILEEND,NITEMS)
      IF (FILEEND.and.nitems.eq.0) EXIT

      ! skip empty lines
      IF (NITEMS.EQ.0) CYCLE

      ! Read in the keyword for this line
      CALL READA(WORD,CASESET=1)

      ! Skip any empty lines or any comment lines
      IF (WORD(1:1).EQ.'#') CYCLE

      SELECT CASE(WORD) ! pick which keyword
        CASE('ACTION')
          CALL READA(ACTION, CASESET=1)
        CASE('BUFLEN')
          CALL READF(BUFLEN)
        CASE('D')
          CALL READF(D(1))
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(D(I))
          END DO
        CASE('DOMLEN')
          CALL READF(DOMLEN)
        CASE('DELT')
          CALL READF(DELT)
        CASE('DOFUSION')
          CALL READO(DOFUSION)
        CASE('KFUSE')
          CALL READI(I)
          DO J = 1,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(KFUSE(I,J))
          END DO
        CASE('KREV')
          CALL READF(KREV(1))
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(KREV(I))
          END DO
        CASE('KSTOP')
          CALL READF(KSTOP(1))
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(KSTOP(I))
          END DO
        CASE('NPART')
          CALL READI(NPART(1))
          NPART = NPART(1)
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READI(NPART(I))
          END DO
        CASE('NPROT')
          CALL READI(NPROT)
        CASE('NREG')
          CALL READI(NREG)
        CASE('NSTEPS')
          CALL READI(NSTEPS)
        CASE('NTRIALS')
          CALL READI(NTRIALS)
        CASE('NTYPE')
          CALL READI(NTYPE)
        CASE('PFUSE')
          CALL READI(I)
          DO J = 1,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(PFUSE(I,J))
          END DO
        CASE('PRINTEVERY')
          CALL READI(PRINTEVERY)
        CASE('PRAD')   
          CALL READF(PRAD(1))
          PRAD = PRAD(1)
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(PRAD(I))
          END DO
        CASE('PROTCONV')
          CALL READI(I)
          DO J = 1,MIN(MAXNPROT,NITEMS-1)
            CALL READF(PROTCONV(I,J))
          END DO
        CASE('REGPOS')
          CALL READF(REGPOS(1))
          DO I = 2,MIN(MAXNREG,NITEMS-1)
            CALL READF(REGPOS(I))
          END DO
        CASE('RNGSEED')
          CALL READI(RNGSEED)
        CASE('SNAPFILE')
          CALL READA(SNAPFILE)  
        CASE('SNAPSHOTS')
          DUMPSNAPSHOTS = .TRUE.
          IF (NITEMS.GT.1) CALL READI(SNAPEVERY)
          IF (NITEMS.GT.2) CALL READA(SNAPFILE)
          IF (NITEMS.GT.3) CALL READO(APPENDSNAPSHOTS)
        CASE('STARTPOS')
          CALL READF(STARTPOS(1))
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(STARTPOS(I))
          END DO
        CASE('STARTPROT')
          CALL READF(STARTPROT(1))
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(STARTPROT(I))
          END DO
        CASE('STARTSTATE')
          CALL READI(STARTSTATE(1))
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READI(STARTSTATE(I))
          END DO
        CASE('USEPFUSE')
          CALL READO(USEPFUSE)
        CASE('VEL')
          CALL READF(VEL(1))
          VEL = VEL(1)
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(VEL(I))
          END DO
        CASE('VERBOSE')
          CALL READO(VERBOSE)
        CASE DEFAULT
          print*, 'ERROR: unidentified keyword ', TRIM(WORD), " Will ignore."
      END SELECT
    ENDDO
    CLOSE(PF)
  ENDDO

  ! -----------------
  ! check validity of some values, raise errors or adjust as necessary
  ! -----------------  

  ! ----------- fix file names -----------
  CALL REPLACESUBSTR(SNAPFILE,'*',TRIM(ADJUSTL(ARG)))
  ! ---------------------------

  ! Initiate random number generator 
  IF (RNGSEED.EQ.0) THEN
    ! use the current time of day in milliseconds
    CALL DATE_AND_TIME(VALUES=TIMEVAL)
    SEED = TIMEVAL(5)*3600*1000 + TIMEVAL(6)*60*1000 + TIMEVAL(7)*1000 + TIMEVAL(8)
  ELSEIF (RNGSEED.EQ.-1) THEN
    ! use the last 5 characters in the command-line argument
    SEED = STRING2NUM(TRIM(ADJUSTL(ARG)))    
  ELSEIF (RNGSEED.EQ.-2) THEN
    ! use the last 4 characters in the command-line argument 
    ! and additionally the millisecond time 
    CALL DATE_AND_TIME(VALUES=TIMEVAL)
    SEED = STRING2NUM(TRIM(ADJUSTL(ARG)),TIMEVAL(8))
  ELSE
    ! use this seed directly
    SEED = RNGSEED
  ENDIF

  print*, 'Initiating Mersenne twister random number generator with seed:', SEED
  CALL SGRND(SEED)

  print*, '------------Parameter values : -------------------'
  print*, 'ACTION: ', TRIM(ADJUSTL(ACTION))
  IF (DUMPSNAPSHOTS) THEN
    PRINT*, 'Dumping snapshot every', SNAPEVERY,'steps. In file:', TRIM(ADJUSTL(SNAPFILE))
  ENDIF
  
  print*, 'number of particle types:', NTYPE
  print*, 'velocities:', VEL

  IF(DOMLEN.GT.0D0) PRINT*, 'domain length: ',DOMLEN
  
  PRINT*, '----------------------------------------------------'
	
END SUBROUTINE READKEY
