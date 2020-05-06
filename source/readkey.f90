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
  INTEGER :: I
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
  VEL = 1D0 ! velocity of particles
  KREV = 1D0 ! reversal rate (per second)
  KJUMP = 1D0 ! jump rate (per second)
  KBIND = 1D0 ! bind rate (per second)
  CRAD = 0.25D0 ! capture radius (microns)
  D = 0.015 ! diffusion coefficient

  ! microtubule properties
  NMT = 5 ! number of microtubules
  MTDYN = .FALSE.
  VGMT = 1D0 !0.172um/s
  VSMT = -1D0 !0.635um/s
  KRMT = 1D0 !0.008/s
  KCMT = 1D0 !0.044/s

  ! domain info
  DOMRAD = 1D0 ! domain radius (micron)
  DOMLEN = 20D0 ! domain length (micron)

  ! input/output  
  SNAPFILE = '*.snap.out' ! snapshot file
  DUMPSNAPSHOTS = .FALSE. ! periodically dump chain snapshots
  APPENDSNAPSHOTS = .FALSE. ! append snapshots to file rather than replacing
  SNAPEVERY = 1 ! how often to dump snapshots
  
  
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
        CASE('D')
          CALL READF(D)
        CASE('DOMRAD')
          CALL READF(DOMRAD)
        CASE('DOMLEN')
          CALL READF(DOMLEN)
        CASE('DELT')
          CALL READF(DELT)
        CASE('KBIND')
          CALL READF(KBIND)
        CASE('KCMT')
          CALL READF(KCMT)
        CASE('KJUMP')
          CALL READF(KJUMP)
        CASE('KREV')
          CALL READF(KREV)
        CASE('KRMT')
          CALL READF(KRMT)
        CASE('MTDYN')
          MTDYN = .TRUE.
        CASE('NMT')
          CALL READI(NMT)
        CASE('NPART')
          CALL READI(NPART)
        CASE('NSTEPS')
          CALL READI(NSTEPS)
        CASE('NTRIALS')
          CALL READI(NTRIALS)
        CASE('PRINTEVERY')
          CALL READI(PRINTEVERY)   
        CASE('RNGSEED')
          CALL READI(RNGSEED)
        CASE('SNAPFILE')
          CALL READA(SNAPFILE)  
        CASE('SNAPSHOTS')
          DUMPSNAPSHOTS = .TRUE.
          IF (NITEMS.GT.1) CALL READI(SNAPEVERY)
          IF (NITEMS.GT.2) CALL READA(SNAPFILE)
          IF (NITEMS.GT.3) CALL READO(APPENDSNAPSHOTS)
        CASE('VEL')
          CALL READF(VEL)
        CASE('VERBOSE')
          CALL READO(VERBOSE)
        CASE('VGMT')
          CALL READF(VGMT)
        CASE('VSMT')
          CALL READF(VSMT)
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
  
  print*, 'number of particles:', NPART
  print*, 'number of microtubules:', NMT
  print*, 'vel, kjump, krev:', VEL,KJUMP,KREV
      
  IF(DOMRAD.GT.0D0) PRINT*, 'domain radius: ',DOMRAD
  
  IF(DOMLEN.GT.0D0) PRINT*, 'domain length: ',DOMLEN
  
  PRINT*, '----------------------------------------------------'
	
END SUBROUTINE READKEY
