MODULE PARTICLEUTILS
  ! utilities for defining parameters and states of individual particles, particle sets, and domain information

  TYPE PARTICLE
  ! object containing information about individual particles

    ! Particle ID
    INTEGER :: ID

    ! type of particle
    ! total number of particle types depends on system to be simulated
    INTEGER :: PTYPE

    ! particle radius
    DOUBLE PRECISION :: PRAD

    ! current particle position
    DOUBLE PRECISION :: POS

    ! particle state
    ! can take different values depending on system to be simulated
    ! default: +1 = moving forward; -1 = moving backward; 0 = stationary
    INTEGER :: CURSTATE

    ! walk velocity
    DOUBLE PRECISION :: VEL

    ! diffusion coefficient for diffusive particles
    DOUBLE PRECISION :: D

    ! rate of direction reversal
    DOUBLE PRECISION :: KREV

    ! can the particle fuse?
    LOGICAL :: CANFUSE

    ! number of fusion events
    INTEGER :: NFUSE

    ! rate of fusion
    DOUBLE PRECISION, POINTER :: KFUSE(:)

    ! probability of fusion per contact
    DOUBLE PRECISION,POINTER :: PFUSE(:)

    ! array of protein content
    ! contains amounts of different protein types
    DOUBLE PRECISION, POINTER :: PROT(:)

    ! conversion rates between protein types
    DOUBLE PRECISION, POINTER :: PROTCONV(:,:)

    ! region in which particle resides
    INTEGER :: CURREG

    ! have arrays been allocated?
    LOGICAL :: ARRAYSET = .FALSE.

  END TYPE PARTICLE

  TYPE PARTICLEGROUP
  ! object containing information about a group of mixed particle types

    ! number of active particles in the group
    INTEGER :: NPART 

    ! maximum number of particles allowed
    INTEGER :: MAXNPART

    ! ID of the last added particle
    INTEGER :: LASTID

    ! pointer containing particles within the group
    TYPE(PARTICLE), POINTER :: PLIST(:)

    ! number of particle types
    INTEGER :: NTYPE

    ! number of protein types
    INTEGER :: NPROT

    ! number of particles in regions
    INTEGER :: NPARTREG

    ! arrays have been allocated
    LOGICAL :: ARRAYSET = .FALSE.
      
  END TYPE PARTICLEGROUP

  TYPE DOMAIN
  ! object containing information about the domain
    
    ! length of the domain
    DOUBLE PRECISION :: DOMLEN

    ! number of special regions
    INTEGER :: NREG

    ! location of special regions
    DOUBLE PRECISION, POINTER :: REGPOS(:)

    ! arrays have been allocated
    LOGICAL :: ARRAYSET = .FALSE.

  END TYPE DOMAIN


CONTAINS

  SUBROUTINE OUTPUTSNAPSHOT(PGROUP,FILENAME,INFO,APPEND)
  ! Output a snapshot of current configuration
    ! ------------------
    ! input parameters:
    ! PGROUP: object describing group of particles
    ! DOMP: object describing domain
    ! FILENAME: output file
    ! INFO: additional informational float (eg: time)
    ! APPEND: append to file or rewrite?
    ! -------------
    ! output format:

    USE KEYS, ONLY: NTRIALS,NSTEPS
    IMPLICIT NONE

    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    TYPE(PARTICLE), POINTER :: PLIST(:)
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    DOUBLE PRECISION, INTENT(IN) :: INFO(:)
    LOGICAL, INTENT(IN) :: APPEND

    INTEGER :: PC,NINFO
    CHARACTER(LEN=2) :: X1
    CHARACTER(LEN=16) :: FMT1

    IF (APPEND) THEN
      OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN',POSITION='APPEND')
    ELSE
      OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN')
    END IF

    NINFO = SIZE(INFO)
    WRITE (X1,'(I2)') NINFO
    FMT1 = '(5I8,'//TRIM(X1)//'ES18.4E2)'
    WRITE(99,FMT1) NTRIALS,NSTEPS,PGROUP%NPART,PGROUP%NTYPE,PGROUP%NPROT,INFO

    PLIST=>PGROUP%PLIST
    DO PC = 1,PGROUP%NPART
      WRITE(99,'(3I8,1ES18.4E2)') PC,PLIST(PC)%ID,PLIST(PC)%PTYPE,PLIST(PC)%POS
    END DO

    CLOSE(99)

  END SUBROUTINE OUTPUTSNAPSHOT

  SUBROUTINE MAKEPARTICLE(PARTP,PARTID,PTYPE,NTYPE,NPROT,PRAD,VEL,D,KREV,INITPOS,&
                          INITSTATE,CANFUSE,KFUSE,PFUSE,INITPROT,PROTCONV,CURREG)
  ! make a new particle and assign parameters
    USE MT19937, ONLY: GRND
    IMPLICIT NONE

    TYPE(PARTICLE), POINTER :: PARTP
    INTEGER, INTENT(IN) :: PTYPE,NTYPE,NPROT,PARTID,INITSTATE
    DOUBLE PRECISION, INTENT(IN) :: PRAD,VEL,D,INITPOS,KREV
    LOGICAL, INTENT(IN),OPTIONAL :: CANFUSE
    DOUBLE PRECISION,INTENT(IN), OPTIONAL:: KFUSE(NTYPE),PFUSE(NTYPE)
    DOUBLE PRECISION,INTENT(IN), OPTIONAL:: INITPROT(NPROT),PROTCONV(NPROT,NPROT)
    INTEGER,INTENT(IN), OPTIONAL :: CURREG

    IF(PARTP%ARRAYSET) THEN
      PRINT*, 'ERROR IN MAKEPARTICLE: trying to overwrite particle'
      PRINT*, PARTID
      STOP 1
    END IF
    
    PARTP%ID = PARTID
    PARTP%PTYPE = PTYPE
    PARTP%PRAD = PRAD
    
    IF(INITPOS.LT.0D0) THEN
      PARTP%POS = GRND()
    ELSE
      PARTP%POS = INITPOS
    END IF
    PARTP%CURSTATE = INITSTATE

    PARTP%VEL = VEL
    PARTP%D = D
    
    PARTP%KREV = KREV

    IF(PRESENT(CANFUSE)) THEN
      PARTP%CANFUSE = CANFUSE
    ELSE
      PARTP%CANFUSE = .FALSE.
    END IF
    PARTP%NFUSE = 0

    IF(PRESENT(KFUSE)) THEN
      ALLOCATE(PARTP%KFUSE(NTYPE))
      PARTP%KFUSE = KFUSE
    END IF

    IF(PRESENT(PFUSE)) THEN
      ALLOCATE(PARTP%PFUSE(NTYPE))
      PARTP%PFUSE = PFUSE
    END IF

    ALLOCATE(PARTP%PROT(NPROT))
    IF(PRESENT(INITPROT)) THEN
      PARTP%PROT = INITPROT
    ELSE
      PARTP%PROT = 0D0
    END IF

    ALLOCATE(PARTP%PROTCONV(NPROT,NPROT))
    IF(PRESENT(PROTCONV)) THEN
      PARTP%PROTCONV = PROTCONV
    ELSE 
      PARTP%PROTCONV = 0D0
    END IF

    IF(PRESENT(CURREG)) THEN
      PARTP%CURREG = CURREG
    ELSE
      PARTP%CURREG = 0
    END IF

    PARTP%ARRAYSET= .TRUE.

  END SUBROUTINE MAKEPARTICLE

  SUBROUTINE DELETEPARTICLE(PARTP)
  ! delete particle and deallocate arrays
    IMPLICIT NONE
    TYPE(PARTICLE), POINTER :: PARTP

    PARTP%ID = 0
    PARTP%PTYPE = 0
    PARTP%PRAD = 0
    PARTP%POS = 0D0
    PARTP%CURSTATE = 0
    PARTP%VEL = 0D0
    PARTP%D = 0D0
    PARTP%KREV = 0D0

    IF(PARTP%CANFUSE) THEN
      PARTP%CANFUSE = .FALSE.
      PARTP%NFUSE = 0
      DEALLOCATE(PARTP%KFUSE,PARTP%PFUSE)
    END IF
    PARTP%CURREG = 0

    DEALLOCATE(PARTP%PROT,PARTP%PROTCONV)
    
    PARTP%ARRAYSET = .FALSE.
  END SUBROUTINE DELETEPARTICLE

  SUBROUTINE SETUPPARTICLEGROUP(PGROUP,NPART,NTYPE,NPROT,PLIST,INITPARTREG)
  ! set up a cluster of particles
    USE KEYS, ONLY: MAXNPART
    IMPLICIT NONE
    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    TYPE(PARTICLE), TARGET :: PLIST(NPART)
    INTEGER, INTENT(IN) :: NPART,NTYPE,NPROT
    INTEGER, INTENT(IN), OPTIONAL :: INITPARTREG

    IF(PGROUP%ARRAYSET) THEN
      PRINT*, 'ERROR IN SETUPPARTICLEGROUP: trying to overwrite'
      STOP 1
    END IF

    PGROUP%NPART = NPART
    ALLOCATE(PGROUP%PLIST(MAXNPART))
    PGROUP%PLIST(1:NPART) = PLIST
    PGROUP%MAXNPART = MAXNPART
    PGROUP%NTYPE = NTYPE
    PGROUP%NPROT = NPROT
    IF(PRESENT(INITPARTREG)) THEN
      PGROUP%NPARTREG = INITPARTREG
    ELSE
      PGROUP%NPARTREG = 0
    END IF
    PGROUP%LASTID = NPART
    
    PGROUP%ARRAYSET = .TRUE.
  END SUBROUTINE SETUPPARTICLEGROUP

  SUBROUTINE CLEANUPPARTICLEGROUP(PGROUP)
  ! clean up particlegroup object and deallocate arrays
    IMPLICIT NONE
    TYPE(PARTICLEGROUP), POINTER :: PGROUP

    PGROUP%NPART = 0
    DEALLOCATE(PGROUP%PLIST)
    PGROUP%NTYPE = 0
    PGROUP%NPARTREG = 0
    PGROUP%LASTID = 0
    PGROUP%ARRAYSET = .FALSE.
  END SUBROUTINE CLEANUPPARTICLEGROUP

  SUBROUTINE SETUPDOMAIN(DOM,DOMLEN,NREG,REGPOS)
  ! set up domain object
    USE GENUTIL, ONLY:LINSPACE
    IMPLICIT NONE
    TYPE(DOMAIN), POINTER :: DOM
    DOUBLE PRECISION, INTENT(IN) :: DOMLEN
    INTEGER, INTENT(IN) :: NREG
    DOUBLE PRECISION, INTENT(IN) :: REGPOS(NREG)
    DOUBLE PRECISION :: TMP(NREG+2)

    DOM%DOMLEN = DOMLEN
    
    IF(NREG.GT.0) THEN
      DOM%NREG = NREG
      ALLOCATE(DOM%REGPOS(NREG))
      IF(REGPOS(1).LT.0D0) THEN
        TMP = LINSPACE(0D0,DOMLEN,NREG+2)
        DOM%REGPOS = TMP(2:NREG+1)
      ELSE 
        DOM%REGPOS = REGPOS
      END IF
    END IF
    

    DOM%ARRAYSET = .TRUE.
  END SUBROUTINE SETUPDOMAIN

  SUBROUTINE CLEANUPDOMAIN(DOM)
  ! clean up domain object and deallocate arrays
    IMPLICIT NONE
    TYPE(DOMAIN), POINTER :: DOM

    DOM%DOMLEN = 0D0
    IF(DOM%NREG.GT.0) THEN
      DOM%NREG = 0
      DEALLOCATE(DOM%REGPOS)
    END IF

    DOM%ARRAYSET = .FALSE.
  END SUBROUTINE CLEANUPDOMAIN

END MODULE PARTICLEUTILS
