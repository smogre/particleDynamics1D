MODULE PARTICLEUTILS
  ! utilities for defining parameters and states of individual particles,
  ! particle sets, and domain information

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

    ! Number of different types of possible fusion events
    INTEGER :: NFUSETYPE
    
    ! rate of fusion
    DOUBLE PRECISION, POINTER :: KFUSE(:)

    ! probability of fusion per contact
    DOUBLE PRECISION,POINTER :: PFUSE(:)

    ! number of different protein types
    INTEGER :: NPROT
    
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
    
    ! number of particles in regions
    INTEGER :: NPARTREG

    ! arrays have been allocated
    LOGICAL :: ARRAYSET = .FALSE.      
  END TYPE PARTICLEGROUP


  TYPE PARTICLEPARAMS
  ! object containing the various parameters associated with each type of particle

    INTEGER :: NTYPE ! how many types of particles?

    INTEGER :: MAXNPART ! maximum number of particles allowed

    ! velocity and diffusivity for each type
    DOUBLE PRECISION, POINTER  :: VEL(:), D(:)

    ! production rate
    ! and reversal rate
    DOUBLE PRECISION, POINTER :: KPROD(:), KREV(:)

    ! starting position of each particle type
    DOUBLE PRECISION, POINTER :: STARTPOS(:)

    ! starting state of each particle type
    INTEGER, POINTER :: STARTSTATE(:)

    ! initial protein content of each particle type
    DOUBLE PRECISION, POINTER :: STARTPROT(:,:)

    ! particle radius
    DOUBLE PRECISION , POINTER :: PRAD(:)
    
    ! number of proteins for each type
    INTEGER, POINTER :: NPROT(:)

    ! conversion rates btwn proteins for each type
    DOUBLE PRECISION, POINTER :: PROTCONV(:,:,:)

    ! rate of fusion and/or prob of fusion
    DOUBLE PRECISION, POINTER :: KFUSE(:,:), PFUSE(:,:)
    LOGICAL :: USEPFUSE

    ! arrays allocated?
    LOGICAL :: ARRAYSET = .FALSE.

  END TYPE PARTICLEPARAMS

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

    IMPLICIT NONE

    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    TYPE(PARTICLE), POINTER :: PLIST(:)
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    DOUBLE PRECISION, INTENT(IN) :: INFO(:)
    LOGICAL, INTENT(IN) :: APPEND

    INTEGER :: PC,NINFO
    CHARACTER(LEN=2) :: X1
    CHARACTER(LEN=16) :: FMT1,FMT2

    IF (APPEND) THEN
      OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN',POSITION='APPEND')
    ELSE
      OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN')
    END IF

    ! General info for this snapshot
    NINFO = SIZE(INFO)
    WRITE (X1,'(I2)') NINFO
    FMT1 = '(1I6,'//TRIM(X1)//'ES15.3E2)'
    WRITE(99,FMT1) PGROUP%NPART,INFO

    PLIST=>PGROUP%PLIST
    
    DO PC = 1,PGROUP%NPART
       WRITE (X1,'(I2)') PLIST(PC)%NPROT+1
       FMT2 = '(4I6,'//TRIM(X1)//'ES15.3E2)'
       WRITE(99,FMT2) PC,PLIST(PC)%ID,PLIST(PC)%PTYPE,PLIST(PC)%NPROT,&
                      PLIST(PC)%POS, PLIST(PC)%PROT
    END DO

    CLOSE(99)

  END SUBROUTINE OUTPUTSNAPSHOT

  SUBROUTINE MAKEPARTICLE(PARTP,PARTID,PTYPE,PARAMP,CANFUSE,CURREG)
  ! make a new particle and assign parameters
    USE MT19937, ONLY: GRND
    IMPLICIT NONE

    TYPE(PARTICLE), POINTER :: PARTP
    INTEGER, INTENT(IN) :: PTYPE, PARTID
    TYPE(PARTICLEPARAMS), POINTER :: PARAMP
    LOGICAL, INTENT(IN), OPTIONAL :: CANFUSE
    INTEGER,INTENT(IN), OPTIONAL :: CURREG
    
    INTEGER :: NTYPE, NPROT
    
    IF(PARTP%ARRAYSET) THEN
      PRINT*, 'ERROR IN MAKEPARTICLE: trying to overwrite particle'
      PRINT*, PARTID
      STOP 1
    END IF
    
    PARTP%ID = PARTID
    PARTP%PTYPE = PTYPE
    PARTP%PRAD = PARAMP%PRAD(PTYPE)
    
    IF(PARAMP%STARTPOS(PTYPE).LT.0D0) THEN
      PARTP%POS = GRND()
    ELSE
      PARTP%POS = PARAMP%STARTPOS(PTYPE)
    END IF
    PARTP%CURSTATE = PARAMP%STARTSTATE(PTYPE)

    PARTP%VEL = PARAMP%VEL(PTYPE)
    PARTP%D = PARAMP%D(PTYPE)
    
    PARTP%KREV = PARAMP%KREV(PTYPE)

    IF(PRESENT(CANFUSE)) THEN
      PARTP%CANFUSE = CANFUSE
    ELSE
      PARTP%CANFUSE = .FALSE.
    END IF
    PARTP%NFUSE = 0

    NTYPE = PARAMP%NTYPE

    ! allocate either rates or probs of fusion
    IF (PARAMP%USEPFUSE) THEN
       ALLOCATE(PARTP%PFUSE(NTYPE))
       PARTP%PFUSE = PARAMP%PFUSE(:,PTYPE)
    ELSE
       ALLOCATE(PARTP%KFUSE(NTYPE))
       PARTP%KFUSE = PARAMP%KFUSE(:,PTYPE)
    ENDIF
    

    NPROT = PARAMP%NPROT(PTYPE)
    IF(NPROT.GT.0) THEN
      ALLOCATE(PARTP%PROT(NPROT))
      PARTP%NPROT = NPROT
      PARTP%PROT = PARAMP%STARTPROT(:,PTYPE)
      ALLOCATE(PARTP%PROTCONV(NPROT,NPROT))
      PARTP%PROTCONV = PARAMP%PROTCONV(1:NPROT,1:NPROT,PTYPE)
    END IF
 
    IF(PRESENT(CURREG)) THEN
      PARTP%CURREG = CURREG
    ELSE
      PARTP%CURREG = 0
    END IF

    PARTP%ARRAYSET= .TRUE.

  END SUBROUTINE MAKEPARTICLE

  SUBROUTINE REMOVEPARTICLE(PGROUP,IND)
  ! remove particle from the group
    IMPLICIT NONE
    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    INTEGER, INTENT(IN) :: IND
    TYPE(PARTICLE), POINTER :: PARTP
    INTEGER :: NPARTOLD
    
    PARTP=>PGROUP%PLIST(IND)
    CALL DELETEPARTICLE(PARTP)

    NPARTOLD = PGROUP%NPART
    PGROUP%PLIST(IND:NPARTOLD-1) = PGROUP%PLIST(IND+1:NPARTOLD)

    PGROUP%NPART = NPARTOLD-1
    
  END SUBROUTINE REMOVEPARTICLE
  
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

    IF(PARTP%NPROT.GT.0) DEALLOCATE(PARTP%PROT,PARTP%PROTCONV)
    
    PARTP%ARRAYSET = .FALSE.
  END SUBROUTINE DELETEPARTICLE

  SUBROUTINE SETUPPARTICLEGROUP(PGROUP,NPART,PLIST,INITPARTREG)
  ! set up a cluster of particles
    USE KEYS, ONLY: MAXNPART
    IMPLICIT NONE
    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    TYPE(PARTICLE), TARGET :: PLIST(NPART)
    INTEGER, INTENT(IN) :: NPART
    INTEGER, INTENT(IN), OPTIONAL :: INITPARTREG

    IF(PGROUP%ARRAYSET) THEN
      PRINT*, 'ERROR IN SETUPPARTICLEGROUP: trying to overwrite'
      STOP 1
    END IF

    PGROUP%NPART = NPART
    ALLOCATE(PGROUP%PLIST(MAXNPART))
    PGROUP%PLIST(1:NPART) = PLIST
    PGROUP%MAXNPART = MAXNPART

    
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
    DOM%NREG = 0
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

  SUBROUTINE SETUPPARTICLEPARAMS(PARAMP)
  ! set up and initiate object to track particle parameters
    USE KEYS, ONLY : NTYPE, VEL, DIFFCONST, KPROD, KREV, STARTPOS, STARTSTATE,&
         & PRAD,NPROT,KFUSE,PFUSE,USEPFUSE,PROTCONV,MAXNPROT, STARTPROT,MAXNPART

    
    IMPLICIT NONE
    TYPE(PARTICLEPARAMS), POINTER :: PARAMP
    
    PARAMP%NTYPE = NTYPE
    ALLOCATE(PARAMP%VEL(NTYPE),PARAMP%D(NTYPE),PARAMP%KPROD(NTYPE),&
         & PARAMP%KREV(NTYPE), PARAMP%STARTPOS(NTYPE), PARAMP%STARTSTATE(NTYPE),&
         & PARAMP%PRAD(NTYPE), PARAMP%NPROT(NTYPE))
    ALLOCATE(PARAMP%PROTCONV(MAXNPROT,MAXNPROT,NTYPE))
    ALLOCATE(PARAMP%KFUSE(NTYPE,NTYPE),PARAMP%PFUSE(NTYPE,NTYPE))
    ALLOCATE(PARAMP%STARTPROT(MAXNPROT,NTYPE))

    PARAMP%MAXNPART = MAXNPART
    PARAMP%VEL = VEL(1:NTYPE)
    PARAMP%D = DIFFCONST(1:NTYPE)
    PARAMP%KPROD = KPROD(1:NTYPE)
    PARAMP%KREV = KREV(1:NTYPE)
    PARAMP%STARTPOS = STARTPOS(1:NTYPE)
    PARAMP%STARTSTATE = STARTSTATE(1:NTYPE)
    PARAMP%PRAD = PRAD(1:NTYPE)
    PARAMP%NPROT = NPROT(1:NTYPE)

    PARAMP%PROTCONV = PROTCONV(:,:,1:NTYPE)
    
    PARAMP%USEPFUSE = USEPFUSE
    PARAMP%KFUSE = KFUSE(1:NTYPE,1:NTYPE)
    PARAMP%PFUSE = PFUSE(1:NTYPE,1:NTYPE)

    PARAMP%STARTPROT = STARTPROT(1:NTYPE,:)

    PARAMP%ARRAYSET = .TRUE.    
  END SUBROUTINE SETUPPARTICLEPARAMS

  SUBROUTINE CLEANUPPARTICLEPARAMS(PARAMP)
  ! deallocate object arrays
    IMPLICIT NONE
    TYPE(PARTICLEPARAMS), POINTER :: PARAMP

    DEALLOCATE(PARAMP%VEL, PARAMP%D, PARAMP%KPROD, PARAMP%KREV, &
         & PARAMP%STARTPOS, PARAMP%STARTSTATE, PARAMP%PRAD, PARAMP%NPROT, &
         & PARAMP%PROTCONV, PARAMP%KFUSE, PARAMP%PFUSE, PARAMP%STARTPROT)

    PARAMP%ARRAYSET = .FALSE.
    
  END SUBROUTINE CLEANUPPARTICLEPARAMS
  
END MODULE PARTICLEUTILS
