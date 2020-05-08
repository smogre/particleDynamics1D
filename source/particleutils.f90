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
    DOUBLE PRECISION, POINTER :: POS

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
    LOGICAL :: ARRAYSET

  END TYPE PARTICLE

  TYPE PARTICLEGROUP
  ! object containing information about a group of mixed particle types

    ! number of active particles in the group
    INTEGER :: NPART 

    ! maximum number of particles allowed
    INTEGER :: MAXNPART

    ! pointer containing particles within the group
    TYPE(PARTICLE), POINTER :: PLIST(:)

    ! number of particle types
    INTEGER :: NTYPE

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
    DOUBLE PRECISION :: NREG

    ! location of special regions
    DOUBLE PRECISION, POINTER :: REGPOS(:)

  END TYPE DOMAIN


CONTAINS

  SUBROUTINE OUTPUTSNAPSHOT(PARTP,FILENAME,INFO,APPEND)
  ! Output a snapshot of current particle configuration
  ! ------------------
  ! input parameters:
  ! PARTP: object describing group of particles
  ! FILENAME: output file
  ! INFO: additional informational float (eg: time)
  ! APPEND: append to file or rewrite?
  ! -------------
  ! output format:
  ! 1st line -> number of particles.
  ! Successive lines -> one line per particle
  ! for each particle, lists INDEX, 3 coordinates of position,
  ! followed by orientation in matrix form (x axis, y axis, z axis of particle frame)
    USE KEYS, ONLY : NSTEPS,NTRIALS,NPART,NMT,SNAPEVERY
    IMPLICIT NONE

    TYPE(PARTICLEGROUP), POINTER :: PARTP
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    DOUBLE PRECISION, INTENT(IN) :: INFO(:)
    LOGICAL, INTENT(IN) :: APPEND

    INTEGER :: PC,MC,NINFO,NSAVE
    CHARACTER(LEN=2) :: X1
    CHARACTER(LEN=19) :: FMT1

    IF (APPEND) THEN
      OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN',POSITION='APPEND')
   ELSE
      OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN')
   END IF

    NSAVE = 1+NSTEPS/SNAPEVERY
    NINFO = SIZE(INFO)
    WRITE (X1,'(I2)') NINFO
    FMT1 = '(4I8,'//TRIM(X1)//'ES18.4E2)'
    WRITE(99,FMT1) NTRIALS,NSAVE,NPART,NMT,INFO

    ! DO MC = 1,PARTP%NMT
    !   WRITE(99,'(I3,I3,3ES18.4E2)') MC, PARTP%MTSTATE(MC), PARTP%MTPOS(:,MC)
    ! END DO

    ! DO PC = 1,PARTP%NPART
    !   WRITE(99,'(I8,3ES18.4E2)') PC, PARTP%PARTPOS(:,PC)
    ! END DO

    CLOSE(99)

  END SUBROUTINE OUTPUTSNAPSHOT

  SUBROUTINE MAKEPARTICLE(PARTP,PARTID,PTYPE,NTYPE,NPROT,PRAD,VEL,D,KREV,INITPOS,&
                          INITSTATE,CANFUSE,KFUSE,PFUSE,INITPROT,PROTCONV,CURREG)
  ! make a new particle and assign appropriate parameters
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
      STOP 1
    END IF

    PARTP%ID = PARTID
    PARTP%PTYPE = PTYPE
    PARTP%PRAD = PRAD
    
    PARTP%POS = INITPOS
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

    PARTP%CANFUSE = .FALSE.
    PARTP%NFUSE = 0
    DEALLOCATE(PARTP%KFUSE,PARTP%PFUSE)
    DEALLOCATE(PARTP%PROT,PARTP%PROTCONV)
    PARTP%CURREG = 0
    
    PARTP%ARRAYSET = .FALSE.
  END SUBROUTINE DELETEPARTICLE

  SUBROUTINE SETUPPARTICLEGROUP(PGROUP,NPART,NTYPE,PLIST,KPROD,KSWITCH,INITPARTREG)
  ! set up a cluster of particles
    USE KEYS, ONLY: MAXNPART
    IMPLICIT NONE
    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    TYPE(PARTICLE), TARGET :: PLIST(NPART)
    INTEGER, INTENT(IN) :: NPART,NTYPE
    DOUBLE PRECISION, INTENT(IN) :: KPROD(NTYPE), KSWITCH(NTYPE)
    INTEGER, INTENT(IN), OPTIONAL :: INITPARTREG
    INTEGER :: PC

    IF(PGROUP%ARRAYSET) THEN
      PRINT*, 'ERROR IN SETUPPARTICLEGROUP: trying to overwrite'
      STOP 1
    END IF

    PGROUP%NPART = NPART
    ALLOCATE(PGROUP%PLIST(NPART))
    PGROUP%PLIST => PLIST
    PGROUP%MAXNPART = MAXNPART
    PGROUP%NTYPE = NTYPE
  
  END SUBROUTINE SETUPPARTICLEGROUP


END MODULE PARTICLEUTILS
