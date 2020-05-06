MODULE PARTICLEUTILS
  ! utilities for defining parameters and states of individual particles, particle sets, and domain information

  TYPE PARTICLE
  ! object containing information about individual particles

    ! type of particle
    ! total number of particle types depends on system to be simulated
    INTEGER :: PTYPE

    ! current particle position
    DOUBLE PRECISION, POINTER :: POS

    ! particle state
    ! can take different values depending on system to be simulated
    ! default: +1 = moving forward; -1 = moving backward; 0 = stationary
    INTEGER :: CURSTATE

    ! walk velocity
    DOUBLE PRECISION :: VEL

    ! can the particle fuse?
    LOGICAL :: CANFUSE


    !MOVE KFUSE AND PFUSE TO MAIN PROGRAM BODY!
    ! rate of fusion
    DOUBLE PRECISION :: KFUSE

    ! probability of fusion per contact
    DOUBLE PRECISION :: PFUSE
    

    ! number of fusion events
    INTEGER :: NFUSE

    ! array of protein content
    ! contains amounts of different protein types
    DOUBLE PRECISION, POINTER :: PROT(:)

    ! conversion rates between protein types
    DOUBLE PRECISION, POINTER :: PROTCONV(:,:)

    ! Particle ID
    INTEGER :: ID

    ! starting position for the particle
    ! represents where the particle is produced
    DOUBLE PRECISION :: STARTPOS

    ! region in which particle resides
    INTEGER :: CURREG

    ! have arrays been allocated?
    LOGICAL :: ARRAYSET

  END TYPE PARTICLE

  TYPE PARTICLEGROUP
  ! object containing information about a group of mixed particle types

    ! number of particles in the group
    INTEGER :: NPART 

    ! maximum number of particles allowed
    INTEGER :: MAXNPART

    ! pointer containing particles within the group
    TYPE(PARTICLE), POINTER :: PLIST(:)

    ! number of particle types
    INTEGER :: NTYPE

    ! rates of production for each particle type
    DOUBLE PRECISION, POINTER :: KPROD(:)

    ! rates of direction reversal for each particle type
    DOUBLE PRECISION, POINTER :: KSWITCH(:)

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

  SUBROUTINE MAKEPARTICLE(PARTP,PTYPE,NPROT,VEL,INITPOS,INITPROT,INITSTATE,PARTID)
  ! make a new particle
    IMPLICIT NONE

    TYPE(PARTICLE), POINTER :: PARTP
    INTEGER, INTENT(IN) :: PTYPE,NPROT,PARTID,INITSTATE
    DOUBLE PRECISION, INTENT(IN) :: VEL,INITPOS,INITPROT(:)
    
    IF(PARTP%ARRAYSET) THEN
      PRINT*, 'ERROR IN MAKEPARTICLE: trying to overwrite particle'
      STOP 1
    END IF

    PARTP%PTYPE = PTYPE
    ALLOCATE(PARTP%PROT(NPROT))

    PARTP%ARRAYSET= .TRUE.

    PARTP%VEL = VEL
    PARTP%POS = INITPOS
    PARTP%CURSTATE = INITSTATE
    PARTP%ID = PARTID
    PARTP%PROT = INITPROT

  END SUBROUTINE MAKEPARTICLE

  SUBROUTINE DELETEPARTICLE(PARTP)
    ! deallocate arrays
    IMPLICIT NONE
    TYPE(PARTICLE), POINTER :: PARTP

    DEALLOCATE(PARTP%PROT,PARTP%PROTCONV)
    PARTP%PTYPE = 0
    PARTP%VEL = 0D0
    PARTP%ARRAYSET = .FALSE.
  END SUBROUTINE DELETEPARTICLE

END MODULE PARTICLEUTILS
