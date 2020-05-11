MODULE PARTDYNAMICS
! subroutines to simulate transport and dynamics of 1d particles
! VERSION BY: SAURABH MOGRE 
! DATE: MAY 2020

USE MT19937, ONLY: GRND,RNORM,GRND
  
IMPLICIT NONE

CONTAINS

  !INSERT PARTICLE DYNAMICS CODE HERE
  SUBROUTINE TRANSPORTSIM(PGROUP,DOMP,DOFUSION)
  !simulates the transport of phagosomes and lysosomes in a 1D domain
    USE KEYS, ONLY: PRINTEVERY,SNAPEVERY,OUTEVERY,NSTEPS,DELT,SNAPFILE,&
                    KPROD,PRAD,VEL,D,KREV,STARTPOS,STARTSTATE,OUTFILE,&
                    USEPFUSE,PFUSE,KFUSE
    USE PARTICLEUTILS, ONLY: PARTICLEGROUP,PARTICLE,DOMAIN,OUTPUTSNAPSHOT,&
                            MAKEPARTICLE
    USE FUSIONUTILS, ONLY: CHECKFUSION,FUSEPART                  
    USE MT19937, ONLY: RNORM,GRND
    
    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    TYPE(PARTICLE), POINTER :: PLIST(:)
    TYPE(PARTICLE), POINTER :: PARTP
    TYPE(DOMAIN), POINTER :: DOMP

    DOUBLE PRECISION, ALLOCATABLE :: INFO(:)
    INTEGER :: TC,PC,NC,MC,PID
    LOGICAL :: ADDPART,CHECKFUSE
    DOUBLE PRECISION :: SEP,PF

    LOGICAL, INTENT(IN) :: DOFUSION

    IF(PRINTEVERY.NE.0) THEN
      PRINT'(A5,I1,A4,I8)', 'Step ',0,' of ',NSTEPS
    END IF

    ! initial snapshot
    IF(SNAPEVERY.GT.0) THEN
      ALLOCATE(INFO(3))
      INFO = (/0D0,DELT,DOMP%DOMLEN/)
      CALL OUTPUTSNAPSHOT(PGROUP,SNAPFILE,INFO,.FALSE.)
    END IF

    DO TC = 1,NSTEPS
      IF((PRINTEVERY.NE.0).AND.MOD(TC,PRINTEVERY).EQ.0) THEN
      ! print information on simulation status
          PRINT'(A5,I8,A4,I8)', 'Step ',TC,' of ',NSTEPS
      END IF

      ! go through each particle
      PLIST=>PGROUP%PLIST
      DO PC = 1,PGROUP%NPART
        IF(PLIST(PC)%ID.EQ.0) CONTINUE
        ! move particle in the current direction
        PLIST(PC)%POS = MAX(MIN(PLIST(PC)%POS+PLIST(PC)%CURSTATE*PLIST(PC)%VEL&
                        *DELT,DOMP%DOMLEN),0D0)

        ! check for fusion
        IF(DOFUSION.AND.PLIST(PC)%CANFUSE) THEN
          DO MC = 1,PGROUP%NPART
            IF(PLIST(MC)%ID.EQ.0) CONTINUE
            CALL CHECKFUSION(CHECKFUSE,PLIST(PC)%POS,PLIST(PC)%POS,&
                            PLIST(PC)%PTYPE,PLIST(PC)%PTYPE,PLIST(PC)%PRAD,&
                            PLIST(PC)%PRAD)
            IF(CHECKFUSE) THEN
              CALL FUSEPART(PGROUP)
            END IF
          END DO
        END IF

        ! check if particle changes direction
        IF(GRND().LE.1-EXP(-PLIST(PC)%KREV*DELT)) THEN
          PLIST(PC)%CURSTATE = -PLIST(PC)%CURSTATE
        END IF

      END DO

      ! particle production
      IF(PGROUP%NPART.LT.PGROUP%MAXNPART) THEN
        DO NC = 1,PGROUP%NTYPE
          ADDPART = .FALSE.
          ! check if particle is produced
          IF(GRND().LE.(1-EXP(-KPROD(NC)*DELT))) THEN
            ! PRINT*, 'TESTX: ADDING PARTICLE'
            PID = PGROUP%LASTID+1
            ! check if particle can be added to existing list
            DO PC = 1,PGROUP%NPART
              IF(PLIST(PC)%ID.EQ.0) THEN
                PARTP=>PLIST(PC)
                CALL MAKEPARTICLE(PARTP,PID,NC,PGROUP%NTYPE,PGROUP%NPROT,&
                                PRAD(NC),VEL(NC),D(NC),KREV(NC),STARTPOS(NC),&
                                STARTSTATE(NC))
                !PLIST(PC) = PARTP
                ADDPART = .TRUE.
                EXIT
              END IF
            END DO
            ! append particle if it wasn't added
            IF(.NOT.ADDPART) THEN
               PARTP=>PLIST(PGROUP%NPART+1)
               CALL MAKEPARTICLE(PARTP,PID,NC,PGROUP%NTYPE,PGROUP%NPROT,&
                                PRAD(NC),VEL(NC),D(NC),KREV(NC),STARTPOS(NC),&
                                STARTSTATE(NC))
            END IF
            PGROUP%NPART = PGROUP%NPART+1
            PGROUP%LASTID = PID
          END IF
        END DO
      END IF

      ! output snapshot to file
      IF(MOD(TC,SNAPEVERY).EQ.0) THEN
        INFO(1) = TC
        CALL OUTPUTSNAPSHOT(PGROUP,SNAPFILE,INFO,.TRUE.)
      END IF

      ! output protein flux
      IF(MOD(TC,OUTEVERY).EQ.0) THEN
        CALL WRITEFLUX(PGROUP,OUTFILE,INFO,.TRUE.)
      END IF
    END DO

  END SUBROUTINE TRANSPORTSIM
 
END MODULE PARTDYNAMICS