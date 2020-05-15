MODULE PARTDYNAMICS
! subroutines to simulate transport and dynamics of 1d particles
! VERSION BY: SAURABH MOGRE 
! DATE: MAY 2020

USE MT19937, ONLY: GRND,RNORM,GRND
  
IMPLICIT NONE

CONTAINS

  !INSERT PARTICLE DYNAMICS CODE HERE
  SUBROUTINE TRANSPORTSIM(PGROUP,DOMP,PARAMP,DOFUSION)
  !simulates the transport of phagosomes and lysosomes in a 1D domain
    USE KEYS, ONLY: PRINTEVERY,SNAPEVERY,OUTEVERY,NSTEPS,DELT,SNAPFILE,OUTFILE
    USE PARTICLEUTILS, ONLY: PARTICLEGROUP,PARTICLE,DOMAIN,OUTPUTSNAPSHOT,&
                            MAKEPARTICLE,PARTICLEPARAMS,DELETEPARTICLE,&
                            REMOVEPARTICLE
    USE FUSIONUTILS, ONLY: CHECKFUSION,FUSEPART                  
    USE MT19937, ONLY: RNORM,GRND
    
    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    TYPE(PARTICLE), POINTER :: PLIST(:)
    TYPE(PARTICLE), POINTER :: PARTP
    TYPE(PARTICLEPARAMS), POINTER :: PARAMP
    TYPE(DOMAIN), POINTER :: DOMP

    DOUBLE PRECISION, ALLOCATABLE :: INFO(:)
    INTEGER :: TC,PC,NC,MC,CC,PID,NPROC,TOTNPART
    LOGICAL :: ADDPART,CHECKFUSE
    DOUBLE PRECISION :: SEP,PF

    LOGICAL, INTENT(IN) :: DOFUSION

    IF(PRINTEVERY.NE.0) THEN
      PRINT'(A5,I1,A4,I8)', 'Step ',0,' of ',NSTEPS
    END IF

    ! initial snapshot
    IF(SNAPEVERY.GT.0) THEN
      ALLOCATE(INFO(2))
      INFO = (/0D0,DOMP%DOMLEN/)
      CALL OUTPUTSNAPSHOT(PGROUP,SNAPFILE,INFO,.FALSE.)
    END IF

    ALLOCATE(PLIST(PARAMP%MAXNPART))

    DO TC = 1,NSTEPS
      
      ! print information on simulation status
      IF((PRINTEVERY.NE.0).AND.MOD(TC,PRINTEVERY).EQ.0) THEN
        PRINT'(A5,I8,A4,I8)', 'Step ',TC,' of ',NSTEPS
      END IF

      ! move particles, change direction, and delete out of bounds
      PC = 1
      NPROC = 0
      TOTNPART = SUM(PGROUP%NPART)
      DO WHILE(NPROC.LT.TOTNPART.AND.PC.LT.PGROUP%MAXNPART &
              .AND.SUM(PGROUP%NPART).GT.0)
        PARTP=>PGROUP%PLIST(PC)
        
        ! move particles
        PARTP%POS = PARTP%POS+PARTP%CURSTATE*PARTP%VEL*DELT

        ! check if particle changes direction
        IF(GRND().LE.(1-EXP(-PARTP%KREV*DELT))) THEN
          PARTP%CURSTATE = -PARTP%CURSTATE
        END IF

        ! delete particles that move out of bounds
        IF((PARTP%POS.GT.DOMP%DOMLEN).OR.(PARTP%POS.LT.0D0)) THEN
          CALL REMOVEPARTICLE(PGROUP,PC)
        ELSE
          PC = PC+1
        END IF
        
        NPROC = NPROC + 1
      END DO

      ! particle production      
      DO NC = 1,PARAMP%NTYPE ! go through each type of particle
        ! check if particle is produced
        IF(GRND().LE.(1-EXP(-PARAMP%KPROD(NC)*DELT))) THEN
          ! check if you have space for a new particle
          IF(SUM(PGROUP%NPART).GE.PGROUP%MAXNPART) THEN
            PRINT*, 'ERROR: too many particles added.', PGROUP%NPART, PGROUP%MAXNPART
            STOP 1
          END IF

          PARTP=>PGROUP%PLIST(SUM(PGROUP%NPART)+1)
          PID = PGROUP%LASTID+1
          CALL MAKEPARTICLE(PARTP,PID,NC,PARAMP)
          PGROUP%NPART(NC)  = PGROUP%NPART(NC) + 1
          PGROUP%LASTID = PID                       
        END IF
      END DO

      ! output snapshot to file
      IF(MOD(TC,SNAPEVERY).EQ.0) THEN
        INFO(1) = TC*DELT
        CALL OUTPUTSNAPSHOT(PGROUP,SNAPFILE,INFO,.TRUE.)
      END IF

    END DO

    PRINT*, 'DYNAMICS DONE'

  END SUBROUTINE TRANSPORTSIM
 
END MODULE PARTDYNAMICS
