MODULE PARTDYNAMICS
! subroutines to simulate transport and dynamics of 1d particles
! VERSION BY: SAURABH MOGRE 
! DATE: MAY 2020

USE MT19937, ONLY: GRND,RNORM,GRND
  
IMPLICIT NONE

CONTAINS

  !INSERT PARTICLE DYNAMICS CODE HERE
  SUBROUTINE TRANSPORTSIM(PGROUP,DOMP,PARAMP,DOFUSION,APPENDTOFILE)
  !simulates the transport of phagosomes and lysosomes in a 1D domain
    USE KEYS, ONLY: PRINTEVERY,SNAPEVERY,OUTEVERY,NSTEPS,DELT,SNAPFILE,OUTFILE,&
                    MAXNPROT
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
    
    LOGICAL, INTENT(IN) :: DOFUSION,APPENDTOFILE

    DOUBLE PRECISION, ALLOCATABLE :: INFO(:)
    INTEGER :: TC,PC,NC,MC,CC,PID,NPROC,TOTNPART
    LOGICAL :: ADDPART,CHECKFUSE
    DOUBLE PRECISION :: SEP,PF,PROTOUT(MAXNPROT)
    CHARACTER(LEN=2) :: X1
    DOUBLE PRECISION,POINTER :: FUSELIST(:)

    IF(PRINTEVERY.NE.0) THEN
      PRINT'(A5,I1,A4,I8)', 'Step ',0,' of ',NSTEPS
    END IF

    ! initial snapshot
    IF(SNAPEVERY.GT.0) THEN
      ALLOCATE(INFO(2))
      INFO = (/0D0,DOMP%DOMLEN/)
      CALL OUTPUTSNAPSHOT(PGROUP,SNAPFILE,INFO,APPENDTOFILE)
    END IF

    ! initial output
    IF(OUTEVERY.GT.0) THEN
      IF(APPENDTOFILE) THEN
        OPEN(UNIT=44,FILE=OUTFILE,STATUS='UNKNOWN',POSITION='APPEND')
      ELSE
        OPEN(UNIT=44,FILE=OUTFILE,STATUS='UNKNOWN')
      END IF
      WRITE(44,'(I8)') MAXNPROT
      CLOSE(44)
    END IF

    ALLOCATE(PLIST(PARAMP%MAXNPART))
    PROTOUT = 0D0

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

        ! delete particles that move out of bounds
        IF((PARTP%POS.GT.DOMP%DOMLEN).OR.(PARTP%POS.LT.0D0)) THEN
          IF(PARTP%POS.LT.0D0) THEN
            PROTOUT(1:PARTP%NPROT) = PROTOUT(1:PARTP%NPROT)+PARTP%PROT
          END IF
          CALL REMOVEPARTICLE(PGROUP,PC)
          NPROC = NPROC + 1
          CONTINUE
        END IF

        ! check if particle changes direction
        IF(GRND().LE.(1-EXP(-PARTP%KREV*DELT))) THEN
          PARTP%CURSTATE = -PARTP%CURSTATE
        END IF

        ! check for fusion

        ! IF(DOFUSION) THEN
        !   ! CALL CHECKFUSION(PGROUP,PARAMP,FUSELIST)
        !   ! CALL FUSEPART(PGROUP,PARAMP,FUSELIST)
        !   NPROC = NPROC+1
        ! END IF

        ! protein conversion
        ! ---------------------------------------
        ! change to matrix exponentiation
        ! ---------------------------------------
        PARTP%PROT = PARTP%PROT+MATMUL(PARTP%PROTCONV,PARTP%PROT)*DELT

        NPROC = NPROC+1
        PC = PC+1
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

      ! output protein flux to file
      IF(MOD(TC,OUTEVERY).EQ.0) THEN
        OPEN(UNIT=44,FILE=OUTFILE,STATUS='UNKNOWN',POSITION='APPEND')
        WRITE(X1,'(I2)') MAXNPROT+1
        WRITE(44,'('//TRIM(X1)//'ES15.3E2)') TC*DELT,PROTOUT
        CLOSE(44)
      END IF

    END DO

    PRINT*, 'DYNAMICS DONE'

  END SUBROUTINE TRANSPORTSIM
 
END MODULE PARTDYNAMICS
