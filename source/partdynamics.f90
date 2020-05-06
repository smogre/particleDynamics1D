MODULE PARTDYNAMICS
! subroutines to simulate transport and dynamics of 1d particles
! VERSION BY: SAURABH MOGRE 
! DATE: MAY 2020

USE MT19937, ONLY: GRND,RNORM,GRND
  
IMPLICIT NONE

CONTAINS

  SUBROUTINE TRANSPORTSIM(PARTP,APPENDSNAPSHOTS)

	!------------------------------------------------------
  !simulates particles undergoing multi-modal transport in a 3d cylinder          
  !particles move on microtubules distributed uniformly within the cross section
  !length of microtubules is uniformly distributed                    
	!NPART - NUMBER OF PARTICLES
	!NMT - NUMBER OF MICROTUBULES
	!NSTEP - NUMBER OF STEPS
	!DELT - TIMESTEP
  !V - PARTICLE VELOCITY ON A MICROTUBULE
  !KR - RATE OF REVERSAL OF A MOVING PARTICLE
	!KJ - RATE OF JUMPING TO ANOTHER MICROTUBULE WITHIN REACH
	!CRAD - RADIUS WITHIN WHICH A MICROTUBULE CAN CAPTURE PARTICLES
  !DOMRAD - RADIUS OF THE DOMAIN
  !DOMLEN - LENGTH OF THE DOMAIN
  !PRINTEVERY - HOW OFTEN STATUS IS DISPLAYED WHILE RUNNING SIMS
  !PARTFILE - FILE TO WHICH PARTICLE POSITIONS ARE WRITTEN
  !SNAPSHOTEVERY - FREQUENCY AT WHICH OUTPUTS ARE WRITTEN TO FILE
  !------------------------------------------------------
    USE KEYS, ONLY : NSTEPS,NPART,DELT,DOMRAD,DOMLEN,KBIND,KREV,KJUMP,CRAD,&
                    MTDYN,VGMT,VSMT,KRMT,KCMT,&
                    PRINTEVERY,SNAPFILE,SNAPEVERY
    USE PARTICLEUTIL, ONLY: PARTICLEGROUP,OUTPUTSNAPSHOT
    USE MT19937, ONLY: RNORM,GRND

    TYPE(PARTICLEGROUP), POINTER :: PARTP
    LOGICAL, INTENT(IN) :: APPENDSNAPSHOTS
  	
    DOUBLE PRECISION,ALLOCATABLE :: INFO(:)
    INTEGER :: TC,PC,MC,NJUMP
    
    INTEGER :: MTJUMPS(PARTP%NMT),NEWIND, DIRVEC(2)
    INTEGER,ALLOCATABLE :: JUMPINDS(:)
    DOUBLE PRECISION :: VMT(3), NEWPOS(3)
    
    IF(PRINTEVERY.NE.0) THEN
      PRINT'(A5,I1,A4,I8)', 'Step ',0,' of ',NSTEPS
    END IF
    
    ! directions available
    DIRVEC = (/-1,1/)
    ! INITIAL SNAPSHOT
    ALLOCATE(INFO(6))
    INFO = (/0D0,DFLOAT(NSTEPS),DELT,DOMRAD,DOMLEN,CRAD/)

    ! array of microtubule velocities
    VMT = (/VGMT,VSMT,0D0/)

    CALL OUTPUTSNAPSHOT(PARTP,SNAPFILE,INFO,APPENDSNAPSHOTS)

    DO TC = 1,NSTEPS
      IF((PRINTEVERY.NE.0).AND.MOD(TC,PRINTEVERY).EQ.0) THEN
      ! print information on particle status
        PRINT'(A5,I8,A4,I8)', 'Step ',TC,' of ',NSTEPS
      END IF
      
      ! microtubule dynamics 
      IF(MTDYN) THEN
        DO MC = 1,PARTP%NMT
          SELECT CASE(PARTP%MTSTATE(MC))
          CASE(1)
            IF(GRND().LE.(1-EXP(-KCMT*DELT))) THEN
              PARTP%MTSTATE(MC) = 2
            END IF
          CASE(2)
            IF(GRND().LE.(1-EXP(-KRMT*DELT))) THEN
              PARTP%MTSTATE(MC) = 1
            END IF
          CASE DEFAULT
          END SELECT
          PARTP%MTPOS(3,MC) = MAX(MIN(PARTP%MTPOS(3,MC)+&
                              VMT(PARTP%MTSTATE(MC))*DELT,DOMLEN),0D0)
        END DO
      END IF

      ! particle motion
      DO PC = 1,NPART
        ! check if particle is on a microtubule
        IF(PARTP%MTINDS(PC).NE.0) THEN
          IF(PARTP%PARTPOS(3,PC).GT.PARTP%MTPOS(3,PARTP%MTINDS(PC))) THEN
            PARTP%MTINDS(PC) = 0
          ELSE
            NJUMP = 0
            IF(PARTP%NMT.GT.1) THEN
              MTJUMPS = 0
              ! check which jumps are possible
              DO MC = 1,PARTP%NMT
                IF((((PARTP%PARTPOS(1,PC)-PARTP%MTPOS(1,MC))**2+&
                  (PARTP%PARTPOS(2,PC)-PARTP%MTPOS(2,MC))**2).LE.CRAD**2).AND.&
                            (PARTP%PARTPOS(3,PC).LE.PARTP%MTPOS(3,MC)).AND.&
                                                (PARTP%PARTPOS(3,PC).GE.0D0)) THEN
                  MTJUMPS(MC) = MC
                  NJUMP = NJUMP+1
                END IF
              END DO
            END IF

            ! check if event occurs
            IF(GRND().LE.(1-EXP(-(KREV+NJUMP*KJUMP)*DELT))) THEN
              ! determine whether event is reversal or jump
              IF(GRND().LE.(KREV/(KREV+NJUMP*KJUMP))) THEN
                !reverse direction
                PARTP%CURDIR(PC) = -PARTP%CURDIR(PC)
              ELSE
                !choose microtubule to jump to
                !update xy position and mtinds
                ALLOCATE(JUMPINDS(NJUMP))
                JUMPINDS = PACK(MTJUMPS,MTJUMPS.GT.0)
                NEWIND = JUMPINDS(FLOOR(1+GRND()*NJUMP))
                PARTP%PARTPOS(1,PC) = PARTP%MTPOS(1,NEWIND)
                PARTP%PARTPOS(2,PC) = PARTP%MTPOS(2,NEWIND)
                PARTP%MTINDS(PC) = NEWIND
                PARTP%CURDIR(PC) = DIRVEC(FLOOR(1+2*GRND()))
                DEALLOCATE(JUMPINDS)
              END IF
            END IF

            ! move particle
            PARTP%PARTPOS(3,PC) = MAX(MIN(PARTP%PARTPOS(3,PC)+PARTP%CURDIR(PC)*PARTP%VEL*DELT,&
                                  PARTP%MTPOS(3,PARTP%MTINDS(PC))),0D0)
          END IF
        END IF
        
        IF(PARTP%MTINDS(PC).EQ.0) THEN
          ! execute Brownian step
          CALL BDPROPAGATE(PARTP%PARTPOS(:,PC),PARTP%D,DOMRAD,DOMLEN,DELT,NEWPOS)
          PARTP%PARTPOS(:,PC) = NEWPOS

          NJUMP = 0
          MTJUMPS = 0
          
          ! check which microtubules the particle can bind to
          DO MC = 1,PARTP%NMT
            IF((((PARTP%PARTPOS(1,PC)-PARTP%MTPOS(1,MC))**2+&
                (PARTP%PARTPOS(2,PC)-PARTP%MTPOS(2,MC))**2).LE.CRAD**2)) THEN
              IF((PARTP%PARTPOS(3,PC).GE.0D0).AND.&
              (PARTP%PARTPOS(3,PC).LE.(PARTP%MTPOS(3,MC)+CRAD))) THEN
                MTJUMPS(MC) = MC
                NJUMP = NJUMP + 1
              END IF
            END IF
          END DO

          ! check if binding event occurs
          IF(GRND().LE.(1-EXP(-KBIND*DELT)).AND.NJUMP.GT.0) THEN
            ALLOCATE(JUMPINDS(NJUMP))
            JUMPINDS = PACK(MTJUMPS,MTJUMPS.GT.0)
            NEWIND = JUMPINDS(FLOOR(1+GRND()*NJUMP))
            PARTP%PARTPOS(1,PC) = PARTP%MTPOS(1,NEWIND)
            PARTP%PARTPOS(2,PC) = PARTP%MTPOS(2,NEWIND)
            PARTP%MTINDS(PC) = NEWIND
            PARTP%CURDIR(PC) = DIRVEC(FLOOR(1+2*GRND()))
            DEALLOCATE(JUMPINDS)
          END IF

        END IF

      END DO



      IF(MOD(TC,SNAPEVERY).EQ.0) THEN
      ! output snapshot to file
        INFO(1) = TC
        CALL OUTPUTSNAPSHOT(PARTP,SNAPFILE,INFO,.TRUE.)
      END IF
    END DO

  END SUBROUTINE TRANSPORTSIM

  SUBROUTINE BDPROPAGATE(PREVPOS,D,DOMRAD,DOMLEN,DELT,NEWPOS)
    ! do a brownian dynamics step to propagate the particle over a time interval
    ! R0, Z0 are initial radial and axial positions
    ! DDT = D*dt (diffusion coefficient times time interval)
    ! RNEW, ZNEQ are final radial and axial positions
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: PREVPOS(3)
    DOUBLE PRECISION, INTENT(IN) :: D,DOMRAD,DOMLEN,DELT
    DOUBLE PRECISION, INTENT(OUT) :: NEWPOS(3)
    DOUBLE PRECISION :: STEPXYZ(3), COEF
    INTEGER :: I
    
    COEF = SQRT(2*D*DELT);
 
    DO I = 1,3
       STEPXYZ(I) = COEF*RNORM()
    END DO
    
    NEWPOS = PREVPOS+STEPXYZ
    
    IF((NEWPOS(3).GT.DOMLEN).OR.(NEWPOS(3).LT.0D0)&
                                .OR.(NEWPOS(1)**2+NEWPOS(2)**2).GT.DOMRAD**2) THEN
      CALL REFLECTPART(PREVPOS,DOMRAD,DOMLEN,NEWPOS)                            
    END IF

    
  END SUBROUTINE BDPROPAGATE

  SUBROUTINE REFLECTPART(POS0,DOMRAD,DOMLEN,POS1)
    !SUBROUTINE TO REFLECT PARTICLES OFF THE BOUNDARY
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: POS0(3),DOMRAD,DOMLEN
    DOUBLE PRECISION, INTENT(OUT) :: POS1(3)
    DOUBLE PRECISION :: X1,Y1,X2,Y2,X0,Y0,A,B,C,C1,M,DISC
      
    DO WHILE ((POS1(3).GT.DOMLEN).OR.(POS1(3).LT.0D0))
      IF(POS1(3).GT.DOMLEN) POS1(3) = 2*DOMLEN-POS1(3)
      IF(POS1(3).LT.0D0) POS1(3) = -POS1(3)
    END DO
    
        
    DO WHILE (SQRT((POS1(1)**2+POS1(2)**2)).GT.DOMRAD)
      
      !FIND CO-ORDINATES OF THE INTERSECTION POINT
      X1 = POS0(1)
      Y1 = POS0(2)
      X2 = POS1(1)
      Y2 = POS1(2)
      M = (Y2-Y1)/(X2-X1)
      C1 = Y1-M*X1
        
      A = 1+M**2
      B = 2*M*C1
      C = C1**2-DOMRAD**2
      
      IF(B**2-4*A*C .LT.0) THEN
        PRINT*, 'TESTXA:', B**2-4*A*C 
        PRINT*, 'TESTXA:',POS0,POS1
        PRINT*, 'TESTXA:',SQRT((POS0(1)**2+POS0(2)**2)),SQRT((POS1(1)**2+POS1(2)**2))
        STOP 1
      END IF
      DISC = SQRT(B**2-4*A*C)
      
        
      IF(((-B-DISC)/(2*A).GT.MIN(X1,X2)).AND.((-B-DISC)/(2*A).LT.MAX(X1,X2))) THEN
        X0 = (-B-DISC)/(2*A)
      ELSE
        X0 = (-B+DISC)/(2*A)
      END IF
        
      Y0 = M*X0+C1
        
      !GET REFLECTION OFF THE TANGENT AT INTERSECTION POINT
      POS1(1) = X2-2*X0*(X0*X2+Y0*Y2-DOMRAD**2)/DOMRAD**2
      POS1(2) = Y2-2*Y0*(X0*X2+Y0*Y2-DOMRAD**2)/DOMRAD**2
        
    END DO
  END SUBROUTINE REFLECTPART
 
END MODULE PARTDYNAMICS