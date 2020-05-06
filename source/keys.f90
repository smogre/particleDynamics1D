MODULE KEYS
  ! keyword parameters that are globally used in many different places in the code
  IMPLICIT NONE

  ! ----------------------------------------------------------------------------
  ! general program control
  ! ----------------------------------------------------------------------------
  CHARACTER*100 :: ACTION
  INTEGER :: RNGSEED
  LOGICAL :: VERBOSE
  
  ! ----------------------------------------------------------------------------
  ! simulation info
  ! ----------------------------------------------------------------------------
  INTEGER:: NTRIALS

  ! ----------------------------------------------------------------------------
  ! particle properties
  ! ----------------------------------------------------------------------------
  INTEGER :: NPART
  DOUBLE PRECISION :: VEL
  DOUBLE PRECISION :: KREV
  DOUBLE PRECISION :: KJUMP
  DOUBLE PRECISION :: CRAD
  DOUBLE PRECISION :: D
  DOUBLE PRECISION :: KBIND

  ! ----------------------------------------------------------------------------
  ! microtubule properties
  ! ----------------------------------------------------------------------------
  INTEGER :: NMT
  LOGICAL :: MTDYN
  DOUBLE PRECISION :: VGMT,VSMT !growth and shrink velocity
  DOUBLE PRECISION :: KRMT,KCMT !rescue and catastrophe rate
  
  ! ----------------------------------------------------------------------------
  ! domain info
  ! ----------------------------------------------------------------------------
  DOUBLE PRECISION :: DOMRAD
  DOUBLE PRECISION :: DOMLEN
  
  ! ----------------------------------------------------------------------------
  ! input/output
  ! ----------------------------------------------------------------------------
  CHARACTER*100 :: SNAPFILE
  LOGICAL :: DUMPSNAPSHOTS, APPENDSNAPSHOTS
  INTEGER :: SNAPEVERY
  
  ! ----------------------------------------------------------------------------
  ! dynamics
  ! ----------------------------------------------------------------------------
  DOUBLE PRECISION :: DELT
  INTEGER :: NSTEPS, PRINTEVERY
	
END MODULE KEYS