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
  INTEGER, PARAMETER :: MAXNPART = 1500 ! maximum number of particles allowed
  INTEGER, PARAMETER :: MAXNTYPE = 2 ! maximum number of particle types

  INTEGER :: NTYPE ! number of particle types
  INTEGER :: NPART(MAXNTYPE) ! number of active particles in a group
  DOUBLE PRECISION :: VEL(MAXNTYPE) ! velocity of each particle type
  DOUBLE PRECISION :: PRAD(MAXNTYPE) ! radius of each particle type
  DOUBLE PRECISION :: BUFLEN ! buffer length to determine contact

  DOUBLE PRECISION :: KREV(MAXNTYPE) ! direction reversal rates for particle types
  DOUBLE PRECISION :: KSTOP(MAXNTYPE) ! stopping rates for particles
  DOUBLE PRECISION :: DIFFCONST(MAXNTYPE) ! diffusion coefficient for diffusive particles
  INTEGER :: STARTSTATE(MAXNTYPE) !state in which particles start

  ! ----------------------------------------------------------------------------
  ! protein dynamics
  ! ----------------------------------------------------------------------------
  INTEGER, PARAMETER :: MAXNPROT = 2 ! maximum number of protein types
  INTEGER :: NPROT(MAXNTYPE) ! number of protein types in each particle type
  DOUBLE PRECISION :: PROTCONV(MAXNPROT,MAXNPROT,MAXNTYPE) ! conversion rate between proteins

  ! ----------------------------------------------------------------------------
  ! fusion parameters
  ! ----------------------------------------------------------------------------
  DOUBLE PRECISION :: KFUSE(MAXNTYPE,MAXNTYPE) !matrix for fusion rates between particle types
  DOUBLE PRECISION :: PFUSE(MAXNTYPE,MAXNTYPE) !matrix for fusion probabilities
  LOGICAL :: DOFUSION ! implement fusion?
  LOGICAL :: USEPFUSE ! use constant fusion probability?
  INTEGER :: MAXNFUSE(MAXNTYPE)

  ! ----------------------------------------------------------------------------
  ! particle production
  ! ----------------------------------------------------------------------------
  DOUBLE PRECISION :: KPROD(MAXNTYPE) ! rate of production
  DOUBLE PRECISION :: STARTPOS(MAXNTYPE) ! starting position
  DOUBLE PRECISION :: STARTPROT(MAXNTYPE,MAXNPROT) ! amount of protein to start with
 
  ! ----------------------------------------------------------------------------
  ! domain info
  ! ----------------------------------------------------------------------------
  DOUBLE PRECISION :: DOMLEN ! length of domain
  INTEGER :: NREG ! number of regions
  INTEGER, PARAMETER :: MAXNREG = 1500 ! maximum number of regions allowed
  DOUBLE PRECISION :: REGPOS(MAXNREG) ! positions of regions
  
  ! ----------------------------------------------------------------------------
  ! input/output
  ! ----------------------------------------------------------------------------
  CHARACTER*100 :: SNAPFILE, OUTFILE, NFUSEFILE,NTOUCHFILE
  LOGICAL :: DUMPSNAPSHOTS, APPENDSNAPSHOTS,CALCNTOUCH
  INTEGER :: SNAPEVERY, OUTEVERY, FUSEDISTEVERY
  
  ! ----------------------------------------------------------------------------
  ! dynamics
  ! ----------------------------------------------------------------------------
  DOUBLE PRECISION :: DELT
  INTEGER :: NSTEPS, PRINTEVERY
	
END MODULE KEYS