!======MODULE CODE FOR ABfluid======

module m

    real(8), parameter :: kB = 1.0
    real(8), parameter :: FPE = 1.0
    real(8), parameter :: ECharge = 1.0
    real(8), parameter :: eV = 1.0
    integer, parameter :: Dim = 3
    integer, dimension(0:4-1), parameter :: BondOrdData = (/ 1, 0, 0, 1 /)
    integer, dimension(0:3-1), parameter :: BondOrdStart = (/ 0, 2, 4 /)
    integer, dimension(0:2-1), parameter :: BondOrdShift = (/ 0, 0 /)
    integer, parameter :: BondOrdLimit = 4
    integer, parameter :: NAID = 2
    integer, parameter :: NSID = 2
    integer, parameter :: NMID = 2
    integer, parameter :: NAtom = 6750
    integer, parameter :: NMol = 6750
    integer, dimension(0:2-1), parameter :: NDOFMID = (/ 3, 3 /)
    integer, dimension(0:2-1), parameter :: AtomsPerMol = (/ 1, 1 /)
    integer, parameter :: MaxAtomsPerMol = 1
    integer, parameter :: MaxRigidAtoms = 0
    logical, dimension(0:6750-1), parameter :: MolIsRigid = (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 &
      & , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
    real(8), dimension(0:1-1), parameter :: RBondLengthSq = (/ -1.0 /)
    integer, dimension(0:1-1, 0:2-1), parameter :: RBondInd = reshape( (/ -1,-1 /) , (/1, 2/) )
    integer, dimension(0:3-1), parameter :: RBondRange = (/ 0, 0, 0 /)
    integer, parameter :: NTerm = 4
    integer, parameter :: NDParam = 19
    integer, parameter :: NDDParam = 117
    integer, parameter :: P0_ArgHistNBin = 10000
    integer, parameter :: P1_ArgHistNBin = 10000
    integer, parameter :: P2_ArgHistNBin = 10000
    integer, parameter :: P3_ArgHistNBin = 10000
    integer, parameter :: M0_NVal = 1
    integer, parameter :: M1_NVal = 1
    integer, parameter :: M2_NVal = 1
    integer, parameter :: M3_NVal = 1
    integer, parameter :: M4_NVal = 1
    integer, parameter :: M5_NVal = 1
    integer, parameter :: M6_NVal = 1
    integer, parameter :: M7_NVal = 1
    integer, parameter :: M8_NVal = 1
    integer, parameter :: M9_NVal = 1
    integer, parameter :: M10_NVal = 1
    integer, parameter :: M11_NVal = 1
    integer, parameter :: M12_NVal = 1
    integer, parameter :: M13_NVal = 19
    integer, parameter :: M14_NVal = 117
    integer, parameter :: M15_NVal = 19
    integer, parameter :: M16_NVal = 117
    integer, parameter :: M17_NVal = 361
    integer, parameter :: M18_NVal = 1
    integer, parameter :: M19_NVal = 1
    integer, parameter :: M20_NVal = 1
    integer, parameter :: NMeasure = 21
    integer, parameter :: VV_NH_N = 2
    integer, parameter :: NMCMoves = 4
    real(8), dimension(0:3-1) :: BoxL
    integer, dimension(0:2-1) :: AIDCount
    integer, dimension(0:6751-1) :: MolRange
    integer, dimension(0:6750-1) :: MolID
    integer :: NDOF
    real(8), dimension(0:2-1, 0:1-1, 0:3-1) :: COMPos
    real(8), dimension(0:6750-1, 0:3-1) :: Pos
    real(8), dimension(0:6750-1, 0:3-1) :: Vel
    real(8), dimension(0:6750-1, 0:3-1) :: Force
    real(8), dimension(0:6750-1) :: Mass
    real(8), dimension(0:6750-1) :: iMass
    real(8), dimension(0:6750-1) :: sqrtMass
    integer, dimension(0:6750-1) :: MInd
    integer :: NActiveMol
    integer :: NInactiveMol
    integer, dimension(0:2-1) :: NActiveMID
    integer, dimension(0:2-1) :: NInactiveMID
    integer, dimension(0:6750-1) :: AID
    integer, dimension(0:6750-1) :: SID
    integer, dimension(0:6750-1) :: MID
    integer, dimension(0:6750-1) :: MolActive
    real(8) :: PEnergy
    real(8) :: KEnergy
    real(8) :: TEnergy
    real(8) :: Virial
    real(8) :: TempSet
    real(8) :: PresSet
    integer, dimension(0:2-1) :: MuSet
    integer :: TargetAtom
    integer :: TargetMol
    real(8), dimension(0:4-1) :: Terms
    real(8), dimension(0:4-1) :: Cut
    real(8), dimension(0:4-1) :: CutSq
    real(8) :: OldPEnergy
    real(8), dimension(0:4-1) :: OldTerms
    real(8) :: FluctE
    real(8) :: FluctE0
    real(8) :: FluctA0
    real(8) :: FluctBeta
    real(8) :: FluctTerm
    real(8), dimension(0:19-1) :: Param
    real(8), dimension(0:19-1) :: DUParam
    real(8), dimension(0:19-1) :: DWParam
    real(8), dimension(0:117-1) :: DDUParam
    real(8), dimension(0:117-1) :: DDWParam
    real(8), dimension(0:13-1) :: P0_UShift
    real(8) :: P0_ArgWeightSumHist
    real(8) :: P0_ArgWeightSumStats
    real(8), dimension(0:1-1) :: P0_ArgMin
    real(8), dimension(0:1-1) :: P0_ArgMax
    real(8), dimension(0:1-1) :: P0_ArgCount
    real(8), dimension(0:1-1) :: P0_ArgSum
    real(8), dimension(0:1-1) :: P0_ArgSumSq
    real(8), dimension(0:1-1) :: P0_ArgReportMin
    real(8), dimension(0:1-1) :: P0_ArgReportMax
    real(8), dimension(0:1-1) :: P0_ArgReportDelta
    real(8), dimension(0:1-1) :: P0_ArgHistMin
    real(8), dimension(0:1-1) :: P0_ArgHistMax
    real(8), dimension(0:1-1) :: P0_ArgHistBinw
    real(8), dimension(0:1-1) :: P0_ArgHistiBinw
    real(8), dimension(0:1-1, 0:10000-1) :: P0_ArgHist
    real(8), dimension(0:13-1) :: P1_UShift
    real(8) :: P1_ArgWeightSumHist
    real(8) :: P1_ArgWeightSumStats
    real(8), dimension(0:1-1) :: P1_ArgMin
    real(8), dimension(0:1-1) :: P1_ArgMax
    real(8), dimension(0:1-1) :: P1_ArgCount
    real(8), dimension(0:1-1) :: P1_ArgSum
    real(8), dimension(0:1-1) :: P1_ArgSumSq
    real(8), dimension(0:1-1) :: P1_ArgReportMin
    real(8), dimension(0:1-1) :: P1_ArgReportMax
    real(8), dimension(0:1-1) :: P1_ArgReportDelta
    real(8), dimension(0:1-1) :: P1_ArgHistMin
    real(8), dimension(0:1-1) :: P1_ArgHistMax
    real(8), dimension(0:1-1) :: P1_ArgHistBinw
    real(8), dimension(0:1-1) :: P1_ArgHistiBinw
    real(8), dimension(0:1-1, 0:10000-1) :: P1_ArgHist
    real(8), dimension(0:13-1) :: P2_UShift
    real(8) :: P2_ArgWeightSumHist
    real(8) :: P2_ArgWeightSumStats
    real(8), dimension(0:1-1) :: P2_ArgMin
    real(8), dimension(0:1-1) :: P2_ArgMax
    real(8), dimension(0:1-1) :: P2_ArgCount
    real(8), dimension(0:1-1) :: P2_ArgSum
    real(8), dimension(0:1-1) :: P2_ArgSumSq
    real(8), dimension(0:1-1) :: P2_ArgReportMin
    real(8), dimension(0:1-1) :: P2_ArgReportMax
    real(8), dimension(0:1-1) :: P2_ArgReportDelta
    real(8), dimension(0:1-1) :: P2_ArgHistMin
    real(8), dimension(0:1-1) :: P2_ArgHistMax
    real(8), dimension(0:1-1) :: P2_ArgHistBinw
    real(8), dimension(0:1-1) :: P2_ArgHistiBinw
    real(8), dimension(0:1-1, 0:10000-1) :: P2_ArgHist
    real(8) :: P3_PlaneLoc
    integer :: P3_PlaneAxis
    real(8) :: P3_ArgWeightSumHist
    real(8) :: P3_ArgWeightSumStats
    real(8), dimension(0:1-1) :: P3_ArgMin
    real(8), dimension(0:1-1) :: P3_ArgMax
    real(8), dimension(0:1-1) :: P3_ArgCount
    real(8), dimension(0:1-1) :: P3_ArgSum
    real(8), dimension(0:1-1) :: P3_ArgSumSq
    real(8), dimension(0:1-1) :: P3_ArgReportMin
    real(8), dimension(0:1-1) :: P3_ArgReportMax
    real(8), dimension(0:1-1) :: P3_ArgReportDelta
    real(8), dimension(0:1-1) :: P3_ArgHistMin
    real(8), dimension(0:1-1) :: P3_ArgHistMax
    real(8), dimension(0:1-1) :: P3_ArgHistBinw
    real(8), dimension(0:1-1) :: P3_ArgHistiBinw
    real(8), dimension(0:1-1, 0:10000-1) :: P3_ArgHist
    integer :: M0_StepFreq
    integer :: M0_CycleFreq
    logical :: M0_Active
    real(8), dimension(0:1-1) :: M0_Val
    real(8), dimension(0:1-1) :: M0_ValSum
    real(8), dimension(0:1-1) :: M0_ValSumSq
    real(8) :: M0_Count
    integer :: M1_StepFreq
    integer :: M1_CycleFreq
    logical :: M1_Active
    real(8), dimension(0:1-1) :: M1_Val
    real(8), dimension(0:1-1) :: M1_ValSum
    real(8), dimension(0:1-1) :: M1_ValSumSq
    real(8) :: M1_Count
    integer :: M2_StepFreq
    integer :: M2_CycleFreq
    logical :: M2_Active
    real(8), dimension(0:1-1) :: M2_Val
    real(8), dimension(0:1-1) :: M2_ValSum
    real(8), dimension(0:1-1) :: M2_ValSumSq
    real(8) :: M2_Count
    integer :: M3_StepFreq
    integer :: M3_CycleFreq
    logical :: M3_Active
    real(8), dimension(0:1-1) :: M3_Val
    real(8), dimension(0:1-1) :: M3_ValSum
    real(8), dimension(0:1-1) :: M3_ValSumSq
    real(8) :: M3_Count
    integer :: M4_StepFreq
    integer :: M4_CycleFreq
    logical :: M4_Active
    real(8), dimension(0:1-1) :: M4_Val
    real(8), dimension(0:1-1) :: M4_ValSum
    real(8), dimension(0:1-1) :: M4_ValSumSq
    real(8) :: M4_Count
    integer :: M5_StepFreq
    integer :: M5_CycleFreq
    logical :: M5_Active
    real(8), dimension(0:1-1) :: M5_Val
    real(8), dimension(0:1-1) :: M5_ValSum
    real(8), dimension(0:1-1) :: M5_ValSumSq
    real(8) :: M5_Count
    integer :: M6_StepFreq
    integer :: M6_CycleFreq
    logical :: M6_Active
    real(8), dimension(0:1-1) :: M6_Val
    real(8), dimension(0:1-1) :: M6_ValSum
    real(8), dimension(0:1-1) :: M6_ValSumSq
    real(8) :: M6_Count
    integer :: M7_StepFreq
    integer :: M7_CycleFreq
    logical :: M7_Active
    real(8), dimension(0:1-1) :: M7_Val
    real(8), dimension(0:1-1) :: M7_ValSum
    real(8), dimension(0:1-1) :: M7_ValSumSq
    real(8) :: M7_Count
    integer :: M8_StepFreq
    integer :: M8_CycleFreq
    logical :: M8_Active
    real(8), dimension(0:1-1) :: M8_Val
    real(8), dimension(0:1-1) :: M8_ValSum
    real(8), dimension(0:1-1) :: M8_ValSumSq
    real(8) :: M8_Count
    integer :: M9_StepFreq
    integer :: M9_CycleFreq
    logical :: M9_Active
    real(8), dimension(0:1-1) :: M9_Val
    real(8), dimension(0:1-1) :: M9_ValSum
    real(8), dimension(0:1-1) :: M9_ValSumSq
    real(8) :: M9_Count
    integer :: M10_StepFreq
    integer :: M10_CycleFreq
    logical :: M10_Active
    real(8), dimension(0:1-1) :: M10_Val
    real(8), dimension(0:1-1) :: M10_ValSum
    real(8), dimension(0:1-1) :: M10_ValSumSq
    real(8) :: M10_Count
    integer :: M11_StepFreq
    integer :: M11_CycleFreq
    logical :: M11_Active
    real(8), dimension(0:1-1) :: M11_Val
    real(8), dimension(0:1-1) :: M11_ValSum
    real(8), dimension(0:1-1) :: M11_ValSumSq
    real(8) :: M11_Count
    integer :: M12_StepFreq
    integer :: M12_CycleFreq
    logical :: M12_Active
    real(8), dimension(0:1-1) :: M12_Val
    real(8), dimension(0:1-1) :: M12_ValSum
    real(8), dimension(0:1-1) :: M12_ValSumSq
    real(8) :: M12_Count
    integer :: M13_StepFreq
    integer :: M13_CycleFreq
    logical :: M13_Active
    real(8), dimension(0:19-1) :: M13_Val
    real(8), dimension(0:19-1) :: M13_ValSum
    real(8), dimension(0:19-1, 0:19-1) :: M13_ValSumSq
    real(8) :: M13_Count
    integer :: M14_StepFreq
    integer :: M14_CycleFreq
    logical :: M14_Active
    real(8), dimension(0:117-1) :: M14_Val
    real(8), dimension(0:117-1) :: M14_ValSum
    real(8), dimension(0:117-1) :: M14_ValSumSq
    real(8) :: M14_Count
    integer :: M15_StepFreq
    integer :: M15_CycleFreq
    logical :: M15_Active
    real(8), dimension(0:19-1) :: M15_Val
    real(8), dimension(0:19-1) :: M15_ValSum
    real(8), dimension(0:19-1) :: M15_ValSumSq
    real(8) :: M15_Count
    integer :: M16_StepFreq
    integer :: M16_CycleFreq
    logical :: M16_Active
    real(8), dimension(0:117-1) :: M16_Val
    real(8), dimension(0:117-1) :: M16_ValSum
    real(8), dimension(0:117-1) :: M16_ValSumSq
    real(8) :: M16_Count
    integer :: M17_StepFreq
    integer :: M17_CycleFreq
    logical :: M17_Active
    real(8), dimension(0:361-1) :: M17_Val
    real(8), dimension(0:361-1) :: M17_ValSum
    real(8), dimension(0:361-1) :: M17_ValSumSq
    real(8) :: M17_Count
    integer :: M18_StepFreq
    integer :: M18_CycleFreq
    logical :: M18_Active
    real(8), dimension(0:1-1) :: M18_Val
    real(8), dimension(0:1-1) :: M18_ValSum
    real(8), dimension(0:1-1) :: M18_ValSumSq
    real(8) :: M18_Count
    integer :: M19_StepFreq
    integer :: M19_CycleFreq
    logical :: M19_Active
    real(8), dimension(0:1-1) :: M19_Val
    real(8), dimension(0:1-1) :: M19_ValSum
    real(8), dimension(0:1-1) :: M19_ValSumSq
    real(8) :: M19_Count
    integer :: M20_StepFreq
    integer :: M20_CycleFreq
    logical :: M20_Active
    real(8), dimension(0:1-1) :: M20_Val
    real(8), dimension(0:1-1) :: M20_ValSum
    real(8), dimension(0:1-1) :: M20_ValSumSq
    real(8) :: M20_Count
    real(8) :: VVQ_TimeStep
    integer :: VVQ_RattleMaxIterPerAtom
    real(8) :: VVQ_RattleTol1
    real(8) :: VVQ_RattleTol2
    real(8) :: VVQ_RattleTol3
    real(8) :: VV_TimeStep
    integer :: VV_RattleMaxIterPerAtom
    real(8) :: VV_RattleTol1
    real(8) :: VV_RattleTol2
    real(8) :: VV_RattleTol3
    integer :: VV_Thermostat
    integer :: VV_AndersenStep
    integer :: VV_AndersenStepFreq
    real(8) :: VV_AndersenCollisionFreq
    real(8), dimension(0:2-1) :: VV_Glogs
    real(8), dimension(0:2-1) :: VV_Vlogs
    real(8), dimension(0:2-1) :: VV_Xlogs
    real(8), dimension(0:2-1) :: VV_QMass
    real(8) :: VV_LangevinGamma
    integer :: VV_Barostat
    integer :: VV_BarostatStepFreq
    real(8) :: VV_BarostatDeltaV
    logical, dimension(0:3-1) :: VV_BarostatUseAxis
    logical :: VV_BarostatDoIsotropic
    real(8) :: VV_BarostatMaxVol
    real(8) :: VV_BarostatMinVol
    real(8) :: VV_BarostatNAtt
    real(8) :: VV_BarostatNAcc
    integer :: VV_RemoveCOMStep
    integer :: VV_RemoveCOMStepFreq
    real(8) :: VV_TEnergySum
    real(8) :: VV_TEnergySqSum
    real(8) :: MC0_NAtt
    real(8) :: MC0_NAcc
    real(8) :: MC0_NAtt2
    real(8) :: MC0_NAcc2
    real(8) :: MC0_Delta
    real(8) :: MC0_Delta2
    real(8) :: MC1_NAtt
    real(8) :: MC1_NAcc
    integer :: MC1_NInsertMols
    integer :: MC1_NDeleteMols
    integer, dimension(0:6750-1) :: MC1_MolMoves
    real(8), dimension(0:2-1, 0:6751-1) :: MC1_BoltzWeights
    integer, dimension(0:2-1) :: MC1_NAccMID
    integer, dimension(0:2-1) :: MC1_NAttMID
    real(8), dimension(0:6751-1, 0:3-1) :: MC1_TM
    real(8) :: MC1_MFactor
    real(8) :: MC2_NAtt
    real(8) :: MC2_NAcc
    logical, dimension(0:3-1) :: MC2_UseAxis
    logical :: MC2_DoIsotropic
    real(8) :: MC2_Delta
    real(8) :: MC2_MaxVol
    real(8) :: MC2_MinVol
    real(8) :: MC3_NAtt
    real(8) :: MC3_NAcc
    integer :: MC3_MID1
    integer :: MC3_MID2
    real(8), dimension(0:3376-1) :: MC3_BoltzWeights
    real(8), dimension(0:3376-1, 0:3-1) :: MC3_TM
    real(8) :: MC3_MFactor
    real(8), dimension(0:4-1) :: MCMoveProbs
    real(8), dimension(0:6750-1, 0:3-1) :: OldPos


contains

subroutine updateactive()
    implicit none
    integer :: i
    integer :: j
    integer :: Start
    integer :: Stop
    integer :: ThisMolID
    integer :: ThisAID

    NDOF = 0
    AIDCount = 0
    NActiveMol = 0
    NInactiveMol = 0
    NActiveMID = 0
    NInactiveMID = 0
    do i = 0, NMol-1
        Start = MolRange(i)
        Stop = MolRange(i+1) - 1
        if (MolActive(i) == 1) then
            NActiveMol = NActiveMol + 1
            ThisMolID = MolID(i)
            NActiveMID(ThisMolID) = NActiveMID(ThisMolID) + 1
            NDOF = NDOF + NDOFMID(ThisMolID)
            do j = Start, Stop
                ThisAID = AID(j)
                AIDCount(ThisAID) = AIDCount(ThisAID) + 1
            enddo
        elseif (MolActive(i) == -1) then
            NInactiveMol = NInactiveMol + 1
            ThisMolID = MolID(i)
            NInactiveMID(ThisMolID) = NInactiveMId(ThisMolID) + 1
        endif
    enddo

end subroutine


subroutine hidemol(MolInd)
    implicit none
    integer, intent(in) :: MolInd
    integer :: i
    integer :: ThisMolID
    integer :: ThisAID

    if (MolActive(MolInd) == 1) then
        ThisMolID = MID(MolInd)
        MolActive(MolInd) = -1
        NActiveMol = NActiveMol - 1
        NInactiveMol = NInactiveMol + 1
        NActiveMID(ThisMolID) = NActiveMID(ThisMolID) - 1
        NInactiveMID(ThisMolID) = NInactiveMID(ThisMolID) + 1
        NDOF = NDOF - NDOFMID(ThisMolID)
        do i = MolRange(MolInd), MolRange(MolInd + 1) - 1
            ThisAID = AID(i)
            AIDCount(ThisAID) = AIDCount(ThisAID) - 1
        enddo
    endif

end subroutine


subroutine showmol(MolInd)
    implicit none
    integer, intent(in) :: MolInd
    integer :: i
    integer :: ThisMolID
    integer :: ThisAID

    if (MolActive(MolInd) == -1) then
        ThisMolID = MID(MolInd)
        MolActive(MolInd) = 1
        NActiveMol = NActiveMol + 1
        NInactiveMol = NInactiveMol - 1
        NActiveMID(ThisMolID) = NActiveMID(ThisMolID) + 1
        NInactiveMID(ThisMolID) = NInactiveMID(ThisMolID) - 1
        NDOF = NDOF + NDOFMID(ThisMolID)
        do i = MolRange(MolInd), MolRange(MolInd+1) - 1
            ThisAID = AID(i)
            AIDCount(ThisAID) = AIDCount(ThisAID) + 1
        enddo
    endif

end subroutine


subroutine calcenergyforces(Mode, CalcForce, CalcVirial, CalcDUParam, CalcDWParam, CalcFluct)
    implicit none
    integer, intent(in) :: Mode
    logical, intent(in) :: CalcForce
    logical, intent(in) :: CalcVirial
    logical, intent(in) :: CalcDUParam
    logical, intent(in) :: CalcDWParam
    logical, intent(in) :: CalcFluct
    real(8), dimension(0:Dim-1) :: iBoxL
    logical :: DoMinImage
    integer :: i
    integer :: j
    integer :: ForceDoMInd
    integer :: ForceDoAtom
    integer :: AIDi
    integer :: AIDj
    integer :: MIndi
    integer :: MIndj
    integer :: istart
    integer :: istop
    integer :: jstart
    integer :: LoopMode
    real(8), dimension(0:Dim-1) :: Posi
    real(8), dimension(0:Dim-1) :: Posj
    real(8), dimension(0:Dim-1) :: PosiMinImage
    real(8), dimension(0:Dim-1) :: rij
    real(8) :: dijsq
    real(8) :: dij
    real(8), parameter :: pi = 3.141592653589d0
    real(8) :: Scale
    real(8), dimension(0:Dim-1) :: Forcei
    real(8) :: ThisU
    real(8) :: ThisW
    real(8) :: idist2
    real(8) :: idist6
    real(8) :: idist12
    real(8) :: val1
    real(8) :: val2
    real(8) :: Sig
    real(8) :: Eps
    real(8) :: val3
    real(8) :: val4
    real(8) :: val5
    real(8) :: val6
    real(8) :: val7
    real(8) :: PlaneDist
    real(8) :: factor

    LoopMode = Mode
    if (Mode == 0) then
        !zero initial quantities
        PEnergy = 0.d0
        Terms = 0.d0
        if (CalcVirial) then
            Virial = 0.d0
        else
            Virial = 0.d0
        endif
        if (CalcForce) Force = 0.d0
        if (CalcDUParam .or. CalcFluct) then
            DUParam = 0.d0
            DDUParam = 0.d0
        endif
        if (CalcDWParam) then
            DWParam = 0.d0
            DDWParam = 0.d0
        endif
        if (CalcFluct) then
            FluctE0 = 0.
            FluctA0 = 0.
        endif
        Scale = 1.d0
    else
        if (CalcDUParam .or. CalcDWParam .or. CalcForce .or. CalcVirial .or. CalcFluct) then
            print *, "Can only use Calcvirial, CalcForce, CalcDUParam, CalcDWParam, CalcFluct for Mode=0."
            stop
        endif
        if (Mode > 0) then
            Scale = 1.d0
        elseif (Mode < 0) then
            Scale = -1.d0
        else
            print *, "Invalid value of Mode in calcenergyforces."
            stop
        endif
    endif

    DoMinImage = any(BoxL > 0.d0)
    iBoxL = 1.d0 / max(1.d-300, BoxL)

    ForceDoMInd = -1
    ForceDoAtom = -1
    if (LoopMode == 0) then
        !normal full atom loop
        istart = 0
        istop = NAtom - 1
    elseif (LoopMode == 1) then
        !single atom interactions for adding TargetAtom
        istart = TargetAtom
        istop = TargetAtom
        ForceDoAtom = TargetAtom
    elseif (LoopMode == -1) then
        !single atom interactions for deleting TargetAtom
        istart = TargetAtom
        istop = TargetAtom
        ForceDoAtom = TargetAtom
    elseif (LoopMode == 2) then
        !molecule interactions fo adding TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
        ForceDoMInd = TargetMol
    elseif (LoopMode == -2) then
        !molecule interactions fo deleting TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
        ForceDoMInd = TargetMol
    else
        print *, "Illegal value for LoopMode."
    endif

    !loop over i
    do i = istart, istop

        MIndi = MInd(i)

        if (MolActive(MIndi) < 0 .and. LoopMode == 0) cycle

        Posi = Pos(i,:)
        PosiMinImage = Posi
        if (DoMinImage) PosiMinImage = Posi - BoxL * dnint(Posi * iBoxL)
        AIDi = AID(i)

        if (AIDi==0) then
            !potential Uext
            !This section calculates the "explicit energy"
            PlaneDist = (POSIMINIMAGE(0) - P3_PlaneLoc)
            if (BOXL(0)> 0.d0) PlaneDist = PlaneDist - BOXL(0) * dnint(PlaneDist * IBOXL(0))
            factor = 2*PI*Param(18)/BOXL(0)
            !THISU = Param(17) * SIN(2*PI*(PlaneDist-P3_PlaneLoc)*Param(18)/BOXL(0))
            THISU = Param(17) * SIN(factor*(PlaneDist-P3_PlaneLoc))
            !THISU = Param(17) * PlaneDist
            PEnergy = PEnergy + ThisU * Scale
            Terms(3) = Terms(3) + ThisU * Scale
            if (CalcDUParam) then
            !    DUParam(17) = DUParam(17) + SIN(2*PI*(PlaneDist-P3_PlaneLoc)*Param(18)/BOXL(0))
                DUParam(17) = DUParam(17) + SIN(factor*(PlaneDist-P3_PlaneLoc))
            endif
            if (CalcForce) then
                FORCEI = 0.d0
            !    FORCEI(0) = -Param(17) * 2*PI*Param(18)/BoxL(0) * COS(2*PI*(PlaneDist-P3_PlaneLoc)*Param(18)/BOXL(0))
                FORCEI(0) = -Param(17) * factor * COS(factor*(PlaneDist-P3_PlaneLoc))
            !    ForceI(0) = - Param(17)
            Force(i,:) = Force(i,:) + Forcei
            endif
        end if

        if (LoopMode == 0) then
            jstart = i+1
        else
            jstart = 0
        endif

        !loop over j
        do j = jstart, NAtom - 1

            !check to see if same atom
            if (i==j) cycle

            MIndj = MInd(j)

            if (LoopMode == 2 .or. LoopMode == -2) then
                !!check to see if we need to skip because of double counting
                if (MIndj == MIndi .and. j < i) cycle
            endif

            if (MolActive(MIndj) < 0 .and. MIndj /= ForceDoMInd .and. j /= ForceDoAtom) cycle

            Posj = Pos(j,:)
            AIDj = AID(j)

            rij = Posj - Posi
            if (DoMinImage) rij = rij - BoxL * dnint(rij * iBoxL)
            dijsq = dot_product(rij, rij)
            dij = -1.d0

            if (dijsq < CutSq(0)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0) .and. (AIDj==0)) .or. ((AIDj==0) .and. (AIDi==0))) then
                    !potential AA
                    Sig = Param(3)
                    Eps = Param(2)
                    idist2 = Sig**2 / DIJSQ
                    idist6 = idist2 * idist2 * idist2
                    idist12 = idist6 * idist6
                    val1 = 4.d0 * (idist12 - idist6)
                    val2 = 24.d0 * idist6 - 48.d0 * idist12
                    val3 = (DIJ - Param(6))
                    val4 = val3*val3
                    val5 = exp(-Param(5) * val4)
                    THISU = Eps * val1 + Param(4) * val5 + P0_UShift(0)
                    PEnergy = PEnergy + ThisU * Scale
                    Terms(0) = Terms(0) + ThisU * Scale
                    if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                        THISW = Eps * val2 - 2.d0 * Param(4) * Param(5) * DIJ * val3 * val5
                        Virial = Virial + ThisW * Scale
                    endif
                    if (CALCFORCE) then
                        FORCEI = (RIJ * THISW / DIJSQ) * SCALE
                        FORCEI = FORCEI * SCALE
                        FORCE(I,:) = FORCE(I,:) + FORCEI
                        FORCE(J,:) = FORCE(J,:) - FORCEI
                    endif
                    if (CALCDUPARAM .or. CALCFLUCT) then
                        DUParam(2) = DUParam(2) + val1 + P0_UShift(1)
                        val6 = 1.d0 / Sig
                        DUParam(3) = DUParam(3) - Eps * val2 * val6 + P0_UShift(2)
                        DDUParam(10) = DDUParam(10) + Param(2) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + P0_UShift(3)
                        DDUParam(5) = DDUParam(5) - val2 * val6 + P0_UShift(4)
                        DUParam(4) = DUParam(4) + val5 + P0_UShift(5)
                        val6 = -val4 * val5
                        DUParam(5) = DUParam(5) + Param(4) * val6 + P0_UShift(6)
                        DDUParam(21) = DDUParam(21) + val6 + P0_UShift(7)
                        val6 = 2.d0 * Param(5) * val3 * val5
                        DUParam(6) = DUParam(6) + Param(4) * val6 + P0_UShift(8)
                        DDUParam(26) = DDUParam(26) + val6 + P0_UShift(9)
                        DDUParam(22) = DDUParam(22) + Param(4) * val4*val4 * val5 + P0_UShift(10)
                        val6 = 2.d0*Param(4) * val5
                        DDUParam(27) = DDUParam(27) + val6 * val3 * (1.d0 - val4*Param(5)) + P0_UShift(11)
                        DDUParam(28) = DDUParam(28) + val6 * Param(5) * (2.d0 * Param(5) * val4 - 1.d0) + P0_UShift(12)
                    endif
                    if (CALCDWPARAM) then
                        DWParam(2) = DWParam(2) + val2
                        val6 = 1.d0 / Sig
                        val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                        DWParam(3) = DWParam(3) + Eps *  val7
                        DDWParam(10) = DDWParam(10) + Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6
                        DDWParam(5) = DDWParam(5) + val7
                        val7 = 2.d0 * DIJ * val3 * val5
                        DWParam(4) = DWParam(4) - val7 * Param(5)
                        val6 = val4 * Param(5)
                        DWParam(5) = DWParam(5) + Param(4) * val7 * (val6 - 1.d0)
                        DDWParam(22) = DDWParam(22) - Param(4) * val7 * val4 * (val6 - 2.d0)
                        DDWParam(21) = DDWParam(21) + val7 * (val6 - 1.d0)
                        val7 = 2.d0 * DIJ * val5
                        DWParam(6) = DWParam(6) - val7 * Param(4) * Param(5) * (2.d0 * Param(5) * val4 - 1.d0)
                        DDWParam(23) = DDWParam(23) + Param(4) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6)
                        DDWParam(28) = DDWParam(28) - 2.d0 * val7 * Param(4) * Param(5)*Param(5) * val3 * (2.d0 * val6 - 3.d0)
                        DDWParam(18) = DDWParam(18) - val7 * Param(5) * (2.d0 * val6 - 1.d0)
                    endif
                end if
            end if

            if (dijsq < CutSq(1)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1) .and. (AIDj==1)) .or. ((AIDj==1) .and. (AIDi==1))) then
                    !potential BB
                    Sig = Param(8)
                    Eps = Param(7)
                    idist2 = Sig**2 / DIJSQ
                    idist6 = idist2 * idist2 * idist2
                    idist12 = idist6 * idist6
                    val1 = 4.d0 * (idist12 - idist6)
                    val2 = 24.d0 * idist6 - 48.d0 * idist12
                    val3 = (DIJ - Param(11))
                    val4 = val3*val3
                    val5 = exp(-Param(10) * val4)
                    THISU = Eps * val1 + Param(9) * val5 + P1_UShift(0)
                    PEnergy = PEnergy + ThisU * Scale
                    Terms(1) = Terms(1) + ThisU * Scale
                    if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                        THISW = Eps * val2 - 2.d0 * Param(9) * Param(10) * DIJ * val3 * val5
                        Virial = Virial + ThisW * Scale
                    endif
                    if (CALCFORCE) then
                        FORCEI = (RIJ * THISW / DIJSQ) * SCALE
                        FORCEI = FORCEI * SCALE
                        FORCE(I,:) = FORCE(I,:) + FORCEI
                        FORCE(J,:) = FORCE(J,:) - FORCEI
                    endif
                    if (CALCDUPARAM .or. CALCFLUCT) then
                        DUParam(7) = DUParam(7) + val1 + P1_UShift(1)
                        val6 = 1.d0 / Sig
                        DUParam(8) = DUParam(8) - Eps * val2 * val6 + P1_UShift(2)
                        DDUParam(35) = DDUParam(35) + Param(7) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + P1_UShift(3)
                        DDUParam(30) = DDUParam(30) - val2 * val6 + P1_UShift(4)
                        DUParam(9) = DUParam(9) + val5 + P1_UShift(5)
                        val6 = -val4 * val5
                        DUParam(10) = DUParam(10) + Param(9) * val6 + P1_UShift(6)
                        DDUParam(46) = DDUParam(46) + val6 + P1_UShift(7)
                        val6 = 2.d0 * Param(10) * val3 * val5
                        DUParam(11) = DUParam(11) + Param(9) * val6 + P1_UShift(8)
                        DDUParam(51) = DDUParam(51) + val6 + P1_UShift(9)
                        DDUParam(47) = DDUParam(47) + Param(9) * val4*val4 * val5 + P1_UShift(10)
                        val6 = 2.d0*Param(9) * val5
                        DDUParam(52) = DDUParam(52) + val6 * val3 * (1.d0 - val4*Param(10)) + P1_UShift(11)
                        DDUParam(53) = DDUParam(53) + val6 * Param(10) * (2.d0 * Param(10) * val4 - 1.d0) + P1_UShift(12)
                    endif
                    if (CALCDWPARAM) then
                        DWParam(7) = DWParam(7) + val2
                        val6 = 1.d0 / Sig
                        val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                        DWParam(8) = DWParam(8) + Eps *  val7
                        DDWParam(35) = DDWParam(35) + Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6
                        DDWParam(30) = DDWParam(30) + val7
                        val7 = 2.d0 * DIJ * val3 * val5
                        DWParam(9) = DWParam(9) - val7 * Param(10)
                        val6 = val4 * Param(10)
                        DWParam(10) = DWParam(10) + Param(9) * val7 * (val6 - 1.d0)
                        DDWParam(47) = DDWParam(47) - Param(9) * val7 * val4 * (val6 - 2.d0)
                        DDWParam(46) = DDWParam(46) + val7 * (val6 - 1.d0)
                        val7 = 2.d0 * DIJ * val5
                        DWParam(11) = DWParam(11) - val7 * Param(9) * Param(10) * (2.d0 * Param(10) * val4 - 1.d0)
                        DDWParam(48) = DDWParam(48) + Param(9) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6)
                        DDWParam(53) = DDWParam(53) - 2.d0 * val7 * Param(9) * Param(10)*Param(10) * val3 * (2.d0 * val6 - 3.d0)
                        DDWParam(43) = DDWParam(43) - val7 * Param(10) * (2.d0 * val6 - 1.d0)
                    endif
                end if
            end if

            if (dijsq < CutSq(2)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0) .and. (AIDj==1)) .or. ((AIDj==0) .and. (AIDi==1))) then
                    !potential AB
                    Sig = Param(13)
                    Eps = Param(12)
                    idist2 = Sig**2 / DIJSQ
                    idist6 = idist2 * idist2 * idist2
                    idist12 = idist6 * idist6
                    val1 = 4.d0 * (idist12 - idist6)
                    val2 = 24.d0 * idist6 - 48.d0 * idist12
                    val3 = (DIJ - Param(16))
                    val4 = val3*val3
                    val5 = exp(-Param(15) * val4)
                    THISU = Eps * val1 + Param(14) * val5 + P2_UShift(0)
                    PEnergy = PEnergy + ThisU * Scale
                    Terms(2) = Terms(2) + ThisU * Scale
                    if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
                        THISW = Eps * val2 - 2.d0 * Param(14) * Param(15) * DIJ * val3 * val5
                        Virial = Virial + ThisW * Scale
                    endif
                    if (CALCFORCE) then
                        FORCEI = (RIJ * THISW / DIJSQ) * SCALE
                        FORCEI = FORCEI * SCALE
                        FORCE(I,:) = FORCE(I,:) + FORCEI
                        FORCE(J,:) = FORCE(J,:) - FORCEI
                    endif
                    if (CALCDUPARAM .or. CALCFLUCT) then
                        DUParam(12) = DUParam(12) + val1 + P2_UShift(1)
                        val6 = 1.d0 / Sig
                        DUParam(13) = DUParam(13) - Eps * val2 * val6 + P2_UShift(2)
                        DDUParam(60) = DDUParam(60) + Param(12) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + P2_UShift(3)
                        DDUParam(55) = DDUParam(55) - val2 * val6 + P2_UShift(4)
                        DUParam(14) = DUParam(14) + val5 + P2_UShift(5)
                        val6 = -val4 * val5
                        DUParam(15) = DUParam(15) + Param(14) * val6 + P2_UShift(6)
                        DDUParam(71) = DDUParam(71) + val6 + P2_UShift(7)
                        val6 = 2.d0 * Param(15) * val3 * val5
                        DUParam(16) = DUParam(16) + Param(14) * val6 + P2_UShift(8)
                        DDUParam(76) = DDUParam(76) + val6 + P2_UShift(9)
                        DDUParam(72) = DDUParam(72) + Param(14) * val4*val4 * val5 + P2_UShift(10)
                        val6 = 2.d0*Param(14) * val5
                        DDUParam(77) = DDUParam(77) + val6 * val3 * (1.d0 - val4*Param(15)) + P2_UShift(11)
                        DDUParam(78) = DDUParam(78) + val6 * Param(15) * (2.d0 * Param(15) * val4 - 1.d0) + P2_UShift(12)
                    endif
                    if (CALCDWPARAM) then
                        DWParam(12) = DWParam(12) + val2
                        val6 = 1.d0 / Sig
                        val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                        DWParam(13) = DWParam(13) + Eps *  val7
                        DDWParam(60) = DDWParam(60) + Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6
                        DDWParam(55) = DDWParam(55) + val7
                        val7 = 2.d0 * DIJ * val3 * val5
                        DWParam(14) = DWParam(14) - val7 * Param(15)
                        val6 = val4 * Param(15)
                        DWParam(15) = DWParam(15) + Param(14) * val7 * (val6 - 1.d0)
                        DDWParam(72) = DDWParam(72) - Param(14) * val7 * val4 * (val6 - 2.d0)
                        DDWParam(71) = DDWParam(71) + val7 * (val6 - 1.d0)
                        val7 = 2.d0 * DIJ * val5
                        DWParam(16) = DWParam(16) - val7 * Param(14) * Param(15) * (2.d0 * Param(15) * val4 - 1.d0)
                        DDWParam(73) = DDWParam(73) + Param(14) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6)
                        DDWParam(78) = DDWParam(78) - 2.d0 * val7 * Param(14) * Param(15)*Param(15) * val3 * (2.d0 * val6 - 3.d0)
                        DDWParam(68) = DDWParam(68) - val7 * Param(15) * (2.d0 * val6 - 1.d0)
                    endif
                end if
            end if

        !end of loop j
        enddo

    !end of loop i
    enddo

end subroutine


subroutine saveenergystate(Mode)
    implicit none
    integer, intent(in) :: Mode
    integer :: istart
    integer :: istop

    if (Mode == 0) then
        !save for all atoms
        istart = 0
        istop = NAtom - 1
    elseif (Mode == 1) then
        !save for single atom TargetAtom
        istart = TargetAtom
        istop = TargetAtom
    elseif (Mode == 2) then
        !save for single molecule TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
    else
        print *, "Illegal value for Mode in SaveEnergyState."
    endif
    OldPEnergy = PEnergy
    OldTerms = Terms

end subroutine


subroutine revertenergystate(Mode)
    implicit none
    integer, intent(in) :: Mode
    integer :: istart
    integer :: istop

    if (Mode == 0) then
        !revert for all atoms
        istart = 0
        istop = NAtom - 1
    elseif (Mode == 1) then
        !revert for single atom TargetAtom
        istart = TargetAtom
        istop = TargetAtom
    elseif (Mode == 2) then
        !revert for single molecule TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
    else
        print *, "Illegal value for Mode in RevertEnergyState."
    endif
    PEnergy = OldPEnergy
    Terms = OldTerms

end subroutine


subroutine calcargstats(Weight)
    implicit none
    real(8), intent(in) :: Weight
    real(8), dimension(0:Dim-1) :: iBoxL
    logical :: DoMinImage
    integer :: i
    integer :: j
    integer :: ForceDoMInd
    integer :: ForceDoAtom
    integer :: AIDi
    integer :: AIDj
    integer :: MIndi
    integer :: MIndj
    integer :: istart
    integer :: istop
    integer :: jstart
    integer :: LoopMode
    real(8), dimension(0:Dim-1) :: Posi
    real(8), dimension(0:Dim-1) :: Posj
    real(8), dimension(0:Dim-1) :: PosiMinImage
    real(8), dimension(0:Dim-1) :: rij
    real(8) :: dijsq
    real(8) :: dij
    real(8), parameter :: pi = 3.141592653589d0
    real(8) :: Scale
    integer :: ArgType
    real(8) :: ArgVal
    real(8) :: PlaneDist
    real(8) :: factor

    Scale = 1.d0
    LoopMode = 0

    !potential AA
    P0_ArgWeightSumStats = P0_ArgWeightSumStats + Weight
    !potential BB
    P1_ArgWeightSumStats = P1_ArgWeightSumStats + Weight
    !potential AB
    P2_ArgWeightSumStats = P2_ArgWeightSumStats + Weight
    !potential Uext
    P3_ArgWeightSumStats = P3_ArgWeightSumStats + Weight

    DoMinImage = any(BoxL > 0.d0)
    iBoxL = 1.d0 / max(1.d-300, BoxL)

    ForceDoMInd = -1
    ForceDoAtom = -1
    if (LoopMode == 0) then
        !normal full atom loop
        istart = 0
        istop = NAtom - 1
    elseif (LoopMode == 1) then
        !single atom interactions for adding TargetAtom
        istart = TargetAtom
        istop = TargetAtom
        ForceDoAtom = TargetAtom
    elseif (LoopMode == -1) then
        !single atom interactions for deleting TargetAtom
        istart = TargetAtom
        istop = TargetAtom
        ForceDoAtom = TargetAtom
    elseif (LoopMode == 2) then
        !molecule interactions fo adding TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
        ForceDoMInd = TargetMol
    elseif (LoopMode == -2) then
        !molecule interactions fo deleting TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
        ForceDoMInd = TargetMol
    else
        print *, "Illegal value for LoopMode."
    endif

    !loop over i
    do i = istart, istop

        MIndi = MInd(i)

        if (MolActive(MIndi) < 0 .and. LoopMode == 0) cycle

        Posi = Pos(i,:)
        PosiMinImage = Posi
        if (DoMinImage) PosiMinImage = Posi - BoxL * dnint(Posi * iBoxL)
        AIDi = AID(i)

        if (AIDi==0) then
            !potential Uext
            PlaneDist = (POSIMINIMAGE(0) - P3_PlaneLoc)
            ARGVAL = PlaneDist
            ARGTYPE = 0
            !factor = 2*PI*Param(18)/BOXL(0)
            if (Weight > 0.d0) then
                P3_ArgMin(ArgType) = min(P3_ArgMin(ArgType), ArgVal)
                P3_ArgMax(ArgType) = max(P3_ArgMax(ArgType), ArgVal)
                P3_ArgCount(ArgType) = P3_ArgCount(ArgType) + Weight
                P3_ArgSum(ArgType) = P3_ArgSum(ArgType) + ArgVal * Weight
                P3_ArgSumSq(ArgType) = P3_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
            endif
        end if

        if (LoopMode == 0) then
            jstart = i+1
        else
            jstart = 0
        endif

        !loop over j
        do j = jstart, NAtom - 1

            !check to see if same atom
            if (i==j) cycle

            MIndj = MInd(j)

            if (LoopMode == 2 .or. LoopMode == -2) then
                !!check to see if we need to skip because of double counting
                if (MIndj == MIndi .and. j < i) cycle
            endif

            if (MolActive(MIndj) < 0 .and. MIndj /= ForceDoMInd .and. j /= ForceDoAtom) cycle

            Posj = Pos(j,:)
            AIDj = AID(j)

            rij = Posj - Posi
            if (DoMinImage) rij = rij - BoxL * dnint(rij * iBoxL)
            dijsq = dot_product(rij, rij)
            dij = -1.d0

            if (dijsq < CutSq(0)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0) .and. (AIDj==0)) .or. ((AIDj==0) .and. (AIDi==0))) then
                    !potential AA
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(0)) cycle
                    if (Weight > 0.d0) then
                        P0_ArgMin(ArgType) = min(P0_ArgMin(ArgType), ArgVal)
                        P0_ArgMax(ArgType) = max(P0_ArgMax(ArgType), ArgVal)
                        P0_ArgCount(ArgType) = P0_ArgCount(ArgType) + Weight
                        P0_ArgSum(ArgType) = P0_ArgSum(ArgType) + ArgVal * Weight
                        P0_ArgSumSq(ArgType) = P0_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(1)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1) .and. (AIDj==1)) .or. ((AIDj==1) .and. (AIDi==1))) then
                    !potential BB
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(1)) cycle
                    if (Weight > 0.d0) then
                        P1_ArgMin(ArgType) = min(P1_ArgMin(ArgType), ArgVal)
                        P1_ArgMax(ArgType) = max(P1_ArgMax(ArgType), ArgVal)
                        P1_ArgCount(ArgType) = P1_ArgCount(ArgType) + Weight
                        P1_ArgSum(ArgType) = P1_ArgSum(ArgType) + ArgVal * Weight
                        P1_ArgSumSq(ArgType) = P1_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(2)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0) .and. (AIDj==1)) .or. ((AIDj==0) .and. (AIDi==1))) then
                    !potential AB
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(2)) cycle
                    if (Weight > 0.d0) then
                        P2_ArgMin(ArgType) = min(P2_ArgMin(ArgType), ArgVal)
                        P2_ArgMax(ArgType) = max(P2_ArgMax(ArgType), ArgVal)
                        P2_ArgCount(ArgType) = P2_ArgCount(ArgType) + Weight
                        P2_ArgSum(ArgType) = P2_ArgSum(ArgType) + ArgVal * Weight
                        P2_ArgSumSq(ArgType) = P2_ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                    endif
                end if
            end if

        !end of loop j
        enddo

    !end of loop i
    enddo

end subroutine


subroutine calcarghist(Weight)
    implicit none
    real(8), intent(in) :: Weight
    real(8), dimension(0:Dim-1) :: iBoxL
    logical :: DoMinImage
    integer :: i
    integer :: j
    integer :: ForceDoMInd
    integer :: ForceDoAtom
    integer :: AIDi
    integer :: AIDj
    integer :: MIndi
    integer :: MIndj
    integer :: istart
    integer :: istop
    integer :: jstart
    integer :: LoopMode
    real(8), dimension(0:Dim-1) :: Posi
    real(8), dimension(0:Dim-1) :: Posj
    real(8), dimension(0:Dim-1) :: PosiMinImage
    real(8), dimension(0:Dim-1) :: rij
    real(8) :: dijsq
    real(8) :: dij
    real(8), parameter :: pi = 3.141592653589d0
    real(8) :: Scale
    integer :: ArgType
    real(8) :: ArgVal
    integer :: BinInd
    real(8) :: PlaneDist
    real(8) :: factor

    Scale = 1.d0
    LoopMode = 0

    !potential AA
    P0_ArgWeightSumHist = P0_ArgWeightSumHist + Weight
    !potential BB
    P1_ArgWeightSumHist = P1_ArgWeightSumHist + Weight
    !potential AB
    P2_ArgWeightSumHist = P2_ArgWeightSumHist + Weight
    !potential Uext
    P3_ArgWeightSumHist = P3_ArgWeightSumHist + Weight

    DoMinImage = any(BoxL > 0.d0)
    iBoxL = 1.d0 / max(1.d-300, BoxL)

    ForceDoMInd = -1
    ForceDoAtom = -1
    if (LoopMode == 0) then
        !normal full atom loop
        istart = 0
        istop = NAtom - 1
    elseif (LoopMode == 1) then
        !single atom interactions for adding TargetAtom
        istart = TargetAtom
        istop = TargetAtom
        ForceDoAtom = TargetAtom
    elseif (LoopMode == -1) then
        !single atom interactions for deleting TargetAtom
        istart = TargetAtom
        istop = TargetAtom
        ForceDoAtom = TargetAtom
    elseif (LoopMode == 2) then
        !molecule interactions fo adding TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
        ForceDoMInd = TargetMol
    elseif (LoopMode == -2) then
        !molecule interactions fo deleting TargetMol
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
        ForceDoMInd = TargetMol
    else
        print *, "Illegal value for LoopMode."
    endif

    !loop over i
    do i = istart, istop

        MIndi = MInd(i)

        if (MolActive(MIndi) < 0 .and. LoopMode == 0) cycle

        Posi = Pos(i,:)
        PosiMinImage = Posi
        if (DoMinImage) PosiMinImage = Posi - BoxL * dnint(Posi * iBoxL)
        AIDi = AID(i)

        if (AIDi==0) then
            !potential Uext
            PlaneDist = (POSIMINIMAGE(0) - P3_PlaneLoc)
            ARGVAL = PlaneDist
            ARGTYPE = 0
            !factor = 2*PI*Param(18)/BOXL(0)
            if (10000 > 0) then
                BinInd = int((ArgVal - P3_ArgHistMin(ArgType)) * P3_ArgHistiBinw(ArgType))
                if (BinInd >= 0 .and. BinInd < 10000) P3_ArgHist(ArgType,BinInd) = P3_ArgHist(ArgType,BinInd) + Weight
            endif
        end if

        if (LoopMode == 0) then
            jstart = i+1
        else
            jstart = 0
        endif

        !loop over j
        do j = jstart, NAtom - 1

            !check to see if same atom
            if (i==j) cycle

            MIndj = MInd(j)

            if (LoopMode == 2 .or. LoopMode == -2) then
                !!check to see if we need to skip because of double counting
                if (MIndj == MIndi .and. j < i) cycle
            endif

            if (MolActive(MIndj) < 0 .and. MIndj /= ForceDoMInd .and. j /= ForceDoAtom) cycle

            Posj = Pos(j,:)
            AIDj = AID(j)

            rij = Posj - Posi
            if (DoMinImage) rij = rij - BoxL * dnint(rij * iBoxL)
            dijsq = dot_product(rij, rij)
            dij = -1.d0

            if (dijsq < CutSq(0)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0) .and. (AIDj==0)) .or. ((AIDj==0) .and. (AIDi==0))) then
                    !potential AA
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(0)) cycle
                    if (10000 > 0) then
                        BinInd = int((ArgVal - P0_ArgHistMin(ArgType)) * P0_ArgHistiBinw(ArgType))
                        if (BinInd >= 0 .and. BinInd < 10000) P0_ArgHist(ArgType,BinInd) = P0_ArgHist(ArgType,BinInd) + Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(1)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==1) .and. (AIDj==1)) .or. ((AIDj==1) .and. (AIDi==1))) then
                    !potential BB
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(1)) cycle
                    if (10000 > 0) then
                        BinInd = int((ArgVal - P1_ArgHistMin(ArgType)) * P1_ArgHistiBinw(ArgType))
                        if (BinInd >= 0 .and. BinInd < 10000) P1_ArgHist(ArgType,BinInd) = P1_ArgHist(ArgType,BinInd) + Weight
                    endif
                end if
            end if

            if (dijsq < CutSq(2)) then
                if (dij < 0.d0) dij = sqrt(dijsq)
                if (((AIDi==0) .and. (AIDj==1)) .or. ((AIDj==0) .and. (AIDi==1))) then
                    !potential AB
                    ARGVAL = DIJ
                    ARGTYPE = 0
                    if (ArgVal*ArgVal > CutSq(2)) cycle
                    if (10000 > 0) then
                        BinInd = int((ArgVal - P2_ArgHistMin(ArgType)) * P2_ArgHistiBinw(ArgType))
                        if (BinInd >= 0 .and. BinInd < 10000) P2_ArgHist(ArgType,BinInd) = P2_ArgHist(ArgType,BinInd) + Weight
                    endif
                end if
            end if

        !end of loop j
        enddo

    !end of loop i
    enddo

end subroutine


subroutine calcargeval(CalcVirial, CalcDUParam, CalcDWParam, CalcFluct)
    implicit none
    logical, intent(in) :: CalcVirial
    logical, intent(in) :: CalcDUParam
    logical, intent(in) :: CalcDWParam
    logical, intent(in) :: CalcFluct
    integer :: i
    real(8) :: dijsq
    real(8) :: dij
    real(8), parameter :: pi = 3.141592653589d0
    real(8) :: Scale
    real(8) :: ThisU
    real(8) :: ThisW
    integer :: ArgType
    real(8) :: ArgVal
    real(8) :: ThisHist
    real(8) :: idist2
    real(8) :: idist6
    real(8) :: idist12
    real(8) :: val1
    real(8) :: val2
    real(8) :: Sig
    real(8) :: Eps
    real(8) :: val3
    real(8) :: val4
    real(8) :: val5
    real(8) :: val6
    real(8) :: val7
    real(8) :: PlaneDist
    real(8) :: factor

    !compute initial quantities
    PEnergy = 0.d0
    Virial = 0.d0
    if (CalcDUParam) then
        DUParam = 0.d0
        DDUParam = 0.d0
    endif
    if (CalcDWParam) then
        DWParam = 0.d0
        DDWParam = 0.d0
    endif
    Terms = 0.d0
    Scale = 1.d0

    !potential AA
    do ArgType = 0, 1 - 1
        do i = 0, 10000 - 1
            ThisHist = P0_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P0_ArgHistMin(ArgType) + P0_ArgHistBinw(ArgType) * (0.5d0 + i)
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            Sig = Param(3)
            Eps = Param(2)
            idist2 = Sig**2 / DIJSQ
            idist6 = idist2 * idist2 * idist2
            idist12 = idist6 * idist6
            val1 = 4.d0 * (idist12 - idist6)
            val2 = 24.d0 * idist6 - 48.d0 * idist12
            val3 = (DIJ - Param(6))
            val4 = val3*val3
            val5 = exp(-Param(5) * val4)
            THISU = Eps * val1 + Param(4) * val5 + P0_UShift(0)
            PEnergy = PEnergy + ThisU * ThisHist
            Terms(0) = Terms(0) + ThisU * ThisHist
            if (CALCVIRIAL .or. CALCDWPARAM) then
                THISW = Eps * val2 - 2.d0 * Param(4) * Param(5) * DIJ * val3 * val5
                Virial = Virial + ThisW * ThisHist
            endif
            if (CALCDUPARAM .or. CALCFLUCT) then
                DUParam(2) = DUParam(2) + (ThisHist) * (val1 + P0_UShift(1))
                val6 = 1.d0 / Sig
                DUParam(3) = DUParam(3) + (ThisHist) * (-Eps * val2 * val6 + P0_UShift(2))
                DDUParam(10) = DDUParam(10) + (ThisHist) * (Param(2) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + &
                  & P0_UShift(3))
                DDUParam(5) = DDUParam(5) + (ThisHist) * (-val2 * val6 + P0_UShift(4))
                DUParam(4) = DUParam(4) + (ThisHist) * (val5 + P0_UShift(5))
                val6 = -val4 * val5
                DUParam(5) = DUParam(5) + (ThisHist) * (Param(4) * val6 + P0_UShift(6))
                DDUParam(21) = DDUParam(21) + (ThisHist) * (val6 + P0_UShift(7))
                val6 = 2.d0 * Param(5) * val3 * val5
                DUParam(6) = DUParam(6) + (ThisHist) * (Param(4) * val6 + P0_UShift(8))
                DDUParam(26) = DDUParam(26) + (ThisHist) * (val6 + P0_UShift(9))
                DDUParam(22) = DDUParam(22) + (ThisHist) * (Param(4) * val4*val4 * val5 + P0_UShift(10))
                val6 = 2.d0*Param(4) * val5
                DDUParam(27) = DDUParam(27) + (ThisHist) * (val6 * val3 * (1.d0 - val4*Param(5)) + P0_UShift(11))
                DDUParam(28) = DDUParam(28) + (ThisHist) * (val6 * Param(5) * (2.d0 * Param(5) * val4 - 1.d0) + P0_UShift(12))
            endif
            if (CALCDWPARAM) then
                DWParam(2) = DWParam(2) + (ThisHist) * (val2)
                val6 = 1.d0 / Sig
                val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                DWParam(3) = DWParam(3) + (ThisHist) * (Eps *  val7)
                DDWParam(10) = DDWParam(10) + (ThisHist) * (Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6)
                DDWParam(5) = DDWParam(5) + (ThisHist) * (val7)
                val7 = 2.d0 * DIJ * val3 * val5
                DWParam(4) = DWParam(4) + (ThisHist) * (-val7 * Param(5))
                val6 = val4 * Param(5)
                DWParam(5) = DWParam(5) + (ThisHist) * (Param(4) * val7 * (val6 - 1.d0))
                DDWParam(22) = DDWParam(22) + (ThisHist) * (-Param(4) * val7 * val4 * (val6 - 2.d0))
                DDWParam(21) = DDWParam(21) + (ThisHist) * (val7 * (val6 - 1.d0))
                val7 = 2.d0 * DIJ * val5
                DWParam(6) = DWParam(6) + (ThisHist) * (-val7 * Param(4) * Param(5) * (2.d0 * Param(5) * val4 - 1.d0))
                DDWParam(23) = DDWParam(23) + (ThisHist) * (Param(4) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6))
                DDWParam(28) = DDWParam(28) + (ThisHist) * (-2.d0 * val7 * Param(4) * Param(5)*Param(5) * val3 * (2.d0 * val6 - &
                  & 3.d0))
                DDWParam(18) = DDWParam(18) + (ThisHist) * (-val7 * Param(5) * (2.d0 * val6 - 1.d0))
            endif
        enddo
    enddo

    !potential BB
    do ArgType = 0, 1 - 1
        do i = 0, 10000 - 1
            ThisHist = P1_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P1_ArgHistMin(ArgType) + P1_ArgHistBinw(ArgType) * (0.5d0 + i)
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            Sig = Param(8)
            Eps = Param(7)
            idist2 = Sig**2 / DIJSQ
            idist6 = idist2 * idist2 * idist2
            idist12 = idist6 * idist6
            val1 = 4.d0 * (idist12 - idist6)
            val2 = 24.d0 * idist6 - 48.d0 * idist12
            val3 = (DIJ - Param(11))
            val4 = val3*val3
            val5 = exp(-Param(10) * val4)
            THISU = Eps * val1 + Param(9) * val5 + P1_UShift(0)
            PEnergy = PEnergy + ThisU * ThisHist
            Terms(1) = Terms(1) + ThisU * ThisHist
            if (CALCVIRIAL .or. CALCDWPARAM) then
                THISW = Eps * val2 - 2.d0 * Param(9) * Param(10) * DIJ * val3 * val5
                Virial = Virial + ThisW * ThisHist
            endif
            if (CALCDUPARAM .or. CALCFLUCT) then
                DUParam(7) = DUParam(7) + (ThisHist) * (val1 + P1_UShift(1))
                val6 = 1.d0 / Sig
                DUParam(8) = DUParam(8) + (ThisHist) * (-Eps * val2 * val6 + P1_UShift(2))
                DDUParam(35) = DDUParam(35) + (ThisHist) * (Param(7) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + &
                  & P1_UShift(3))
                DDUParam(30) = DDUParam(30) + (ThisHist) * (-val2 * val6 + P1_UShift(4))
                DUParam(9) = DUParam(9) + (ThisHist) * (val5 + P1_UShift(5))
                val6 = -val4 * val5
                DUParam(10) = DUParam(10) + (ThisHist) * (Param(9) * val6 + P1_UShift(6))
                DDUParam(46) = DDUParam(46) + (ThisHist) * (val6 + P1_UShift(7))
                val6 = 2.d0 * Param(10) * val3 * val5
                DUParam(11) = DUParam(11) + (ThisHist) * (Param(9) * val6 + P1_UShift(8))
                DDUParam(51) = DDUParam(51) + (ThisHist) * (val6 + P1_UShift(9))
                DDUParam(47) = DDUParam(47) + (ThisHist) * (Param(9) * val4*val4 * val5 + P1_UShift(10))
                val6 = 2.d0*Param(9) * val5
                DDUParam(52) = DDUParam(52) + (ThisHist) * (val6 * val3 * (1.d0 - val4*Param(10)) + P1_UShift(11))
                DDUParam(53) = DDUParam(53) + (ThisHist) * (val6 * Param(10) * (2.d0 * Param(10) * val4 - 1.d0) + P1_UShift(12))
            endif
            if (CALCDWPARAM) then
                DWParam(7) = DWParam(7) + (ThisHist) * (val2)
                val6 = 1.d0 / Sig
                val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                DWParam(8) = DWParam(8) + (ThisHist) * (Eps *  val7)
                DDWParam(35) = DDWParam(35) + (ThisHist) * (Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6)
                DDWParam(30) = DDWParam(30) + (ThisHist) * (val7)
                val7 = 2.d0 * DIJ * val3 * val5
                DWParam(9) = DWParam(9) + (ThisHist) * (-val7 * Param(10))
                val6 = val4 * Param(10)
                DWParam(10) = DWParam(10) + (ThisHist) * (Param(9) * val7 * (val6 - 1.d0))
                DDWParam(47) = DDWParam(47) + (ThisHist) * (-Param(9) * val7 * val4 * (val6 - 2.d0))
                DDWParam(46) = DDWParam(46) + (ThisHist) * (val7 * (val6 - 1.d0))
                val7 = 2.d0 * DIJ * val5
                DWParam(11) = DWParam(11) + (ThisHist) * (-val7 * Param(9) * Param(10) * (2.d0 * Param(10) * val4 - 1.d0))
                DDWParam(48) = DDWParam(48) + (ThisHist) * (Param(9) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6))
                DDWParam(53) = DDWParam(53) + (ThisHist) * (-2.d0 * val7 * Param(9) * Param(10)*Param(10) * val3 * (2.d0 * val6 &
                  & - 3.d0))
                DDWParam(43) = DDWParam(43) + (ThisHist) * (-val7 * Param(10) * (2.d0 * val6 - 1.d0))
            endif
        enddo
    enddo

    !potential AB
    do ArgType = 0, 1 - 1
        do i = 0, 10000 - 1
            ThisHist = P2_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P2_ArgHistMin(ArgType) + P2_ArgHistBinw(ArgType) * (0.5d0 + i)
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            Sig = Param(13)
            Eps = Param(12)
            idist2 = Sig**2 / DIJSQ
            idist6 = idist2 * idist2 * idist2
            idist12 = idist6 * idist6
            val1 = 4.d0 * (idist12 - idist6)
            val2 = 24.d0 * idist6 - 48.d0 * idist12
            val3 = (DIJ - Param(16))
            val4 = val3*val3
            val5 = exp(-Param(15) * val4)
            THISU = Eps * val1 + Param(14) * val5 + P2_UShift(0)
            PEnergy = PEnergy + ThisU * ThisHist
            Terms(2) = Terms(2) + ThisU * ThisHist
            if (CALCVIRIAL .or. CALCDWPARAM) then
                THISW = Eps * val2 - 2.d0 * Param(14) * Param(15) * DIJ * val3 * val5
                Virial = Virial + ThisW * ThisHist
            endif
            if (CALCDUPARAM .or. CALCFLUCT) then
                DUParam(12) = DUParam(12) + (ThisHist) * (val1 + P2_UShift(1))
                val6 = 1.d0 / Sig
                DUParam(13) = DUParam(13) + (ThisHist) * (-Eps * val2 * val6 + P2_UShift(2))
                DDUParam(60) = DDUParam(60) + (ThisHist) * (Param(12) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + &
                  & P2_UShift(3))
                DDUParam(55) = DDUParam(55) + (ThisHist) * (-val2 * val6 + P2_UShift(4))
                DUParam(14) = DUParam(14) + (ThisHist) * (val5 + P2_UShift(5))
                val6 = -val4 * val5
                DUParam(15) = DUParam(15) + (ThisHist) * (Param(14) * val6 + P2_UShift(6))
                DDUParam(71) = DDUParam(71) + (ThisHist) * (val6 + P2_UShift(7))
                val6 = 2.d0 * Param(15) * val3 * val5
                DUParam(16) = DUParam(16) + (ThisHist) * (Param(14) * val6 + P2_UShift(8))
                DDUParam(76) = DDUParam(76) + (ThisHist) * (val6 + P2_UShift(9))
                DDUParam(72) = DDUParam(72) + (ThisHist) * (Param(14) * val4*val4 * val5 + P2_UShift(10))
                val6 = 2.d0*Param(14) * val5
                DDUParam(77) = DDUParam(77) + (ThisHist) * (val6 * val3 * (1.d0 - val4*Param(15)) + P2_UShift(11))
                DDUParam(78) = DDUParam(78) + (ThisHist) * (val6 * Param(15) * (2.d0 * Param(15) * val4 - 1.d0) + P2_UShift(12))
            endif
            if (CALCDWPARAM) then
                DWParam(12) = DWParam(12) + (ThisHist) * (val2)
                val6 = 1.d0 / Sig
                val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6
                DWParam(13) = DWParam(13) + (ThisHist) * (Eps *  val7)
                DDWParam(60) = DDWParam(60) + (ThisHist) * (Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6)
                DDWParam(55) = DDWParam(55) + (ThisHist) * (val7)
                val7 = 2.d0 * DIJ * val3 * val5
                DWParam(14) = DWParam(14) + (ThisHist) * (-val7 * Param(15))
                val6 = val4 * Param(15)
                DWParam(15) = DWParam(15) + (ThisHist) * (Param(14) * val7 * (val6 - 1.d0))
                DDWParam(72) = DDWParam(72) + (ThisHist) * (-Param(14) * val7 * val4 * (val6 - 2.d0))
                DDWParam(71) = DDWParam(71) + (ThisHist) * (val7 * (val6 - 1.d0))
                val7 = 2.d0 * DIJ * val5
                DWParam(16) = DWParam(16) + (ThisHist) * (-val7 * Param(14) * Param(15) * (2.d0 * Param(15) * val4 - 1.d0))
                DDWParam(73) = DDWParam(73) + (ThisHist) * (Param(14) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6))
                DDWParam(78) = DDWParam(78) + (ThisHist) * (-2.d0 * val7 * Param(14) * Param(15)*Param(15) * val3 * (2.d0 * &
                  & val6 - 3.d0))
                DDWParam(68) = DDWParam(68) + (ThisHist) * (-val7 * Param(15) * (2.d0 * val6 - 1.d0))
            endif
        enddo
    enddo

    !potential Uext
    do ArgType = 0, 1 - 1
        do i = 0, 10000 - 1
            ThisHist = P3_ArgHist(ArgType, i)
            if (ThisHist == 0.d0) cycle
            ArgVal = P3_ArgHistMin(ArgType) + P3_ArgHistBinw(ArgType) * (0.5d0 + i)
            !This section is for histogram evaluation
            PlaneDist = ARGVAL
            factor = 2*PI*Param(18)/BOXL(0)
            !THISU = Param(17) * SIN(2*PI*(PlaneDist-P3_PlaneLoc)*Param(18)/BOXL(0))
            THISU = Param(17) * SIN(factor*(PlaneDist-P3_PlaneLoc))
            !THISU = Param(17) * PlaneDist
            if (CalcDUParam) then
            !    DUParam(17) = DUParam(17) + (ThisHist) * (SIN(2*PI*(PlaneDist-P3_PlaneLoc)*Param(18)/BOXL(0)))
                DUParam(17) = DUParam(17) + (ThisHist) * (SIN(factor*(PlaneDist-P3_PlaneLoc)))
            endif
            PEnergy = PEnergy + ThisU * ThisHist
            Terms(3) = Terms(3) + ThisU * ThisHist
        enddo
    enddo

end subroutine


subroutine calcmeasures(MeasureAll, StepNum, CycleNum, Weight, ConservesMomentum)
    implicit none
    logical, intent(in) :: MeasureAll
    integer, intent(in) :: StepNum
    integer, intent(in) :: CycleNum
    real(8), intent(in) :: Weight
    logical, intent(in) :: ConservesMomentum
    integer :: i
    integer :: j
    logical, dimension(0:NMeasure-1) :: UseMeasure
    integer :: v1
    integer :: v2
    real(8) :: Val0

    !measure KEnergy
    if (M0_Active) then
        if (MeasureAll) then
            UseMeasure(0) = .true.
            M0_Val = 0.d0
        elseif (M0_StepFreq > 0 .and. mod(StepNum, M0_StepFreq)==0) then
            UseMeasure(0) = .true.
            M0_Val = 0.d0
        elseif (M0_CycleFreq > 0 .and. mod(CycleNum, M0_CycleFreq)==0) then
            UseMeasure(0) = .true.
            M0_Val = 0.d0
        else
            UseMeasure(0) = .false.
        endif
    else
        UseMeasure(0) = .false.
    endif
    !measure PEnergy
    if (M1_Active) then
        if (MeasureAll) then
            UseMeasure(1) = .true.
            M1_Val = 0.d0
        elseif (M1_StepFreq > 0 .and. mod(StepNum, M1_StepFreq)==0) then
            UseMeasure(1) = .true.
            M1_Val = 0.d0
        elseif (M1_CycleFreq > 0 .and. mod(CycleNum, M1_CycleFreq)==0) then
            UseMeasure(1) = .true.
            M1_Val = 0.d0
        else
            UseMeasure(1) = .false.
        endif
    else
        UseMeasure(1) = .false.
    endif
    !measure TEnergy
    if (M2_Active) then
        if (MeasureAll) then
            UseMeasure(2) = .true.
            M2_Val = 0.d0
        elseif (M2_StepFreq > 0 .and. mod(StepNum, M2_StepFreq)==0) then
            UseMeasure(2) = .true.
            M2_Val = 0.d0
        elseif (M2_CycleFreq > 0 .and. mod(CycleNum, M2_CycleFreq)==0) then
            UseMeasure(2) = .true.
            M2_Val = 0.d0
        else
            UseMeasure(2) = .false.
        endif
    else
        UseMeasure(2) = .false.
    endif
    !measure KTemp
    if (M3_Active) then
        if (MeasureAll) then
            UseMeasure(3) = .true.
            M3_Val = 0.d0
        elseif (M3_StepFreq > 0 .and. mod(StepNum, M3_StepFreq)==0) then
            UseMeasure(3) = .true.
            M3_Val = 0.d0
        elseif (M3_CycleFreq > 0 .and. mod(CycleNum, M3_CycleFreq)==0) then
            UseMeasure(3) = .true.
            M3_Val = 0.d0
        else
            UseMeasure(3) = .false.
        endif
    else
        UseMeasure(3) = .false.
    endif
    !measure Vol
    if (M4_Active) then
        if (MeasureAll) then
            UseMeasure(4) = .true.
            M4_Val = 0.d0
        elseif (M4_StepFreq > 0 .and. mod(StepNum, M4_StepFreq)==0) then
            UseMeasure(4) = .true.
            M4_Val = 0.d0
        elseif (M4_CycleFreq > 0 .and. mod(CycleNum, M4_CycleFreq)==0) then
            UseMeasure(4) = .true.
            M4_Val = 0.d0
        else
            UseMeasure(4) = .false.
        endif
    else
        UseMeasure(4) = .false.
    endif
    !measure Pressure
    if (M5_Active) then
        if (MeasureAll) then
            UseMeasure(5) = .true.
            M5_Val = 0.d0
        elseif (M5_StepFreq > 0 .and. mod(StepNum, M5_StepFreq)==0) then
            UseMeasure(5) = .true.
            M5_Val = 0.d0
        elseif (M5_CycleFreq > 0 .and. mod(CycleNum, M5_CycleFreq)==0) then
            UseMeasure(5) = .true.
            M5_Val = 0.d0
        else
            UseMeasure(5) = .false.
        endif
    else
        UseMeasure(5) = .false.
    endif
    !measure Virial
    if (M6_Active) then
        if (MeasureAll) then
            UseMeasure(6) = .true.
            M6_Val = 0.d0
        elseif (M6_StepFreq > 0 .and. mod(StepNum, M6_StepFreq)==0) then
            UseMeasure(6) = .true.
            M6_Val = 0.d0
        elseif (M6_CycleFreq > 0 .and. mod(CycleNum, M6_CycleFreq)==0) then
            UseMeasure(6) = .true.
            M6_Val = 0.d0
        else
            UseMeasure(6) = .false.
        endif
    else
        UseMeasure(6) = .false.
    endif
    !measure N_A
    if (M7_Active) then
        if (MeasureAll) then
            UseMeasure(7) = .true.
            M7_Val = 0.d0
        elseif (M7_StepFreq > 0 .and. mod(StepNum, M7_StepFreq)==0) then
            UseMeasure(7) = .true.
            M7_Val = 0.d0
        elseif (M7_CycleFreq > 0 .and. mod(CycleNum, M7_CycleFreq)==0) then
            UseMeasure(7) = .true.
            M7_Val = 0.d0
        else
            UseMeasure(7) = .false.
        endif
    else
        UseMeasure(7) = .false.
    endif
    !measure N_B
    if (M8_Active) then
        if (MeasureAll) then
            UseMeasure(8) = .true.
            M8_Val = 0.d0
        elseif (M8_StepFreq > 0 .and. mod(StepNum, M8_StepFreq)==0) then
            UseMeasure(8) = .true.
            M8_Val = 0.d0
        elseif (M8_CycleFreq > 0 .and. mod(CycleNum, M8_CycleFreq)==0) then
            UseMeasure(8) = .true.
            M8_Val = 0.d0
        else
            UseMeasure(8) = .false.
        endif
    else
        UseMeasure(8) = .false.
    endif
    !measure r_A
    if (M9_Active) then
        if (MeasureAll) then
            UseMeasure(9) = .true.
            M9_Val = 0.d0
        elseif (M9_StepFreq > 0 .and. mod(StepNum, M9_StepFreq)==0) then
            UseMeasure(9) = .true.
            M9_Val = 0.d0
        elseif (M9_CycleFreq > 0 .and. mod(CycleNum, M9_CycleFreq)==0) then
            UseMeasure(9) = .true.
            M9_Val = 0.d0
        else
            UseMeasure(9) = .false.
        endif
    else
        UseMeasure(9) = .false.
    endif
    !measure r_B
    if (M10_Active) then
        if (MeasureAll) then
            UseMeasure(10) = .true.
            M10_Val = 0.d0
        elseif (M10_StepFreq > 0 .and. mod(StepNum, M10_StepFreq)==0) then
            UseMeasure(10) = .true.
            M10_Val = 0.d0
        elseif (M10_CycleFreq > 0 .and. mod(CycleNum, M10_CycleFreq)==0) then
            UseMeasure(10) = .true.
            M10_Val = 0.d0
        else
            UseMeasure(10) = .false.
        endif
    else
        UseMeasure(10) = .false.
    endif
    !measure x_A
    if (M11_Active) then
        if (MeasureAll) then
            UseMeasure(11) = .true.
            M11_Val = 0.d0
        elseif (M11_StepFreq > 0 .and. mod(StepNum, M11_StepFreq)==0) then
            UseMeasure(11) = .true.
            M11_Val = 0.d0
        elseif (M11_CycleFreq > 0 .and. mod(CycleNum, M11_CycleFreq)==0) then
            UseMeasure(11) = .true.
            M11_Val = 0.d0
        else
            UseMeasure(11) = .false.
        endif
    else
        UseMeasure(11) = .false.
    endif
    !measure x_B
    if (M12_Active) then
        if (MeasureAll) then
            UseMeasure(12) = .true.
            M12_Val = 0.d0
        elseif (M12_StepFreq > 0 .and. mod(StepNum, M12_StepFreq)==0) then
            UseMeasure(12) = .true.
            M12_Val = 0.d0
        elseif (M12_CycleFreq > 0 .and. mod(CycleNum, M12_CycleFreq)==0) then
            UseMeasure(12) = .true.
            M12_Val = 0.d0
        else
            UseMeasure(12) = .false.
        endif
    else
        UseMeasure(12) = .false.
    endif
    !measure DUParam
    if (M13_Active) then
        if (MeasureAll) then
            UseMeasure(13) = .true.
            M13_Val = 0.d0
        elseif (M13_StepFreq > 0 .and. mod(StepNum, M13_StepFreq)==0) then
            UseMeasure(13) = .true.
            M13_Val = 0.d0
        elseif (M13_CycleFreq > 0 .and. mod(CycleNum, M13_CycleFreq)==0) then
            UseMeasure(13) = .true.
            M13_Val = 0.d0
        else
            UseMeasure(13) = .false.
        endif
    else
        UseMeasure(13) = .false.
    endif
    !measure DDUParam
    if (M14_Active) then
        if (MeasureAll) then
            UseMeasure(14) = .true.
            M14_Val = 0.d0
        elseif (M14_StepFreq > 0 .and. mod(StepNum, M14_StepFreq)==0) then
            UseMeasure(14) = .true.
            M14_Val = 0.d0
        elseif (M14_CycleFreq > 0 .and. mod(CycleNum, M14_CycleFreq)==0) then
            UseMeasure(14) = .true.
            M14_Val = 0.d0
        else
            UseMeasure(14) = .false.
        endif
    else
        UseMeasure(14) = .false.
    endif
    !measure DWParam
    if (M15_Active) then
        if (MeasureAll) then
            UseMeasure(15) = .true.
            M15_Val = 0.d0
        elseif (M15_StepFreq > 0 .and. mod(StepNum, M15_StepFreq)==0) then
            UseMeasure(15) = .true.
            M15_Val = 0.d0
        elseif (M15_CycleFreq > 0 .and. mod(CycleNum, M15_CycleFreq)==0) then
            UseMeasure(15) = .true.
            M15_Val = 0.d0
        else
            UseMeasure(15) = .false.
        endif
    else
        UseMeasure(15) = .false.
    endif
    !measure DDWParam
    if (M16_Active) then
        if (MeasureAll) then
            UseMeasure(16) = .true.
            M16_Val = 0.d0
        elseif (M16_StepFreq > 0 .and. mod(StepNum, M16_StepFreq)==0) then
            UseMeasure(16) = .true.
            M16_Val = 0.d0
        elseif (M16_CycleFreq > 0 .and. mod(CycleNum, M16_CycleFreq)==0) then
            UseMeasure(16) = .true.
            M16_Val = 0.d0
        else
            UseMeasure(16) = .false.
        endif
    else
        UseMeasure(16) = .false.
    endif
    !measure DUParamDWParam
    if (M17_Active) then
        if (MeasureAll) then
            UseMeasure(17) = .true.
            M17_Val = 0.d0
        elseif (M17_StepFreq > 0 .and. mod(StepNum, M17_StepFreq)==0) then
            UseMeasure(17) = .true.
            M17_Val = 0.d0
        elseif (M17_CycleFreq > 0 .and. mod(CycleNum, M17_CycleFreq)==0) then
            UseMeasure(17) = .true.
            M17_Val = 0.d0
        else
            UseMeasure(17) = .false.
        endif
    else
        UseMeasure(17) = .false.
    endif
    !measure FluctTerm
    if (M18_Active) then
        if (MeasureAll) then
            UseMeasure(18) = .true.
            M18_Val = 0.d0
        elseif (M18_StepFreq > 0 .and. mod(StepNum, M18_StepFreq)==0) then
            UseMeasure(18) = .true.
            M18_Val = 0.d0
        elseif (M18_CycleFreq > 0 .and. mod(CycleNum, M18_CycleFreq)==0) then
            UseMeasure(18) = .true.
            M18_Val = 0.d0
        else
            UseMeasure(18) = .false.
        endif
    else
        UseMeasure(18) = .false.
    endif
    !measure FluctE0
    if (M19_Active) then
        if (MeasureAll) then
            UseMeasure(19) = .true.
            M19_Val = 0.d0
        elseif (M19_StepFreq > 0 .and. mod(StepNum, M19_StepFreq)==0) then
            UseMeasure(19) = .true.
            M19_Val = 0.d0
        elseif (M19_CycleFreq > 0 .and. mod(CycleNum, M19_CycleFreq)==0) then
            UseMeasure(19) = .true.
            M19_Val = 0.d0
        else
            UseMeasure(19) = .false.
        endif
    else
        UseMeasure(19) = .false.
    endif
    !measure FluctA0
    if (M20_Active) then
        if (MeasureAll) then
            UseMeasure(20) = .true.
            M20_Val = 0.d0
        elseif (M20_StepFreq > 0 .and. mod(StepNum, M20_StepFreq)==0) then
            UseMeasure(20) = .true.
            M20_Val = 0.d0
        elseif (M20_CycleFreq > 0 .and. mod(CycleNum, M20_CycleFreq)==0) then
            UseMeasure(20) = .true.
            M20_Val = 0.d0
        else
            UseMeasure(20) = .false.
        endif
    else
        UseMeasure(20) = .false.
    endif

    if (UseMeasure(0)) then
        !measure KEnergy
        M0_Val = KEnergy
        Val0 = M0_Val(0)
        M0_ValSum(0) = M0_ValSum(0) + Val0*Weight
        M0_ValSumSq(0) = M0_ValSumSq(0) + Val0*Val0*Weight
        M0_Count = M0_Count + Weight
    end if

    if (UseMeasure(1)) then
        !measure PEnergy
        M1_Val = PEnergy
        Val0 = M1_Val(0)
        M1_ValSum(0) = M1_ValSum(0) + Val0*Weight
        M1_ValSumSq(0) = M1_ValSumSq(0) + Val0*Val0*Weight
        M1_Count = M1_Count + Weight
    end if

    if (UseMeasure(2)) then
        !measure TEnergy
        M2_Val = TEnergy
        Val0 = M2_Val(0)
        M2_ValSum(0) = M2_ValSum(0) + Val0*Weight
        M2_ValSumSq(0) = M2_ValSumSq(0) + Val0*Val0*Weight
        M2_Count = M2_Count + Weight
    end if

    if (UseMeasure(3)) then
        !measure KTemp
        if (merge(NDOF-Dim, NDOF, ConservesMomentum) <= 0) then
            M3_Val = 0.
        else
            M3_Val = 2.d0 * KEnergy / (kB * merge(NDOF-Dim, NDOF, ConservesMomentum))
        endif
        Val0 = M3_Val(0)
        M3_ValSum(0) = M3_ValSum(0) + Val0*Weight
        M3_ValSumSq(0) = M3_ValSumSq(0) + Val0*Val0*Weight
        M3_Count = M3_Count + Weight
    end if

    if (UseMeasure(4)) then
        !measure Vol
        M4_Val = product(BoxL)
        Val0 = M4_Val(0)
        M4_ValSum(0) = M4_ValSum(0) + Val0*Weight
        M4_ValSumSq(0) = M4_ValSumSq(0) + Val0*Val0*Weight
        M4_Count = M4_Count + Weight
    end if

    if (UseMeasure(5)) then
        !measure Pressure
        M5_Val = (kB*NDOF*TempSet - Virial) / (Dim * product(BoxL))
        Val0 = M5_Val(0)
        M5_ValSum(0) = M5_ValSum(0) + Val0*Weight
        M5_ValSumSq(0) = M5_ValSumSq(0) + Val0*Val0*Weight
        M5_Count = M5_Count + Weight
    end if

    if (UseMeasure(6)) then
        !measure Virial
        M6_Val = Virial
        Val0 = M6_Val(0)
        M6_ValSum(0) = M6_ValSum(0) + Val0*Weight
        M6_ValSumSq(0) = M6_ValSumSq(0) + Val0*Val0*Weight
        M6_Count = M6_Count + Weight
    end if

    if (UseMeasure(7)) then
        !measure N_A
        M7_Val = float(NActiveMID(0))
        Val0 = M7_Val(0)
        M7_ValSum(0) = M7_ValSum(0) + Val0*Weight
        M7_ValSumSq(0) = M7_ValSumSq(0) + Val0*Val0*Weight
        M7_Count = M7_Count + Weight
    end if

    if (UseMeasure(8)) then
        !measure N_B
        M8_Val = float(NActiveMID(1))
        Val0 = M8_Val(0)
        M8_ValSum(0) = M8_ValSum(0) + Val0*Weight
        M8_ValSumSq(0) = M8_ValSumSq(0) + Val0*Val0*Weight
        M8_Count = M8_Count + Weight
    end if

    if (UseMeasure(9)) then
        !measure r_A
        M9_Val = float(NActiveMID(0)) / product(BoxL)
        Val0 = M9_Val(0)
        M9_ValSum(0) = M9_ValSum(0) + Val0*Weight
        M9_ValSumSq(0) = M9_ValSumSq(0) + Val0*Val0*Weight
        M9_Count = M9_Count + Weight
    end if

    if (UseMeasure(10)) then
        !measure r_B
        M10_Val = float(NActiveMID(1)) / product(BoxL)
        Val0 = M10_Val(0)
        M10_ValSum(0) = M10_ValSum(0) + Val0*Weight
        M10_ValSumSq(0) = M10_ValSumSq(0) + Val0*Val0*Weight
        M10_Count = M10_Count + Weight
    end if

    if (UseMeasure(11)) then
        !measure x_A
        M11_Val = float(NActiveMID(0)) / sum(NActiveMID)
        Val0 = M11_Val(0)
        M11_ValSum(0) = M11_ValSum(0) + Val0*Weight
        M11_ValSumSq(0) = M11_ValSumSq(0) + Val0*Val0*Weight
        M11_Count = M11_Count + Weight
    end if

    if (UseMeasure(12)) then
        !measure x_B
        M12_Val = float(NActiveMID(1)) / sum(NActiveMID)
        Val0 = M12_Val(0)
        M12_ValSum(0) = M12_ValSum(0) + Val0*Weight
        M12_ValSumSq(0) = M12_ValSumSq(0) + Val0*Val0*Weight
        M12_Count = M12_Count + Weight
    end if

    if (UseMeasure(13)) then
        !measure DUParam
        M13_Val = DUParam
        M13_ValSum = M13_ValSum + (M13_Val)*Weight
        do v1 = 0, 19 - 1
            do v2 = 0, 19 - 1
                M13_ValSumSq(v1,v2) = M13_ValSumSq(v1,v2) + M13_Val(v1)*M13_Val(v2)*Weight
            enddo
        enddo
        M13_Count = M13_Count + Weight
    end if

    if (UseMeasure(14)) then
        !measure DDUParam
        M14_Val = DDUParam
        M14_ValSum = M14_ValSum + (M14_Val)*Weight
        M14_ValSumSq = M14_ValSumSq + (M14_Val)*(M14_Val)*Weight
        M14_Count = M14_Count + Weight
    end if

    if (UseMeasure(15)) then
        !measure DWParam
        M15_Val = DWParam
        M15_ValSum = M15_ValSum + (M15_Val)*Weight
        M15_ValSumSq = M15_ValSumSq + (M15_Val)*(M15_Val)*Weight
        M15_Count = M15_Count + Weight
    end if

    if (UseMeasure(16)) then
        !measure DDWParam
        M16_Val = DDWParam
        M16_ValSum = M16_ValSum + (M16_Val)*Weight
        M16_ValSumSq = M16_ValSumSq + (M16_Val)*(M16_Val)*Weight
        M16_Count = M16_Count + Weight
    end if

    if (UseMeasure(17)) then
        !measure DUParamDWParam
        do i = 0, NDParam - 1
            do j = 0, NDParam - 1
                M17_Val(i*NDParam + j) = DUParam(i) * DWParam(j)
            enddo
        enddo
        M17_ValSum = M17_ValSum + (M17_Val)*Weight
        M17_ValSumSq = M17_ValSumSq + (M17_Val)*(M17_Val)*Weight
        M17_Count = M17_Count + Weight
    end if

    if (UseMeasure(18)) then
        !measure FluctTerm
        M18_Val = FluctTerm
        Val0 = M18_Val(0)
        M18_ValSum(0) = M18_ValSum(0) + Val0*Weight
        M18_ValSumSq(0) = M18_ValSumSq(0) + Val0*Val0*Weight
        M18_Count = M18_Count + Weight
    end if

    if (UseMeasure(19)) then
        !measure FluctE0
        M19_Val = FluctE0
        Val0 = M19_Val(0)
        M19_ValSum(0) = M19_ValSum(0) + Val0*Weight
        M19_ValSumSq(0) = M19_ValSumSq(0) + Val0*Val0*Weight
        M19_Count = M19_Count + Weight
    end if

    if (UseMeasure(20)) then
        !measure FluctA0
        M20_Val = FluctA0
        Val0 = M20_Val(0)
        M20_ValSum(0) = M20_ValSum(0) + Val0*Weight
        M20_ValSumSq(0) = M20_ValSumSq(0) + Val0*Val0*Weight
        M20_Count = M20_Count + Weight
    end if

end subroutine


subroutine vvquench(NSteps, CalcVirial, CalcDUParam, CalcDWParam, CalcFluct)
    implicit none
    integer, intent(in) :: NSteps
    logical, intent(in) :: CalcVirial
    logical, intent(in) :: CalcDUParam
    logical, intent(in) :: CalcDWParam
    logical, intent(in) :: CalcFluct
    real(8) :: dtsq2
    real(8) :: dt2
    real(8) :: idt
    real(8), dimension(0:Dim-1) :: Accel
    integer :: i
    integer :: Step
    logical :: ThisCalcVirial
    logical :: ThisCalcDUParam
    logical :: ThisCalcDWParam

    !velocity verlet quench 1
    dtsq2 = VVQ_TimeStep*VVQ_TimeStep*0.5
    dt2 = 0.5*VVQ_TimeStep
    idt = 1./VVQ_TimeStep

    do Step = 0, NSteps-1
        KEnergy = 0.d0
        TEnergy = KEnergy + PEnergy
        ThisCalcVirial= (CalcVirial .and. Step == NSteps - 1)
        ThisCalcDUParam = (CalcDUParam .and. Step == NSteps - 1)
        ThisCalcDWParam = (CalcDWParam .and. Step == NSteps - 1)

        do i = 0, NAtom-1
            if (.not. MolActive(MInd(i))==1) cycle
            Accel = Force(i,:) * iMass(i)
            Vel(i,:) = 0.d0
            Pos(i,:) = Pos(i,:) + dtsq2*Accel
        enddo

        !no rattle for this system

        call calcenergyforces(0, .true., ThisCalcVirial, ThisCalcDUParam, ThisCalcDWParam, CalcFluct)
    enddo

end subroutine


subroutine vvupdatekenergy()
    implicit none



end subroutine


subroutine vvintegrate(NSteps, CalcVirial, CalcDUParam, CalcDWParam, CalcFluct)
    implicit none
    integer, intent(in) :: NSteps
    logical, intent(in) :: CalcVirial
    logical, intent(in) :: CalcDUParam
    logical, intent(in) :: CalcDWParam
    logical, intent(in) :: CalcFluct
    real(8) :: dtsq2
    real(8) :: dt2
    real(8) :: idt
    real(8), dimension(0:Dim-1) :: Accel
    integer :: i
    integer :: j
    integer :: istart
    integer :: istop
    integer :: ind
    integer :: Step
    integer :: m
    real(8) :: NH_wdti0
    real(8) :: NH_wdti1
    real(8) :: NH_wdti2
    real(8) :: NH_wdti3
    real(8) :: NH_kT
    real(8) :: NH_NkT
    real(8) :: NH_scale
    real(8) :: AA
    real(8) :: NH_akin
    integer :: inos
    real(8) :: rn
    real(8), dimension(0:Dim-1) :: ranvec
    real(8) :: dtfreq
    real(8) :: sqrtkt
    real(8) :: langevin1
    real(8) :: langevin2
    real(8), parameter :: randomvelclip = -1.d0
    logical :: ThisCalcVirial
    logical :: ThisCalcDUParam
    logical :: ThisCalcDWParam
    integer :: NActiveAxes
    real(8) :: Beta
    real(8), dimension(0:Dim-1) :: DeltaPos
    real(8), dimension(0:Dim-1) :: CentroidPos
    real(8) :: DeltaVol
    real(8), dimension(0:Dim-1) :: ScaleFactor
    real(8) :: OldVol
    real(8) :: NewVol
    real(8) :: OldE
    real(8) :: NewE
    real(8) :: r
    real(8) :: lnP
    real(8), dimension(0:Dim-1) :: OldBoxL
    logical :: Acc
    logical, dimension(0:Dim-1) :: AxisMask
    real(8) :: Scale

    dtsq2 = VV_TimeStep*VV_TimeStep*0.5
    dt2 = 0.5*VV_TimeStep
    idt = 1./VV_TimeStep

    do Step = 0, NSteps-1

        ThisCalcVirial = (CalcVirial .and. Step == NSteps - 1)
        ThisCalcDUParam = (CalcDUParam .and. Step == NSteps - 1)
        ThisCalcDWParam = (CalcDWParam .and. Step == NSteps - 1)

        if (VV_Barostat == 1 .and. mod(Step, VV_BarostatStepFreq)==0) then
            Beta = 1. / (TempSet * kB)
            !update the attempt
            VV_BarostatNAtt = VV_BarostatNAtt + 1.

            NActiveAxes = count(VV_BarostatUseAxis)
            OldBoxL = BoxL
            OldVol = product(BoxL)
            OldE = PEnergy
            call saveenergystate(0)

            !check if we need to scale independently
            if (.not. VV_BarostatDoIsotropic) then
                !find a random active axis
                call ran2int(NActiveAxes, i)
                j = -1
                do ind = 0, Dim - 1
                    if (VV_BarostatUseAxis(ind)) j = j + 1
                    if (j == i) exit
                enddo
                AxisMask = .false.
                AxisMask(ind) = .true.
            else
                AxisMask = VV_BarostatUseAxis
            endif

            !choose a random volume change
            call ran2(r)
            DeltaVol = VV_BarostatDeltaV * (2.d0 * r - 1.d0)
            NewVol = OldVol + DeltaVol

            if (NewVol > 0. .and. NewVol >= VV_BarostatMinVol .and. NewVol <= VV_BarostatMaxVol) then

                ScaleFactor = merge((NewVol / OldVol)**(1.d0 / dble(count(AxisMask))), 1.d0, AxisMask)
                BoxL = BoxL * ScaleFactor
                OldPos = Pos

                !now scale the molecule centers of mass
                do m = 0, NMol - 1
                    !skip frozen and inactive mols
                    if (MolActive(m) < 1) cycle

                    !find the current centroid
                    istart = MolRange(m)
                    istop = MolRange(m+1) - 1
                    CentroidPos = sum(Pos(istart:istop,:), dim=1) / AtomsPerMol(MolID(m))

                    !find displacement
                    DeltaPos = CentroidPos * (ScaleFactor - 1.d0)

                    !update atom positions
                    do i = istart, istop
                        Pos(i,:) = Pos(i,:) + DeltaPos
                    enddo
                enddo

                !update energy
                call calcenergyforces(0, .false., .false., .false., .false., CalcFluct)
                NewE = PEnergy

                lnP = Beta * (OldE - NewE) - Beta * PresSet * DeltaVol
                lnP = lnP + NActiveMol * log(NewVol / OldVol)
                if (lnP >= 0) then
                    Acc = .true.
                else
                    call ran2(r)
                    Acc = (exp(lnP) > r)
                endif

                if (Acc) then
                    VV_BarostatNAcc = VV_BarostatNAcc + 1.
                else
                    BoxL = OldBoxL
                    Pos = OldPos
                    call revertenergystate(0)
                endif
            endif
        endif

        !Andersen thermostats
        if (VV_Thermostat == 1) then
            !do an andersen massive collision update
            if (mod(VV_AndersenStep, VV_AndersenStepFreq)==0) then
                VV_AndersenStep = 0
                sqrtkt = sqrt(kB * TempSet)
                do i = 0, NAtom-1
                    if (.not. MolActive(MInd(i))==1) cycle
                    call ran2normarray(Dim, ranvec)
                    !clip the limits on the random variate for stability
                    if (randomvelclip > 0.d0) then
                        ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
                    endif
                    ranvec = ranvec * sqrtkt / sqrtMass(i)
                    Vel(i,:) = ranvec
                enddo
            endif
            VV_AndersenStep = VV_AndersenStep + 1
        endif

        if (VV_Thermostat == 2) then
            !do an andersen particle collision update
            dtfreq = VV_TimeStep * VV_AndersenCollisionFreq
            sqrtkt = sqrt(kB * TempSet)
            do m = 0, NMol-1
                if (.not. MolActive(m)==1) cycle
                call ran2(rn)
                if (rn < dtfreq) then
                    do i = MolRange(m), MolRange(m+1)-1
                        call ran2normarray(Dim, ranvec)
                        !clip the limits on the random variate for stability
                        if (randomvelclip > 0.d0) then
                            ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
                        endif
                        ranvec = ranvec * sqrtkt / sqrtMass(i)
                        Vel(i,:) = ranvec
                    enddo
                endif
            enddo
        endif

        !remove the center of mass
        if (VV_RemoveCOMStepFreq > 0) then
            if (mod(VV_RemoveCOMStep, VV_RemoveCOMStepFreq)==0) then
                VV_RemoveCOMStep = 0
                ranvec = 0.d0
                do i = 0, NAtom-1
                    if (.not. MolActive(MInd(i))==1) cycle
                    ranvec = ranvec + Mass(i) * Vel(i,:)
                enddo
                ranvec = ranvec / NAtom
                do i = 0, NAtom-1
                    if (.not. MolActive(MInd(i))==1) cycle
                    Vel(i,:) = Vel(i,:) - ranvec * iMass(i)
                enddo
            endif
            VV_RemoveCOMStep = VV_RemoveCOMStep + 1
        endif

        !NOSE-HOOVER ROUTINES
        if (VV_Thermostat == 3) then
            !update kinetic energy
            KEnergy = 0.d0
            do i = 0, NAtom-1
                KEnergy = KEnergy + dot_product(Vel(i,:), Vel(i,:)) * Mass(i)
            enddo
            KEnergy = KEnergy * 0.5d0
            TEnergy = KEnergy + PEnergy
            !set frequently used variables
            NH_wdti0 = VV_TimeStep
            NH_wdti1 = VV_TimeStep * 0.5d0
            NH_wdti2 = VV_TimeStep * 0.25d0
            NH_wdti3 = VV_TimeStep * 0.125
            NH_kT = TempSet * kB
            NH_NkT = NH_kT * dble(NDOF - Dim)
            NH_scale = 1.D0
            !get kinetic energy
            NH_akin = 2.d0 * KEnergy
            !update the forces
            VV_Glogs(0) = (NH_akin - NH_NkT) / VV_QMass(0)
            !update the VV_Thermostat velocities
            VV_Vlogs(VV_NH_N-1) = VV_Vlogs(VV_NH_N-1) + VV_Glogs(VV_NH_N-1) * NH_wdti2
            do inos = 1, VV_NH_N - 1
                AA = exp( -NH_wdti3 * VV_Vlogs(VV_NH_N-inos) )
                VV_Vlogs(VV_NH_N-inos-1) = VV_Vlogs(VV_NH_N-inos-1)*AA*AA + NH_wdti2*VV_Glogs(VV_NH_N-inos-1)*AA
            enddo
            !update the particle velocities
            AA = exp( -NH_wdti1 * VV_Vlogs(0) )
            NH_scale = NH_scale * AA
            !update the forces
            VV_Glogs(0) = (NH_scale*NH_scale*NH_akin - NH_NkT) / VV_QMass(0)
            !update the VV_Thermostat positions
            do inos = 0, VV_NH_N - 1
                VV_Xlogs(inos) = VV_Xlogs(inos) + VV_Vlogs(inos) * NH_wdti1
            enddo
            !update the VV_Thermostat velocities
            do inos = 1, VV_NH_N-1
                AA = exp( -NH_wdti3 * VV_Vlogs(inos) )
                VV_Vlogs(inos-1) = VV_Vlogs(inos-1)*AA*AA + NH_wdti2*VV_Glogs(inos-1)*AA
                VV_Glogs(inos) = (VV_QMass(inos-1)*VV_Vlogs(inos-1)*VV_Vlogs(inos-1)-NH_kT) / VV_QMass(inos)
            enddo
            VV_Vlogs(VV_NH_N-1) = VV_Vlogs(VV_NH_N-1) + NH_wdti2*VV_Glogs(VV_NH_N-1)
            !update the particle velocities
            if (NH_scale > 0.) then
                Vel = Vel * NH_scale
            endif
        endif

        !velocity verlet integration part 1
        if (VV_Thermostat == 4) then
            !langevin VV_Thermostat
            !based on Bussi and Parrinello, Physical Review E 75, 056707, 2007
            langevin1 = exp(-0.5d0 * VV_LangevinGamma * VV_TimeStep)
            langevin2 = sqrt((1.d0 - langevin1**2) * (kB * TempSet))
            do i = 0, NAtom-1
                if (.not. MolActive(MInd(i))==1) cycle
                Accel = Force(i,:) * iMass(i)
                call ran2normarray(Dim, ranvec)
                !clip the limits on the random variate for stability
                if (randomvelclip > 0.d0) then
                    ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
                endif
                Vel(i,:) = langevin1*Vel(i,:) + langevin2*SqrtMass(i)*iMass(i)*ranvec
                Pos(i,:) = Pos(i,:) + VV_TimeStep*Vel(i,:) + dtsq2*Accel
                Vel(i,:) = Vel(i,:) + dt2*Accel
            enddo
        else
            !normal constant energy dynamics
            do i = 0, NAtom-1
                if (.not. MolActive(MInd(i))==1) cycle
                Accel = Force(i,:) * iMass(i)
                Pos(i,:) = Pos(i,:) + VV_TimeStep*Vel(i,:) + dtsq2*Accel
                Vel(i,:) = Vel(i,:) + dt2*Accel
            enddo
        endif

        !no rattle for this system

        call calcenergyforces(0, .true., ThisCalcVirial, ThisCalcDUParam, ThisCalcDWParam, CalcFluct)

        !velocity verlet integration part 2

        if (VV_Thermostat == 4) then
            !langevin VV_Thermostat
            langevin1 = exp(-0.5d0 * VV_LangevinGamma * VV_TimeStep)
            langevin2 = sqrt((1.d0 - langevin1**2) * (kB * TempSet))
            do i = 0, NAtom-1
                if (.not. MolActive(MInd(i))==1) cycle
                call ran2normarray(Dim, ranvec)
                !clip the limits on the random variate for stability
                if (randomvelclip > 0.d0) then
                    ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
                endif
                Accel = Force(i,:) * iMass(i)
                Vel(i,:) = Vel(i,:) + dt2*Accel
                Vel(i,:) = langevin1*Vel(i,:) + langevin2*SqrtMass(i)*iMass(i)*ranvec
            enddo
        else
            !normal constant energy dynamics
            do i = 0, NAtom-1
                if (.not. MolActive(MInd(i))==1) cycle
                Accel = Force(i,:) * iMass(i)
                Vel(i,:) = Vel(i,:) + dt2*Accel
            enddo
        endif

        !no rattle for this system

        !NOSE-HOOVER ROUTINES
        if (VV_Thermostat == 3) then
            !update kinetic energy
            KEnergy = 0.d0
            do i = 0, NAtom-1
                KEnergy = KEnergy + dot_product(Vel(i,:), Vel(i,:)) * Mass(i)
            enddo
            KEnergy = KEnergy * 0.5d0
            TEnergy = KEnergy + PEnergy
            !set frequently used variables
            NH_wdti0 = VV_TimeStep
            NH_wdti1 = VV_TimeStep * 0.5d0
            NH_wdti2 = VV_TimeStep * 0.25d0
            NH_wdti3 = VV_TimeStep * 0.125
            NH_kT = TempSet * kB
            NH_NkT = NH_kT * dble(NDOF - Dim)
            NH_scale = 1.D0
            !get kinetic energy
            NH_akin = 2.d0 * KEnergy
            !update the forces
            VV_Glogs(0) = (NH_akin - NH_NkT) / VV_QMass(0)
            !update the VV_Thermostat velocities
            VV_Vlogs(VV_NH_N-1) = VV_Vlogs(VV_NH_N-1) + VV_Glogs(VV_NH_N-1) * NH_wdti2
            do inos = 1, VV_NH_N - 1
                AA = exp( -NH_wdti3 * VV_Vlogs(VV_NH_N-inos) )
                VV_Vlogs(VV_NH_N-inos-1) = VV_Vlogs(VV_NH_N-inos-1)*AA*AA + NH_wdti2*VV_Glogs(VV_NH_N-inos-1)*AA
            enddo
            !update the particle velocities
            AA = exp( -NH_wdti1 * VV_Vlogs(0) )
            NH_scale = NH_scale * AA
            !update the forces
            VV_Glogs(0) = (NH_scale*NH_scale*NH_akin - NH_NkT) / VV_QMass(0)
            !update the VV_Thermostat positions
            do inos = 0, VV_NH_N - 1
                VV_Xlogs(inos) = VV_Xlogs(inos) + VV_Vlogs(inos) * NH_wdti1
            enddo
            !update the VV_Thermostat velocities
            do inos = 1, VV_NH_N-1
                AA = exp( -NH_wdti3 * VV_Vlogs(inos) )
                VV_Vlogs(inos-1) = VV_Vlogs(inos-1)*AA*AA + NH_wdti2*VV_Glogs(inos-1)*AA
                VV_Glogs(inos) = (VV_QMass(inos-1)*VV_Vlogs(inos-1)*VV_Vlogs(inos-1)-NH_kT) / VV_QMass(inos)
            enddo
            VV_Vlogs(VV_NH_N-1) = VV_Vlogs(VV_NH_N-1) + NH_wdti2*VV_Glogs(VV_NH_N-1)
            !update the particle velocities
            if (NH_scale > 0.) then
                Vel = Vel * NH_scale
            endif
        endif

        !update kinetic energy
        KEnergy = 0.d0
        do i = 0, NAtom-1
            if (.not. MolActive(MInd(i))==1) cycle
            KEnergy = KEnergy + dot_product(Vel(i,:), Vel(i,:)) * Mass(i)
        enddo
        KEnergy = KEnergy * 0.5d0
        TEnergy = KEnergy + PEnergy

        !update running sums for total energy
        VV_TEnergySum = VV_TEnergySum + TEnergy
        VV_TEnergySqSum = VV_TEnergySqSum + TEnergy*TEnergy

    enddo

end subroutine


subroutine montecarlocycle(NSteps, CalcForce, CalcVirial, CalcDUParam, CalcDWParam, CalcFluct)
    implicit none
    integer, intent(in) :: NSteps
    logical, intent(in) :: CalcForce
    logical, intent(in) :: CalcVirial
    logical, intent(in) :: CalcDUParam
    logical, intent(in) :: CalcDWParam
    logical, intent(in) :: CalcFluct
    integer :: StepNum
    real(8) :: mcmoverannum
    real(8) :: Beta
    real(8), dimension(0:Dim-1) :: iBoxL
    integer :: i
    integer :: istart
    integer :: istop
    integer :: m
    integer :: ActiveInd
    integer :: ind
    integer :: ThisAtomsPerMol
    real(8), dimension(0:Dim-1) :: OldAtomPos
    real(8), dimension(0:Dim-1) :: dPos
    real(8), dimension(0:Dim-1) :: CentroidPos
    real(8) :: OldE
    real(8) :: NewE
    real(8) :: r
    real(8) :: lnP
    real(8), dimension(0:Dim-1,0:Dim-1) :: RotMat
    logical :: Acc
    logical :: Rotate
    integer :: ThisMID
    real(8), dimension(0:Dim-1) :: NewPos
    real(8) :: Pacc0
    real(8) :: ProbIns
    real(8) :: ProbDel
    logical :: DoInsert
    integer :: j
    integer :: NActiveAxes
    real(8), dimension(0:Dim-1) :: DeltaPos
    real(8) :: DeltaVol
    real(8), dimension(0:Dim-1) :: ScaleFactor
    real(8) :: OldVol
    real(8) :: NewVol
    real(8), dimension(0:Dim-1) :: OldBoxL
    logical, dimension(0:Dim-1) :: AxisMask
    integer :: DeleteMID
    integer :: DeleteMol
    integer :: DeleteAtomsPerMol
    integer :: InsertMol
    integer :: InsertMID
    integer :: InsertAtomsPerMol
    integer :: Delta1
    integer :: N1
    real(8), dimension(0:Dim-1) :: DeleteCentroidPos
    real(8), dimension(0:Dim-1) :: InsertCentroidPos
    real(8) :: Scale

    Beta = 1.d0 / (kB * TempSet)
    iBoxL = 1.d0 / BoxL
    KEnergy = 0.d0
    TEnergy = 0.d0

    do StepNum = 0, NSteps - 1

        !draw a random number to decide which move to do
        call ran2(mcmoverannum)

        if (mcmoverannum <= MCMoveProbs(0)) then

            !check if active
            if (NActiveMol == 0) then
                MC0_NAtt = MC0_NAtt + 1
                return
            endif

            !pick a random molecule that's active
            call ran2int(NActiveMol, ind)
            ActiveInd = -1
            TargetMol = -1
            do while (ActiveInd < ind)
                TargetMol = TargetMol + 1
                if (MolActive(TargetMol) == 1) ActiveInd = ActiveInd + 1
            enddo
            ThisAtomsPerMol = AtomsPerMol(MolID(TargetMol))

            !check if this is a rigid molecule or not
            if (MolIsRigid(TargetMol)) then

                !rigid molecule; find start and stop atoms
                istart = MolRange(TargetMol)
                istop = MolRange(TargetMol+1) - 1

                !save the old information
                OldE = PEnergy
                OldPos(istart:istop,:) = Pos(istart:istop,:)
                call saveenergystate(2)

                !calculate energy without molecule
                call calcenergyforces(-2, .false., .false., .false., .false., CalcFluct)

                !decide whether to move or rotate
                call ran2(r)
                Rotate = (r > 0.5d0)

                if (Rotate) then

                    MC0_NAtt2 = MC0_NAtt2 + 1
                    !rotate about a random axis
                    call RandomRotMatxyz(MC0_Delta2, RotMat)
                    CentroidPos = sum(Pos(istart:istop,:), dim=1) / ThisAtomsPerMol
                    do i = istart, istop
                        Pos(i,:) = matmul(RotMat, Pos(i,:) - CentroidPos) + CentroidPos
                    enddo

                else

                    MC0_NAtt = MC0_NAtt + 1
                    !translation
                    call ran2array(Dim, dPos)
                    dPos = MC0_Delta * (2.d0 * dPos - 1.d0)
                    do i = istart, istop
                        Pos(i,:) = Pos(i,:) + dPos
                    enddo

                endif

                !calculate energy with molecule
                call calcenergyforces(+2, .false., .false., .false., .false., CalcFluct)
                NewE = PEnergy

                lnP = Beta * (OldE - NewE)
                if (lnP >= 0) then
                    Acc = .true.
                else
                    call ran2(r)
                    Acc = (exp(lnP) > r)
                endif
                if (Acc) then
                    if (Rotate) then
                        MC0_NAcc2 = MC0_NAcc2 + 1
                    else
                        MC0_NAcc = MC0_NAcc + 1
                    endif
                else
                    !revert to old state
                    Pos(istart:istop,:) = OldPos(istart:istop,:)
                    call revertenergystate(2)
                endif

            else

                !this is a non-rigid molecule; pick a random atom to displace
                call ran2int(ThisAtomsPerMol, i)
                TargetAtom = i + MolRange(TargetMol)

                !save the old information
                OldE = PEnergy
                OldAtomPos = Pos(TargetAtom,:)
                call saveenergystate(1)

                !calculate energy without atom
                call calcenergyforces(-1, .false., .false., .false., .false., CalcFluct)

                MC0_NAtt = MC0_NAtt + 1

                !random displacement
                call ran2array(Dim, dPos)
                dPos = MC0_Delta * (2.d0 * dPos - 1.d0)
                Pos(TargetAtom,:) = Pos(TargetAtom,:) + dPos

                !calculate energy with atom
                call calcenergyforces(+1, .false., .false., .false., .false., CalcFluct)
                NewE = PEnergy

                lnP = Beta * (OldE - NewE)
                if (lnP >= 0) then
                    Acc = .true.
                else
                    call ran2(r)
                    Acc = (exp(lnP) > r)
                endif
                if (Acc) then
                    MC0_NAcc = MC0_NAcc + 1
                else
                    !revert to old state
                    Pos(TargetAtom,:) = OldAtomPos
                    call revertenergystate(1)
                endif

            endif

        elseif (mcmoverannum <= MCMoveProbs(1)) then

            !choose whether to insert or delete
            call ran2(r)
            DoInsert = (r > 0.5d0)

            !update the attempt
            MC1_NAtt = MC1_NAtt + 1

            !do an insertion
            if (DoInsert) then

                !make sure there are available molecules to insert
                if (MC1_NInsertMols > 0) then

                    call ran2int(MC1_NInsertMols, ind)
                    TargetMol = MC1_MolMoves(ind)
                    ThisMID = MolID(TargetMol)
                    ThisAtomsPerMol = AtomsPerMol(ThisMID)
                    istart = MolRange(TargetMol)
                    istop = MolRange(TargetMol+1) - 1

                    MC1_NAttMID(ThisMID) = MC1_NAttMID(ThisMID) + 1
                    OldE = PEnergy
                    call saveenergystate(2)

                    !pick a random location
                    call ran2array(Dim, NewPos)
                    NewPos = (NewPos - 0.5d0) * BoxL
                    if (ThisAtomsPerMol > 1) then
                        !pick a random rotation
                        CentroidPos = sum(Pos(istart:istop,:), dim=1) / ThisAtomsPerMol
                        call RandomRotMat3D(3.1415926535897931D0, RotMat)
                        do i = istart, istop
                            Pos(i,:) = matmul(RotMat, Pos(i,:) - CentroidPos) + NewPos
                        enddo
                    else
                        Pos(istart,:) = NewPos
                    endif

                    !update energy for adding this molecule
                    call calcenergyforces(+2, .false., .false., .false., .false., CalcFluct)
                    NewE = PEnergy

                    lnP = Beta * (OldE - NewE) + Beta * MuSet(ThisMID)
                    ProbIns = dble(NInactiveMID(ThisMID)) / dble(MC1_NInsertMols) / product(BoxL)
                    ProbDel = 1.d0 / dble(MC1_NDeleteMols + 1)
                    lnP = lnP + log(ProbDel / ProbIns)
                    Pacc0 = exp(min(0.d0, lnP))
                    lnP = lnP + MC1_BoltzWeights(ThisMID, NActiveMID(ThisMID)+1) - MC1_BoltzWeights(ThisMID, NActiveMID(ThisMID))

                    !update transition matrix
                    MC1_TM(NActiveMol, 1) = MC1_TM(NActiveMol, 1) + (1.d0 - Pacc0)
                    MC1_TM(NActiveMol, 2) = MC1_TM(NActiveMol, 2) + Pacc0

                    if (lnP >= 0) then
                        Acc = .true.
                    else
                        call ran2(r)
                        Acc = (exp(lnP) > r)
                    endif

                    if (Acc) then
                        MC1_NAcc = MC1_NAcc + 1
                        NActiveMol = NActiveMol + 1
                        NActiveMID(ThisMID) = NActiveMID(ThisMID) + 1
                        NInactiveMol = NInactiveMol - 1
                        NInactiveMID(ThisMID) = NInactiveMID(ThisMID) - 1
                        MolActive(TargetMol) = 1
                        MC1_NAccMID(ThisMID) = MC1_NAccMID(ThisMID) + 1
                        !swap the indices of the inserted atom the last possible inserted
                        MC1_MolMoves(ind) = MC1_MolMoves(MC1_NInsertMols - 1)
                        MC1_MolMoves(MC1_NInsertMols - 1) = TargetMol
                        MC1_NInsertMols = MC1_NInsertMols - 1
                        MC1_NDeleteMols = MC1_NDeleteMols + 1
                    else
                        MolActive(TargetMol) = -1
                        call revertenergystate(2)
                    endif

                else

                    !update transition matrix
                    MC1_TM(NActiveMol, 1) = MC1_TM(NActiveMol, 1) + 1.d0

                endif

            !do deletion
            else

                !make sure there are available molecules to delete
                if (MC1_NDeleteMols > 0) then

                    call ran2int(MC1_NDeleteMols, ind)
                    ind = ind + MC1_NInsertMols
                    TargetMol = MC1_MolMoves(ind)
                    ThisMID = MolID(TargetMol)
                    istart = MolRange(TargetMol)
                    istop = MolRange(TargetMol+1) - 1

                    MC1_NAttMID(ThisMID) = MC1_NAttMID(ThisMID) + 1
                    OldE = PEnergy
                    call saveenergystate(2)

                    !update energy for deleting this molecule
                    call calcenergyforces(-2, .false., .false., .false., .false., CalcFluct)
                    NewE = PEnergy

                    lnP = Beta * (OldE - NewE) - Beta * MuSet(ThisMID)
                    ProbDel = 1.d0 / dble(MC1_NDeleteMols)
                    ProbIns = dble(NInactiveMID(ThisMID) + 1) / dble(MC1_NInsertMols+1) / product(BoxL)
                    lnP = lnP + log(ProbIns / ProbDel)
                    Pacc0 = exp(min(0.d0, lnP))
                    lnP = lnP + MC1_BoltzWeights(ThisMID, NActiveMID(ThisMID)-1) - MC1_BoltzWeights(ThisMID, NActiveMID(ThisMID))

                    !update transition matrix
                    MC1_TM(NActiveMol, 1) = MC1_TM(NActiveMol, 1) + (1.d0 - Pacc0)
                    MC1_TM(NActiveMol, 0) = MC1_TM(NActiveMol, 0) + Pacc0

                    if (lnP >= 0) then
                        Acc = .true.
                    else
                        call ran2(r)
                        Acc = (exp(lnP) > r)
                    endif

                    if (Acc) then
                        MC1_NAcc = MC1_NAcc + 1
                        NActiveMol = NActiveMol - 1
                        NActiveMID(ThisMID) = NActiveMID(ThisMID) - 1
                        NInactiveMol = NInactiveMol + 1
                        NInactiveMID(ThisMID) = NInactiveMID(ThisMID) + 1
                        MolActive(TargetMol) = -1
                        MC1_NAccMID(ThisMID) = MC1_NAccMID(ThisMID) + 1
                        !swap the indices of the deleted atom and the first possible deleted
                        MC1_MolMoves(ind) = MC1_MolMoves(MC1_NInsertMols)
                        MC1_MolMoves(MC1_NInsertMols) = TargetMol
                        MC1_NInsertMols = MC1_NInsertMols + 1
                        MC1_NDeleteMols = MC1_NDeleteMols - 1
                    else
                        call revertenergystate(2)
                    endif

                else

                    !update transition matrix
                    MC1_TM(NActiveMol, 1) = MC1_TM(NActiveMol, 1) + 1.d0

                endif

            endif

            !update weights
            do i = 0, NMID-1
                MC1_BoltzWeights(i, NActiveMID(i)) = MC1_BoltzWeights(i, NActiveMID(i)) + MC1_MFactor
            enddo

        elseif (mcmoverannum <= MCMoveProbs(2)) then

            !update the attempt
            MC2_NAtt = MC2_NAtt + 1

            NActiveAxes = count(MC2_UseAxis)
            OldBoxL = BoxL
            OldVol = product(BoxL)
            OldE = PEnergy
            call saveenergystate(0)

            !check if we need to scale independently
            if (.not. MC2_DoIsotropic) then
                !find a random active axis
                call ran2int(NActiveAxes, i)
                j = -1
                do ind = 0, Dim - 1
                    if (MC2_UseAxis(ind)) j = j + 1
                    if (j == i) exit
                enddo
                AxisMask = .false.
                AxisMask(ind) = .true.
            else
                AxisMask = MC2_UseAxis
            endif

            !choose a random volume change
            call ran2(r)
            DeltaVol = MC2_Delta * (2.d0 * r - 1.d0)
            NewVol = OldVol + DeltaVol

            if (NewVol > 0. .and. NewVol >= MC2_MinVol .and. NewVol <= MC2_MaxVol) then

                ScaleFactor = merge((NewVol / OldVol)**(1.d0 / dble(count(AxisMask))), 1.d0, AxisMask)
                BoxL = BoxL * ScaleFactor
                OldPos = Pos

                !now scale the molecule centers of mass
                do m = 0, NMol - 1
                    !skip frozen and inactive mols
                    if (MolActive(m) < 1) cycle

                    !find the current centroid
                    istart = MolRange(m)
                    istop = MolRange(m+1) - 1
                    CentroidPos = sum(Pos(istart:istop,:), dim=1) / AtomsPerMol(MolID(m))

                    !find displacement
                    DeltaPos = CentroidPos * (ScaleFactor - 1.d0)

                    !update atom positions
                    do i = istart, istop
                        Pos(i,:) = Pos(i,:) + DeltaPos
                    enddo
                enddo

                !update energy
                call calcenergyforces(0, .false., .false., .false., .false., CalcFluct)
                NewE = PEnergy

                lnP = Beta * (OldE - NewE) - Beta * PresSet * DeltaVol
                lnP = lnP + NActiveMol * log(NewVol / OldVol)
                if (lnP >= 0) then
                    Acc = .true.
                else
                    call ran2(r)
                    Acc = (exp(lnP) > r)
                endif

                if (Acc) then
                    MC2_NAcc = MC2_NAcc + 1
                else
                    BoxL = OldBoxL
                    Pos = OldPos
                    call revertenergystate(0)
                endif
            endif

        elseif (mcmoverannum <= MCMoveProbs(3)) then

            !update the attempt
            MC3_NAtt = MC3_NAtt + 1

            !pick a random molecule to mutate that's active and type MC3_MID1 or MC3_MID2
            call ran2int(NActiveMID(MC3_MID1) + NActiveMID(MC3_MID2), ind)
            ActiveInd = -1
            DeleteMol = -1
            do while (ActiveInd < ind)
                DeleteMol = DeleteMol + 1
                DeleteMID = MolID(DeleteMol)
                if (MolActive(DeleteMol) == 1 .and. (DeleteMID == MC3_MID1 .or. DeleteMID == MC3_MID2)) then
                    ActiveInd = ActiveInd + 1
                endif
            enddo
            DeleteAtomsPerMol = AtomsPerMol(DeleteMID)

            !find the molecule type to replace with
            if (DeleteMID == MC3_MID1 .and. NInactiveMID(MC3_MID2) > 0) then
                InsertMID = MC3_MID2
                Delta1 = -1
            elseif (DeleteMID == MC3_MID2 .and. NInactiveMID(MC3_MID1) > 0) then
                InsertMID = MC3_MID1
                Delta1 = 1
            else
                InsertMID = -1
            endif

            !order parameter for weights and MC3_TM
            N1 = NActiveMID(MC3_MID1)

            if (InsertMID >= 0) then

                !find a replacement molecule
                call ran2int(NInactiveMID(InsertMID), ind)
                ActiveInd = -1
                InsertMol = -1
                do while (ActiveInd < ind)
                    InsertMol = InsertMol + 1
                    if (MolActive(InsertMol) == -1 .and. MolID(InsertMol) == InsertMID) then
                        ActiveInd = ActiveInd + 1
                    endif
                enddo
                InsertAtomsPerMol = AtomsPerMol(InsertMID)

                !save current energy
                OldE = PEnergy
                call saveenergystate(0)

                !get the center of mass of the original molecule
                istart = MolRange(DeleteMol)
                istop = MolRange(DeleteMol+1) - 1
                DeleteCentroidPos = sum(Pos(istart:istop,:), dim=1) / DeleteAtomsPerMol

                !update energy for deleting the original
                TargetMol = DeleteMol
                MolActive(DeleteMol) = -1
                call calcenergyforces(-2, .false., .false., .false., .false., CalcFluct)

                !now add the new molecule at the same location
                istart = MolRange(InsertMol)
                istop = MolRange(InsertMol+1) - 1
                InsertCentroidPos = sum(Pos(istart:istop,:), dim=1) / InsertAtomsPerMol

                if (InsertAtomsPerMol > 1) then
                    !pick a random rotation
                    call RandomRotMat3D(3.1415926535897931D0, RotMat)
                    do i = istart, istop
                        Pos(i,:) = matmul(RotMat, Pos(i,:) - InsertCentroidPos) + DeleteCentroidPos
                    enddo
                else
                    Pos(istart,:) = DeleteCentroidPos
                endif

                !update energy for adding this molecule
                TargetMol = InsertMol
                MolActive(InsertMol) = 1
                call calcenergyforces(+2, .false., .false., .false., .false., CalcFluct)

                !get the new energy
                NewE = PEnergy

                lnP = Beta * (OldE - NewE) + Beta * (MuSet(InsertMID) - MuSet(DeleteMID))
                Pacc0 = exp(min(0.d0, lnP))
                lnP = lnP + MC3_BoltzWeights(N1 + Delta1) - MC3_BoltzWeights(N1)

                !update transition matrix
                MC3_TM(N1, 1) = MC3_TM(N1, 1) + (1.d0 - Pacc0)
                MC3_TM(N1, 1+Delta1) = MC3_TM(N1, 1+Delta1) + Pacc0

                if (lnP >= 0) then
                    Acc = .true.
                else
                    call ran2(r)
                    Acc = (exp(lnP) > r)
                endif

                if (Acc) then
                    MC3_NAcc = MC3_NAcc + 1
                    NActiveMID(MC3_MID1) = NActiveMID(MC3_MID1) + Delta1
                    NActiveMID(MC3_MID2) = NActiveMID(MC3_MID2) - Delta1
                    NInactiveMID(MC3_MID1) = NInactiveMID(MC3_MID1) - Delta1
                    NInactiveMID(MC3_MID2) = NInactiveMID(MC3_MID2) + Delta1
                    N1 = N1 + Delta1
                    MC3_BoltzWeights(N1) = MC3_BoltzWeights(N1) + MC3_MFactor
                else
                    MolActive(InsertMol) = -1
                    MolActive(DeleteMol) = 1
                    MC3_BoltzWeights(N1) = MC3_BoltzWeights(N1) + MC3_MFactor
                    call revertenergystate(0)
                endif

            else

                !update transition matrix
                MC3_TM(N1, 1) = MC3_TM(N1, 1) + 1.d0

            endif

        end if

    enddo

    call calcenergyforces(0, CalcForce, CalcVirial, CalcDUParam, CalcDWParam, CalcFluct)

end subroutine





subroutine ran2(r)
    implicit none
    real(8), intent(out) :: r
    integer, parameter :: NTAB = 32
    integer, parameter :: IM1=2147483563, IM2=2147483399, IMM1=2147483562
    integer, parameter :: IA1=40014, IA2=40692, IQ1=53668, IQ2=52774
    integer, parameter :: IR1=12211, IR2=3791, NDIV=1+IMM1/NTAB  
    real(8), parameter :: AM=1.d0/IM1, EPS=1.2d-7, RNMX=1.d0-EPS
    integer, dimension(3+32) :: state 
    COMMON /ran2data/ state
    integer :: j, k

    if (state(1).le.0) then
        state(1)=max(-state(1),1)
        state(2)=state(1)
        do j=NTAB+8,1,-1
             k=state(1)/IQ1
             state(1)=IA1*(state(1)-k*IQ1)-k*IR1
             if (state(1).lt.0) state(1)=state(1)+IM1
             if (j.le.NTAB) state(3+j)=state(1)
        enddo
        state(3)=state(4)
    endif
    k=state(1)/IQ1
    state(1)=IA1*(state(1)-k*IQ1)-k*IR1
    if (state(1).lt.0) state(1)=state(1)+IM1
    k=state(2)/IQ2
    state(2)=IA2*(state(2)-k*IQ2)-k*IR2
    if (state(2).lt.0) state(2)=state(2)+IM2
    j=1+state(3)/NDIV
    state(3)=state(3+j)-state(2)
    state(3+j)=state(1)
    if(state(3).lt.1)state(3)=state(3)+IMM1
    r=min(AM*state(3),RNMX)
end subroutine

subroutine ran2seed(seedval)
    implicit none
    integer, intent(in) :: seedval
    integer, dimension(3+32) :: state 
    COMMON /ran2data/ state
    state(1) = abs(seedval)
end subroutine

subroutine ran2array(n, rarray)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(out) :: rarray
    integer :: i
    real(8) :: r
    do i = 1, n
        call ran2(r)
        rarray(i) = r
    enddo
end subroutine

subroutine ran2int(n, i)
    !returns an integer on [0,n)
    implicit none
    integer, intent(in) :: n
    integer, intent(out) :: i
    real(8) :: r
    call ran2(r)
    i = int(real(n) * r)
    i = min(n-1, max(0, i))
end subroutine

subroutine ran2norm(r)
    implicit none
    real(8), intent(out) :: r
    real(8) :: r1, r2, rsq
    real(8), save :: rsaved = 0.1
    logical, save :: hassaved = .false.
    if (hassaved) then
        r = rsaved
        hassaved = .false.
    else
        rsq = 2.d0
        do while (rsq == 0.d0 .or. rsq >= 1.d0)
            call ran2(r1)
            call ran2(r2)
            r1 = 2.d0 * r1 - 1.d0
            r2 = 2.d0 * r2 - 1.d0
            rsq = r1*r1 + r2*r2
        enddo
        rsq = sqrt(-2.d0 * log(rsq) / rsq)
        r = r1 * rsq
        rsaved = r2 * rsq
        hassaved = .true.
    endif
end subroutine

subroutine ran2normarray(n, rarray)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(out) :: rarray
    integer :: i
    real(8) :: r
    do i = 1, n
        call ran2norm(r)
        rarray(i) = r
    enddo
end subroutine


subroutine erfc(x, value)
    real(8), intent(in) :: x
    real(8), intent(out) :: value
    !Returns the complementary error function erfc(x) with fractional 
    !error everywhere less than 1:2^10-7
    real(8) :: t, z
    z = abs(x)
    t = 1./(1.+0.5*z)
    value = t*exp(-z*z - 1.26551223 + t*(1.00002368 + t*(.37409196 + &
          & t*(.09678418 + t*(-.18628806 + t*(.27886807 + t*(-1.13520398 + &
          & t*(1.48851587 + t*(-.82215223 + t*.17087277)))))))))
    if (x < 0.) value = 2. - value
end subroutine erfc

subroutine erf(x, value)
    real(8), intent(in) :: x
    real(8), intent(out) :: value
    !Returns the error function erf(x) with fractional 
    !error everywhere less than 1:2^10-7
    real(8) :: t, z
    z = abs(x)
    t = 1./(1.+0.5*z)
    value = t*exp(-z*z - 1.26551223 + t*(1.00002368 + t*(.37409196 + &
          & t*(.09678418 + t*(-.18628806 + t*(.27886807 + t*(-1.13520398 + &
          & t*(1.48851587 + t*(-.82215223 + t*.17087277)))))))))
    if (x < 0.) value = 2. - value
    value = 1. - value
end subroutine erf


integer function GetPairIndFromij(i, j)
    integer, intent(in) :: i, j
    if (i > j) then
        GetPairIndFromij = i * (i + 1) / 2 + j
    else
        GetPairIndFromij = j * (j + 1) / 2 + i
    endif      
end function

subroutine GetijFromPairInd(ind, i, j)
    integer, intent(in) :: ind
    integer, intent(out) :: i, j
    real(8) :: r
    r = dble(ind) + 0.001d0
    r = sqrt(1.d0 + 8.d0 * r) * 0.5d0 - 0.5d0
    i = int(r)
    j = ind - i*(i+1)/2
end subroutine

subroutine RandomRotMat3D(MaxAng, RotMat)
    real(8), intent(in) :: MaxAng
    real(8), dimension(3,3) :: RotMat
    real(8) :: theta, sphi, cphi
    real(8) :: Ang, q0, q1, q2, q3
    real(8), dimension(3) :: Vec
    real(8), parameter :: pi = 3.1415926535897931D0
    !get a random rotation angle
    call ran2(Ang)
    Ang = (2.*Ang - 1) * MaxAng
    !get a random vector about which to rotate
    call ran2(cphi)
    cphi = 2.*cphi-1.
    cphi = max(-1.,min(1.,cphi))
    sphi = sqrt(1.-cphi*cphi)
    call ran2(theta)
    theta = 2.d0 * pi * theta
    Vec = (/ cos(theta)*sphi, sin(theta)*sphi, cphi /) 
    !make intermediate variables
    q0 = cos(0.5*Ang)
    Vec = sin(0.5*Ang) * Vec
    q1 = Vec(1)
    q2 = Vec(2)
    q3 = Vec(3)
    !assemble the rotation matrix
    RotMat(1,1) = q0*q0 + q1*q1 - q2*q2 - q3*q3
    RotMat(1,2) = 2.*(q1*q2 + q0*q3)
    RotMat(1,3) = 2.*(q1*q3 - q0*q2)
    RotMat(2,1) = 2.*(q1*q2 - q0*q3)
    RotMat(2,2) = q0*q0 - q1*q1 + q2*q2 - q3*q3
    RotMat(2,3) = 2.*(q2*q3 + q0*q1)
    RotMat(3,1) = 2.*(q1*q3 + q0*q2)
    RotMat(3,2) = 2.*(q2*q3 - q0*q1)
    RotMat(3,3) = q0*q0 - q1*q1 - q2*q2 + q3*q3
end subroutine

subroutine RandomRotMatxyz(MaxAng, RotMat)
    real(8), intent(in) :: MaxAng
    real(8), dimension(3,3) :: RotMat
    real(8) :: theta, costheta, sintheta, axis
    !get a random rotation angle
    call ran2(theta)
    theta = (2.d0 * theta - 1.d0) * MaxAng
    sintheta = sin(theta)
    costheta = cos(theta)
    !get a random axis about which to rotate
    call ran2(axis)
    RotMat = 0.d0
    if (axis < 0.33333333333d0) then
        RotMat(1,1) = costheta
        RotMat(1,2) = -sintheta
        RotMat(2,1) = sintheta
        RotMat(2,2) = costheta
        RotMat(3,3) = 1.d0
    elseif (axis < 0.66666666667d0) then
        RotMat(2,2) = costheta
        RotMat(2,3) = -sintheta
        RotMat(3,2) = sintheta
        RotMat(3,3) = costheta
        RotMat(1,1) = 1.d0    
    else
        RotMat(1,1) = costheta
        RotMat(1,3) = -sintheta
        RotMat(3,1) = sintheta
        RotMat(3,3) = costheta
        RotMat(2,2) = 1.d0
    endif
end subroutine



function modulehash()
    character(len=40) :: modulehash
    modulehash = 'e0d99b4785cb33d14dbf8e79f8c7a8ce0debf1c2'
end function


end module

