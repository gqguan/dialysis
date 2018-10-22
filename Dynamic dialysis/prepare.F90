!     =============================================================70
!     Declare the common variables
!     Original by Guan, Guoqiang (GGQ) at 19-10-2015
!     Revision 2017/08/15 GGQ 

      MODULE ComParams
      
      IMPLICIT NONE
      
      REAL*8 pi, zero, half, one, two, three, four
      PARAMETER (pi = 3.1415926D0, zero = 0.0D0, half = 0.5d0)
      PARAMETER (one = 1.0d0, two = 2.0d0, three = 3.0d0, four = 4.0d0)
      integer NumComp, ETOH, H2O, Tubeside, Shellside
      parameter (NumComp = 2, ETOH = 1, H2O = 2, Tubeside = 1, Shellside = 2)
      
      REAL*8 :: OP_Temperature = 298.15d0 ! operating temperature [K]
      REAL*8 :: OP_Pressure = 1.01325d5 ! operating pressure [Pa]
            
!     Common variables declaration
      TYPE GeomParams
!           ID: inner diameters [m]
!           OD: outer diameters [m]
!           Length: effective membrane length [m]
!           DE: equivalent diameter [m]
!           DH: hydraulic diameter ( = 4CSA/PERIM) [m]
!           CSA: cross-sectional area of the conduit [m2]
!           CIRCUM: wetted perimeter [m]
            REAL*8 :: ID, OD, Length
            REAL*8 :: DE, DH, CSA, CIRCUM ! Derived properties by supplied calculators
      END TYPE GeomParams
      TYPE MemProps
            REAL*8 :: Thickness ! Membrane thickness [m]
            REAL*8 :: Permeability(NumComp) ! Membrane permeability of each component [kg/m/s]
      END TYPE MemProps
      TYPE MemMods
            TYPE (GeomParams) :: FlowChannel(2)
            TYPE (MemProps) :: Membrane
      END TYPE MemMods
      TYPE Bases
            REAL*8 :: Mole, Mass, Volume
      END TYPE Bases
      TYPE MatProps
            CHARACTER(12) :: Formula ! Molecular formula [12-letter-length string]
            REAL*8 :: MW ! Molecular weight [kg/mol]
            REAL*8 :: Density ! Density [kg/m3]
            TYPE (Bases) :: FracConc ! Fraction-based concentrations, related to the multi-components mixture
            TYPE (Bases) :: CorrCoeff ! Correction coefficients, related to the multi-components mixture
      END TYPE MatProps
      TYPE StrmProps
            TYPE (MatProps) :: Component(NumComp)
            REAL*8 :: Density ! Density [kg/m3]
            TYPE (Bases) :: Flowrate ! Flowrate
            REAL*8 :: Velocity ! Could be derived with known CSA and Volumetric Flow [m/s]
      END TYPE StrmProps    

      TYPE (MemMods) :: DD_MOD
      TYPE (StrmProps) :: Influent(2), Effluent(2)
      
      TYPE (MatProps) :: Chemical(NumComp)
      DATA Chemical(ETOH)%Formula, Chemical(ETOH)%MW, Chemical(ETOH)%Density /"ETOH", 0.04606844d0, 879.d0/
      DATA Chemical(H2O)%Formula,  Chemical(H2O)%MW,  Chemical(H2O)%Density  /"H2O" , 0.01801527d0, 998.d0/
     
!     Common area - ComIVars, to store the global integer variables
!                   COM_OPT - running option
      INTEGER COM_OPT
      COMMON /ComIVars/ COM_OPT
      
      CONTAINS
      
      subroutine Initiate()
            type (MemProps) :: Prod_68035
            real*8 :: VolFrac(2) ! Alcohol concentration in volumetric fraction (v/v)
            integer :: iside, icomp
!           Initiate the membrane properties for dynamic dialysis
            Prod_68035%Thickness          = 25.d-6
            Prod_68035%Permeability(ETOH) = 1.3d-4 ! <<<=== check [kg/m-s]
            Prod_68035%Permeability(H2O)  = 2.0d-4 ! <<<=== check [kg/m-s]
!           Initiate the geometrics of the membrane module
            DD_MOD%Membrane                      = Prod_68035
            DD_MOD%FlowChannel(Tubeside)%ID      = 2.20d-2
            DD_MOD%FlowChannel(Tubeside)%OD      = DD_MOD%FlowChannel(Tubeside)%ID+two*DD_MOD%Membrane%Thickness
            DD_MOD%FlowChannel(Tubeside)%CSA     = pi/four*DD_MOD%FlowChannel(Tubeside)%ID**two
            DD_MOD%FlowChannel(Tubeside)%CIRCUM  = pi*DD_MOD%FlowChannel(Tubeside)%ID
            DD_MOD%FlowChannel(Tubeside)%DE      = dsqrt(four/pi*DD_MOD%FlowChannel(Tubeside)%CSA)
            DD_MOD%FlowChannel(Tubeside)%DH      = four*DD_MOD%FlowChannel(Tubeside)%CSA/DD_MOD%FlowChannel(Tubeside)%CIRCUM
            DD_MOD%FlowChannel(Tubeside)%Length  = 2.00d-1
            DD_MOD%FlowChannel(Shellside)%ID     = 2.68d-2
            DD_MOD%FlowChannel(Shellside)%OD     = 3.00d-2
            DD_MOD%FlowChannel(Shellside)%CSA    = pi/four*(DD_MOD%FlowChannel(Shellside)%ID**two-DD_MOD%FlowChannel(Tubeside)%OD**two)
            DD_MOD%FlowChannel(Shellside)%CIRCUM = pi*DD_MOD%FlowChannel(Tubeside)%OD
            DD_MOD%FlowChannel(Shellside)%DE     = dsqrt(four/pi*DD_MOD%FlowChannel(Shellside)%CSA)
            DD_MOD%FlowChannel(Shellside)%DH     = four*DD_MOD%FlowChannel(Shellside)%CSA/DD_MOD%FlowChannel(Shellside)%CIRCUM
            DD_MOD%FlowChannel(shellside)%Length = DD_MOD%FlowChannel(Tubeside)%Length
!           Initiate the lumen-side influent
            VolFrac(Tubeside)  = 0.00d0 ! lumen-side volumetric fraction of methanol
            Influent(Tubeside)%Component(ETOH)                 = Chemical(ETOH)
            Influent(Tubeside)%Component(ETOH)%FracConc%Volume = VolFrac(Tubeside)
            Influent(Tubeside)%Component(H2O)                  = Chemical(H2O)
            Influent(Tubeside)%Component(H2O)%FracConc%Volume  = one-VolFrac(Tubeside)
!           Get the initial mixing properties
            do icomp = 1, NumComp
                  Influent(Tubeside)%Component(icomp)%CorrCoeff%Volume = one
                  Influent(Tubeside)%Component(icomp)%CorrCoeff%Mass   = one
                  Influent(Tubeside)%Component(icomp)%CorrCoeff%Mole   = one
            end do
            Influent(Tubeside)%Component%FracConc%Mass         = ConvertFrac(NumComp, Influent(Tubeside), "Vol2Mass")
            Influent(Tubeside)%Component%FracConc%Mole         = ConvertFrac(NumComp, Influent(Tubeside), "Vol2Mole")
            Influent(Tubeside)%Density = zero
            do icomp = 1, NumComp
                  Influent(Tubeside)%Density = Influent(Tubeside)%Density + Influent(Tubeside)%Component(icomp)%Density * Influent(Tubeside)%Component(icomp)%FracConc%Mass
            end do
            Influent(Tubeside)%Flowrate%Volume                 = 5.d-6/60.d0 ! 5 mLPM
            Influent(Tubeside)%Velocity = Influent(Tubeside)%Flowrate%Volume/DD_MOD%FlowChannel(Tubeside)%CSA
!      	Initiate shell-side influent
            VolFrac(Shellside)  = 0.25d0 ! shell-side volumetric fraction of methanol
            Influent(Shellside)%Component(ETOH)                 = Chemical(ETOH)
            Influent(Shellside)%Component(ETOH)%FracConc%Volume = VolFrac(Shellside)
            Influent(Shellside)%Component(H2O)                  = Chemical(H2O)
            Influent(Shellside)%Component(H2O)%FracConc%Volume  = one-VolFrac(Shellside)
!           Get the initial mixing properties
            do icomp = 1, NumComp
                  Influent(Shellside)%Component(icomp)%CorrCoeff%Volume = one
                  Influent(Shellside)%Component(icomp)%CorrCoeff%Mass   = one
                  Influent(Shellside)%Component(icomp)%CorrCoeff%Mole   = one
            end do
            Influent(Shellside)%Component%FracConc%Mass         = ConvertFrac(NumComp, Influent(Shellside), "Vol2Mass")
            Influent(Shellside)%Component%FracConc%Mole         = ConvertFrac(NumComp, Influent(Shellside), "Vol2Mole")
            Influent(Shellside)%Density = zero
            do icomp = 1, NumComp
                  Influent(Shellside)%Density = Influent(Shellside)%Density + Influent(Shellside)%Component(icomp)%Density * Influent(Shellside)%Component(icomp)%FracConc%Mass
            end do
            Influent(Shellside)%Flowrate%Volume                 = 5.d-6/60.d0 ! 5 mLPM
            Influent(Shellside)%Velocity = Influent(Shellside)%Flowrate%Volume/DD_MOD%FlowChannel(Shellside)%CSA
      end subroutine Initiate
      
      function ConvertFrac(ncomp, mixture, str_opt)
!     Convert the volumetric fraction into mass or mole fraction
            integer, intent(in) :: ncomp ! number of components
            type (StrmProps), intent(in) :: mixture
            character*(*), intent(in) :: str_opt ! string option
            real*8, dimension(ncomp) :: CompQuantity, ConvertFrac
            integer :: icomp
!           Get the mixing propeties according to the given composition <<<=== Not done yet
            select case(str_opt)
                  case("Vol2Mass")
                        do icomp = 1, ncomp
                              CompQuantity(icomp) = mixture%Component(icomp)%FracConc%Volume / mixture%Component(icomp)%Density * mixture%Component(icomp)%CorrCoeff%Volume
                        end do
                        do icomp = 1, ncomp
                              ConvertFrac(icomp) = CompQuantity(icomp) / SUM(CompQuantity)
                        end do
                  case("Vol2Mole")
                        do icomp = 1, ncomp
                              CompQuantity(icomp) = mixture%Component(icomp)%FracConc%Volume / mixture%Component(icomp)%Density * mixture%Component(icomp)%CorrCoeff%Volume / mixture%Component(icomp)%MW
                        end do
                        do icomp = 1, ncomp
                              ConvertFrac(icomp) = CompQuantity(icomp) / SUM(CompQuantity)
                        end do
                  case default
                        do icomp = 1, ncomp
                              ConvertFrac(icomp) = zero
                        end do
            end select
            return
      end function ConvertFrac
      
      real*8 function AlcoholDensity(T, P, frac, str_opt)
!     Calculate the density of alcohol aqueous solution for given EtOH fraction at specified temperature and pressure
!     The given EtOH fraction can be volumetric, mass or molar, which depends on the string option of "VolFrac", "MassFrac", or "MolFrac"
            real*8, intent(in) :: T, P, frac
            character*(*), intent(in) :: str_opt
            real*8, dimension(NumComp) :: CompMass, MassFrac
            select case(str_opt)
            case("VolFrac") ! simplified calculation w/o considering the partial molar properties
                  CompMass(ETOH) = frac/Chemical(ETOH)%Density
                  CompMass(H2O)  = (one-frac)/Chemical(H2O)%Density
                  MassFrac(ETOH) = CompMass(ETOH)/(CompMass(ETOH)+CompMass(H2O))
                  MassFrac(H2O)  = CompMass(H2O) /(CompMass(ETOH)+CompMass(H2O))
                  AlcoholDensity = Chemical(ETOH)%Density*MassFrac(ETOH)+Chemical(H2O)%Density*MassFrac(H2O)
            case default
                  AlcoholDensity = zero
            end select
            return
      end function AlcoholDensity
      
      END MODULE ComParams