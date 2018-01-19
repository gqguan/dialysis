!     =============================================================70
!     Declare the common variables
!     by Guan, Guoqiang at 19-10-2015

      MODULE ComParams
      
      IMPLICIT NONE
      
      REAL*8 pi, zero, half, one, two, three, four
      PARAMETER (pi = 3.1415926D0, zero = 0.0D0, half = 0.5d0)
      PARAMETER (one = 1.0d0, two = 2.0d0, three = 3.0d0, four = 4.0d0)
      integer NumComp, Tubeside, Shellside
      parameter (NumComp = 2, Tubeside = 1, Shellside = 2)
            
!     Common variables declaration
!     Common area - GeoParams, to store the global variables of
!                              tube-side and shell-side geometric 
!                              parameters
!                   ID: inner diameters
!                   OD: outer diameters
!                   DE: equivalent diameter
!                   DH: hydraulic diameter
!                   CSA: cross section area
      REAL*8 ID(2), OD(2), DE(2), DH(2), CSA(2)
      COMMON /GeoParas/ ID, OD, DE, DH, CSA
      
!     Common area - MemProps, to store the global variables of
!                             membrane properties
!                   Length: membrane length
!                   Thickness: membrane thickness
!                   PA: membrane permeability of component A
!                   PB: membrane permeability of component B
      REAL*8 Length, Thickness, PA, PB
      COMMON /MemProps/ Length, Thickness, PA, PB
      
!     Common area - InitStream, to store the global variables of
!                               tube-side and shell-side streams
!                               at z=0
!                   InitMassConcA: inlet mass concentration of component A
!                   InitMassConcB: inlet mass concentration of component B
!                   InitVelocity: inlet velocity
!                   InitRe: inlet Reynolds number
!                   InitDensity: inlet density
!                   InitMassFracA: inlet mass fraction of component A
!                   InitMassFracB: inlet mass fraction of component B
!                   InitMassFluxA: inlet mass flux of component A
!                   InitMassFluxB: inlet mass flux of component B
!                   InitMassFlux: inlet mass flux of solution
!                   InitViscosity: influent viscosity
!                   DimlessPresGrad: gradient of dimensionless pressure
!                   Pressure_in: inlet pressure
!                   Pressure_out: outlet pressure
!                   SpecVolA: specific mass-based volume of material A
!                   SpecVolB: specific mass-based volume of material B
      REAL*8 InitMassConcA(2), InitMassConcB(2), InitVelocity(2)
      REAL*8 InitRe(2), InitDensity(2), InitMassFracA(2), InitMassFracB(2)
      REAL*8 InitMassFluxA(2), InitMassFluxB(2), InitMassFlux(2), InitViscosity(2)
      REAL*8 DimlessPresGrad(2), Pressure_in(2), Pressure_out(2)
      REAL*8 SpecVolA, SpecVolB
      COMMON /InitStream/ InitMassConcA, InitMassConcB, InitVelocity, InitRe, InitDensity, InitMassFracA, InitMassFracB, InitMassFluxA, InitMassFluxB, InitMassFlux, InitViscosity
      COMMON /FlowDrivingForces/ Pressure_in, Pressure_out, DimlessPresGrad
      COMMON /StreamProps/ SpecVolA, SpecVolB
      
!     Common area - ComIVars, to store the global integer variables
!                   COM_OPT - running option
      INTEGER COM_OPT
      COMMON /ComIVars/ COM_OPT
      
      CONTAINS
      
      subroutine Initiate()
      real*8 :: InitMass(NumComp), MoleWeight(NumComp), SpecVol(NumComp)
      real*8 :: rho(NumComp), omega(NumComp), c(NumComp), x(NumComp)
      integer :: i
      
!     Initiate tube-side influent
      data InitMass /1.d0,9.d0/
      InitVelocity(Tubeside)  = 1.0d-1
      InitViscosity(Tubeside) = 1.0d-3
      Pressure_in(Tubeside)   = 1.5d5
      Pressure_out(Tubeside)  = 1.0d5
      call CalcConc(NumComp, InitMass, MoleWeight, SpecVol, rho, c, omega, x)
      InitMassConcA(Tubeside) = rho(1)
      InitMassConcB(Tubeside) = rho(2)
      InitMassFracA(Tubeside) = omega(1)
      InitMassFracB(Tubeside) = omega(2)
      InitDensity(Tubeside) = InitMassConcA(Tubeside)+InitMassConcB(Tubeside)
      InitRe(Tubeside) = DRE(InitDensity(Tubeside), InitVelocity(Tubeside), ID(Tubeside), InitViscosity(Tubeside))

!     Initiate shell-side influent
      InitMass = (/9.d0,1.d0/)
      InitVelocity(shellside)  = 1.0d-1
      InitViscosity(shellside) = 1.0d-3
      Pressure_in(shellside)   = 1.5d5
      Pressure_out(shellside)  = 1.0d5
      call CalcConc(NumComp, InitMass, MoleWeight, SpecVol, rho, c, omega, x)
      InitMassConcA(Shellside) = rho(1)
      InitMassConcB(Shellside) = rho(2)
      InitMassFracA(Shellside) = omega(1)
      InitMassFracB(Shellside) = omega(2)
      InitDensity(Shellside) = InitMassConcA(Shellside)+InitMassConcB(Shellside)
      InitRe(Shellside) = DRE(InitDensity(Shellside), InitVelocity(Shellside), ID(Shellside), InitViscosity(Shellside))

!     Geometric properties of membrane and module 
      ID(Tubeside)  = 1.0d-2
      ID(Shellside) = 1.5d-2
      Length = 2.0d-1
      Thickness = 4.2d-5
      PA = 1.0d-3
      PB = 1.0d-3
!     Fluid properties: A means ethanol, B means water
      SpecVolA = 1.267d-3
      SpecVolB = 1.000d-3
      data MoleWeight /18.d0, 46.d0/


      contains

            SUBROUTINE CalcConc(n, CompMass, MoleWeight, SpecVol, MassConc, MoleConc, MassFrac, MoleFrac)
            INTEGER, INTENT(IN) :: n
            REAL*8, INTENT(IN) :: CompMass(n), MoleWeight(n), SpecVol(n)
            REAL*8, INTENT(OUT) :: MassConc(n), MassFrac(n), MoleConc(n), MoleFrac(n)
            INTEGER :: i
            real*8 :: volume, CompMole(n)
            volume = zero
            do i = 1, n
                  volume = volume+CompMass(i)*SpecVol(i)
                  CompMole(i) = CompMass(i)/MoleWeight(i)
            end do
            do i = 1, n
                  MassConc(i) = CompMass(i)/volume
                  MoleConc(i) = CompMole(i)/volume
            end do
            do i = 1, n
                  MassFrac(i) = CompMass(i)/sum(CompMass)
                  MoleFrac(i) = CompMole(i)/sum(CompMole)
            end do
            END SUBROUTINE CalcConc

      end subroutine Initiate

      REAL*8 FUNCTION DRE(density, velocity, distance, viscosity)
      REAL*8, INTENT(IN) :: density, velocity, distance, viscosity
      DRE = density*velocity*distance/viscosity
      END FUNCTION DRE

      END MODULE ComParams