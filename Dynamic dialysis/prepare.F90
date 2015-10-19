!     =============================================================70
!     Declare the common variables
!     by Guan, Guoqiang at 19-10-2015

      MODULE ComParams
      
      IMPLICIT NONE
      
      REAL*8 pi, zero, half, one, two, three, four
      PARAMETER (pi = 3.1415926D0, zero = 0.0D0, half = 0.5d0)
      PARAMETER (one = 1.0d0, two = 2.0d0, three = 3.0d0, four = 4.0d0)
            
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
      
!     Common area - StmProps, to store the global variables of
!                             tube-side and shell-side stream's
!                             properties
!                   InitMassConcA: initial mass concentration of component A
!                   InitMassConcB: initial mass concentration of component B
!                   InitVelocity: initial velocity
!                   InitViscosity: initial viscosity

      REAL*8 InitMassConcA(2), InitMassConcB(2), InitVelocity(2), InitViscosity(2)
      COMMON /StmProps/ InitMassConcA, InitMassConcB, InitViscosity
      
!     Common area - ComIVars, to store the global integer variables
!                   COM_OPT - running option
      INTEGER COM_OPT
      COMMON /ComIVars/ COM_OPT
      
      END MODULE ComParams