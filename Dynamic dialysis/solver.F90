!******************************************************************************* 
!                              INTEL CONFIDENTIAL 
!   Copyright(C) 2007-2008 Intel Corporation. All Rights Reserved. 
!   The source code contained  or  described herein and all documents related to 
!   the source code ("Material") are owned by Intel Corporation or its suppliers 
!   or licensors.  Title to the  Material remains with  Intel Corporation or its 
!   suppliers and licensors. The Material contains trade secrets and proprietary 
!   and  confidential  information of  Intel or its suppliers and licensors. The 
!   Material  is  protected  by  worldwide  copyright  and trade secret laws and 
!   treaty  provisions. No part of the Material may be used, copied, reproduced, 
!   modified, published, uploaded, posted, transmitted, distributed or disclosed 
!   in any way without Intel's prior express written permission. 
!   No license  under any  patent, copyright, trade secret or other intellectual 
!   property right is granted to or conferred upon you by disclosure or delivery 
!   of the Materials,  either expressly, by implication, inducement, estoppel or 
!   otherwise.  Any  license  under  such  intellectual property  rights must be 
!   express and approved by Intel in writing.
! 
!*******************************************************************************
!
!   To solve the concentration profile along the tubular dynamic dialysis module
!   
!   by Guoqiang GUAN Dr., 2017/08/15
!
!*******************************************************************************

      SUBROUTINE dd1s_sol()
      
      INTEGER n, ierr, i
! It is higly recommended to declare ipar array of size 128 
! for compatibility with future versions of ODE solvers
      INTEGER kd(2), ipar(128)
      DOUBLE PRECISION z, z_end, h, hm, ep, zr
! As ODE system has size n=6, than the size of dpar array is equal to 
! max{13*n,(7+2*n)*n}=max{78,114}=114. More details on dpar array can be 
! found in the Manual
      DOUBLE PRECISION y(6), dpar(114)
      EXTERNAL rhs_dd1s, jacmat_dd1s
      REAL time_begin, time_end

! global parameter settings suitable for all 6 dodesol routines  
! minimal step size for the methods
         hm = 1.d-12
! relative tolerance. The code cannot guarantee the requested accuracy for ep<1.d-9
         ep = 1.d-6
! absolute tolerance
         zr = 1.d-3

! ****************************** dodesol ********************************
! initialize ipar array with zeros before the first call to dodesol routines
         DO i=1,128
            ipar(i) = 0
         END DO  
         
         z = 0.d0
         h = 1.d-7       
         
! setting size of the system n, end of integration interval z_end, and initial 
! value y at z=0
         CALL dd1s_ic(n,z_end,y)
! display ic         
         write(*,'(A3, F8.5)') 't=', z
         PRINT*, 'Solution of y(i)'
         write(*,'(A5, I1, A4, F8.5)') ('y(', i, ') = ', y(i), i = 1, n)

         CALL CPU_TIME(time_begin)
! universal solver
         CALL dodesol(ipar,n,z,z_end,y,rhs_dd1s,jacmat_dd1s,h,hm,ep,zr,dpar,kd,ierr)

         CALL CPU_TIME(time_end)

         IF(ierr.ne.0) THEN
            PRINT*,'========================'
            PRINT*,'DODESOL FORTRAN example FAILED'
            PRINT*,'dodesol routine exited with error code',ierr  
            STOP 1
         END IF
   
         PRINT*
         PRINT*, 'dodesol results'
         PRINT*
         write(*,'(A11, I1, A12, I1)') 'ipar(2) = ', ipar(2), ', ipar(4) = ', ipar(4)
         write(*,'(A3, F8.5)') 't=', z
         PRINT*, 'Solution of y(i)'
         write(*,'(A5, I1, A4, F8.5)') ('y(', i, ') = ', y(i), i = 1, n)
         PRINT*, '-----------------------------------------------------'
         PRINT*, 'CPU time=',time_end-time_begin,' seconds'
         PRINT*, '====================================================='
         PRINT*
         pause

      RETURN
      END SUBROUTINE


! ********************** Initial condition for odes **************  
      SUBROUTINE dd1s_ic(n,z_end,y)
! The routine initializes the size of the system n, the end of 
! integration interval t_end, and inital data y at t=0.0 
      USE ComParams

      IMPLICIT NONE

      INTEGER n
      DOUBLE PRECISION z_end,y(*)
          
      n = 6 
      z_end = DD_MOD%FlowChannel(Tubeside)%Length

      y(1) = Influent(Tubeside)%Component(ETOH)%FracConc%Mass
      y(2) = Influent(Shellside)%Component(ETOH)%FracConc%Mass
      y(3) = Influent(Tubeside)%Component(H2O)%FracConc%Mass
      y(4) = Influent(Shellside)%Component(H2O)%FracConc%Mass
      y(5) = Influent(Tubeside)%Velocity
      y(6) = Influent(Shellside)%Velocity

      RETURN
      END

! ******************* Right hand side ******************  
      SUBROUTINE rhs_dd1s(n,z,y,f)

      USE ComParams
      
      INTEGER n
      DOUBLE PRECISION z, y(*), f(*)
      
      integer :: icomp, iside
      real*8 :: RH(2), Density(2), JM(NumComp)

!     Initiation
<<<<<<< HEAD
      call Initiate()

!     Influent properties

      AM = PI*ID(1)*Length
      alpha(1) = four*Length/ID(1)
      alpha(2) = four*ID(1)*Length/(ID(2)*ID(2)-ID(1)*ID(1))
      NMA = PA/Thickness*(y(1)*InitMassConcA(1)-y(3)*InitMassConcA(2))
      NMB = PB/Thickness*(y(2)*InitMassConcB(2)-y(4)*InitMassConcB(1))
      InitMassFluxA(1) = InitMassConcA(1)*InitVelocity(1)
      InitMassFluxB(1) = InitMassConcB(1)*InitVelocity(1)
      InitMassFluxA(2) = InitMassConcA(2)*InitVelocity(2)
      InitMassFluxB(2) = InitMassConcB(2)*InitVelocity(2)
      InitMassFlux(1) = InitMassFluxA(1)+InitMassFluxB(1)
      InitMassFlux(2) = InitMassFluxA(2)+InitMassFluxB(2)
      DimlessPresGrad(1) = (Pressure_out(1)-Pressure_in(1))/InitMassFlux(1)/InitVelocity(1)
      DimlessPresGrad(2) = (Pressure_out(2)-Pressure_in(2))/InitMassFlux(2)/InitVelocity(2)

      f(1) = (-alpha(1)*NMA/InitMassFluxA(1)-y(1)*AM*(NMA*SpecVolA-NMB*SpecVolB)/InitVelocity(1))/y(5)
      f(2) = (+alpha(1)*NMB/InitMassFluxB(1)-y(2)*AM*(NMA*SpecVolA-NMB*SpecVolB)/InitVelocity(1))/y(5)
      f(3) = (+alpha(2)*NMA/InitMassFluxA(2)-y(3)*AM*(NMB*SpecVolB-NMA*SpecVolA)/InitVelocity(2))/y(6)
      f(4) = (-alpha(2)*NMB/InitMassFluxB(2)-y(4)*AM*(NMB*SpecVolB-NMA*SpecVolA)/InitVelocity(2))/y(6)
      f(5) = y(7)
      f(6) = y(8)
!     DimlessDensityGrad(1) = InitMassFracA(1)*f(1)+InitMassFracB(1)*f(2)
!     DimlessDensityGrad(2) = InitMassFracA(2)*f(3)+InitMassFracB(2)*f(4)
!     DimlessDensity(1) = InitMassFracA(1)*y(1)+InitMassFracB(1)*y(2)
!     DimlessDensity(2) = InitMassFracA(2)*y(3)+InitMassFracB(2)*y(4)
      f(7) = InitRe(1)/two*(DimlessPresGrad(1)+y(5)*y(5)*(InitMassFracA(1)*f(1)+InitMassFracB(1)*f(2))+two*(InitMassFracA(1)*y(1)+InitMassFracB(1)*y(2))*y(5)*y(7))
      f(8) = InitRe(2)/two*(DimlessPresGrad(2)+y(6)*y(6)*(InitMassFracA(2)*f(3)+InitMassFracB(2)*f(4))+two*(InitMassFracA(2)*y(3)+InitMassFracB(2)*y(4))*y(6)*y(8))

=======
	call Initiate()
	
	do iside = Tubeside, Shellside
	      RH(iside) = DD_MOD%FlowChannel(iside)%DH/four
	end do
      Density(Tubeside)  = y(1)*Chemical(ETOH)%Density+y(3)*Chemical(H2O)%Density
      Density(Shellside) = y(2)*Chemical(ETOH)%Density+y(4)*Chemical(H2O)%Density
      JM(ETOH)           = DD_MOD%Membrane%Permeability(ETOH)/DD_MOD%Membrane%Thickness
      JM(H2O)            = DD_MOD%Membrane%Permeability(H2O)/DD_MOD%Membrane%Thickness
      
!     right-hand side of odes
      f(1) = RH(Tubeside)*(y(1)*JM(H2O)+y(3)*JM(ETOH)) / (Density(Tubeside)*y(5))
      f(2) = -RH(Shellside)*(y(4)*JM(ETOH)+y(2)*JM(H2O)) / (Density(Shellside)*y(6))
      f(3) = -RH(Tubeside)*(y(1)*JM(H2O)+y(3)*JM(ETOH)) / (Density(Tubeside)*y(5))
      f(4) = RH(Shellside)*(y(4)*JM(ETOH)+y(2)*JM(H2O)) / (Density(Shellside)*y(6))
      f(5) = RH(Tubeside)*(JM(H2O)-JM(ETOH))/Density(Tubeside) - y(5)/Density(Tubeside)*(Chemical(ETOH)%Density*f(1)+Chemical(H2O)%Density*f(3))
      f(6) = RH(Shellside)*(JM(H2O)-JM(ETOH))/Density(Shellside) - y(6)/Density(Shellside)*(Chemical(ETOH)%Density*f(2)+Chemical(H2O)%Density*f(4))
       
>>>>>>> develop
      RETURN
      END

! ************* analytical Jacoby matrix **************
      SUBROUTINE jacmat_dd1s(n,t,y,a)
<<<<<<< HEAD
=======
      
      USE ComParams
      
>>>>>>> develop
      IMPLICIT NONE 
      INTEGER n
<<<<<<< HEAD
      DOUBLE PRECISION t,y(*),a(n,*)
      a(1,1) = 0.d0
      a(1,2) = 1.d0
      a(2,1) = -1.d6*(1.d0+2.d0*y(1)*y(2))
      a(2,2) = 1.d6*(1.d0-y(1)* y(1))
=======
      DOUBLE PRECISION t, y(*), a(n,*)
      
>>>>>>> develop
      RETURN
      END

! ************************* End of subroutine *************************
