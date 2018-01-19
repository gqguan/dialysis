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
!*******************************************************************************

      SUBROUTINE dd1s_sol()
      
      USE ComParams
      
      INTEGER n, ierr, i
! It is higly recommended to declare ipar array of size 128 
! for compatibility with future versions of ODE solvers
      INTEGER kd(2), ipar(128)
      DOUBLE PRECISION z, z_end, h, hm, ep, zr
! As ODE system has size n=8, than the size of dpar array is equal to 
! max{13*n,(7+2*n)*n}=max{104,184}=184. More details on dpar array can be 
! found in the Manual
      DOUBLE PRECISION y(8), dpar(184)
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
         PRINT*, 'ipar(2)=',ipar(2),', ipar(4)=',ipar(4)
         PRINT*, 't=',z
         PRINT*, 'Solution','  y1=',y(1),'  y2=',y(2)
         PRINT*, '-----------------------------------------------------'
         PRINT*, 'CPU time=',time_end-time_begin,' seconds'
         PRINT*, '====================================================='
         PRINT*

      RETURN
      END SUBROUTINE


! ********************** Initial condition for odes **************  
      SUBROUTINE dd1s_ic(n,z_end,y)
! The routine initializes the size of the system n, the end of 
! integration interval t_end, and inital data y at t=0.0 
      IMPLICIT NONE
      
      INTEGER n
      DOUBLE PRECISION z_end,y(*)
            
         n = 8 
         z_end = 1.d0

         y(1) = 1.d0
         y(2) = 1.d0
         y(3) = 1.d0
         y(4) = 1.d0
         y(5) = 1.d0
         y(6) = 1.d0
         y(7) = 0.d0
         y(8) = 0.d0

      RETURN
      END

! ******************* Right hand side ******************  
      SUBROUTINE rhs_dd1s(n,z,y,f)

      USE ComParams
      
      INTEGER n
      DOUBLE PRECISION z, y(*), f(*)
      DOUBLE PRECISION AM, alpha(2), NMA, NMB

!     Initiation
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

      RETURN
      END

! ************* analytical Jacoby matrix **************
      SUBROUTINE jacmat_dd1s(n,t,y,a)
      IMPLICIT NONE 
      INTEGER n
      DOUBLE PRECISION t,y(*),a(n,*)
      a(1,1) = 0.d0
      a(1,2) = 1.d0
      a(2,1) = -1.d6*(1.d0+2.d0*y(1)*y(2))
      a(2,2) = 1.d6*(1.d0-y(1)* y(1))
      RETURN
      END

! ************************* End of subroutine *************************
