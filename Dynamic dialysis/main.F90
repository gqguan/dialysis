program dd1s
! Nomenclature:
! DD means the dynamic dialysis
! 1  means the 1-D model
! S  means the steady state
!
! Problem description:
! The steady-state simulation of mass transfer in a tubular membrane
! module will be performed for an extensive application of dynamic dialysis, 
! where a tubular dialysis membrane divides the champers into tube side and 
! shell side respectively. The models could be found in the document of 
! "20170814 Numerical solution of governing equations.docx".
!
! Application example:
! For a typical application, where the ethanol and water binary solution is fed 
! in the tube side and the water is fed in the shell side, 
! the tube-side components will permeate into the shell side through the membrane
! due to the concentration gradient. The shell-side materials will also reversely
! permeate into the tube side. 
! The mass concentration of each component in the influents of both sides, 
! together with the velocities and Reynolds number, should be known in advance.
! A constant pressure gradient along the flowing path was also assumed to
! simplify the equations of the momentum balance.

! Program file structure:
! main.f90 
! |-- prepare.f90    (modules of collecting the variables and parameters)
!     MODULE ComParams
!     |-- SUBROUTINE Initiate()
! |-- solver.f90     (modules of modelling equations and ode solver)
!     SUBROUTINE dd1s_sol()
! |-- supplement.f90 (modules of subroutines and customized functions)
!
! License:
! This code is distributed under the GNU LGPL license.
!
! Authors:
! Guoqiang GUAN Dr. (gqguan@scut.edu.cn)
!
! Revision log:
! Refer to https://github.com/gqguan/dialysis

use ComParams

  call Initiate()
  call dd1s_sol()

end program