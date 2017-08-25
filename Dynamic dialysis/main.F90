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
! For a typical application, where the alcohol aqueous solution is fed 
! in the shell side and the pure water is fed in the lumen side, 
! the shell-side solute will permeate into the lumen side through the membrane
! due to the concentration gradient. At the same time, the lumen-side solvent
! will also reversely permeate into the shell side. 
! The mass concentration of each component, and inlet velocities of both sides
! should be known in advance.
!
! Program structure:
! main.f90 
! |-- prepare.f90    (modules of collecting the variables and parameters)
! |-- solver.f90     (modules of modelling equations and INTEL ode solver)
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