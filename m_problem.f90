module m_problem
   use iso_fortran_env
   implicit none
   ! defines problem-specific but global constants
   integer, parameter :: n_spc = 8, n_reaction = 5
   !     Species Numbering  H2O(1) CO2(2) CO(3) OH(4) O2(5) H2(6) O(7) H(8)
   integer, parameter :: H2O = 1, CO2 = 2, CO = 3, OH = 4, O2 = 5, H2 = 6, O = 7, H = 8
end module m_problem
