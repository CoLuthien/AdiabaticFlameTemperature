module m_specie
   use iso_fortran_env

   implicit none
   type Specie
      integer :: idx, comp(4)
      character(len=20) :: name ! for debug
      real(real64) :: mole_weight, &
                      molar_hof, mass_hof, & ! heat of formation
                      mole_diam, & !molecular diameter
                      chrt_temp, & ! charateristic temperature
                      poly_cp(7, 2), & ! polynomial for cp fit
                      poly_h(7, 2), & !polynomial fitting for enthalpy
                      poly_s(7, 2) ! polynomial coefficients for entropy
   contains
      procedure :: H => enthalpy
      procedure :: S => entropy
      procedure :: G => gibbs
      procedure :: mG => m_gibbs
   end type Specie

   character(len=80), parameter, private :: & ! Ref. USC Mech-II
      f1 = "(*)", &
      f2 = "(3F10.3)", &
      f3 = "(A18,A6,4(A2,I3),A1,D10.0,D10.0,D8.0,5x,I1)", &
      f4 = "(5(F15.8),I5)", &
      f5 = "(5(F15.8),I5)", &
      f6 = "(5(F15.8),I5)", &
      f7 = "(F10.5)" ! Molar weight

contains

   subroutine read_thermo(self, spc_name)
      type(Specie), intent(inout) :: self
      character(len=*), intent(in) :: spc_name
      real(real64) :: temp_rng(3), temp(3), coef(7, 2), garb
      integer :: fd, j, i, formular(4)
      character(len=80) :: THERMO
      character(len=18) :: chem_name
      character(len=6) :: date
      character(len=2) :: symbol(4)
      character :: phase
      fd = abs(fd) ! negative trash value makes problem
      open (fd, file=spc_name)

      !read (fd, f3) temp_rng(1), temp_rng(2), n_coeff, (t_exp(j), j=1, 8), zero_hof

      read (fd, *) THERMO
      read (fd, f2) temp_rng(:)
      read (fd, f3) chem_name(:), date, (symbol(j), formular(j), j=1, 4), phase, temp(:), i
      read (fd, f4) coef(:5, 1), i
      read (fd, f5) coef(6:7, 1), coef(:3, 2), i
      read (fd, f6) coef(4:7, 2), garb, i

      close (fd)
      self%comp(:) = formular(:)
      self%poly_cp(:, :) = coef(:, :)

      self%poly_h(:, :) = self%poly_cp(:, :)
      self%poly_s(:, :) = self%poly_cp(:, :)

      self%poly_h(2, :) = self%poly_cp(2, :)*0.50d0
      self%poly_h(3, :) = self%poly_cp(3, :)/3.0d0
      self%poly_h(4, :) = self%poly_cp(4, :)*0.25d0
      self%poly_h(5, :) = self%poly_cp(5, :)*0.20d0

      self%poly_s(3, :) = self%poly_cp(3, :)*.50d0
      self%poly_s(4, :) = self%poly_cp(4, :)/3.0d0
      self%poly_s(5, :) = self%poly_cp(5, :)*.25d0

   end subroutine read_thermo

   pure elemental function enthalpy(self, temp, idx) result(H)
      class(Specie), intent(in) :: self
      real(real64), intent(in) :: temp
      integer, intent(in) :: idx
      real(real64) :: H
      H = self%poly_h(1, idx) &
          + temp*(self%poly_h(2, idx) &
                  + temp*(self%poly_h(3, idx) &
                          + temp*(self%poly_h(4, idx) &
                                  + temp*(self%poly_h(5, idx))))) &
          + self%poly_h(6, idx)/temp
   end function enthalpy

   pure elemental function entropy(self, temp, idx) result(S)
      class(Specie), intent(in) :: self
      real(real64), intent(in) :: temp
      real(real64) :: S
      integer, intent(in) :: idx
      S = self%poly_s(1, idx)*log(temp) &
          + temp*(self%poly_s(2, idx) &
                  + temp*(self%poly_s(3, idx) &
                          + temp*(self%poly_s(4, idx) &
                                  + temp*(self%poly_s(5, idx))))) &
          + self%poly_s(7, idx)
   end function entropy

   pure elemental function gibbs(self, temp, idx) result(G)
      class(Specie), intent(in) :: self
      integer, intent(in) :: idx
      real(real64), intent(in) :: temp
      real(real64) :: G
      ! G = H - TS
      G = enthalpy(self, temp, idx) - entropy(self, temp, idx)
   end function gibbs

   pure elemental function m_gibbs(self, temp, idx) result(G)
      class(Specie), intent(in) :: self
      integer, intent(in) :: idx
      real(real64), intent(in) :: temp
      real(real64) :: G
      ! -G = TS -H
      G = entropy(self, temp, idx) - enthalpy(self, temp, idx)
   end function m_gibbs

   pure elemental function cp_heat(self, temp, idx) result(C)
      class(Specie), intent(in) :: self
      integer, intent(in) :: idx
      real(real64), intent(in) :: temp
      real(real64) :: C
      C = self%poly_cp(1, idx) &
          + temp*(self%poly_cp(2, idx) &
                  + temp*(self%poly_cp(3, idx) &
                          + temp*(self%poly_cp(4, idx) &
                                  + temp*(self%poly_cp(5, idx)))))

   end function cp_heat

end module m_specie
