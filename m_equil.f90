module m_equil
   use iso_fortran_env
   use m_specie
   use m_problem
   implicit none

contains

   pure subroutine calc_eq_const(spcs, temp, kp)
      class(Specie), intent(in) :: spcs(n_spc)
      real(real64), intent(in) :: temp
      real(real64), intent(out) :: kp(n_reaction)
      real(real64) :: G(n_spc)
      integer :: idx
      idx = merge(1, 2, temp > 1000.d0)
      G(:) = m_gibbs(spcs(:), temp, idx)
      kp(1) = exp(2.d0*(G(CO) - G(CO2)) + G(O2))
      kp(2) = exp(2.d0*(G(H2) - G(H2O)) + G(O2))
      kp(3) = exp(2.d0*G(OH) - (G(O2) + G(H2)))
      kp(4) = exp(2.d0*G(O) - G(O2))
      kp(5) = exp(2.d0*G(H) - G(H2))
   end subroutine calc_eq_const

   ! constant pressure, equilibrium
   ! all the input is non-dimensional value
   subroutine calc_equilibrium(frac, kp, r_carb, r_hydro, p_ref)
      use m_numerics, only: LUBKSB, LUDCMP
      real(real64), intent(in) :: r_carb, r_hydro, p_ref, kp(n_reaction)
      real(real64), intent(inout) :: frac(n_spc)
      real(real64) :: rp_ref, deno, eqn(n_spc), jaco(n_spc, n_spc)
      integer :: indx(n_spc)
      rp_ref = 1.d0/p_ref
      jaco(:, :) = 0.d0
      ! need comment
      !     Species Numbering  H2O(1) CO2(2) CO(3) OH(4) O2(5) H2(6) O(7) H(8)
      deno = frac(H2O) + frac(CO) + frac(OH) + frac(O) + 2.d0*(frac(CO2) + frac(O2))

      eqn(1) = r_carb*deno - (frac(CO2) + frac(CO)) ! C/O ratio
      jaco(1, :) = [r_carb, 2.d0*r_carb - 1.d0, r_carb - 1.d0, r_carb, 2.d0*r_carb, 0.d0, r_carb, 0.d0]

      eqn(2) = r_hydro*deno & ! H/ O ratio
               - (frac(OH) + frac(H) + 2.d0*(frac(H2O) + frac(H2)))
      jaco(2, :) = [r_hydro - 2.d0, 2.d0*r_hydro, r_hydro, r_hydro - 1.d0, 2.d0*r_hydro, -2.d0, r_hydro, -1.d0]

      eqn(3) = sum(frac(:)) - 1.d0
      jaco(3, :) = 1.d0

      eqn(4) = (frac(CO)**2)*frac(O2) & ! C, O atom
               - (kp(1)*rp_ref)*(frac(CO2)**2)
      jaco(4, :) = [0.d0, -2.d0*frac(CO2)*(kp(1)*rp_ref), 2.d0*frac(CO)*frac(O2), 0.d0, frac(CO)**2, 0.d0, 0.d0, 0.d0]

      eqn(5) = (frac(H2)**2)*frac(O2) - (kp(2)*rp_ref)*(frac(H2O)**2) ! H, O atom
      jaco(5, 1) = -2.d0*frac(H2O)*(kp(2)*rp_ref); jaco(5, 6) = 2.d0*frac(O2)*frac(H2)

      eqn(6) = frac(OH)**2 - kp(3)*frac(O2)*frac(H2) ! O, H atom
      jaco(6, 4) = 2.d0*frac(4); jaco(6, 5) = -frac(6)*kp(3); jaco(6, 6) = -frac(5)*kp(3)

      eqn(7) = frac(O)**2 - (kp(4)*rp_ref)*frac(O2) ! O atom
      jaco(7, 5) = -(kp(5)*rp_ref); jaco(7, 7) = 2.d0*frac(O)

      eqn(8) = frac(H)**2 - (kp(5)*rp_ref)*frac(H2) ! H atom
      jaco(8, 6) = -(kp(4)*rp_ref); jaco(8, 8) = 2.d0*frac(H)

      call LUDCMP(jaco, n_spc, n_spc, indx)
      call LUBKSB(jaco, n_spc, n_spc, indx, eqn)

      frac(:) = frac(:) - eqn(:)*0.1
   end subroutine calc_equilibrium

   subroutine fuel_property(fuels, frac, temp, avg_mass, H, n_carb, n_hydro, n_fuel)
      integer, intent(in) :: n_fuel
      class(Specie), intent(in) :: fuels(n_fuel)
      real(real64), intent(in) :: frac(n_fuel), temp
      real(real64), intent(out) :: avg_mass, n_carb, n_hydro, H
      real(real64) :: Hs(n_fuel)
      integer :: idx

      avg_mass = sum(fuels(:)%mole_weight*frac(:)) ! fuel and mole fraction in same order

      idx = merge(1, 2, temp > 1000.d0)
      Hs(:) = enthalpy(fuels(:), temp, idx)
      H = sum(Hs(:)*frac(:))

      ! this considering only CxHy species
      n_carb = sum(fuels(:)%comp(1)*frac(:))
      n_hydro = sum(fuels(:)%comp(2)*frac(:))
   end subroutine fuel_property

   pure function temperature(spcs, T_guess, frac, H_init) result(T)
      class(Specie), intent(in) :: spcs(n_spc)
      real(real64), intent(in) :: T_guess, frac(n_spc), H_init
      real(real64) :: T, H, Cp, delta, f1, f2, f3, f4, f5, shf, &
                      fp1, fp2, fp3, fp4, fp5
      integer :: idx
      T = T_guess

      idx = merge(1, 2, T >= 1000.d0)

      f1 = sum(frac(:)*spcs(:)%poly_h(1, idx))
      f2 = sum(frac(:)*spcs(:)%poly_h(2, idx))
      f3 = sum(frac(:)*spcs(:)%poly_h(3, idx))
      f4 = sum(frac(:)*spcs(:)%poly_h(4, idx))
      f5 = sum(frac(:)*spcs(:)%poly_h(5, idx))

      fp1 = sum(frac(:)*spcs(:)%poly_h(1, idx))
      fp2 = sum(frac(:)*spcs(:)%poly_h(2, idx)*2.)
      fp3 = sum(frac(:)*spcs(:)%poly_h(3, idx)*3.)
      fp4 = sum(frac(:)*spcs(:)%poly_h(4, idx)*4.)
      fp5 = sum(frac(:)*spcs(:)%poly_h(5, idx)*5.)
      shf = sum(frac(:)*spcs(:)%poly_h(6, idx))
      ! H = sum(frac(:) * T * enthalpy(spcs(:), T, idx)) - H_init
      H = t*(f1 + t*(f2 + t*(f3 + t*(f4 + t*f5)))) + shf - H_init
      Cp = fp1 + t*(fp2 + t*(fp3 + t*(fp4 + t*fp5)))

      delta = -(H/Cp)
      T = T + delta

      do while (abs(delta) >= 1e-1)
         ! idx = merge(1, 2, T >= 1000.d0)
         H = t*(f1 + t*(f2 + t*(f3 + t*(f4 + t*f5)))) + shf - H_init
         ! H = sum(frac(:) * enthalpy(spcs(:), T, idx)* T)  - H_init
         Cp = fp1 + t*(fp2 + t*(fp3 + t*(fp4 + t*fp5)))
         delta = -(H/Cp)
         T = T + delta
      end do
   end function temperature

end module m_equil
