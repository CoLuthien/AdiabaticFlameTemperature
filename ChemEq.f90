! Program to solve Genenral CxHy-O2 Equilibrium
! In this program, we are only dealing with these species
! Species Numbering  H2O, CO2, CO, OH, O2, H2, O, H
program ChemEq
   use, intrinsic :: iso_fortran_env
   use m_specie
   use m_problem
   use m_equil
   implicit none
   type(Specie), allocatable :: hydro_carbons(:)
   type(Specie) :: spcs(8)
   real(real64) :: mole_fraction(8), Pinit, Tinit, eps, eps_t, eps_to, &
                   rel_x, rel_t, kp(n_reaction), chwm, H_fuel, H_oxi, H_mix, &
                   rng(3), phi, n_carb, n_hydro, n_oxi, r_co, r_oh, &
                   b, c, d, T_old, T_cur ! will be change.. too tired
   real(real64), allocatable :: fuel_fractions(:)
   integer :: fd, n_fuel, iter, sub_iter, i, j, iter_max
   iter_max = 10000
   fd = abs(fd)
   open (fd, file='EQ.inp')

   read (fd, *) Pinit, Tinit
   read (fd, *) rng(:3)
   read (fd, *) n_fuel
   allocate (hydro_carbons(n_fuel))
   allocate (fuel_fractions(n_fuel))
   read (fd, *) fuel_fractions(:n_fuel)
   read (fd, *) eps, eps_t, eps_to
   read (fd, *) rel_x, rel_t
   read (fd, *) iter, sub_iter
   open (2, file='CXHYAD2.OUT', status='unknown')
   write (2, 1000)

   call read_thermo(hydro_carbons(1), 'ch4.dat')
   call read_thermo(hydro_carbons(2), 'c2h4.dat')
   call read_prop(spcs)

   hydro_carbons(1)%mole_weight = 16.04303d0
   hydro_carbons(2)%mole_weight = 28.05418d0
   spcs(1)%mole_weight = 18.016d0 ! H2O
   spcs(2)%mole_weight = 44.011d0 ! CO2
   spcs(3)%mole_weight = 28.010d0 ! CO
   spcs(4)%mole_weight = 17.007d0 ! OH
   spcs(5)%mole_weight = 31.999d0 ! O2
   spcs(6)%mole_weight = 2.016d0 ! H2
   spcs(7)%mole_weight = 16.000d0 ! O
   spcs(8)%mole_weight = 1.088d0 ! H

   call fuel_property(hydro_carbons, fuel_fractions, Tinit, chwm, H_fuel, n_carb, n_hydro, n_fuel)

   H_oxi = enthalpy(spcs(O2), Tinit, merge(1, 2, Tinit >= 1000.d0))
   print *, Tinit, H_oxi, H_fuel

   do i = 1, 10
      phi = 0.2d0 + 0.2d0*dble(i - 1)
      n_oxi = (n_carb + n_hydro*0.25d0)/phi
      H_mix = (n_oxi/(n_oxi + 1.d0))*(H_oxi*Tinit) &
              + (1.d0/(n_oxi + 1.d0))*(H_fuel*Tinit)

      r_co = n_carb/(2.d0*n_oxi)
      r_oh = n_hydro/(2.d0*n_oxi)

      ! initial guess of mole fraction
      mole_fraction(:) = 1.d-5
      b = n_carb + n_hydro*0.25d0
      c = n_oxi + n_hydro*0.25d0
      d = n_carb + n_hydro*0.5d0
      if (n_oxi >= b) then !if lean condition
         mole_fraction(H2O) = (n_hydro*0.5d0)/c
         mole_fraction(CO2) = n_carb/c
         mole_fraction(O2) = (n_oxi - n_carb - n_hydro*0.25d0)/c
      else ! fuel rich condition
         mole_fraction(H2O) = (n_hydro*0.5d0)/d
         mole_fraction(CO2) = (2.d0*n_oxi - n_carb - n_hydro*0.5d0)/d
         mole_fraction(CO) = (2.d0*n_carb + n_hydro*0.5d0 - 2.d0*n_oxi)/d
      end if

      T_cur = 300.0
      do j = 1, iter_max
         T_old = T_cur
         call calc_eq_const(spcs(:), T_cur, kp)
         call calc_equilibrium(mole_fraction(:), kp(:), r_co, r_oh, Pinit)

         where (mole_fraction <= 0.d0)
            mole_fraction = 0.
         end where

         T_cur = temperature(spcs(:), T_old, mole_fraction(:), H_mix)
         if (abs((T_cur - T_old)) <= eps_t) then
            exit
         end if
      end do

      print *, phi, T_cur
      write (2, 1100) Pinit, phi, T_cur, mole_fraction(:)
      ! print*, mole_fraction
   end do
   close (2)
   print*, STORAGE_SIZE(phi)

1000 format(' p     phi    Temp.(K)      H2O         CO2          CO&
           &          OH          O2          H2          O          H')
1100 format(f4.1, 2x, f5.2, 2x, f8.2, 2x, 2x, 8(F10.5, 2x), i4)

contains

   subroutine read_prop(spcs)
      class(Specie) :: spcs(n_spc)
      call read_thermo(spcs(H2O), 'h2o.dat')
      call read_thermo(spcs(CO2), 'co2.dat')
      call read_thermo(spcs(CO), 'co.dat')
      call read_thermo(spcs(OH), 'oh.dat')
      call read_thermo(spcs(O2), 'o2.dat')
      call read_thermo(spcs(H2), 'h2.dat')
      call read_thermo(spcs(O), 'o.dat')
      call read_thermo(spcs(H), 'h.dat')
   end subroutine read_prop

end program ChemEq
