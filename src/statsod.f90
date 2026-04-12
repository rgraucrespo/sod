!*******************************************************************************
!    Copyright (c) 2014 Ricardo Grau-Crespo, Said Hamad
!
!    This file is part of the SOD package.
!
!    SOD is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SOD is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SOD.  If not, see <http://www.gnu.org/licenses/>.
!
!******************************************************************************

program stats

  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: nconfmax = 1000000, ncolmax = 10, npointsmax = 800, ntempmax = 1000
  real(real64), parameter :: kb = 8.61734e-5_real64, tolprob = 1.0e-12_real64, tolminspec = 1.0e-6_real64
  integer :: m, auxm, mm, ncol, npoints, col, tt, ntt, nsubs, point, ios
  integer :: memin, memax
  real(real64)  :: emin, emax, maxspec
  real(real64), dimension(ntempmax) :: z, e, f, s, t
  real(real64) :: zinf, einf, sinf
  real(real64), dimension(nconfmax) :: ene, erel
  real(real64), dimension(:,:), allocatable :: p
  real(real64), dimension(nconfmax) :: pinf
  integer, dimension(nconfmax) :: omega
  real(real64), dimension(:,:), allocatable :: data
  real(real64), dimension(ncolmax, ntempmax) :: avedata
  real(real64), dimension(:,:), allocatable :: spec
  real(real64), dimension(npointsmax, ntempmax) :: avespec
  real(real64), dimension(ncolmax) :: avedatainf
  real(real64), dimension(npointsmax) :: xspec, avespecinf
  character(len=30) fmtemplist
  character(len=200) :: outsod_line
  character(len=20) :: arg
  logical :: temperatures_exists, data_exists, spectra_exists, quiet

  quiet = .false.
  if (command_argument_count() > 0) then
    call get_command_argument(1, arg)
    if (trim(arg) == '-q') quiet = .true.
  end if

  if (.not. quiet) then
    write (*, *) "============================================================================"
    write (*, *) "         SOD (Site Occupancy Disorder) version 0.71"
    write (*, *) ""
    write (*, *) "         Authors: R. Grau-Crespo and S. Hamad"
    write (*, *) "         Contact: <r.grau-crespo@qmul.ac.uk>"
    write (*, *) "============================================================================"
    write (*, *) ""
  end if
  write (*, *) " > Statistical mechanics analysis..."
  write (*, *) ""

!Input files

  inquire (file="TEMPERATURES", exist=temperatures_exists)
  if (temperatures_exists) then
    write (*, *) " > TEMPERATURES file found."
    open (unit=10, file="TEMPERATURES")
  end if

  open (unit=11, file="OUTSOD")
  open (unit=12, file="ENERGIES")

  inquire (file="DATA", exist=data_exists)
  if (data_exists) then
    write (*, *) " > DATA file found: averaging of observables will be performed."
    open (unit=13, file="DATA")
  end if

  inquire (file="SPECTRA", exist=spectra_exists)
  if (spectra_exists) then
    write (*, *) " > SPECTRA file found: averaging of spectra will be performed."
    open (unit=14, file="SPECTRA")
    open (unit=15, file="XSPEC")
  end if
  write (*, *) ""

!Output files
  open (unit=20, file="probabilities.dat")
  open (unit=21, file="thermodynamics.dat")
  write (*, *) " > Writing probabilities.dat and thermodynamics.dat..."

  if (data_exists) then
    open (unit=22, file="ave_data.dat")
  end if

  if (spectra_exists) then
    open (unit=23, file="ave_spectra.dat")
  end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Read the TEMPERATURES files
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if (temperatures_exists) then
    tt = 1
    do
      read (10, *, iostat=ios) t(tt)
      if (ios /= 0) exit
      tt = tt + 1
    end do
    ntt = tt - 1
    close (10)
  else
    t(1) = 0.0_real64
    t(2) = 300.0_real64
    t(3) = 1000.0_real64
    ntt = 3
  end if

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      Read the OUTSOD file, which is the output from SOD,
  !      giving configuration numbers and degeneracies
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  do
    read (11, '(a)') outsod_line
    if (outsod_line(1:1) /= '#') exit
  end do
  read (outsod_line, *) nsubs
  read (11, *) mm
  do m = 1, mm
    read (11, *) auxm, omega(m)
  end do
  close (11)
  allocate (p(mm, ntt))
  allocate (data(ncolmax, mm))

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      Read ENERGIES
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  do m = 1, mm
    read (12, *) ene(m)
  end do
  close (12)

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      Read the data file, DATA
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if (data_exists) then

    read (13, *) ncol
    do m = 1, mm
      read (13, *) data(1:ncol, m)
    end do
    close (13)
  end if

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      Read the SPECTRA and XSPEC files, and calculate maximum spectrum intensity
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if (spectra_exists) then
    write (*, *) "Reading SPECTRA file..."
    write (*, *)
    read (14, *) npoints
    allocate (spec(npoints, mm))
    do m = 1, mm
      read (14, *) spec(1:npoints, m)
    end do
    close (14)

    maxspec = maxval(spec)

    do point = 1, npoints
      read (15, *) xspec(point)
    end do
    close (15)

  end if

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !       Start writing thermodynamics.dat file
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  write (21, *) "      T/K            E/eV            F/eV         S/(eV/K)"

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !       Start writing ave_data.dat file
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if (data_exists) then
    write (22, *) "       T    Average data"
  end if

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !       Start writing ave_spectra.dat file
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if (spectra_exists) then
    write (23, *) "x   ", t(1:ntt), "    Infinity"
  end if

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !          Get minimum and maximum energies and relative energies with
  !          respect to minimum in order to get better accuracy for the
  !          exponential calculations
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  emin = ene(1)
  emax = ene(1)
  memin = 1
  memax = 1
  do m = 2, mm
    if (emin > ene(m)) then
      emin = ene(m)
      memin = m
    end if
    if (emax < ene(m)) then
      emax = ene(m)
      memax = m
    end if
  end do

  write (20, *) "Configuration with minimum energy: ", memin
  write (20, *) "Configuration with maximum energy: ", memax
  write (20, *)

  do m = 1, mm
    erel(m) = ene(m) - emin
  end do

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !       This starts the loop over all temperature values
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  do tt = 1, ntt

    if (t(tt) == 0.0_real64) then

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!          T=0 K: trivial case — all probability on the lowest-energy config
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      p(1:mm, tt) = 0.0_real64
      p(memin, tt) = 1.0_real64
      e(tt) = emin
      f(tt) = emin
      s(tt) = 0.0_real64

      if (data_exists) then
        avedata(1:ncol, tt) = data(1:ncol, memin)
      end if

      if (spectra_exists) then
        avespec(1:npoints, tt) = spec(1:npoints, memin)
      end if

    else

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!          Calculate the partition function and probabilities
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      z(tt) = 0.0
      do m = 1, mm
        z(tt) = z(tt) + omega(m)*exp(-erel(m)/(kb*t(tt)))
      end do

      do m = 1, mm
        p(m, tt) = omega(m)*exp(-erel(m)/(kb*t(tt)))/z(tt)
      end do

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculate the energy (E) and Helmholtz free energy (F)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      e(tt) = 0.0
      s(tt) = 0.0

      do m = 1, mm
        e(tt) = e(tt) + ene(m)*p(m, tt)
      end do

      f(tt) = emin - kb*t(tt)*log(z(tt))

      s(tt) = (e(tt) - f(tt))/t(tt)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculate the average data
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (data_exists) then

        avedata(1:ncol, tt) = 0.0

        do col = 1, ncol
          do m = 1, mm
            avedata(col, tt) = avedata(col, tt) + data(col, m)*p(m, tt)
          end do
        end do

      end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculate the average spectra
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (spectra_exists) then

        avespec(1:npoints, tt) = 0.0

        do point = 1, npoints
          do m = 1, mm
            avespec(point, tt) = avespec(point, tt) + spec(point, m)*p(m, tt)
          end do
          if (avespec(point, tt)/maxspec < tolminspec) then
            avespec(point, tt) = 0.0
          end if
        end do

      end if

    end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!          Write PROBABILITIES for each temperature T
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    write (20, 100) "Temperature=", t(tt)
100 format(a12, 2x, f10.4)
    write (20, *) "        m  omega(m)   Erel(m)/eV       p(m)        p(m)/omega(m)"
    do m = 1, mm
      write (20, 101) m, omega(m), erel(m), p(m, tt), p(m, tt)/omega(m)
101   format(i10, 2x, i8, 2x, e12.6, 2x, f10.4, 2x, e12.4)
    end do

    write (20, *)
    write (20, *)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!          Write thermodynamic information for each temperature T
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    write (21, 200) t(tt), e(tt), f(tt), s(tt)
200 format(f10.1, 2x, 2(f14.4, 2x), e12.6)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       Writing ave_data for each temperature
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    if (data_exists) then
      write (22, 201) t(tt), avedata(1:ncol, tt)
201   format(f10.1, 2x, 10(f10.4, 2x))
    end if

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !       This ends the loop over all temperature values
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  end do

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !       Calculating results in the limit of infinite temperature
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  zinf = 0.0
  einf = 0.0
  do m = 1, mm
    zinf = zinf + omega(m)
  end do

  do m = 1, mm
    pinf(m) = omega(m)/zinf
  end do

  do m = 1, mm
    einf = einf + ene(m)*pinf(m)
  end do

  sinf = kb*log(zinf)

  write (20, *) "Infinite Temperature Limit"
  write (20, *) "        m  omega(m)   Erel(m)/eV       p(m)      p(m)/omega(m)"
  do m = 1, mm
    write (20, 101) m, omega(m), erel(m), pinf(m), pinf(m)/omega(m)
  end do

  write (21, 300) "Infinite", einf, " - ", sinf
300 format(a10, 2x, f14.4, 8x, a3, 7x, e12.6)

  if (data_exists) then
    avedatainf(1:ncol) = 0.0
    do m = 1, mm
      do col = 1, ncol
        avedatainf(col) = avedatainf(col) + data(col, m)*pinf(m)
      end do
    end do
    write (22, 301) adjustr("Infinite"), avedatainf(1:ncol)
301 format(a10, 2x, 10(f10.4, 2x))
  end if

  if (spectra_exists) then
    avespecinf(1:npoints) = 0.0
    do point = 1, npoints
      do m = 1, mm
        avespecinf(point) = avespecinf(point) + spec(point, m)*pinf(m)
      end do
      if (avespecinf(point)/maxspec < tolminspec) then
        avespecinf(point) = 0.0
      end if
    end do
  end if

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !       Writing ave_spectra for all temperatures
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if (spectra_exists) then
    write (fmtemplist, '(a,i0,a)') '(f10.3,2x,', ntt + 1, '(e12.6,2x))'
    do point = 1, npoints
      write (23, fmtemplist) xspec(point), avespec(point, 1:ntt), avespecinf(point)
    end do
  end if

  close (20)
  close (21)
  if (data_exists) close (22)
  if (spectra_exists) close (23)

  write (*, *) " > Done!"
  write (*, *) ""

end

