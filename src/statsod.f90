!*******************************************************************************
!    Copyright (c) 2014 Ricardo Grau-Crespo and co-authors
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
  use ensemble_io,     only: read_energies_file
  implicit none

  integer, parameter :: nconfmax = 1000000, ncolmax = 10, npointsmax = 800, ntempmax = 1000
  real(real64), parameter :: kb = 8.61734e-5_real64, tolprob = 1.0e-12_real64, tolminspec = 1.0e-6_real64
  integer :: m, auxm, mm, ncol, npoints, col, tt, ntt, point, ios
  integer :: colon_pos, kpos, kpos2
  integer :: memin, memax, n_missing
  logical :: ene_ok
  real(real64)  :: emin, emax, maxspec
  real(real64), dimension(ntempmax) :: z, e, f, s, t
  real(real64) :: zinf, einf, sinf
  real(real64) :: tsampling, sum_omega
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
  character(len=200) :: ensemble_line
  character(len=20) :: arg
  logical :: temperatures_exists, data_exists, spectra_exists, quiet
  logical :: metropolis_sample

  tsampling = -1.0_real64
  metropolis_sample = .false.
  quiet = .false.
  if (command_argument_count() > 0) then
    call get_command_argument(1, arg)
    if (trim(arg) == '-q') quiet = .true.
  end if

  if (.not. quiet) then
    write (*, '(A)') "SOD (Site-Occupancy Disorder) version 0.83 - statsod"
  end if
  write (*, *) " > Statistical mechanics analysis..."
  write (*, *) ""

!Input files

  inquire (file="TEMPERATURES", exist=temperatures_exists)
  if (temperatures_exists) then
    write (*, *) " > TEMPERATURES file found."
    open (unit=10, file="TEMPERATURES")
  end if

  open (unit=11, file="ENSEMBLE")

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
  !      Read the ENSEMBLE file (supports v2 and v3 formats)
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! Read first non-blank line
  do
    read (11, '(a)', iostat=ios) ensemble_line
    if (ios /= 0) then
      write (*, *) "Error: could not read ENSEMBLE."
      stop 1
    end if
    if (len_trim(ensemble_line) > 0) exit
  end do

  if (ensemble_line(1:1) == '#') then
    ! v2: first non-blank line starts with '#'
    do while (ensemble_line(1:1) == '#')
      if (index(ensemble_line, 'Sampling temperature') > 0) then
        colon_pos = index(ensemble_line, ':')
        if (colon_pos > 0 .and. colon_pos < len_trim(ensemble_line)) then
          read (ensemble_line(colon_pos+1:), *, iostat=ios) tsampling
          if (ios /= 0) tsampling = -1.0_real64
        end if
      end if
      read (11, '(a)', iostat=ios) ensemble_line
      if (ios /= 0) then
        write (*, *) "Error: could not read ENSEMBLE."
        stop 1
      end if
    end do
    ! ensemble_line is "N substitutions in M sites" — skip it, read nic next
    read (11, *) mm
    do m = 1, mm
      read (11, *) auxm, omega(m)
    end do
  else
    ! v3: first line is "<Type> ensemble [...]: nic configurations"
    if (index(ensemble_line, 'Metropolis') > 0) then
      kpos  = index(ensemble_line, '(')
      kpos2 = index(ensemble_line, 'K)')
      if (kpos > 0 .and. kpos2 > kpos) then
        read (ensemble_line(kpos+1:kpos2-1), *, iostat=ios) tsampling
        if (ios /= 0) tsampling = -1.0_real64
      end if
    end if
    colon_pos = index(ensemble_line, ':', back=.true.)
    kpos      = index(ensemble_line, 'configurations')
    if (colon_pos > 0 .and. kpos > colon_pos) then
      read (ensemble_line(colon_pos+1:kpos-1), *, iostat=ios) mm
      if (ios /= 0) then
        write (*, *) "Error: could not parse ENSEMBLE configuration count."
        stop 1
      end if
    else
      write (*, *) "Error: could not parse ENSEMBLE header."
      stop 1
    end if
    ! Skip target lines (contain 'sites' and '->') and column-header comment ('#')
    do
      read (11, '(a)', iostat=ios) ensemble_line
      if (ios /= 0) then
        write (*, *) "Error: could not read ENSEMBLE data."
        stop 1
      end if
      if (len_trim(ensemble_line) == 0) cycle
      if (ensemble_line(1:1) == '#') cycle
      if (index(ensemble_line, 'sites') > 0 .and. index(ensemble_line, '->') > 0) cycle
      exit
    end do
    ! ensemble_line holds first data row; read omega(1) from it, then continue
    read (ensemble_line, *, iostat=ios) auxm, omega(1)
    if (ios /= 0) then
      write (*, *) "Error: could not read ENSEMBLE data row."
      stop 1
    end if
    do m = 2, mm
      read (11, *, iostat=ios) auxm, omega(m)
      if (ios /= 0) then
        write (*, *) "Error: could not read ENSEMBLE data row."
        stop 1
      end if
    end do
  end if
  close (11)

  metropolis_sample = (tsampling >= 0.0_real64)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Read the TEMPERATURES file, unless ENSEMBLE is already Metropolis biased
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if (metropolis_sample) then
    if (temperatures_exists) then
      close (10)
      write (*, *) "Warning: Metropolis ENSEMBLE detected; TEMPERATURES is ignored."
    end if
    t(1) = tsampling
    ntt = 1
    write (*, *) " > Metropolis ENSEMBLE detected."
    write (*, *) " > statsod probabilities use Omega/sum(Omega) at Tsampling = ", tsampling, " K."
  else
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
  end if

  allocate (p(mm, ntt))
  allocate (data(ncolmax, mm))

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      Read ENERGIES
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  call read_energies_file("ENERGIES", mm, ene, ene_ok, n_missing)
  if (.not. ene_ok) then
    write (*, *) "Error: could not open or read ENERGIES."
    stop 1
  end if
  if (n_missing > 0) then
    write (*, '(A,I0,A)') "Error: missing energies for ", n_missing, " configuration(s) in ENERGIES."
    stop 1
  end if

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
  if (metropolis_sample) then
    write (21, *) "# Metropolis-sampled ENSEMBLE: probabilities use Omega/sum(Omega); F and S are not evaluated."
  end if

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
    if (metropolis_sample) then
      write (23, *) "x   ", t(1:ntt)
    else
      write (23, *) "x   ", t(1:ntt), "    Infinity"
    end if
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
  if (metropolis_sample) then
    write (20, *) "Metropolis sampling temperature (K): ", tsampling
    write (20, *) "Probabilities use Omega/sum(Omega); ENERGIES are not Boltzmann-weighted again."
    write (20, *)
  end if

  do m = 1, mm
    erel(m) = ene(m) - emin
  end do

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !       This starts the loop over all temperature values
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  do tt = 1, ntt

    if (metropolis_sample) then

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!          Metropolis-sampled ENSEMBLE: the sample already contains energy bias
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      sum_omega = sum(real(omega(1:mm), real64))
      if (sum_omega <= 0.0_real64) then
        write (*, *) "Error: non-positive total Omega in ENSEMBLE."
        stop 1
      end if

      do m = 1, mm
        p(m, tt) = real(omega(m), real64)/sum_omega
      end do

      e(tt) = 0.0_real64
      do m = 1, mm
        e(tt) = e(tt) + ene(m)*p(m, tt)
      end do

      f(tt) = e(tt)
      s(tt) = 0.0_real64

      if (data_exists) then
        avedata(1:ncol, tt) = 0.0_real64
        do col = 1, ncol
          do m = 1, mm
            avedata(col, tt) = avedata(col, tt) + data(col, m)*p(m, tt)
          end do
        end do
      end if

      if (spectra_exists) then
        avespec(1:npoints, tt) = 0.0_real64
        do point = 1, npoints
          do m = 1, mm
            avespec(point, tt) = avespec(point, tt) + spec(point, m)*p(m, tt)
          end do
          if (avespec(point, tt)/maxspec < tolminspec) then
            avespec(point, tt) = 0.0
          end if
        end do
      end if

    else if (t(tt) == 0.0_real64) then

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

  if (.not. metropolis_sample) then
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
301   format(a10, 2x, 10(f10.4, 2x))
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
  end if

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !       Writing ave_spectra for all temperatures
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if (spectra_exists) then
    if (metropolis_sample) then
      write (fmtemplist, '(a,i0,a)') '(f10.3,2x,', ntt, '(e12.6,2x))'
      do point = 1, npoints
        write (23, fmtemplist) xspec(point), avespec(point, 1:ntt)
      end do
    else
      write (fmtemplist, '(a,i0,a)') '(f10.3,2x,', ntt + 1, '(e12.6,2x))'
      do point = 1, npoints
        write (23, fmtemplist) xspec(point), avespec(point, 1:ntt), avespecinf(point)
      end do
    end if
  end if

  close (20)
  close (21)
  if (data_exists) close (22)
  if (spectra_exists) close (23)

  write (*, *) " > Done!"
  write (*, *) ""

end
