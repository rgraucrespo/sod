!v*******************************************************************************
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

program gcstatsod

  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: nconfmax = 1000000, ncolmax = 10, npointsmax = 800, ntempmax = 1000, noutsodsmax = 1000, nbismax = 1000
  real(real64), parameter :: kb = 8.61734e-5_real64, tolprob = 1.0e-12_real64, tolminspec = 1.0e-6_real64, eva3togpa = 160.2176621_real64
  integer :: m, auxm, ncol, col, tt, ntt, nsubs, nbis, npoints, point
  integer :: memin, nsubsemin, iostatus
  real(real64)  :: emin, maxspec
  integer :: nsubsmin, nsubsmax, nconfigmax
  real(real64), dimension(:), allocatable :: z, e, f, s, t, gpot
  real(real64) :: einf, sinf
  real(real64), dimension(:, :, :), allocatable :: data
  real(real64), dimension(ncolmax, ntempmax) :: avedata
  real(real64), dimension(:, :, :), allocatable :: spec
  real(real64), dimension(npointsmax, ntempmax) :: avespec
  real(real64), dimension(ncolmax) :: avedatainf
  real(real64), dimension(:), allocatable :: xspec, avespecinf

  logical :: temperatures_exists, data_exists, ingc_exists, spectra_exists, booldata, boolspec
  real(real64) :: xormuvalue, x, xeq, mu, nx, oldmu
  character(len=2) :: xormu
  character(len=9) :: filenameout
  character(len=15) :: auxstring
  character(len=11) :: filenameene
  character(len=200) :: outsod_line
  character(len=9), dimension(:), allocatable :: filenamedat
  character(len=12), dimension(:), allocatable :: filenamespec
  integer :: nsubsread, npos
  integer, dimension(:, :), allocatable :: omega
  real(real64), dimension(:, :), allocatable :: ene, enemun, enemunrel
  integer, dimension(:), allocatable:: mm
!REAL (kind=8),DIMENSION(0:NOUTSODSMAX,NCONFMAX,NTEMPMAX)  :: p
  real(real64), dimension(:, :, :), allocatable  :: p
  real(real64), dimension(0:noutsodsmax, ntempmax)  :: pn
  integer :: omegasum
  real(real64)  :: naver, eneaver, epsilona, epsilonb, epsilon, e0
  real(real64)  :: a0, a1, a2, b0, b1, b2
  real(real64)  :: momentaa0, momentaa1, momentaa2, momentab0, momentab1, momentab2
  real(real64)  :: a11, a12, a21, a22, bb1, bb2, alpha, beta
  real(real64), dimension(0:noutsodsmax)  :: qn, c, tau
  real(real64), dimension(:, :), allocatable  :: pinf
  real(real64), parameter :: tolmu = 1.0e-10, tolq = 1.0e-10
  real(real64) valpolynom, qtest
  real(real64) valpolynomnew, r1, r2, qa, qb, q, va, vb
  real(real64) lambda, v0, v1, bv, bm0, bm1, bb, vol, bm, voln, xn, eta
  real(real64), dimension(0:noutsodsmax)  :: evsc

  write (*, *) "============================================================================"
  write (*, *) "         SOD (Site Occupancy Disorder) version 0.70"
  write (*, *) ""
  write (*, *) "         Authors: R. Grau-Crespo and S. Hamad"
  write (*, *) "         Contact: <r.grau-crespo@qmul.ac.uk>"
  write (*, *) "============================================================================"
  write (*, *) ""
  write (*, *) " > Grand-canonical statistical analysis..."
  write (*, *) ""

  !Reading the TEMPERATURES files (checking first if it exists)

  inquire (file="TEMPERATURES", exist=temperatures_exists)

  ntt = 2 !this is the value by defect if TEMPERATURES does not exist

  if (temperatures_exists) then
    write (*, *) " > TEMPERATURES file found."
    open (unit=10, file="TEMPERATURES", status='OLD', action='READ', iostat=iostatus)
    ntt = 0
    do
      read (10, '(A)', iostat=iostatus) auxstring
      if (iostatus /= 0) exit
      ntt = ntt + 1
    end do
  else
    write (*, *) " > TEMPERATURES file not found: analysis will be done at 300 K and 1000 K."
  end if

  call random_seed()

!Allocating arrays with tt as index
  allocate (t(1:ntt))
  allocate (z(1:ntt))
  allocate (e(1:ntt))
  allocate (f(1:ntt))
  allocate (s(1:ntt))
  allocate (gpot(1:ntt))

  if (temperatures_exists) then
    rewind (10)
    do tt = 1, ntt
      read (10, *) t(tt)
    end do
    close (10)
  else
    t(1) = 300.0
    t(2) = 1000.0
  end if

  !Reading the INGC input file

  inquire (file="INGC", exist=ingc_exists)
  if (ingc_exists) then
    write (*, *) " > INGC file found."
    open (unit=14, file="INGC")
  else
    write (*, *) "ERROR: INGC file does not exist, but it is needed for grand-canonical analysis."
    stop 1
  end if

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      Read the INGC file
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  read (14, *)
  read (14, *) nsubsmin, nsubsmax
  read (14, *)
  read (14, *) xormu, xormuvalue
  read (14, *)
  read (14, *) lambda
  if (lambda > 0) then
    read (14, *)
    read (14, *) v0, v1, bv
    read (14, *)
    read (14, *) bm0, bm1, bb
  end if
  close (14)

  if (xormu == "x ") then
    x = xormuvalue
  else
    mu = xormuvalue
  end if

  !Allocating arrays including nsubs as index
  allocate (filenamedat(0:nsubsmax))
  allocate (filenamespec(0:nsubsmax))
  allocate (mm(0:nsubsmax))

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      Read the OUTSOD_xx files
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  do nsubs = nsubsmin, nsubsmax
    if (nsubs < 10) then
      write (filenameout, "(a8,I1)") "OUTSOD_0", nsubs
      open (unit=100 + nsubs, file=trim(filenameout), status='old')
    else
      write (filenameout, "(a7,I2)") "OUTSOD_", nsubs
      open (unit=100 + nsubs, file=trim(filenameout), status='old')
    end if
  end do

  do nsubs = nsubsmin, nsubsmax

    do
      read (100 + nsubs, '(a)') outsod_line
      if (outsod_line(1:1) /= '#') exit
    end do
    read (outsod_line, *) nsubsread, auxstring
    if (trim(auxstring) /= 'substitutions') then
      write (*, *) 'ERROR: grand-canonical analysis supports only single-site binary substitutions.'
      write (*, *) '       Multi-species or multi-nary OUTSOD detected. Aborting.'
      stop 1
    end if
    read (outsod_line, *) nsubsread, auxstring, auxstring, npos
    if (nsubsread /= nsubs) then
      write (*, *) 'ERROR: nsubsread.ne.nsubs'
      stop 1
    end if

    read (100 + nsubs, *) mm(nsubs)
    write (*, '(a, i4, a, i10)') "   - Level nsubs =", nsubs, ": ", mm(nsubs), " independent configurations"
  end do

  nconfigmax = maxval(mm)

  !Allocating arrays including also m as index
  allocate (omega(0:nsubsmax, 1:nconfigmax))
  allocate (ene(0:nsubsmax, 1:nconfigmax))
  allocate (enemun(0:nsubsmax, 1:nconfigmax))
  allocate (enemunrel(0:nsubsmax, 1:nconfigmax))
  allocate (p(0:nsubsmax, 1:nconfigmax, 1:ntt))
  allocate (pinf(0:nsubsmax, 1:nconfigmax))

  do nsubs = nsubsmin, nsubsmax

    do m = 1, mm(nsubs)
      read (100 + nsubs, *) auxm, omega(nsubs, m)
    end do
    close (100 + nsubs)

  end do

  if ((xormu == "x ") .and. ((x < real(nsubsmin)/real(npos)) .or. (x > real(nsubsmax)/real(npos)))) then
    write (*, *) "x is out of range, given the number of substitutions considered"
    write (*, *) "x should be between", real(nsubsmin)/real(npos), "and", real(nsubsmax)/real(npos)
    write (*, *)
    stop 1
  end if

  !Reading ENERGIES_xx files

  do nsubs = nsubsmin, nsubsmax
    if (nsubs < 10) then
      write (filenameene, "(a10,I1)") "ENERGIES_0", nsubs
      open (unit=200 + nsubs, file=trim(filenameene), status='old')
    else
      write (filenameene, "(a9,I2)") "ENERGIES_", nsubs
      open (unit=200 + nsubs, file=trim(filenameene), status='old')
    end if
  end do

  do nsubs = nsubsmin, nsubsmax
    do m = 1, mm(nsubs)
      read (200 + nsubs, *) ene(nsubs, m)
    end do
    close (200 + nsubs)
  end do

  !Reading DATA_xx files

  do nsubs = nsubsmin, nsubsmax
    if (nsubs < 10) then
      write (filenamedat(nsubs), "(a6,I1)") "DATA_0", nsubs
    else
      write (filenamedat(nsubs), "(a5,I2)") "DATA_", nsubs
    end if
  end do

  data_exists = .true.
  do nsubs = nsubsmin, nsubsmax
    inquire (file=filenamedat(nsubs), exist=booldata)
    if (booldata) then
      open (unit=300 + nsubs, file=trim(filenamedat(nsubs)), status='old')
      read (300 + nsubs, *) ncol
    end if
    data_exists = (booldata .and. data_exists)
  end do

  if (data_exists) then
    allocate (data(0:nsubsmax, 1:nconfigmax, 1:ncol))
    write (*, *) " > DATA files found: averaging of observables will be performed."
    do nsubs = nsubsmin, nsubsmax
      do m = 1, mm(nsubs)
        read (300 + nsubs, *) data(nsubs, m, 1:ncol)
      end do
      close (300 + nsubs)
    end do
  end if

  !Reading SPECTRA_xx files

  do nsubs = nsubsmin, nsubsmax
    if (nsubs < 10) then
      write (filenamespec(nsubs), "(a9,I1)") "SPECTRA_0", nsubs
    else
      write (filenamespec(nsubs), "(a8,I2)") "SPECTRA_", nsubs
    end if
  end do

  spectra_exists = .true.
  do nsubs = nsubsmin, nsubsmax
    inquire (file=filenamespec(nsubs), exist=boolspec)
    if (boolspec) then
      open (unit=400 + nsubs, file=trim(filenamespec(nsubs)), status='old')
      read (400 + nsubs, *) npoints
    end if
    spectra_exists = (boolspec .and. spectra_exists)
  end do

  if (spectra_exists) then
    write (*, *) " > SPECTRA files found: averaging of spectra will be performed."
    allocate (spec(0:nsubsmax, 1:nconfigmax, 1:npoints))
    allocate (xspec(1:npoints))
    do nsubs = nsubsmin, nsubsmax
      do m = 1, mm(nsubs)
        read (400 + nsubs, *) spec(nsubs, m, 1:npoints)
      end do
      close (400 + nsubs)
    end do
    open (unit=15, file="XSPEC")
    do point = 1, npoints
      read (15, *) xspec(point)
    end do
    maxspec = maxval(spec)
    close (15)
  end if
  write (*, *) ""

  !Opening output files
  open (unit=20, file="probabilities.dat")
  open (unit=21, file="thermodynamics.dat")
  write (*, *) " > Writing probabilities.dat and thermodynamics.dat..."
  if (data_exists) then
    open (unit=22, file="ave_data.dat")
  end if
  if (spectra_exists) then
    open (unit=23, file="ave_spectra.dat")
  end if

  !Starting to write thermodynamics.dat file
  write (21, *) "      T/K            E/eV            F/eV         S/(eV/K)"

  !Starting to write ave_data.dat file
  if (data_exists) then
    write (22, *) "       T    Average data"
  end if

  !Starting to write ave_spectra.dat file
  if (spectra_exists) then
!       WRITE(fmtemplist,'(a4,i1,a15)') "(a6,",Ntt,adjustl("(f13.1,2x),a12)")
!       write(23,fmtemplist) "x   ",     T(1:Ntt),"    Infinity"
!       FIX THIS
    write (23, *) "x   ", t(1:ntt)
  end if

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !       This starts the loop over all temperature values
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  do tt = 1, ntt  ! loop over all temperature values

    write (*, '(a, f10.2, a)') " > Processing T =", t(tt), " K..."

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      If x is specified, get mu
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    if (xormu == "x ") then  ! if for xormu.eq."x "

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !      Calculate omegasum
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      omegasum = 0

      do nsubs = nsubsmin, nsubsmax
        do m = 1, mm(nsubs)
          omegasum = omegasum + omega(nsubs, m)
        end do
      end do

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !      Calculate naver
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      naver = 0.0
      do nsubs = nsubsmin, nsubsmax
        do m = 1, mm(nsubs)
          naver = naver + omega(nsubs, m)*nsubs
        end do
      end do
      naver = naver/omegasum

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !      Calculate eneaver
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      eneaver = 0.0
      do nsubs = nsubsmin, nsubsmax
        do m = 1, mm(nsubs)
          eneaver = eneaver + omega(nsubs, m)*ene(nsubs, m)
        end do
      end do
      eneaver = eneaver/omegasum

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !      Calculate epsilon
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      epsilona = 0.0
      epsilonb = 0.0
      do nsubs = nsubsmin, nsubsmax
        do m = 1, mm(nsubs)
          epsilona = epsilona + omega(nsubs, m)*(nsubs - naver)*(ene(nsubs, m) - eneaver)
          epsilonb = epsilonb + omega(nsubs, m)*(nsubs - naver)*(nsubs - naver)
        end do
      end do
      epsilon = epsilona/epsilonb

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !      Calculate E0
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      e0 = eneaver - epsilon*naver

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !      Calculate EVSC (Energy of Volume-Stress Correction)
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (lambda == 0.0) then
        ! write (*, *) "No Volume-Stress Correction applied."
        do nsubs = nsubsmin, nsubsmax
          evsc(nsubs) = 0.0
        end do
      else
        ! Calculate Volume at this x
        vol = v0*(1 - x) + v1*x + bv*x*(1 - x)
        bm = (bm0*(1 - x) + bm1*x + bb*x*(1 - x))/eva3togpa
        do nsubs = nsubsmin, nsubsmax
          xn = real(nsubs)/real(npos)
          voln = v0*(1 - xn) + v1*xn + bv*xn*(1 - xn)
          eta = (vol/voln)**(2.0d0/3.0d0)
          evsc(nsubs) = (9.0d0/8.0d0)*vol*bm*(eta - 1.0d0)**2
          ! write (*, *) "Energy of Volume-Stress Correction: nsubs=", nsubs, "voln=", voln, "EVSC(nsubs)=", evsc(nsubs)
        end do
      end if

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !      Calculate Qn
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do nsubs = nsubsmin, nsubsmax
        qn(nsubs) = 0.0
        do m = 1, mm(nsubs)
          qn(nsubs) = qn(nsubs) + real(omega(nsubs, m))*exp(-(ene(nsubs, m) + evsc(nsubs) - e0 - nsubs*epsilon)/(kb*t(tt)))
        end do
      end do

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !      Calculate polynomial coefficients Cn
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do nsubs = nsubsmin, nsubsmax
        c(nsubs) = (real(nsubs)/real(npos) - x)*qn(nsubs)
      end do

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !      Calculate q with the bisection method
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write (*, *) "   - Finding chemical potential (bisection method)..."
      q = x/(1 - x)
      mu = epsilon + kb*t(tt)*log(q)
      ! write (*, *) "Initial chemical potential:  mu =", mu, " eV"

      valpolynom = 0.0
      do nsubs = nsubsmin, nsubsmax
        valpolynom = valpolynom + c(nsubs)*(q**nsubs)
      end do

      valpolynomnew = valpolynom
      do while (valpolynom*valpolynomnew > 0)

        call random_number(r1)
        call random_number(r2)
        do while (r1 < tolq)
          call random_number(r1)
        end do
        do while (r2 < tolq)
          call random_number(r2)
        end do

        qtest = q*r1/r2
        valpolynomnew = 0.0
        do nsubs = nsubsmin, nsubsmax
          valpolynomnew = valpolynomnew + c(nsubs)*(qtest**nsubs)
        end do

      end do

      if (q < qtest) then
        qa = q
        qb = qtest
        va = valpolynom
        vb = valpolynomnew
      else
        qa = qtest
        qb = q
        va = valpolynomnew
        vb = valpolynom
      end if

      do nbis = 1, nbismax

        q = (qa + qb)/2
        oldmu = mu
        mu = epsilon + kb*t(tt)*log(q)
        valpolynomnew = 0.0
        do nsubs = nsubsmin, nsubsmax
          valpolynomnew = valpolynomnew + c(nsubs)*(q**nsubs)
        end do

        if (valpolynomnew*va > 0) then
          qa = q
          va = valpolynomnew
        else
          qb = q
          vb = valpolynomnew
        end if

        if ((abs(qa - qb) < tolq) .and. (abs(mu - oldmu) < tolmu)) then
          ! write (*, *) "Convergence achieved:  delta(mu) = ", mu - oldmu, " eV"
          exit
        end if

      end do

      mu = epsilon + kb*t(tt)*log(q)
      write (*, '(a, f12.6, a)') "   - Converged mu =", mu, " eV"

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !          Get minimum value of the grand potential to calculate enemunrel
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do nsubs = nsubsmin, nsubsmax
        do m = 1, mm(nsubs)
          enemun(nsubs, m) = ene(nsubs, m) + evsc(nsubs) - nsubs*mu
        end do
      end do

      emin = enemun(nsubsmin, 1)
      memin = 1
      nsubsemin = nsubsmin
      do nsubs = nsubsmin, nsubsmax
        do m = 1, mm(nsubs)
          if (emin > enemun(nsubs, m)) then
            emin = enemun(nsubs, m)
            memin = m
            nsubsemin = nsubs
          end if
        end do
      end do

      do nsubs = nsubsmin, nsubsmax
        do m = 1, mm(nsubs)
          enemunrel(nsubs, m) = enemun(nsubs, m) - emin
        end do
      end do

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !          Calculate the partition function from mu
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      z(tt) = 0.0
      do nsubs = nsubsmin, nsubsmax
        do m = 1, mm(nsubs)
          z(tt) = z(tt) + omega(nsubs, m)*exp(-enemunrel(nsubs, m)/(kb*t(tt)))
        end do
      end do

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !          Calculate probabilities and x (equilibrium composition) from mu
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      nx = 0
      do nsubs = nsubsmin, nsubsmax
        do m = 1, mm(nsubs)
          p(nsubs, m, tt) = omega(nsubs, m)*exp(-enemunrel(nsubs, m)/(kb*t(tt)))/z(tt)
          nx = nx + nsubs*p(nsubs, m, tt)
        end do
      end do

      xeq = nx/npos
      write (*, '(a, f12.6)') "   - Equilibrium composition x =", xeq

    else !(xormu = 'mu')
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !          Get minimum value of the grand potential to calculate enemunrel
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do nsubs = nsubsmin, nsubsmax
        do m = 1, mm(nsubs)
          enemun(nsubs, m) = ene(nsubs, m) - nsubs*mu
        end do
      end do

      emin = enemun(nsubsmin, 1)
      memin = 1
      nsubsemin = nsubsmin
      do nsubs = nsubsmin, nsubsmax
        do m = 1, mm(nsubs)
          if (emin > enemun(nsubs, m)) then
            emin = enemun(nsubs, m)
            memin = m
            nsubsemin = nsubs
          end if
        end do
      end do

      do nsubs = nsubsmin, nsubsmax
        do m = 1, mm(nsubs)
          enemunrel(nsubs, m) = enemun(nsubs, m) - emin
        end do
      end do
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !          Calculate the partition function from mu
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      z(tt) = 0.0
      do nsubs = nsubsmin, nsubsmax
        do m = 1, mm(nsubs)
          z(tt) = z(tt) + omega(nsubs, m)*exp(-enemunrel(nsubs, m)/(kb*t(tt)))
        end do
      end do

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !          Calculate probabilities and x (equilibrium composition) from mu
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      nx = 0
      do nsubs = nsubsmin, nsubsmax
        do m = 1, mm(nsubs)
          p(nsubs, m, tt) = omega(nsubs, m)*exp(-enemunrel(nsubs, m)/(kb*t(tt)))/z(tt)
          nx = nx + nsubs*p(nsubs, m, tt)
        end do
      end do

      xeq = nx/npos
      write (*, '(a, f12.6)') "   - Equilibrium composition x =", xeq
      x = xeq

    end if ! End if for xormu.eq."x "

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!          Writing  probabilities.txt file
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    write (20, *)
    write (20, *) "____________________________________________________________________________________________"
    write (20, *) "Temperature: T = ", t(tt), "K"
    write (20, *) "Chemical potential: mu = ", mu, "eV"
    write (20, *) "Composition: x = ", xeq
    write (20, *)

    nx = 0
    write (20, *) " n      m     omega(n,m)       E(n,m)          E-n*mu         prob(n,m,T)   prob(n,m,T)/omega    "
    do nsubs = nsubsmin, nsubsmax
      pn(nsubs, tt) = 0.0
      do m = 1, mm(nsubs)
        write (20, 101) nsubs, m, omega(nsubs, m), ene(nsubs, m), enemun(nsubs, m), p(nsubs, m, tt), p(nsubs, m, tt)/omega(nsubs, m)
101     format(i3, 1x, i6, 1x, i8, 5x, 2(4x, f14.5), 2(4x, e12.6))
        pn(nsubs, tt) = pn(nsubs, tt) + p(nsubs, m, tt)
      end do
    end do

    write (20, *) "------------------------------"
    write (20, *) " n      Cumulative prob(n,T) "
    do nsubs = nsubsmin, nsubsmax
      write (20, 102) nsubs, pn(nsubs, tt)
102   format(i6, 5x, f14.5)
    end do
    write (20, *) "------------------------------"

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculate the energy (E), free energy (F), entropy (S)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    e(tt) = 0.0
    s(tt) = 0.0

    do nsubs = nsubsmin, nsubsmax
      do m = 1, mm(nsubs)
        e(tt) = e(tt) + ene(nsubs, m)*p(nsubs, m, tt)
        s(tt) = s(tt) - kb*p(nsubs, m, tt)*log(p(nsubs, m, tt)/omega(nsubs, m))
      end do
    end do

    f(tt) = e(tt) - t(tt)*s(tt)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculate the average data
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    if (data_exists) then

      avedata(1:ncol, tt) = 0.0

      do col = 1, ncol
        do nsubs = nsubsmin, nsubsmax
          do m = 1, mm(nsubs)
            avedata(col, tt) = avedata(col, tt) + data(nsubs, m, col)*p(nsubs, m, tt)
          end do
        end do
      end do

    end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculate the average spectra
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    if (spectra_exists) then

      avespec(1:npoints, tt) = 0.0

      do point = 1, npoints
        do nsubs = nsubsmin, nsubsmax
          do m = 1, mm(nsubs)
            avespec(point, tt) = avespec(point, tt) + spec(nsubs, m, point)*p(nsubs, m, tt)
          end do
        end do
        if (avespec(point, tt)/maxspec < tolminspec) then
          avespec(point, tt) = 0.0
        end if
      end do

    end if

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
  !      Calculate the limit of full disorder
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      Calculate residual momenta
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  a0 = momentaa0(nsubsmax, npos, x)
  a1 = momentaa1(nsubsmax, npos, x)
  a2 = momentaa2(nsubsmax, npos, x)
  b0 = momentab0(nsubsmin, npos, x)
  b1 = momentab1(nsubsmin, npos, x)
  b2 = momentab2(nsubsmin, npos, x)

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      Calculate alpha and beta
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  a11 = 1.0 - b0 - a0
  a12 = x*npos - b1 - a1
  a21 = a12
  a22 = (x*npos*(1.0 + x*(npos - 1))) - b2 - a2
  bb1 = 1.0
  bb2 = x*npos

  alpha = ((bb1*a22) - (a12*bb2))/((a11*a22) - (a12*a21))
  beta = ((a11*bb2) - (a21*bb1))/((a11*a22) - (a12*a21))

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      Calculate tau
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  do nsubs = nsubsmin, nsubsmax
    tau(nsubs) = alpha + beta*nsubs
  end do

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      Calculate renormalised probabilities and equilibrium concentration
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  nx = 0.0
  do nsubs = nsubsmin, nsubsmax
    do m = 1, mm(nsubs)
      pinf(nsubs, m) = tau(nsubs)*omega(nsubs, m)*(x**nsubs)*((1 - x)**(npos - nsubs))
      nx = nx + nsubs*pinf(nsubs, m)
    end do
  end do

  xeq = nx/npos

  write (*, *) " > Processing T = infinity (ideal disorder limit)..."
  write (*, '(a, f12.6)') "   - Composition x =", xeq

  if (xormu == "mu") then
    x = xeq
  end if

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !          Write probabilities in the limit of infinite temperature
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  write (20, *)
  write (20, *) "____________________________________________________________________________________________"
  write (20, *) "Ideal disorder limit"
  write (20, *) "Composition: x = ", xeq
  write (20, *) " n      m     omega(n,m)     E(n,m)                           prob(n,m,T)   prob(n,m,T)/omega    "
  do nsubs = nsubsmin, nsubsmax
    do m = 1, mm(nsubs)
      write (20, 103) nsubs, m, omega(nsubs, m), ene(nsubs, m), pinf(nsubs, m), pinf(nsubs, m)/omega(nsubs, m)
103   format(i3, 1x, i6, 1x, i8, 5x, 4x, f14.5, 16x, 2(4x, e12.6))
    end do
  end do

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !          Calculate and write results in the limit of infinite temperature
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  einf = 0.0
  do nsubs = nsubsmin, nsubsmax
    do m = 1, mm(nsubs)
      einf = einf + ene(nsubs, m)*pinf(nsubs, m)
    end do
  end do

  if (x == 0.0d0 .or. x == 1.0d0) then
    sinf = 0.0d0
  else
    sinf = -kb*(x*log(x) + (1 - x)*log(1 - x))*npos
  end if

  ! Write thermodynamics.dat
  write (21, 300) "Infinite", einf, " - ", sinf
300 format(a10, 2x, f14.4, 8x, a3, 7x, e12.6)

  ! Calculate and write avedatainf
  if (data_exists) then
    avedatainf(1:ncol) = 0.0
    do col = 1, ncol
      do nsubs = nsubsmin, nsubsmax
        do m = 1, mm(nsubs)
          avedatainf(col) = avedatainf(col) + data(nsubs, m, col)*pinf(nsubs, m)
        end do
      end do
    end do
    write (22, 301) adjustr("Infinite"), avedatainf(1:ncol)
301 format(a10, 2x, 10(f10.4, 2x))
  end if

  ! Calculate avespecinf
  if (spectra_exists) then
    allocate (avespecinf(1:npoints))
    avespecinf(1:npoints) = 0.0
    do point = 1, npoints
      do nsubs = nsubsmin, nsubsmax
        do m = 1, mm(nsubs)
          avespecinf(point) = avespecinf(point) + spec(nsubs, m, point)*pinf(nsubs, m)
        end do
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
    do point = 1, npoints
      write (23, 302) xspec(point), avespec(point, 1:ntt), avespecinf(point)
302   format(1(f10.3, 2x), 8(e12.6, 2x))
    end do
  end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  close (20)
  close (21)
  if (data_exists) close (22)
  if (spectra_exists) close (23)

  write (*, *) ""
  write (*, *) " > Grand-canonical analysis completed."
  write (*, *) ""

end

