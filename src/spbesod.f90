!*******************************************************************************
!    Copyright (c) 2018 Ricardo Grau-Crespo, Said Hamad
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

program spbesod
  implicit none

  integer :: p, q, m, m1, m2, m1ref, m2ref, mmin, mmax, mm, mm1, mm2, op, nop, npos, nsubs, aux, irescale
  real(kind=8) :: e0, mu1, mu2, e1ref, e2ref, emin, emax
  logical:: inspbe_exists
  real(kind=8), dimension(:), allocatable:: de1
  real(kind=8), dimension(:, :), allocatable:: de2
  real(kind=8), dimension(:), allocatable:: energies1
  real(kind=8), dimension(:), allocatable:: energies2
  real(kind=8), dimension(:), allocatable:: a1, a2, energies
  integer, dimension(:, :), allocatable:: eqmatrix
  integer, dimension(:), allocatable:: omega, omega1, omega2
  integer, dimension(:), allocatable:: conf1
  integer, dimension(:, :), allocatable:: conf2
  integer, dimension(:, :), allocatable:: conf

!!!!!! Input files

  open (unit=10, file="ENERGIES0")
  open (unit=11, file="ENERGIES1")
  open (unit=12, file="ENERGIES2")
  open (unit=13, file="EQMATRIX")
  open (unit=14, file="OUTSOD")
  open (unit=15, file="OUTSOD1")
  open (unit=16, file="OUTSOD2")

  inquire (file="INSPBE", exist=inspbe_exists)
  if (inspbe_exists) then
    open (unit=17, file="INSPBE")
  else
    open (unit=18, file="INSPBE.tmp")
  end if

!!!!!! Output files

  open (unit=50, file="ENERGIES")
  open (unit=99, file="OUTSPBE")

! DEFINITION OF VARIABLES:
!
! nsubs               numbers of substitutions
! op                  Index for the operators in the supercell
! nop                 Total number of operators in the supercell
! npos                Number of atoms of the target species
! pos                 Index for the positions
! conf                List of all independent configurations for the structure with nsubs substitutions
! conf1               List of all independent configurations for the structure with 1 substitutions
! conf2               List of all independent configurations for the structure with 2 substitutions
!                     All configurations are given as a list of the substituted positions
! m                   Index for the configurations (conf)
! Mn                  Total number of configurations in conf (m=1,Mm)
! m1                  Index for the configurations (conf1) for the structure with 1 defect
! Mm1                 Total number of configurations in conf1 (m1=1,Mm1) for the structure with 1 defect
! m2                  Index for the configurations (conf2) for the structure with 2 defects
! Mm2                 Total number of configurations in conf2 (m2=1,Mm2) for the structure with 2 defects
! mu1                 Scaling parameter for single defect energies
! mu2                 Scaling parameter for pair-wise interactions
! EN                  Energy of the x=1 endmember

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the EQMATRIX file
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  read (13, *) nop, npos
  allocate (eqmatrix(1:nop, 1:npos))

  do op = 1, nop
    read (13, *) eqmatrix(op, 1:npos)
  end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the OUTSOD files:
!              - OUTSOD0 is not read, since we know the energy is E0
!              - OUTSOD1 for the structure with 1   substitution
!              - OUTSOD2 for the structure with 2   substitutions
!              - OUTSOD  for the structure with n>2 substitutions
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!       OUTSOD
!
  read (14, *) nsubs
  read (14, *) mm

  allocate (conf(1:mm, 1:nsubs))
  allocate (omega(1:mm))

  do m = 1, mm
    read (14, *) aux, omega(m), conf(m, 1:nsubs)
  end do

!       OUTSOD1
!
  read (15, *)
  read (15, *) mm1

  allocate (conf1(1:mm1))
  allocate (omega1(1:mm1))

  do m = 1, mm1
    read (15, *) aux, omega1(m), conf1(m)
  end do

!       OUTSOD2
!
  read (16, *)
  read (16, *) mm2

  allocate (conf2(1:mm2, 1:2))
  allocate (omega2(1:mm2))

  do m = 1, mm2
    read (16, *) aux, omega2(m), conf2(m, 1:2)
  end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the files with the energies of 0, 1 and 2 substitutions (ENERGIES0, ENERGIES1 and  ENERGIES2)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  allocate (energies1(1:mm1))
  allocate (energies2(1:mm2))

  read (10, *) e0

  do m = 1, mm1
    read (11, *) energies1(m)
  end do

  do m = 1, mm2
    read (12, *) energies2(m)
  end do

  write (99, *) "-----------------------------------------------------------------"
  write (99, *) "Parameters calculated from data for 0, 1 and 2 substitutions:"
  write (99, *) "-----------------------------------------------------------------"
  write (99, *)
  write (99, *) "E0 (eV)"
  write (99, *) e0

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculating dE1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  allocate (de1(1:npos))

  do m1 = 1, mm1
    do op = 1, nop
      de1(eqmatrix(op, conf1(m1))) = energies1(m1) - e0
    end do
  end do

  write (99, *)
  write (99, *) "dE1 (eV)"
  write (99, *) de1(:)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculating dE2
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  allocate (de2(1:npos, 1:npos))

  do m2 = 1, mm2
    do op = 1, nop
      de2(eqmatrix(op, conf2(m2, 1)), eqmatrix(op, conf2(m2, 2))) = &
        energies2(m2) - e0 - de1(eqmatrix(op, conf2(m2, 1))) - de1(eqmatrix(op, conf2(m2, 2)))
    end do
  end do

  write (99, *)
  write (99, *) "dE2 (eV)"
  do p = 1, npos
    write (99, *) de2(p, 1:npos)
  end do
  write (99, *)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculating configuration-dependent summations
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  allocate (a1(1:mm))
  allocate (a2(1:mm))

  do m = 1, mm
    ! Initializing
    a1(m) = 0.0
    a2(m) = 0.0

    do p = 1, nsubs
      a1(m) = a1(m) + de1(conf(m, p))
    end do

    do p = 1, nsubs - 1
      do q = p + 1, nsubs
        a2(m) = a2(m) + de2(conf(m, p), conf(m, q))
      end do
    end do

  end do

  write (99, *)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the INSPBE file, if it exists, and calculating rescaling parameters mu1 and mu2
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if (inspbe_exists) then

    read (17, *)
    read (17, *) irescale
    read (17, *)
    select case (irescale)
    case (0)
      read (17, *)
      write (*, *) "No rescaling (mu1=mu2=1.0)."
      write (99, *) "No rescaling (mu1=mu2=1.0)."
      mu1 = 1.0
      mu2 = 1.0
    case (1)
      read (17, *) mu1, mu2
      write (*, *) "Rescaling parameters mu1 and mu2 read from INSPBE."
      write (99, *) "Rescaling parameters mu1 and mu2 read from INSPBE."
      write (99, *) mu1, mu2
    case (2)
      read (17, *) m1ref, e1ref
      read (17, *) m2ref, e2ref
      mu1 = (a2(m2ref)*(e1ref - e0) - a2(m1ref)*(e2ref - e0))/(a1(m1ref)*a2(m2ref) - a2(m1ref)*a1(m2ref))
      mu2 = (a1(m1ref)*(e2ref - e0) - a1(m2ref)*(e1ref - e0))/(a1(m1ref)*a2(m2ref) - a2(m1ref)*a1(m2ref))
      write (*, *) "Rescaling parameters calculated using reference energies for two configurations."
      write (99, *) "mu1 = ", mu1, "mu2 = ", mu2
    case default
      mu1 = 0.0
      mu2 = 0.0
      write (*, *) "Invalid case for irescale in INSPBE - no rescaling will be applied (mu1=mu2=1.0)."
    end select

  else

    mu1 = 1.0
    mu2 = 1.0
    write (*, *) "No INSPBE file, therefore no rescaling applied (mu1=mu2=1.0)."
    write (99, *) "No INSPBE file, therefore no rescaling applied (mu1=mu2=1.0)."

  end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculating ENERGIES
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  allocate (energies(1:mm))

  write (99, *)
  write (99, *)
  write (99, *) "Zeroth-, first-, and second-order contribution by configuration for ", nsubs, "substitutions:"
  write (99, *)
  write (99, *) "   m     E0(eV)             dE1(eV)         dE2(eV)      E[total] (eV)"
  write (99, *) " --------------------------------------------------------------------------"

  do m = 1, mm
    energies(m) = e0 + mu1*a1(m) + mu2*a2(m)
    write (99, '(1I5,4F16.6)') m, e0, mu1*a1(m), mu2*a2(m), energies(m)
    write (50, *) energies(m)
  end do

  write (99, *)
  write (99, *)

  mmin = minloc(energies(1:mm), dim=1)
  emin = minval(energies(1:mm))
  mmax = maxloc(energies(1:mm), dim=1)
  emax = maxval(energies(1:mm))

  write (99, *) "Minimum-energy configuration: ", mmin, " with energy: ", emin, " eV."
  write (99, *) "Maximum-energy configuration: ", mmax, " with energy: ", emax, " eV."
  write (99, *)
  write (99, *)

  if (.not. (inspbe_exists)) then
    write (99, *) "To improve agreement with reference set, calculate reference energies for these two configurations,"
    write (99, *) "then input reference energies in INSPBE file. Use INSPBE.tmp as template."
    write (99, *)
    write (18, *) "# irescale case: 0 = no rescaling; 1 = enter mu1 and mu2 manually; 2= enter two reference energies &
                &(e.g. from dft)"
    write (18, *) "2"
    write (18, *) "# If irescale=1, enter one line with mu1, mu2; if irescale=2, enter two lines (m1, E1), (m2, E2)"
    write (18, *) mmin, emin
    write (18, *) mmax, emax
  end if

!!!!!!Deallocating arrays
  deallocate (eqmatrix)
  deallocate (conf)
  deallocate (omega)
  deallocate (conf1)
  deallocate (omega1)
  deallocate (conf2)
  deallocate (omega2)
  deallocate (energies1)
  deallocate (energies2)
  deallocate (de1)
  deallocate (de2)
  deallocate (energies)
  deallocate (a1)
  deallocate (a2)

!!!!!!Reporting the end
  write (*, *) "Done!!!"
  write (*, *)
  write (*, *)

end program spbesod

