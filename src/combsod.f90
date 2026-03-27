!*******************************************************************************
!    Copyright (c) 2022 Ricardo Grau-Crespo, Said Hamad
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

program combsod
  implicit none

  integer, parameter :: nspmax = 10, natmax = 10000, nopmax = 10000, ncellmax = 1000
  real, parameter :: tol0 = 0.001, kb = 8.61734e-5

  integer :: i, j, t, ina, inb, inc, elei, ierr
  integer :: op1, nop1, op, nop, op1new, nop1new, opc
  integer :: sp, nsp, cumnatsp, sptarget
  integer :: at0, nat0, at1, nat1, nat1r, at, nat, at1r, at1i, attmp, att, filer, genxtl, genarc
  integer :: na, nb, nc, nsubs, atini, atfin
  integer :: pos, npos
  integer(kind=8):: ntc, count, nic, equivcount, indcount, combinations, r
  logical :: found, foundnoind, mores
  integer, dimension(:), allocatable:: newconf
  logical, dimension(:), allocatable:: visited
  integer, dimension(2)   :: newshell
  integer, dimension(nopmax)   :: op1good
  real, dimension(nopmax, 3, 3) :: mgroup1, mgroup, mgroup1new
  real, dimension(nopmax, 3)   :: vgroup1, vgroup, vgroup1new
  integer, dimension(:, :), allocatable :: fulleqmatrix, eqmatrixtarget
  integer, dimension(natmax) :: spat0, spat1, spat, spat1r
  real, dimension(natmax, 3) :: coords0, coords1, coords, coords1r
  integer, dimension(natmax)   :: as
  integer, dimension(:), allocatable:: degen
  integer, dimension(nspmax) :: natsp0, natsp1, natsp, ishell
  real, dimension(3) :: coordstemp
  real, dimension(ncellmax, 3)   :: vt
  integer, dimension(:, :), allocatable:: conf, indconf
  real :: a1, b1, c1, alpha, beta, gamma, a, b, c, cc
  real(kind=8) :: prod, x, maxentropy, ientropy, perc
  character, dimension(nspmax) :: symbol*3
  character, dimension(2) :: newsymbol*3
  character :: runtitle*40

! Input files

  open (unit=9, file="INSOD")
  open (unit=12, file="SGO")

! Output files

  open (unit=25, file="coordinates.xyz")
  open (unit=26, file="EQMATRIX")
  open (unit=30, file="OUTSOD")
  open (unit=31, file="SUPERCELL")
  open (unit=43, file="filer")
  open (unit=46, file="OPERATORS")
  open (unit=47, file="cSGO")

!
! DEFINITION OF VARIABLES:
!
! op                  Index for the operators in the supercell
! op1                 Index for the operators in the unit cell
! nop                 Total number of operators in the supercell
! nop1                Total number of operators in the unit cell
! sp                  Index for the species
! nsp                 Total number of species
! sptarget            Number of the species to be substituted
! at0,at1,at          Indexes for the atoms in the asymmetric unit, unit cell and supercell
! nat0,nat1,nat       Total numbers of atoms in the asymmetric unit, unit cell and supercell
! at1r,nat1r            Idem for the atoms in the redundant cell (with repeated positions)
! atini,atfin           Initial and final atom indexes of the species to be substituted
! pos                 Index for atomic positions of the target species
! npos                      Number of atoms of the target species
! conf                      List of all configurations (each configuration is a list of the substituted positions)
! count                     Index for the configurations (conf)
! ntc                 Total number of configurations in conf (count=1,ntc)
! indconf             List of independent configurations
! indcount            Index for the independent configurations
! nic                 Total number of independent configurations in indconf (indcount=1,nic)
! equivconf           Temporary list containing the equivalent configurations at every step of the algorithm
! equivcount          Index for the equivalent configurations (equivconf)
! tol0                      General tolerance
! tol1                      Tolerance used for correcting the x-FLOOR(x) function
!
!
!

  write (*, *) "**************************************************************************** "
  write (*, *) "         SOD (Site Occupancy Disorder) version 0.61  "
  write (*, *) " "
  write (*, *) "         Authors: R. Grau-Crespo and S. Hamad                                   "
  write (*, *) " "
  write (*, *) "         Contact:  <r.grau-crespo@qmul.ac.uk> "
  write (*, *) "**************************************************************************** "
  write (*, *) " "
  write (*, *) " "
  write (*, *) " "
  write (*, *) "Reading input files..."
  write (*, *) " "

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the input file with the uc space group information: SGO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  read (12, *)
  read (12, *) op1
  do while (op1 > 0)
    do i = 1, 3
      read (12, *) (mgroup1(op1, i, j), j=1, 3), vgroup1(op1, i)
    end do
    nop1 = op1
    read (12, *) op1
  end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the INSOD file
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  read (9, *)

  read (9, *) runtitle
  read (9, *)
  read (9, *)
  read (9, *) a1, b1, c1, alpha, beta, gamma
  read (9, *)
  read (9, *)
  read (9, *) nsp
  read (9, *)
  read (9, *)
  read (9, *) (symbol(sp), sp=1, nsp)
  read (9, *)
  read (9, *)
  read (9, *) (natsp0(sp), sp=1, nsp)
  read (9, *)
  read (9, *)

  nat0 = 0
  do sp = 1, nsp
    nat0 = nat0 + natsp0(sp)
  end do

  do at0 = 1, nat0
    read (9, *) (coords0(at0, i), i=1, 3)
  end do
  read (9, *)
  read (9, *)
  read (9, *) na, nb, nc
  read (9, *)
  read (9, *)
  read (9, *) sptarget
  read (9, *)
  read (9, *)
  read (9, *) nsubs
  if (nsubs < 0) then
    write (*, *) "Error: number of substitutions must be >= 0"
    stop
  end if
  read (9, *)
  read (9, *)
  read (9, *)
  read (9, *) (newsymbol(i), i=1, 2)
  read (9, *)
  read (9, *)
  read (9, *) filer

  if (filer < 10) then
    read (9, *)        ! blank line after filer
    read (9, *)
    read (9, *)
    read (9, *) (ishell(sp), sp=1, nsp)
    read (9, *)
    read (9, *) (newshell(i), i=1, 2)
    read (9, *)
    read (9, *) genxtl, genarc
  end if

!cccccccccccccccccccccccccccccccccccc
! Generating spat0 array
!cccccccccccccccccccccccccccccccccccc

  do at0 = 1, natsp0(1)
    spat0(at0) = 1
    cumnatsp = natsp0(1)
  end do
  do sp = 2, nsp
    do at0 = cumnatsp + 1, cumnatsp + natsp0(sp)
      spat0(at0) = sp
    end do
    cumnatsp = cumnatsp + natsp0(sp)
  end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      First create the redundant unit cell from the asymmetric unit
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  write (*, *) ""
  write (*, *) "Generating the supercell..."
  write (*, *) ""

  at1r = 0
  do at0 = 1, nat0
    do op1 = 1, nop1
      at1r = at1r + 1
      coords1r(at1r, 1:3) = matmul(mgroup1(op1, 1:3, 1:3), coords0(at0, 1:3)) + vgroup1(op1, 1:3)
      coords1r(at1r, 1) = cc(coords1r(at1r, 1))
      coords1r(at1r, 2) = cc(coords1r(at1r, 2))
      coords1r(at1r, 3) = cc(coords1r(at1r, 3))
      spat1r(at1r) = spat0(at0)
    end do
  end do
  nat1r = at1r

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Get rid of the redundant atoms, to create the unit cell
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  coords1(1, :) = coords1r(1, :)
  at1r = 1
  at1 = 1
  coords1(1, :) = coords1r(1, :)
  spat1(1) = spat1r(1)
  do at1r = 2, nat1r
    found = .false.
    do at1i = 1, at1
      prod = dot_product(coords1r(at1r, :) - coords1(at1i, :), &
                         coords1r(at1r, :) - coords1(at1i, :))
      if (prod <= tol0) found = .true.
    end do
    if (.not. found) then
      at1 = at1 + 1
      coords1(at1, :) = coords1r(at1r, :)
      spat1(at1) = spat1r(at1r)
    end if
  end do
  nat1 = at1

  do sp = 1, nsp
    natsp1(sp) = 0
  end do
  do at1 = 1, nat1
    natsp1(spat1(at1)) = natsp1(spat1(at1)) + 1
  end do

  natsp(:) = na*nb*nc*natsp1(:)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Create the traslation vectors of the supercell
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  t = 0
  do ina = 0, na - 1
    do inb = 0, nb - 1
      do inc = 0, nc - 1
        t = t + 1
        vt(t, 1) = real(ina)/real(na)
        vt(t, 2) = real(inb)/real(nb)
        vt(t, 3) = real(inc)/real(nc)
      end do
    end do
  end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Create the traslation vectors of the supercell
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  do op1 = 1, nop1
    op1good(op1) = 1
  end do

  if (.not. (na == nb .and. na == nc)) then
    if (na == nb) then
      do op1 = 1, nop1
        do i = 1, 3
          if ((mgroup1(op1, 1, 3) /= 0) .or. (mgroup1(op1, 3, 1) /= 0) .or. &
              (mgroup1(op1, 2, 3) /= 0) .or. (mgroup1(op1, 3, 2) /= 0)) then
            op1good(op1) = 0
          end if
        end do
      end do
    else
      if (na == nc) then
        do op1 = 1, nop1
          do i = 1, 3
            if ((mgroup1(op1, 1, 2) /= 0) .or. (mgroup1(op1, 2, 1) /= 0) .or. &
                (mgroup1(op1, 3, 2) /= 0) .or. (mgroup1(op1, 2, 3) /= 0)) then
              op1good(op1) = 0
            end if
          end do
        end do
      end if
      if (nb == nc) then
        do op1 = 1, nop1
          do i = 1, 3
            if ((mgroup1(op1, 2, 1) /= 0) .or. (mgroup1(op1, 1, 2) /= 0) .or. &
                (mgroup1(op1, 3, 1) /= 0) .or. (mgroup1(op1, 1, 3) /= 0)) then
              op1good(op1) = 0
            end if
          end do
        end do
      end if
      if ((nb /= nc) .and. (na /= nc)) then
        do op1 = 1, nop1
          do i = 1, 3
            if ((mgroup1(op1, 2, 1) /= 0) .or. (mgroup1(op1, 1, 2) /= 0) .or. &
                (mgroup1(op1, 3, 1) /= 0) .or. (mgroup1(op1, 1, 3) /= 0) .or. &
                (mgroup1(op1, 3, 2) /= 0) .or. (mgroup1(op1, 2, 3) /= 0)) then
              op1good(op1) = 0
            end if
          end do
        end do
      end if
    end if
  end if

  op1new = 0
  do op1 = 1, nop1
    if (op1good(op1) == 1) then
      op1new = op1new + 1
      mgroup1new(op1new, :, :) = mgroup1(op1, :, :)
      vgroup1new(op1new, :) = vgroup1(op1, :)
    end if
  end do
  nop1new = op1new

  nop1 = nop1new
  do op1 = 1, nop1
    mgroup1(op1, :, :) = mgroup1new(op1, :, :)
    vgroup1(op1, :) = vgroup1new(op1, :)
  end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Generate the supercell coordinates, by applying these vectors
!      to the unit cell. Also calculate the supercell parameters.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  a = na*a1
  b = nb*b1
  c = nc*c1

  at = 0
  do at1 = 1, nat1
    do t = 1, na*nb*nc
      at = at + 1
      coords(at, 1) = vt(t, 1) + coords1(at1, 1)/na
      coords(at, 2) = vt(t, 2) + coords1(at1, 2)/nb
      coords(at, 3) = vt(t, 3) + coords1(at1, 3)/nc
      spat(at) = spat1(at1)
    end do
  end do
  nat = at

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Write a file with the fractional coordinates of the supercell
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  write (31, 111) a, b, c, alpha, beta, gamma
111 format(6(f10.4, 2x))
  write (31, *) natsp(1:nsp)
  do at = 1, nat
    write (31, 110) symbol(spat(at)), coords(at, 1), coords(at, 2), coords(at, 3)
110 format(a3, 2x, 3(f11.7, 2x))
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculate the supercell operators
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  op = 0
  do op1 = 1, nop1
    do t = 1, na*nb*nc
      op = op + 1
      mgroup(op, :, :) = mgroup1(op1, :, :)
      vgroup(op, 1) = vgroup1(op1, 1)/na + vt(t, 1)
      vgroup(op, 2) = vgroup1(op1, 2)/nb + vt(t, 2)
      vgroup(op, 3) = vgroup1(op1, 3)/nc + vt(t, 3)
    end do
  end do

  nop = op

  allocate (fulleqmatrix(1:nop, 1:nat), stat=ierr)
  if (ierr /= 0) then
    write (*, *) "Error: problem too large for SOD - insufficient memory (fulleqmatrix)"
    stop
  end if

  write (46, *) nop
  do op = 1, nop
    write (46, *) "Operator number ", op
    do i = 1, 3
      write (46, *) mgroup(op, i, 1:3), vgroup(op, i)
    end do
  end do

  write (*, *) "       Composition of the original supercell (parent structure):"
  do sp = 1, nsp
    write (*, *) "                                                         ", symbol(sp), natsp(sp)
  end do
  write (*, *) " "
  write (*, *) "       Number of symmetry operators in the supercell:       ", nop
  write (*, *) " "
  write (*, *) "       Composition of the substituted supercell:"
  do sp = 1, nsp
    if (sp == sptarget) then
      write (*, *) "                                                         ", newsymbol(1), nsubs
      write (*, *) "                                                         ", newsymbol(2), natsp(sp) - nsubs
    else
      write (*, *) "                                                         ", symbol(sp), natsp(sp)
    end if
  end do

  write (*, *) ""
  write (*, *) ""
  write (*, *) "Generating the complete configurational space..."
  write (*, *) " "

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Create Full Equivalence Matrix
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  coordstemp(:) = matmul(mgroup(1, :, :), coords(1, :)) + vgroup(1, :)

  do op = 1, nop
    do at = 1, nat
      coordstemp(:) = matmul(mgroup(op, :, :), coords(at, :)) + vgroup(op, :)
      coordstemp(1) = cc(coordstemp(1))
      coordstemp(2) = cc(coordstemp(2))
      coordstemp(3) = cc(coordstemp(3))
      found = .false.
      attmp = 0
      do while ((.not. found) .and. (attmp < nat))
        attmp = attmp + 1
        if ((abs(coordstemp(1) - coords(attmp, 1)) < tol0) .and. &
            (abs(coordstemp(2) - coords(attmp, 2)) < tol0) .and. &
            (abs(coordstemp(3) - coords(attmp, 3)) < tol0)) &
          found = .true.
      end do
      if ((attmp == nat) .and. (.not. found)) then
        write (*, *) "Error!!! Operator", &
          op, "applied on atom", at, "does not produce another atom in the list!!!"
        attmp = 0
      end if

      fulleqmatrix(op, at) = attmp

    end do
  end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       Calculate the initial and final nat of the target species
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if (sptarget == 1) then
    atini = 1

  else
    atini = 1
    do sp = 1, sptarget - 1
      atini = atini + natsp(sp)
    end do
  end if

  atfin = atini + natsp(sptarget) - 1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Obtain the Equivalence Matrix for the target species from the Full Equivalence Matrix
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  npos = natsp(sptarget)
  allocate (eqmatrixtarget(1:nop, 1:npos), stat=ierr)
  if (ierr /= 0) then
    write (*, *) "Error: problem too large for SOD - insufficient memory (eqmatrixtarget)"
    stop
  end if
  do op = 1, nop
    att = 0
    do at = atini, atfin
      att = att + 1
      eqmatrixtarget(op, att) = fulleqmatrix(op, at) - atini + 1
    end do
  end do

  deallocate (fulleqmatrix)

  write (26, *) nop, npos
  do op = 1, nop
    write (26, *) (eqmatrixtarget(op, pos), pos=1, npos)
  end do

  if (nsubs > npos) then
    write (*, *) "Error: number of substitutions (", nsubs, ") exceeds number of available sites (", npos, ")"
    stop
  end if

! Handle trivial cases: 0 substitutions or all sites substituted (1 configuration each)
  if (nsubs == 0 .or. nsubs == npos) then
    write (*, *) " "
    write (*, *) "       Number of inequivalent configurations:                1"
    write (*, *) " "
    write (30, *) nsubs, " substitutions in ", npos, "sites"
    write (30, *) 1, " configurations"
    if (nsubs == 0) then
      write (30, '(i6, 1x, i6)') 1, 1
    else
      write (30, '(i6, 1x, i6, 30(1x, i4))') 1, 1, (pos, pos=1, npos)
    end if
    write (43, *) filer
    deallocate (eqmatrixtarget)
    write (*, *) "Done!!!"
    write (*, *) ""
    stop
  end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Generate the list of configurations
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ntc = combinations(nsubs, npos)
  if (ntc <= 0) then
    write (*, *) "Error: the total number of configurations exceeds the integer range."
    write (*, *) "       The problem is too large for SOD. Reduce supercell or substitutions."
    stop
  end if
  maxentropy = kb*log(real(ntc))
  write (*, *) " "
  write (*, *) "       Total number of configurations in the supercell:     ", ntc
  write (*, *) " "
  write (*, *) "       Maximum entropy for this composition and supercell:", maxentropy, " eV/K"
  write (*, *) " "

  x = real(nsubs)/real(npos)
  write (*, *) "       Fraction of substituted sites:                 x = ", x
  write (*, *) " "
  if (x == 0.0 .or. x == 1.0) then
    ientropy = 0.0
  else
    ientropy = -npos*kb*(x*log(x) + (1 - x)*log(1 - x))
  end if
  write (*, *) "       Ideal entropy (per cell) for this composition:     ", ientropy, " eV/K"
  write (*, *) " "
!        ntcmax=exp(ientropy/kB)

!!!!!!!!Allocating array sizes

  allocate (degen(1:ntc), stat=ierr)
  if (ierr /= 0) then
    write (*, *) "Error: problem too large for SOD - insufficient memory"
    write (*, *) "       Total configurations: ", ntc, " - reduce supercell size or number of substitutions"
    stop
  end if
  allocate (newconf(1:nsubs), stat=ierr)
  allocate (conf(1:ntc, 1:nsubs), stat=ierr)
  if (ierr /= 0) then
    write (*, *) "Error: problem too large for SOD - insufficient memory"
    write (*, *) "       Total configurations: ", ntc, " - reduce supercell size or number of substitutions"
    stop
  end if
  allocate (indconf(1:ntc, 1:nsubs), stat=ierr)
  if (ierr /= 0) then
    write (*, *) "Error: problem too large for SOD - insufficient memory"
    write (*, *) "       Total configurations: ", ntc, " - reduce supercell size or number of substitutions"
    stop
  end if
  allocate (visited(1:ntc), stat=ierr)
  if (ierr /= 0) then
    write (*, *) "Error: problem too large for SOD - insufficient memory"
    write (*, *) "       Total configurations: ", ntc, " - reduce supercell size or number of substitutions"
    stop
  end if
  visited(:) = .false.

  mores = .false.
  as(:) = 1
  count = 1
  call ksubset(npos, nsubs, as, mores)
  conf(1, 1:nsubs) = as(1:nsubs)
  do while (mores)
    call ksubset(npos, nsubs, as, mores)
    count = count + 1
    conf(count, 1:nsubs) = as(1:nsubs)
  end do

  if (count /= ntc) then
    write (*, *) "Error in KSUBSET subroutine"
    stop
  end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Generate the list of independent configurations
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  write (*, *) " "
  write (*, *) "Finding the inequivalent configurations..."
  write (*, *) " "
  write (*, *) "       Found    Completion "
  write (*, *) "       =====    ========== "

! With this loop we get a function that returns the configuration equivalent, by the operator op,
! to the configuration count.

! For the first configuration (count=1)
  indconf(:, :) = 0
  newconf(:) = 0
  degen(:) = 1

  count = 1
  equivcount = 1
  indcount = 1
  write (47, *) "List of operators for configuration: ", indcount

  indconf(indcount, 1:nsubs) = conf(count, 1:nsubs)
  opc = 1
  write (47, *) opc
  do i = 1, 3
    write (47, *) (mgroup(1, i, j), j=1, 3), vgroup(1, i)
  end do
  visited(count) = .true.
  do op = 2, nop
    do i = 1, nsubs
      elei = conf(count, i)                 ! elei is the element i of the configuration number count
      newconf(i) = eqmatrixtarget(op, elei) ! newconf is the configuration obtained by applying the
      ! operator op to the configuration number conf(count,i)
    end do
    call bubble(newconf, nsubs)
    r = combinations(nsubs, npos)
    do i = 1, nsubs
      if (npos - newconf(i) >= nsubs - i + 1) r = r - combinations(nsubs - i + 1, npos - newconf(i))
    end do
    found = .false.
    if (all(newconf(1:nsubs) == conf(count, 1:nsubs))) then
      found = .true.
      opc = opc + 1
      write (47, *) opc
      do i = 1, 3
        write (47, *) (mgroup(op, i, j), j=1, 3), vgroup(op, i)
      end do
    else
      if (visited(r)) found = .true.
    end if
    if (.not. found) then
      equivcount = equivcount + 1
      visited(r) = .true.
      degen(indcount) = degen(indcount) + 1
    end if
  end do
  write (47, *) 0

! For the rest of configurations (count>1)
  perc = 0.0
  do while (equivcount < ntc)
    count = count + 1
    foundnoind = visited(count)
    if (.not. foundnoind) then
      indcount = indcount + 1
      write (47, *) "List of operators for configuration: ", indcount
      indconf(indcount, 1:nsubs) = conf(count, 1:nsubs)
      equivcount = equivcount + 1
      visited(count) = .true.
      opc = 1
      write (47, *) opc
      do i = 1, 3
        write (47, *) (mgroup(1, i, j), j=1, 3), vgroup(1, i)
      end do
      op_loop: do op = 2, nop
        do i = 1, nsubs
          elei = conf(count, i)
          newconf(i) = eqmatrixtarget(op, elei)
        end do
        call bubble(newconf, nsubs)
        r = combinations(nsubs, npos)
        do i = 1, nsubs
          if (npos - newconf(i) >= nsubs - i + 1) r = r - combinations(nsubs - i + 1, npos - newconf(i))
        end do
        found = .false.
        if (all(newconf(1:nsubs) == conf(count, 1:nsubs))) then
          found = .true.
          opc = opc + 1
          write (47, *) opc
          do i = 1, 3
            write (47, *) (mgroup(op, i, j), j=1, 3), vgroup(op, i)
          end do
        else
          if (visited(r)) cycle op_loop
        end if
        if (.not. found) then
          equivcount = equivcount + 1
          visited(r) = .true.
          degen(indcount) = degen(indcount) + 1
        end if
      end do op_loop
      write (47, *) 0
      if ((100.0d0*equivcount/ntc - perc > 5.0d0) .or. (equivcount == ntc)) then
        perc = 100.0d0*equivcount/ntc
        write (*, '(4x,i6,7x,f5.1,a2)') indcount, perc, "% "
      end if
    end if
  end do

  nic = indcount
  deallocate (eqmatrixtarget)
  deallocate (conf)
  deallocate (visited)

! Trim indconf and degen to the actual number of independent configurations
  block
    integer, dimension(:, :), allocatable :: tmp2d
    integer, dimension(:), allocatable :: tmp1d
    allocate (tmp2d(1:nic, 1:nsubs))
    tmp2d = indconf(1:nic, 1:nsubs)
    call move_alloc(tmp2d, indconf)
    allocate (tmp1d(1:nic))
    tmp1d = degen(1:nic)
    call move_alloc(tmp1d, degen)
  end block

  write (*, *) " "
  write (*, *) "       Number of inequivalent configurations:               ", nic
  write (*, *) " "

  write (30, *) nsubs, " substitutions in ", npos, "sites"
  write (30, *) nic, " configurations"
  do indcount = 1, nic
    write (30, 330) indcount, degen(indcount), indconf(indcount, 1:nsubs)
330 format(i6, 1x, i6, 30(1x, i4))
  end do

!!!!!!!!End of search for independent configurations


!!!!!!! Write FILER to file, to be read by the shell script
  write (43, *) filer
! Calculation input file generation is handled by genersod,
! called automatically by sod_comb.sh when FILER is not 0.

!!!!!!!Deallocating arrays
  deallocate (newconf)
  deallocate (degen)
  deallocate (indconf)

!!!!!!!Reporting the end
  write (*, *) "Done!!!"
  write (*, *) ""
  write (*, *) ""

end program combsod

