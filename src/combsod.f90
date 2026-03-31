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
  integer :: at0, nat0, at1, nat1, nat1r, at, nat, at1r, at1i, attmp, att, atsp, filer
  integer :: na, nb, nc, nsubs, nsubs_min, nsubs_max, nsubs_loop, atini, atfin
  integer :: pos, npos
  integer(kind=8):: ntc, count, nic, equivcount, indcount, combinations, r
  logical :: found, foundnoind, mores
  integer, dimension(:), allocatable:: newconf
  logical, dimension(:), allocatable:: visited
  character(len=256) :: line_buffer
  character(len=10) :: nxx_dir
  integer(kind=8) :: t_start_total, t_start_level, t_now, clock_rate, clock_max
  real(kind=8) :: elapsed_level, elapsed_total
  real(kind=8), dimension(:), allocatable :: elapsed_per_level
  integer, dimension(:), allocatable :: nsubs_level_list
  integer(kind=8), allocatable :: binom(:,:)

! Stage 2: recursive enumeration variables
  integer :: nsubs_prev, ncand, ncand_max, ic, jpar, jj, k, j_remove, nic_prev, idummy, npos_check
  integer(kind=8) :: cand
  logical :: going_upward, use_recursive
  character(len=20) :: prev_outsod, prev_dir
  integer, dimension(:, :), allocatable :: indconf_prev
  integer, dimension(:), allocatable :: degen_prev


  integer, dimension(nopmax)   :: op1good
  real, dimension(nopmax, 3, 3) :: mgroup1, mgroup, mgroup1new
  real, dimension(nopmax, 3)   :: vgroup1, vgroup, vgroup1new
  integer, dimension(:, :), allocatable :: fulleqmatrix, eqmatrixtarget
  integer, dimension(natmax) :: spat0, spat1, spat, spat1r
  real, dimension(natmax, 3) :: coords0, coords1, coords, coords1r
  integer, dimension(natmax)   :: as
  integer, dimension(:), allocatable:: degen
  integer, dimension(nspmax) :: natsp0, natsp1, natsp
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

! Output files written once (not per substitution level)

  open (unit=26, file="EQMATRIX")
  open (unit=31, file="supercell.cif")
  open (unit=43, file="filer")
  open (unit=46, file="OPERATORS")

! Note: OUTSOD (unit 30) and cSGO (unit 47) are opened per level inside the loop below

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
! nsubs               Current number of substitutions (loop variable)
! nsubs_min           Minimum number of substitutions requested
! nsubs_max           Maximum number of substitutions requested
! conf                      List of all configurations (direct) or candidate list (recursive)
! count                     Index for the configurations (conf), used in direct mode
! cand                      Index for candidates, used in recursive mode
! ntc                 Total number of configurations C(npos,nsubs), used in direct mode
! ncand               Number of candidates generated from parent level (recursive mode)
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
  write (*, *) "         SOD (Site Occupancy Disorder) version 0.62  "
  write (*, *) " "
  write (*, *) "         Authors: R. Grau-Crespo and S. Hamad                                   "
  write (*, *) " "
  write (*, *) "         Contact:  <r.grau-crespo@qmul.ac.uk> "
  write (*, *) "**************************************************************************** "
  write (*, *) " "
  write (*, *) " "
  write (*, *) " "
  call system_clock(t_start_total, clock_rate, clock_max)

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
! Read nsubs_min and nsubs_max from the same line.
! If only one value is present (old-style INSOD), both are set equal to it.
  read (9, '(A)') line_buffer
  read (line_buffer, *, IOSTAT=ierr) nsubs_min, nsubs_max
  if (ierr /= 0) then
    read (line_buffer, *) nsubs_min
    nsubs_max = nsubs_min
  end if
  if (nsubs_min < 0) then
    write (*, *) "Error: number of substitutions must be >= 0"
    stop
  end if
  if (nsubs_max < nsubs_min) then
    write (*, *) "Error: nsubs_max must be >= nsubs_min"
    stop
  end if
  read (9, *)
  read (9, *)
  read (9, *)
  read (9, *) (newsymbol(i), i=1, 2)
  read (9, *)
  read (9, *)
  read (9, *) filer

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
! Write supercell.cif with P1 symmetry; atom order matches OUTSOD indices
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  write (31, '(a)') "data_supercell"
  write (31, '(a)') " "
  write (31, '(a, f12.6)') "_cell_length_a     ", a
  write (31, '(a, f12.6)') "_cell_length_b     ", b
  write (31, '(a, f12.6)') "_cell_length_c     ", c
  write (31, '(a, f12.6)') "_cell_angle_alpha  ", alpha
  write (31, '(a, f12.6)') "_cell_angle_beta   ", beta
  write (31, '(a, f12.6)') "_cell_angle_gamma  ", gamma
  write (31, '(a)') " "
  write (31, '(a)') "_symmetry_space_group_name_H-M  'P 1'"
  write (31, '(a)') "_symmetry_Int_Tables_number      1"
  write (31, '(a)') " "
  write (31, '(a)') "loop_"
  write (31, '(a)') "_symmetry_equiv_pos_as_xyz"
  write (31, '(a)') "'x, y, z'"
  write (31, '(a)') " "
  write (31, '(a)') "loop_"
  write (31, '(a)') "_atom_site_label"
  write (31, '(a)') "_atom_site_type_symbol"
  write (31, '(a)') "_atom_site_fract_x"
  write (31, '(a)') "_atom_site_fract_y"
  write (31, '(a)') "_atom_site_fract_z"
  atsp = 0
  do at = 1, nat
    if (at == 1 .or. spat(at) /= spat(at - 1)) atsp = 0
    atsp = atsp + 1
    write (31, '(a, i0, 2x, a, 3(2x, f11.7))') trim(symbol(spat(at))), atsp, &
          trim(symbol(spat(at))), coords(at, 1), coords(at, 2), coords(at, 3)
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

  if (nsubs_max > npos) then
    write (*, *) "Error: nsubs_max (", nsubs_max, ") exceeds number of available sites (", npos, ")"
    stop
  end if

! Direction selection: choose the end that minimises the starting candidate count
  going_upward = (nsubs_min <= npos - nsubs_max)
  if (going_upward) then
    write (*, *) "Direction: UPWARD (direct start at nsubs =", nsubs_min, ")"
  else
    write (*, *) "Direction: DOWNWARD (direct start at nsubs =", nsubs_max, ")"
  end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Loop over substitution levels
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  allocate (elapsed_per_level(0:nsubs_max - nsubs_min))
  allocate (nsubs_level_list(0:nsubs_max - nsubs_min))
  elapsed_per_level(:) = 0.0d0
  nsubs_level_list(:) = 0

  do nsubs_loop = 0, nsubs_max - nsubs_min
    if (going_upward) then
      nsubs = nsubs_min + nsubs_loop
    else
      nsubs = nsubs_max - nsubs_loop
    end if

! Create output directory nXX/ and open per-level output files

    call system_clock(t_start_level)
    write (nxx_dir, '("n", i2.2)') nsubs
    call execute_command_line("mkdir -p " // trim(nxx_dir), wait=.true.)
    open (unit=30, file=trim(nxx_dir) // "/OUTSOD")
    open (unit=47, file=trim(nxx_dir) // "/cSGO")

    write (*, *) " "
    write (*, *) "=== Substitution level: nsubs =", nsubs, " ==="
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
        write (30, '(i6, 1x, i6, *(1x, i4))') 1, 1, (pos, pos=1, npos)
      end if
      close (30)
      close (47)
      cycle
    end if

! Precompute binomial coefficient table: binom(k,n) = C(n,k), k=0..nsubs, n=0..npos
    allocate (binom(0:nsubs, 0:npos))
    binom(0, :) = 1_8
    do k = 1, nsubs
      binom(k, 0:k - 1) = 0_8
      do j = k, npos
        binom(k, j) = binom(k - 1, j - 1) + binom(k, j - 1)
      end do
    end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Check for adjacent OUTSOD to decide direct vs recursive path
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    if (going_upward) then
      write (prev_dir, '("n", i2.2)') nsubs - 1
    else
      write (prev_dir, '("n", i2.2)') nsubs + 1
    end if
    prev_outsod = trim(prev_dir) // "/OUTSOD"

    open (unit=50, file=trim(prev_outsod), status='old', IOSTAT=ierr)
    use_recursive = (ierr == 0)

    if (use_recursive) then
! Read header line 1: "N substitutions in M sites"
      read (50, '(A)') line_buffer
      read (line_buffer, *) nsubs_prev
      k = index(line_buffer, ' in ')
      read (line_buffer(k + 4:), *) npos_check
      if (npos_check /= npos) then
        write (*, *) "Warning: npos in adjacent OUTSOD (", npos_check, &
                     ") does not match current npos (", npos, "). Falling back to direct."
        use_recursive = .false.
        close (50)
      end if
    end if

    if (use_recursive) then
! Read header line 2: "K configurations"
      read (50, '(A)') line_buffer
      read (line_buffer, *) nic_prev
      allocate (indconf_prev(1:nic_prev, 1:nsubs_prev))
      allocate (degen_prev(1:nic_prev))
      if (nsubs_prev == 0) then
        do ic = 1, nic_prev
          read (50, *) idummy, degen_prev(ic)
        end do
      else
        do ic = 1, nic_prev
          read (50, *) idummy, degen_prev(ic), (indconf_prev(ic, j), j=1, nsubs_prev)
        end do
      end if
      close (50)
      write (*, *) "Using recursive generation from ", trim(prev_outsod), &
                   " (", nic_prev, " parents)"
    else
      write (*, *) "Generating the complete configurational space..."
    end if
    write (*, *) " "

    if (use_recursive) then

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Recursive path: generate candidate list from parent level
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (going_upward) then
! Upward: add one substitution to each parent config
! Each parent (nsubs-1 sites) generates npos-(nsubs-1) children
        ncand_max = nic_prev * (npos - nsubs_prev)
        allocate (conf(1:ncand_max, 1:nsubs), stat=ierr)
        if (ierr /= 0) then
          write (*, *) "Error: insufficient memory for candidate list"
          stop
        end if
        ncand = 0
        do ic = 1, nic_prev
          jpar = 1
          do pos = 1, npos
            if (jpar <= nsubs_prev .and. indconf_prev(ic, jpar) == pos) then
              jpar = jpar + 1  ! pos is in the parent; skip it
            else
! pos is not in the parent; create child by inserting pos into the sorted parent
! Sites indconf_prev(ic,1:jpar-1) are < pos; sites indconf_prev(ic,jpar:nsubs_prev) are > pos
              ncand = ncand + 1
              conf(ncand, 1:jpar - 1) = indconf_prev(ic, 1:jpar - 1)
              conf(ncand, jpar) = pos
              conf(ncand, jpar + 1:nsubs) = indconf_prev(ic, jpar:nsubs_prev)
            end if
          end do
        end do

      else
! Downward: remove one substitution from each parent config
! Each parent (nsubs+1 sites) generates nsubs+1 children
        ncand_max = nic_prev * nsubs_prev
        allocate (conf(1:ncand_max, 1:nsubs), stat=ierr)
        if (ierr /= 0) then
          write (*, *) "Error: insufficient memory for candidate list"
          stop
        end if
        ncand = 0
        do ic = 1, nic_prev
          do j_remove = 1, nsubs_prev
            ncand = ncand + 1
! Child = parent with site at position j_remove removed; result is already sorted
            k = 0
            do jj = 1, nsubs_prev
              if (jj /= j_remove) then
                k = k + 1
                conf(ncand, k) = indconf_prev(ic, jj)
              end if
            end do
          end do
        end do

      end if

      deallocate (indconf_prev, degen_prev)

      write (*, *) "       Number of candidate configurations:            ", ncand
      write (*, *) " "

! Allocate arrays for symmetry reduction of candidate list
! visited is still indexed by combinatorial rank in C(npos,nsubs) space
      ntc = binom(nsubs, npos)
      if (ntc <= 0) then
        write (*, *) "Error: total configuration count overflows integer range."
        write (*, *) "       The problem is too large for SOD. Reduce supercell or substitutions."
        stop
      end if
      allocate (visited(1:ntc), stat=ierr)
      if (ierr /= 0) then
        write (*, *) "Error: insufficient memory for visited array (size", ntc, ")"
        write (*, *) "       Try reducing supercell or substitutions."
        stop
      end if
      visited(:) = .false.
      allocate (newconf(1:nsubs), stat=ierr)
      allocate (indconf(1:ncand, 1:nsubs), stat=ierr)
      if (ierr /= 0) then
        write (*, *) "Error: insufficient memory for indconf"
        stop
      end if
      allocate (degen(1:ncand), stat=ierr)
      if (ierr /= 0) then
        write (*, *) "Error: insufficient memory for degen"
        stop
      end if
      indconf(:, :) = 0
      newconf(:) = 0
      degen(:) = 1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Symmetry reduction of candidate list
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      write (*, *) "Finding the inequivalent configurations (recursive)..."
      write (*, *) " "
      write (*, *) "       Found    Completion "
      write (*, *) "       =====    ========== "

      indcount = 0
      perc = 0.0d0
      do cand = 1, int(ncand, kind=8)
! Compute combinatorial rank of this candidate
        r = binom(nsubs, npos)
        do i = 1, nsubs
          if (npos - conf(cand, i) >= nsubs - i + 1) &
            r = r - binom(nsubs - i + 1, npos - conf(cand, i))
        end do
        if (visited(r)) cycle  ! already accounted for as equivalent of a known config

! New independent configuration found
        indcount = indcount + 1
        indconf(indcount, 1:nsubs) = conf(cand, 1:nsubs)
        visited(r) = .true.
        write (47, *) "List of operators for configuration: ", indcount
        opc = 1
        write (47, *) opc
        do i = 1, 3
          write (47, *) (mgroup(1, i, j), j=1, 3), vgroup(1, i)
        end do

! Apply all symmetry operators to find equivalents and count degeneracy
        do op = 2, nop
          do i = 1, nsubs
            elei = conf(cand, i)
            newconf(i) = eqmatrixtarget(op, elei)
          end do
          call bubble(newconf, nsubs)
          r = binom(nsubs, npos)
          do i = 1, nsubs
            if (npos - newconf(i) >= nsubs - i + 1) &
              r = r - binom(nsubs - i + 1, npos - newconf(i))
          end do
          found = .false.
          if (all(newconf(1:nsubs) == conf(cand, 1:nsubs))) then
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
            visited(r) = .true.
            degen(indcount) = degen(indcount) + 1
          end if
        end do
        write (47, *) 0
        if ((100.0d0*cand/ncand - perc > 5.0d0) .or. (cand == ncand)) then
          perc = 100.0d0*cand/ncand
          write (*, '(4x,i6,7x,f5.1,a2)') indcount, perc, "% "
        end if
      end do
      nic = indcount
      deallocate (conf)
      deallocate (visited)

    else

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Direct path: enumerate C(npos,nsubs) configurations on-the-fly
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ntc = binom(nsubs, npos)
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

!!!!!!!!Allocating array sizes
! conf(ntc, nsubs) is NOT allocated; configurations are generated one at a time via ksubset.
! indconf is allocated at ntc (safe upper bound on the number of independent configs).

      allocate (degen(1:ntc), stat=ierr)
      if (ierr /= 0) then
        write (*, *) "Error: problem too large for SOD - insufficient memory"
        write (*, *) "       Total configurations: ", ntc, " - reduce supercell size or number of substitutions"
        stop
      end if
      allocate (newconf(1:nsubs), stat=ierr)
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

! Initialise ksubset: first call (mores=.false.) sets as to the first configuration.
      mores = .false.
      call ksubset(npos, nsubs, as, mores)
      ! as(1:nsubs) now holds configuration rank 1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Generate the list of independent configurations
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      write (*, *) " "
      write (*, *) "Finding the inequivalent configurations..."
      write (*, *) " "
      write (*, *) "       Found    Completion "
      write (*, *) "       =====    ========== "

      indconf(:, :) = 0
      newconf(:) = 0
      degen(:) = 1

! First configuration (count=1, always new — visited(1) is false)
      count = 1
      equivcount = 1
      indcount = 1
      write (47, *) "List of operators for configuration: ", indcount

      indconf(indcount, 1:nsubs) = as(1:nsubs)
      opc = 1
      write (47, *) opc
      do i = 1, 3
        write (47, *) (mgroup(1, i, j), j=1, 3), vgroup(1, i)
      end do
      visited(count) = .true.
      do op = 2, nop
        do i = 1, nsubs
          elei = as(i)
          newconf(i) = eqmatrixtarget(op, elei)
        end do
        call bubble(newconf, nsubs)
        r = binom(nsubs, npos)
        do i = 1, nsubs
          if (npos - newconf(i) >= nsubs - i + 1) r = r - binom(nsubs - i + 1, npos - newconf(i))
        end do
        found = .false.
        if (all(newconf(1:nsubs) == as(1:nsubs))) then
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

! Advance ksubset to next configuration
      call ksubset(npos, nsubs, as, mores)

! Remaining configurations (count > 1): advance ksubset one step per iteration
      perc = 0.0
      do while (equivcount < ntc)
        count = count + 1
        foundnoind = visited(count)
        if (.not. foundnoind) then
          indcount = indcount + 1
          write (47, *) "List of operators for configuration: ", indcount
          indconf(indcount, 1:nsubs) = as(1:nsubs)
          equivcount = equivcount + 1
          visited(count) = .true.
          opc = 1
          write (47, *) opc
          do i = 1, 3
            write (47, *) (mgroup(1, i, j), j=1, 3), vgroup(1, i)
          end do
          op_loop: do op = 2, nop
            do i = 1, nsubs
              elei = as(i)
              newconf(i) = eqmatrixtarget(op, elei)
            end do
            call bubble(newconf, nsubs)
            r = binom(nsubs, npos)
            do i = 1, nsubs
              if (npos - newconf(i) >= nsubs - i + 1) r = r - binom(nsubs - i + 1, npos - newconf(i))
            end do
            found = .false.
            if (all(newconf(1:nsubs) == as(1:nsubs))) then
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
        call ksubset(npos, nsubs, as, mores)
      end do

      nic = indcount
      deallocate (visited)

    end if  ! use_recursive / direct

!!!!!!!!End of search for independent configurations

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
    end do
330 format(i6, 1x, i6, *(1x, i4))

!!!!!!!Deallocating arrays
    deallocate (newconf)
    deallocate (degen)
    deallocate (indconf)

    deallocate (binom)

    call system_clock(t_now)
    elapsed_level = real(t_now - t_start_level, kind=8) / real(clock_rate, kind=8)
    elapsed_per_level(nsubs_loop) = elapsed_level
    nsubs_level_list(nsubs_loop) = nsubs

    close (30)
    close (47)

  end do  ! nsubs_loop

!!!!!!! Write FILER to file, to be read by the shell script
  write (43, *) filer
! Calculation input file generation is handled by genersod,
! called automatically by sod_comb.sh when FILER is not -1.

  deallocate (eqmatrixtarget)

!!!!!!!Reporting the end
  call system_clock(t_now)
  elapsed_total = real(t_now - t_start_total, kind=8) / real(clock_rate, kind=8)
  write (*, *) " "
  write (*, *) "Timing summary:"
  write (*, *) "       nsubs    Wall time (s)"
  write (*, *) "       -----    -------------"
  do nsubs_loop = 0, nsubs_max - nsubs_min
    write (*, '(7x, i5, 4x, f12.2)') nsubs_level_list(nsubs_loop), elapsed_per_level(nsubs_loop)
  end do
  write (*, *) "       -----    -------------"
  write (*, '(a, f12.2)') "       Total         ", elapsed_total
  write (*, *) " "
  deallocate (elapsed_per_level, nsubs_level_list)
  write (*, *) "Done!!!"
  write (*, *) ""
  write (*, *) ""

end program combsod
