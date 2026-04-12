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
  use iso_fortran_env, only: int64, real64
  implicit none

  integer, parameter :: nspmax = 10, natmax = 10000, nopmax = 10000, ncellmax = 1000
  integer, parameter :: ntargetmax = 5, nkmax = 5
  real(real64), parameter :: tol0 = 0.001_real64
  real(real64), parameter :: kb = 8.61734e-5_real64

  integer :: i, j, t, ina, inb, inc, ierr
  integer :: op1, nop1, op, nop, op1new, nop1new, opc
  integer :: sp, nsp, cumnatsp
  integer :: ntarget, sptarget(ntargetmax)
  integer :: nk(ntargetmax)
  integer :: at0, nat0, at1, nat1, nat1r, at, nat, at1r, at1i, attmp, att, atsp, filer
  integer :: na, nb, nc, nsubs, nsubs_min, nsubs_max, nsubs_loop, atini, atfin
  integer :: nsubs_t(ntargetmax, nkmax), npos_t(ntargetmax), atini_t(ntargetmax), atfin_t(ntargetmax)
  integer :: npos_max, t_A, tt, offset, nsubs_tot, nsubs_max_t
  integer :: pos, npos
  integer(int64):: ntc, count, nic, equivcount, indcount, r
  logical :: found, foundnoind, mores
  integer, dimension(:), allocatable:: newconf
  logical, dimension(:), allocatable:: visited
  character(len=256) :: line_buffer
  character(len=10) :: nxx_dir
  integer(int64) :: t_start_total, t_start_level, t_now, clock_rate, clock_max
  integer(int64) :: ntc_A, ntc_B, ntc_joint
  integer(int64) :: ntc_t(ntargetmax), r_t(ntargetmax), new_r, new_r_tt
  logical :: mores_t(ntargetmax)
  integer(int64) :: r_A, r_B, new_r_A, new_r_B
  integer, dimension(:), allocatable :: as_A, as_B, newconf_A, newconf_B, new_as_B
  integer, dimension(:,:), allocatable :: as_t, newconf_t
  integer, dimension(:), allocatable :: tmp_conf
  integer(int64), allocatable :: binom_A(:,:), binom_B(:,:)
  integer(int64), allocatable :: binom_t(:,:,:)
  logical :: mores_A, mores_B
  integer(int64) :: ntc_C
  integer(int64) :: r_C, new_r_C
  integer, dimension(:), allocatable :: as_C, newconf_C, new_as_C
  integer, dimension(:), allocatable :: as_notB_tmp, new_notB_tmp
  integer(int64), allocatable :: binom_C(:,:)
  logical :: mores_C
  integer :: n2rem
! Phase 4: multi-target multi-nary
  integer :: nsubs_tot_t(ntargetmax), nrem2_t(ntargetmax), nsubs_pos_max
  integer(int64) :: ntc_pos_t(ntargetmax), ntc_col1_t(ntargetmax), ntc_col2_t(ntargetmax)
  integer(int64) :: r_pos, r_col1, r_col2, new_r_pos, new_r_col1, new_r_col2
  logical :: mores_pos_t(ntargetmax), mores_col1_t(ntargetmax), mores_col2_t(ntargetmax)
  integer, allocatable :: as_pos_t(:,:), as_col1_t(:,:), as_col2_t(:,:), notcol1_t(:,:)
  integer, allocatable :: newpos_t(:,:), newval_col1_t(:,:), new_col1_t(:,:)
  integer, allocatable :: new_notcol1_t(:,:), newval_col2_t(:,:), new_col2_t(:,:)
  integer(int64), allocatable :: binom_pos_t(:,:,:), binom_col_t(:,:,:)
  real(real64) :: elapsed_level, elapsed_total
  real(real64), dimension(:), allocatable :: elapsed_per_level
  integer, dimension(:), allocatable :: nsubs_level_list
  integer(int64), allocatable :: binom(:,:)

! Stage 2: recursive enumeration variables
  integer :: nsubs_prev, ncand, ncand_max, ic, jpar, jj, k, j_remove, nic_prev, idummy, npos_check
  integer(int64) :: cand
  logical :: going_upward, use_recursive
  character(len=20) :: prev_outsod, prev_dir
  integer, dimension(:, :), allocatable :: indconf_prev
  integer, dimension(:), allocatable :: degen_prev


  integer, dimension(nopmax)   :: op1good
  real(real64), dimension(nopmax, 3, 3) :: mgroup1, mgroup, mgroup1new
  real(real64), dimension(nopmax, 3)   :: vgroup1, vgroup, vgroup1new
  integer, dimension(:, :), allocatable :: fulleqmatrix
  integer, dimension(:, :, :), allocatable :: eqmt
  integer, dimension(natmax) :: spat0, spat1, spat, spat1r
  real(real64), dimension(natmax, 3) :: coords0, coords1, coords, coords1r
  integer, dimension(natmax)   :: as
  integer, dimension(:), allocatable:: degen
  integer, dimension(nspmax) :: natsp0, natsp1, natsp
  real(real64), dimension(3) :: coordstemp
  real(real64), dimension(ncellmax, 3)   :: vt
  integer, dimension(:, :), allocatable:: conf, indconf
  real(real64) :: a1, b1, c1, alpha, beta, gamma, a, b, c, cc
  real(real64) :: prod, x, maxentropy, ientropy, perc
  character(len=3), dimension(nspmax) :: symbol
  character(len=5), dimension(ntargetmax, nkmax+1) :: newsymbol
  character(len=40) :: runtitle

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

  write (*, *) "============================================================================"
  write (*, *) "         SOD (Site Occupancy Disorder) version 0.71"
  write (*, *) ""
  write (*, *) "         Authors: R. Grau-Crespo and S. Hamad"
  write (*, *) "         Contact: <r.grau-crespo@qmul.ac.uk>"
  write (*, *) "============================================================================"
  write (*, *) ""
  call system_clock(t_start_total, clock_rate, clock_max)

  write (*, *) " > Reading input files..."
  write (*, *) ""

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
  if (nat0 > natmax) then
    write (*, '(a,i0,a,i0,a)') " Error: nat0 = ", nat0, " exceeds natmax = ", natmax, ". Increase natmax in combsod.f90."
    stop 1
  end if

  do at0 = 1, nat0
    read (9, *) (coords0(at0, i), i=1, 3)
  end do
  read (9, *)
  read (9, *)
  read (9, *) na, nb, nc
  read (9, *)
  read (9, *)
! Read sptarget line: 1..ntargetmax species indices
  read (9, '(A)') line_buffer
  ntarget = 0
  do t = 1, ntargetmax
    read (line_buffer, *, IOSTAT=ierr) (sptarget(i), i=1,t)
    if (ierr /= 0) exit
    ntarget = t
  end do
  if (ntarget == 0) then
    write (*, *) "Error: could not parse sptarget from: ", trim(line_buffer)
    stop 1
  end if
  read (9, *)
  read (9, *)
  read (9, *)
  read (9, *)
! Read nsubs: one line per target site.
! ntarget==1: single integer, space-separated integers (multi-nary), or colon range X:Y.
! ntarget==2: first line for target 1, second line for target 2.
  read (9, '(A)') line_buffer
  nk(:) = 1
  if (ntarget >= 2) then
!   Multi-target: space-separated integers per line; 1 integer = binary, 2+ = multi-nary.
    if (index(trim(line_buffer), '/') > 0) then
      write (*, *) "Error: '/' separator in nsubs is no longer supported."
      write (*, *) "  Use one line per target site instead (e.g. first line: 2, second line: 1)."
      stop 1
    end if
!   Parse target 1 (line_buffer already read above)
    nk(1) = 0
    do j = 1, nkmax
      read (line_buffer, *, IOSTAT=ierr) (nsubs_t(1,i), i=1,j)
      if (ierr /= 0) exit
      nk(1) = j
    end do
    if (nk(1) == 0) then
      write (*, *) "Error: could not parse nsubs for target 1 from: ", trim(line_buffer)
      stop 1
    end if
    if (nk(1) > 3) then
      write (*, *) "Error: more than 3 species per site (k=", nk(1), ") not yet supported."
      stop 1
    end if
    do j = 1, nk(1)
      if (nsubs_t(1,j) < 0) then
        write (*, *) "Error: number of substitutions must be >= 0"
        stop 1
      end if
    end do
!   Parse targets 2..ntarget
    do t = 2, ntarget
      read (9, '(A)') line_buffer
      nk(t) = 0
      do j = 1, nkmax
        read (line_buffer, *, IOSTAT=ierr) (nsubs_t(t,i), i=1,j)
        if (ierr /= 0) exit
        nk(t) = j
      end do
      if (nk(t) == 0) then
        write (*, '(a,i0,a,a)') "Error: could not parse nsubs for target ", t, " from: ", trim(line_buffer)
        stop 1
      end if
      if (nk(t) > 3) then
        write (*, *) "Error: more than 3 species per site (k=", nk(t), ") not yet supported."
        stop 1
      end if
      do j = 1, nk(t)
        if (nsubs_t(t,j) < 0) then
          write (*, *) "Error: number of substitutions must be >= 0"
          stop 1
        end if
      end do
    end do
    nsubs_min = nsubs_t(1,1)
    nsubs_max = nsubs_t(1,1)
  else if (index(trim(line_buffer), ':') > 0) then
!   Colon range notation: nsubs_min:nsubs_max (single target only)
    j = index(trim(line_buffer), ':')
    read (line_buffer(1:j-1), *, IOSTAT=ierr) nsubs_min
    if (ierr /= 0) then
      write (*, *) "Error: could not parse nsubs_min from colon notation: ", trim(line_buffer)
      stop 1
    end if
    read (line_buffer(j+1:), *, IOSTAT=ierr) nsubs_max
    if (ierr /= 0) then
      write (*, *) "Error: could not parse nsubs_max from colon notation: ", trim(line_buffer)
      stop 1
    end if
    if (nsubs_min < 0) then
      write (*, *) "Error: number of substitutions must be >= 0"
      stop 1
    end if
    if (nsubs_max < nsubs_min) then
      write (*, *) "Error: nsubs_max must be >= nsubs_min"
      stop 1
    end if
  else
!   Single integer or space-separated integers (binary, or multi-nary; single target)
    nk(1) = 0
    do j = 1, nkmax
      read (line_buffer, *, IOSTAT=ierr) (nsubs_t(1,i), i=1,j)
      if (ierr /= 0) exit
      nk(1) = j
    end do
    if (nk(1) == 0) then
      write (*, *) "Error: could not parse nsubs line: ", trim(line_buffer)
      stop 1
    end if
    if (nk(1) > 3) then
      write (*, *) "Error: more than 3 species per site (k=", nk(1), ") not yet supported."
      write (*, *) "  Supported: k=1 (binary), k=2 and k=3 (multi-nary)."
      stop 1
    end if
    do j = 1, nk(1)
      if (nsubs_t(1,j) < 0) then
        write (*, *) "Error: number of substitutions must be >= 0"
        stop 1
      end if
    end do
    nsubs_min = nsubs_t(1,1)
    nsubs_max = nsubs_t(1,1)
  end if
  read (9, *)
  read (9, *)
  read (9, *)
  read (9, *)
  read (9, *) (newsymbol(1, j), j=1, nk(1)+1)
  do t = 2, ntarget
    read (9, *) (newsymbol(t, j), j=1, nk(t)+1)
  end do
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
    stop 1
  end if

  write (46, *) nop
  do op = 1, nop
    write (46, *) "Operator number ", op
    do i = 1, 3
      write (46, *) mgroup(op, i, 1:3), vgroup(op, i)
    end do
  end do

  write (*, *) " > Composition of the parent supercell:"
  do sp = 1, nsp
    write (*, '(a, a, i10)') "    - ", symbol(sp), natsp(sp)
  end do
  write (*, *) ""
  write (*, '(a, i10)') " > Number of symmetry operators in the supercell: ", nop
  write (*, *) ""

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

! Compute atini/atfin and eqmt for all target species

  do t_A = 1, ntarget
    atini_t(t_A) = 1
    do sp = 1, sptarget(t_A) - 1
      atini_t(t_A) = atini_t(t_A) + natsp(sp)
    end do
    atfin_t(t_A) = atini_t(t_A) + natsp(sptarget(t_A)) - 1
    npos_t(t_A) = natsp(sptarget(t_A))
  end do

! Scalar aliases for ntarget==1 code paths
  atini = atini_t(1)
  atfin = atfin_t(1)
  npos  = npos_t(1)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Obtain the Equivalence Matrix for the target species from the Full Equivalence Matrix
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  npos_max = maxval(npos_t(1:ntarget))
  allocate (eqmt(1:ntarget, 1:npos_max, 1:nop), stat=ierr)
  if (ierr /= 0) then
    write (*, *) "Error: problem too large for SOD - insufficient memory (eqmt)"
    stop 1
  end if
  do t = 1, ntarget
    do op = 1, nop
      att = 0
      do at = atini_t(t), atfin_t(t)
        att = att + 1
        eqmt(t, att, op) = fulleqmatrix(op, at) - atini_t(t) + 1
      end do
    end do
    write (26, *) nop, npos_t(t)
    do op = 1, nop
      write (26, *) (eqmt(t, pos, op), pos=1, npos_t(t))
    end do
  end do

  deallocate (fulleqmatrix)

  if (ntarget == 1 .and. nk(1) >= 2) then
!   Multi-nary: check total substitutions against available sites
    if (sum(nsubs_t(1, 1:nk(1))) > npos_t(1)) then
      write (*, *) "Error: total substitutions (", sum(nsubs_t(1, 1:nk(1))), &
        ") exceed available sites (", npos_t(1), ")"
      stop 1
    end if
  else if (nsubs_max > npos_t(1)) then
    write (*, *) "Error: nsubs_max (", nsubs_max, ") exceeds number of available sites (", npos_t(1), ")"
    stop 1
  end if
  do t = 2, ntarget
    if (sum(nsubs_t(t, 1:nk(t))) > npos_t(t)) then
      write (*, '(a,i0,a,i0,a,i0,a)') "Error: total nsubs for target ", t, &
        " (", sum(nsubs_t(t, 1:nk(t))), ") exceeds available sites (", npos_t(t), ")"
      stop 1
    end if
  end do

! Direction selection: only for binary substitution (single new-species range)
  if (ntarget == 1 .and. nk(1) == 1) then
    going_upward = (nsubs_min <= npos - nsubs_max)
    if (nsubs_min < nsubs_max) then
      if (going_upward) then
        write (*, '(a,i0,a)') " > Direction: UPWARD (starting at nsubs = ", nsubs_min, ")"
      else
        write (*, '(a,i0,a)') " > Direction: DOWNWARD (starting at nsubs = ", nsubs_max, ")"
      end if
    end if
  end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Loop over substitution levels
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if (ntarget == 1 .and. nk(1) == 1) then

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
    write (30, '(a)') "# SOD OUTSOD format version 2"
    open (unit=47, file=trim(nxx_dir) // "/cSGO")

    write (*, *) "----------------------------------------------------------------------------"
    write (*, '(a,i0)') " Substitution level: nsubs = ", nsubs
    write (*, *) "----------------------------------------------------------------------------"
    write (*, *) ""
    write (*, *) " > Composition of the substituted supercell:"
    do sp = 1, nsp
      if (ntarget == 1 .and. sp == sptarget(1)) then
        write (*, '(a, a, i10)') "    - ", newsymbol(1,1), nsubs
        write (*, '(a, a, i10)') "    - ", newsymbol(1,2), natsp(sp) - nsubs
      else
        write (*, '(a, a, i10)') "    - ", symbol(sp), natsp(sp)
      end if
    end do
    write (*, *) ""
    write (*, *) ""

! Handle trivial cases: 0 substitutions or all sites substituted (1 configuration each)

    if (nsubs == 0 .or. nsubs == npos) then
      write (*, *) " "
      write (*, '(a)') "       Number of inequivalent configurations: 1"
      write (*, *) " "
      write (30, *) nsubs, " substitutions in ", npos, "sites"
      write (30, *) 1, " configurations"
      if (nsubs == 0) then
        write (30, '(i6, 1x, i6)') 1, 1
      else
        write (30, '(i6, 1x, i6, *(1x, i4))') 1, 1, (pos + atini - 1, pos=1, npos)
      end if
      close (30)
      close (47)
      cycle
    end if

! Precompute binomial coefficient table: binom(k,n) = C(n,k), k=0..nsubs, n=0..npos
    allocate (binom(0:nsubs, 0:npos))
    binom(0, :) = 1_int64
    do k = 1, nsubs
      binom(k, 0:k - 1) = 0_int64
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
! Read header line 1: "N substitutions in M sites" (skip any leading # comment lines)
      do
        read (50, '(A)') line_buffer
        if (line_buffer(1:1) /= '#') exit
      end do
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
        ! OUTSOD stores global atom indices; convert to local (1..npos) for the expansion
        indconf_prev(1:nic_prev, 1:nsubs_prev) = indconf_prev(1:nic_prev, 1:nsubs_prev) - atini + 1
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
          stop 1
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
          stop 1
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
        stop 1
      end if
      allocate (visited(1:ntc), stat=ierr)
      if (ierr /= 0) then
        write (*, *) "Error: insufficient memory for visited array (size", ntc, ")"
        write (*, *) "       Try reducing supercell or substitutions."
        stop 1
      end if
      visited(:) = .false.
      allocate (newconf(1:nsubs), stat=ierr)
      allocate (indconf(1:ncand, 1:nsubs), stat=ierr)
      if (ierr /= 0) then
        write (*, *) "Error: insufficient memory for indconf"
        stop 1
      end if
      allocate (degen(1:ncand), stat=ierr)
      if (ierr /= 0) then
        write (*, *) "Error: insufficient memory for degen"
        stop 1
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
      do cand = 1, int(ncand, int64)
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
          do concurrent (i = 1:nsubs)
            newconf(i) = eqmt(1, conf(cand, i), op)
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
        stop 1
      end if
      maxentropy = kb * log(real(ntc, real64))
      x = real(nsubs, real64) / real(npos, real64)
      if (x > 0.d0 .and. x < 1.d0) then
        ientropy = -kb * real(npos, real64) * (x * log(x) + (1.d0 - x) * log(1.d0 - x))
      else
        ientropy = 0.d0
      end if
      if (ientropy > 0.d0) then
        perc = maxentropy / ientropy * 100.d0
      else
        perc = 0.d0
      end if
      write (*, *) " "
      write (*, '(a, i0)')       " > Total configurations in supercell:           ", ntc
      write (*, '(a, f12.6, a)') " > Maximum entropy (supercell ensemble):        ", maxentropy * 1000.d0, " meV/K"
      write (*, '(a, f12.6)')    " > Fraction of substituted sites (x):           ", x
      write (*, '(a, f12.6, a)') " > Ideal mixing entropy (per cell):             ", ientropy * 1000.d0, " meV/K"
      write (*, '(a, f7.2, a)')  " > Supercell captures                           ", perc, "% of ideal entropy"
      write (*, *) ""

!!!!!!!!Allocating array sizes
! conf(ntc, nsubs) is NOT allocated; configurations are generated one at a time via ksubset.
! indconf is allocated at ntc (safe upper bound on the number of independent configs).

      allocate (degen(1:ntc), stat=ierr)
      if (ierr /= 0) then
        write (*, *) "Error: problem too large for SOD - insufficient memory"
        write (*, *) "       Total configurations: ", ntc, " - reduce supercell size or number of substitutions"
        stop 1
      end if
      allocate (newconf(1:nsubs), stat=ierr)
      allocate (indconf(1:ntc, 1:nsubs), stat=ierr)
      if (ierr /= 0) then
        write (*, *) "Error: problem too large for SOD - insufficient memory"
        write (*, *) "       Total configurations: ", ntc, " - reduce supercell size or number of substitutions"
        stop 1
      end if
      allocate (visited(1:ntc), stat=ierr)
      if (ierr /= 0) then
        write (*, *) "Error: problem too large for SOD - insufficient memory"
        write (*, *) "       Total configurations: ", ntc, " - reduce supercell size or number of substitutions"
        stop 1
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
        do concurrent (i = 1:nsubs)
          newconf(i) = eqmt(1, as(i), op)
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
            do concurrent (i = 1:nsubs)
              newconf(i) = eqmt(1, as(i), op)
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
            write (*, '(3x,i10,10x,f6.2)') indcount, perc
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

    write (*, '(a,i0)') "  > Number of inequivalent configurations: ", nic
    write (*, *) ""

    write (30, *) nsubs, " substitutions in ", npos, "sites"
    write (30, *) nic, " configurations"
    do indcount = 1, nic
      write (30, 330) indcount, degen(indcount), indconf(indcount, 1:nsubs) + atini - 1
    end do
330 format(i6, 1x, i6, *(1x, i4))

!!!!!!!Deallocating arrays
    deallocate (newconf)
    deallocate (degen)
    deallocate (indconf)

    deallocate (binom)

    call system_clock(t_now)
    elapsed_level = real(t_now - t_start_level, real64) / real(clock_rate, real64)
    elapsed_per_level(nsubs_loop) = elapsed_level
    nsubs_level_list(nsubs_loop) = nsubs

    close (30)
    close (47)

  end do  ! nsubs_loop

  else if (ntarget == 1 .and. nk(1) == 2) then

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Multi-nary substitution on one site (k=2), direct enumeration
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! nsubs = total substitutions on site 1
  nsubs = nsubs_t(1,1) + nsubs_t(1,2)

! Build directory name: n01_04/ etc.
  write (nxx_dir, '("n", i2.2, "_", i2.2)') nsubs_t(1,1), nsubs_t(1,2)
  call execute_command_line("mkdir -p " // trim(nxx_dir), wait=.true.)
  open (unit=30, file=trim(nxx_dir) // "/OUTSOD")
  write (30, '(a)') "# SOD OUTSOD format version 2"
  open (unit=47, file=trim(nxx_dir) // "/cSGO")

  write (*, *) " "
  write (*, '(a,i4,a,i4,a)') " === Multi-nary substitution: nsubs = [", &
    nsubs_t(1,1), ",", nsubs_t(1,2), "] ==="
  write (*, *) " "
  write (*, *) "       Composition of the substituted supercell:"
  do sp = 1, nsp
    if (sp == sptarget(1)) then
      write (*, *) "                                                         ", newsymbol(1,1), nsubs_t(1,1)
      write (*, *) "                                                         ", newsymbol(1,2), nsubs_t(1,2)
      write (*, *) "                                                         ", newsymbol(1,3), natsp(sp) - nsubs
    else
      write (*, *) "                                                         ", symbol(sp), natsp(sp)
    end if
  end do
  write (*, *) ""

! Binomial tables:
!   binom_A(k, n) = C(n,k) for k=0..nsubs (combined selection from npos)
!   binom_B(k, n) = C(n,k) for k=0..nsubs_t(1,1) (colouring: sp1 from combined)
  allocate (binom_A(0:nsubs, 0:npos))
  binom_A(0, :) = 1_int64
  do k = 1, nsubs
    binom_A(k, 0:k-1) = 0_int64
    do j = k, npos
      binom_A(k, j) = binom_A(k-1, j-1) + binom_A(k, j-1)
    end do
  end do

  allocate (binom_B(0:nsubs_t(1,1), 0:nsubs))
  binom_B(0, :) = 1_int64
  do k = 1, nsubs_t(1,1)
    binom_B(k, 0:k-1) = 0_int64
    do j = k, nsubs
      binom_B(k, j) = binom_B(k-1, j-1) + binom_B(k, j-1)
    end do
  end do

  ntc_A = binom_A(nsubs, npos)                      ! C(npos, nsubs_tot)
  ntc_B = binom_B(nsubs_t(1,1), nsubs)              ! C(nsubs_tot, n1) = multinomial_size
  ntc_joint = ntc_A * ntc_B

  write (*, '(a,i0,a,i0,a,i0)') "       Combined-position space: C(", npos, ", ", nsubs, ") = ", ntc_A
  write (*, '(a,i0,a,i0,a,i0)') "       sp1 colouring space:      C(", nsubs, ", ", nsubs_t(1,1), ") = ", ntc_B
  write (*, '(a,i0,a)')         "       sp2 fills remaining ", nsubs_t(1,2), " positions"
  write (*, '(a,i0)')           "       Total joint configurations: ", ntc_joint
  maxentropy = kb * log(real(ntc_joint, real64))
  ientropy = 0.d0
  x = real(nsubs_t(1,1), kind=8) / real(npos, real64)
  if (x > 0.d0) ientropy = ientropy - kb * x * log(x)
  x = real(nsubs_t(1,2), kind=8) / real(npos, real64)
  if (x > 0.d0) ientropy = ientropy - kb * x * log(x)
  x = real(npos - nsubs, real64) / real(npos, real64)
  if (x > 0.d0) ientropy = ientropy - kb * x * log(x)
  ientropy = ientropy * real(npos, real64)
  if (ientropy > 0.d0) then
    perc = maxentropy / ientropy * 100.d0
  else
    perc = 0.d0
  end if
  write (*, '(a,f12.6,a)') "       Maximum entropy (supercell ensemble):    ", maxentropy * 1000.d0, " meV/K"
  write (*, '(a,f12.6,a)') "       Ideal mixing entropy (per cell):          ", ientropy * 1000.d0, " meV/K"
  write (*, '(a,f7.2,a)')  "       Supercell captures                        ", perc, "% of ideal entropy"
  write (*, *) " "
  allocate (visited(1:ntc_joint), stat=ierr)
  if (ierr /= 0) then
    write (*, *) "Error: insufficient memory for visited array (size", ntc_joint, ")"
    stop 1
  end if
  visited(:) = .false.

  allocate (as_A(1:nsubs), stat=ierr)
  allocate (as_B(1:nsubs_t(1,1)), stat=ierr)
  allocate (newconf_A(1:nsubs), stat=ierr)
  allocate (newconf_B(1:nsubs_t(1,1)), stat=ierr)
  allocate (new_as_B(1:nsubs_t(1,1)), stat=ierr)
  allocate (indconf(1:ntc_joint, 1:nsubs), stat=ierr)
  if (ierr /= 0) then
    write (*, *) "Error: insufficient memory for indconf"
    stop 1
  end if
  allocate (degen(1:ntc_joint), stat=ierr)
  indconf(:,:) = 0
  degen(:) = 1

  write (*, *) "Finding the inequivalent configurations..."
  write (*, *) " "
  write (*, *) "       Found    Completion "
  write (*, *) "       =====    ========== "

  indcount = 0
  perc = 0.0d0
  equivcount = 0

! Outer loop: enumerate combined positions (nsubs from npos)
! Initialize as_A = [1, 2, ..., nsubs]; no ksubset (avoids shared save-variable conflict)
  do i = 1, nsubs
    as_A(i) = i
  end do
  mores_A = (nsubs < npos)

  do   ! outer combined loop
!   Rank of combined selection: combinatorial number system
    r_A = binom_A(nsubs, npos)
    do i = 1, nsubs
      if (npos - as_A(i) >= nsubs - i + 1) &
        r_A = r_A - binom_A(nsubs - i + 1, npos - as_A(i))
    end do

!   Inner loop: enumerate sp1 colouring (nsubs_t(1,1) from nsubs)
!   Initialize as_B = [1, 2, ..., nsubs_t(1,1)]
    do i = 1, nsubs_t(1,1)
      as_B(i) = i
    end do
    mores_B = (nsubs_t(1,1) < nsubs)

    do   ! inner colouring loop
!     Rank of colouring
      r_B = binom_B(nsubs_t(1,1), nsubs)
      do i = 1, nsubs_t(1,1)
        if (nsubs - as_B(i) >= nsubs_t(1,1) - i + 1) &
          r_B = r_B - binom_B(nsubs_t(1,1) - i + 1, nsubs - as_B(i))
      end do

!     Joint rank (1-based)
      r = (r_A - 1_int64) * ntc_B + r_B

      if (.not. visited(r)) then
        equivcount = equivcount + 1
        indcount = indcount + 1

!       Store sp1 positions, then sp2 positions (as relative indices 1..npos)
        do i = 1, nsubs_t(1,1)
          indconf(indcount, i) = as_A(as_B(i))
        end do
        k = nsubs_t(1,1) + 1
        do i = 1, nsubs
          found = .false.
          do jj = 1, nsubs_t(1,1)
            if (as_B(jj) == i) then
              found = .true.
              exit
            end if
          end do
          if (.not. found) then
            indconf(indcount, k) = as_A(i)
            k = k + 1
          end if
        end do

        visited(r) = .true.
        write (47, *) "List of operators for configuration: ", indcount
        opc = 1
        write (47, *) opc
        do i = 1, 3
          write (47, *) (mgroup(1, i, j), j=1, 3), vgroup(1, i)
        end do

!       Apply symmetry operators
        do op = 2, nop
!         Image of all combined positions
          do concurrent (i = 1:nsubs)
            newconf_A(i) = eqmt(1, as_A(i), op)
          end do
          call bubble(newconf_A, nsubs)

!         Image of sp1 positions (values, then find indices in new combined)
          do concurrent (i = 1:nsubs_t(1,1))
            newconf_B(i) = eqmt(1, as_A(as_B(i)), op)
          end do
          call bubble(newconf_B, nsubs_t(1,1))

!         Find indices of sp1 in new combined (two-pointer: both sorted)
          jj = 1
          do i = 1, nsubs_t(1,1)
            do while (newconf_A(jj) /= newconf_B(i))
              jj = jj + 1
            end do
            new_as_B(i) = jj
            jj = jj + 1
          end do
          call bubble(new_as_B, nsubs_t(1,1))

!         Compute new ranks
          new_r_A = binom_A(nsubs, npos)
          do i = 1, nsubs
            if (npos - newconf_A(i) >= nsubs - i + 1) &
              new_r_A = new_r_A - binom_A(nsubs - i + 1, npos - newconf_A(i))
          end do
          new_r_B = binom_B(nsubs_t(1,1), nsubs)
          do i = 1, nsubs_t(1,1)
            if (nsubs - new_as_B(i) >= nsubs_t(1,1) - i + 1) &
              new_r_B = new_r_B - binom_B(nsubs_t(1,1) - i + 1, nsubs - new_as_B(i))
          end do

          new_r_A = (new_r_A - 1_int64) * ntc_B + new_r_B

          if (all(newconf_A(1:nsubs) == as_A(1:nsubs)) .and. &
              all(new_as_B(1:nsubs_t(1,1)) == as_B(1:nsubs_t(1,1)))) then
            opc = opc + 1
            write (47, *) opc
            do i = 1, 3
              write (47, *) (mgroup(op, i, j), j=1, 3), vgroup(op, i)
            end do
          else if (.not. visited(new_r_A)) then
            visited(new_r_A) = .true.
            degen(indcount) = degen(indcount) + 1
            equivcount = equivcount + 1
          end if
        end do  ! op loop
        write (47, *) 0

        if ((100.0d0*equivcount/ntc_joint - perc > 5.0d0) .or. (equivcount == ntc_joint)) then
          perc = 100.0d0*equivcount/ntc_joint
          write (*, '(4x,i6,7x,f5.1,a2)') indcount, perc, "% "
        end if
      end if  ! not visited

      if (.not. mores_B) exit
!     Advance as_B to next nsubs_t(1,1)-subset of {1..nsubs} (inline, no shared state)
      do i = nsubs_t(1,1), 1, -1
        if (as_B(i) < nsubs - nsubs_t(1,1) + i) then
          as_B(i) = as_B(i) + 1
          do j_remove = i+1, nsubs_t(1,1)
            as_B(j_remove) = as_B(j_remove-1) + 1
          end do
          mores_B = (as_B(1) /= nsubs - nsubs_t(1,1) + 1)
          exit
        end if
      end do
    end do  ! inner colouring loop

    if (.not. mores_A) exit
!   Advance as_A to next nsubs-subset of {1..npos} (inline, no shared state)
    do i = nsubs, 1, -1
      if (as_A(i) < npos - nsubs + i) then
        as_A(i) = as_A(i) + 1
        do j_remove = i+1, nsubs
          as_A(j_remove) = as_A(j_remove-1) + 1
        end do
        mores_A = (as_A(1) /= npos - nsubs + 1)
        exit
      end if
    end do
  end do  ! outer combined loop

  nic = indcount

! Trim arrays to actual size
  block
    integer, dimension(:,:), allocatable :: tmp2d
    integer, dimension(:), allocatable :: tmp1d
    allocate (tmp2d(1:nic, 1:nsubs))
    tmp2d = indconf(1:nic, 1:nsubs)
    call move_alloc(tmp2d, indconf)
    allocate (tmp1d(1:nic))
    tmp1d = degen(1:nic)
    call move_alloc(tmp1d, degen)
  end block

  write (*, *) " "
  write (*, '(a,i0)') "       Number of inequivalent configurations: ", nic
  write (*, *) " "

! Write OUTSOD header and data (absolute positions)
  write (30, *) nsubs_t(1,1), nsubs_t(1,2), " substitutions in ", npos, "sites"
  write (30, *) nic, " configurations"
  do indcount = 1, nic
    write (30, 331) indcount, degen(indcount), indconf(indcount, 1:nsubs) + atini - 1
  end do

  deallocate (as_A, as_B, newconf_A, newconf_B, new_as_B)
  deallocate (degen)
  deallocate (indconf)
  deallocate (binom_A, binom_B)
  deallocate (visited)

  close (30)
  close (47)

  else if (ntarget == 1 .and. nk(1) == 3) then

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Multi-nary substitution on one site (k=3), direct enumeration
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! nsubs = total substitutions on site 1
  nsubs = nsubs_t(1,1) + nsubs_t(1,2) + nsubs_t(1,3)
  n2rem = nsubs_t(1,2) + nsubs_t(1,3)   ! positions after removing sp1

! Build directory name: n01_02_03/ etc.
  write (nxx_dir, '("n", i2.2, "_", i2.2, "_", i2.2)') nsubs_t(1,1), nsubs_t(1,2), nsubs_t(1,3)
  call execute_command_line("mkdir -p " // trim(nxx_dir), wait=.true.)
  open (unit=30, file=trim(nxx_dir) // "/OUTSOD")
  write (30, '(a)') "# SOD OUTSOD format version 2"
  open (unit=47, file=trim(nxx_dir) // "/cSGO")

  write (*, *) " "
  write (*, '(a,i4,a,i4,a,i4,a)') " === Multi-nary substitution: nsubs = [", &
    nsubs_t(1,1), ",", nsubs_t(1,2), ",", nsubs_t(1,3), "] ==="
  write (*, *) " "
  write (*, *) "       Composition of the substituted supercell:"
  do sp = 1, nsp
    if (sp == sptarget(1)) then
      write (*, *) "                                                         ", newsymbol(1,1), nsubs_t(1,1)
      write (*, *) "                                                         ", newsymbol(1,2), nsubs_t(1,2)
      write (*, *) "                                                         ", newsymbol(1,3), nsubs_t(1,3)
      write (*, *) "                                                         ", newsymbol(1,4), natsp(sp) - nsubs
    else
      write (*, *) "                                                         ", symbol(sp), natsp(sp)
    end if
  end do
  write (*, *) ""

! Binomial tables:
!   binom_A(k, n) = C(n,k) for k=0..nsubs (combined selection from npos)
!   binom_B(k, n) = C(n,k) for k=0..nsubs_t(1,1) (sp1 colouring from nsubs)
!   binom_C(k, n) = C(n,k) for k=0..nsubs_t(1,2) (sp2 colouring from n2rem)
  allocate (binom_A(0:nsubs, 0:npos))
  binom_A(0, :) = 1_int64
  do k = 1, nsubs
    binom_A(k, 0:k-1) = 0_int64
    do j = k, npos
      binom_A(k, j) = binom_A(k-1, j-1) + binom_A(k, j-1)
    end do
  end do

  allocate (binom_B(0:nsubs_t(1,1), 0:nsubs))
  binom_B(0, :) = 1_int64
  do k = 1, nsubs_t(1,1)
    binom_B(k, 0:k-1) = 0_int64
    do j = k, nsubs
      binom_B(k, j) = binom_B(k-1, j-1) + binom_B(k, j-1)
    end do
  end do

  allocate (binom_C(0:nsubs_t(1,2), 0:n2rem))
  binom_C(0, :) = 1_int64
  do k = 1, nsubs_t(1,2)
    binom_C(k, 0:k-1) = 0_int64
    do j = k, n2rem
      binom_C(k, j) = binom_C(k-1, j-1) + binom_C(k, j-1)
    end do
  end do

  ntc_A = binom_A(nsubs, npos)                       ! C(npos, nsubs_tot)
  ntc_B = binom_B(nsubs_t(1,1), nsubs)               ! C(nsubs_tot, n1)
  ntc_C = binom_C(nsubs_t(1,2), n2rem)               ! C(n2rem, n2)
  ntc_joint = ntc_A * ntc_B * ntc_C

  write (*, '(a,i0,a,i0,a,i0)') "       Combined-position space: C(", npos, ", ", nsubs, ") = ", ntc_A
  write (*, '(a,i0,a,i0,a,i0)') "       sp1 colouring space:      C(", nsubs, ", ", nsubs_t(1,1), ") = ", ntc_B
  write (*, '(a,i0,a,i0,a,i0)') "       sp2 colouring space:      C(", n2rem, ", ", nsubs_t(1,2), ") = ", ntc_C
  write (*, '(a,i0,a)')         "       sp3 fills remaining ", nsubs_t(1,3), " positions"
  write (*, '(a,i0)')           "       Total joint configurations: ", ntc_joint
  maxentropy = kb * log(real(ntc_joint, real64))
  ientropy = 0.d0
  x = real(nsubs_t(1,1), kind=8) / real(npos, real64)
  if (x > 0.d0) ientropy = ientropy - kb * x * log(x)
  x = real(nsubs_t(1,2), kind=8) / real(npos, real64)
  if (x > 0.d0) ientropy = ientropy - kb * x * log(x)
  x = real(nsubs_t(1,3), kind=8) / real(npos, real64)
  if (x > 0.d0) ientropy = ientropy - kb * x * log(x)
  x = real(npos - nsubs, real64) / real(npos, real64)
  if (x > 0.d0) ientropy = ientropy - kb * x * log(x)
  ientropy = ientropy * real(npos, real64)
  if (ientropy > 0.d0) then
    perc = maxentropy / ientropy * 100.d0
  else
    perc = 0.d0
  end if
  write (*, '(a,f12.6,a)') "       Maximum entropy (supercell ensemble):    ", maxentropy * 1000.d0, " meV/K"
  write (*, '(a,f12.6,a)') "       Ideal mixing entropy (per cell):          ", ientropy * 1000.d0, " meV/K"
  write (*, '(a,f7.2,a)')  "       Supercell captures                        ", perc, "% of ideal entropy"
  write (*, *) " "

  allocate (visited(1:ntc_joint), stat=ierr)
  if (ierr /= 0) then
    write (*, *) "Error: insufficient memory for visited array (size", ntc_joint, ")"
    stop 1
  end if
  visited(:) = .false.

  allocate (as_A(1:nsubs), stat=ierr)
  allocate (as_B(1:nsubs_t(1,1)), stat=ierr)
  allocate (as_C(1:nsubs_t(1,2)), stat=ierr)
  allocate (newconf_A(1:nsubs), stat=ierr)
  allocate (newconf_B(1:nsubs_t(1,1)), stat=ierr)
  allocate (newconf_C(1:nsubs_t(1,2)), stat=ierr)
  allocate (new_as_B(1:nsubs_t(1,1)), stat=ierr)
  allocate (new_as_C(1:nsubs_t(1,2)), stat=ierr)
  allocate (as_notB_tmp(1:n2rem), stat=ierr)
  allocate (new_notB_tmp(1:n2rem), stat=ierr)
  allocate (indconf(1:ntc_joint, 1:nsubs), stat=ierr)
  if (ierr /= 0) then
    write (*, *) "Error: insufficient memory for indconf"
    stop 1
  end if
  allocate (degen(1:ntc_joint), stat=ierr)
  indconf(:,:) = 0
  degen(:) = 1

  write (*, *) "Finding the inequivalent configurations..."
  write (*, *) " "
  write (*, *) "       Found    Completion "
  write (*, *) "       =====    ========== "

  indcount = 0
  perc = 0.0d0
  equivcount = 0

! Outer loop: enumerate combined positions (nsubs from npos)
  do i = 1, nsubs
    as_A(i) = i
  end do
  mores_A = (nsubs < npos)

  do   ! outer combined loop
!   Rank of combined selection: combinatorial number system
    r_A = binom_A(nsubs, npos)
    do i = 1, nsubs
      if (npos - as_A(i) >= nsubs - i + 1) &
        r_A = r_A - binom_A(nsubs - i + 1, npos - as_A(i))
    end do

!   Middle loop: enumerate sp1 colouring (nsubs_t(1,1) from nsubs)
    do i = 1, nsubs_t(1,1)
      as_B(i) = i
    end do
    mores_B = (nsubs_t(1,1) < nsubs)

    do   ! middle colouring loop (sp1)
!     Rank of sp1 colouring
      r_B = binom_B(nsubs_t(1,1), nsubs)
      do i = 1, nsubs_t(1,1)
        if (nsubs - as_B(i) >= nsubs_t(1,1) - i + 1) &
          r_B = r_B - binom_B(nsubs_t(1,1) - i + 1, nsubs - as_B(i))
      end do

!     Build as_notB_tmp: indices in {1..nsubs} not in as_B (complement)
      jj = 0
      j = 1
      do i = 1, nsubs
        if (j <= nsubs_t(1,1) .and. as_B(j) == i) then
          j = j + 1
        else
          jj = jj + 1
          as_notB_tmp(jj) = i
        end if
      end do

!     Inner loop: enumerate sp2 colouring (nsubs_t(1,2) from n2rem)
      do i = 1, nsubs_t(1,2)
        as_C(i) = i
      end do
      mores_C = (nsubs_t(1,2) < n2rem)

      do   ! inner colouring loop (sp2)
!       Rank of sp2 colouring
        r_C = binom_C(nsubs_t(1,2), n2rem)
        do i = 1, nsubs_t(1,2)
          if (n2rem - as_C(i) >= nsubs_t(1,2) - i + 1) &
            r_C = r_C - binom_C(nsubs_t(1,2) - i + 1, n2rem - as_C(i))
        end do

!       Joint rank (1-based)
        r = (r_A - 1_int64) * ntc_B * ntc_C + (r_B - 1_int64) * ntc_C + r_C

        if (.not. visited(r)) then
          equivcount = equivcount + 1
          indcount = indcount + 1

!         Store sp1 positions: as_A(as_B(i))
          do i = 1, nsubs_t(1,1)
            indconf(indcount, i) = as_A(as_B(i))
          end do
!         Store sp2 positions: as_A(as_notB_tmp(as_C(i)))
          do i = 1, nsubs_t(1,2)
            indconf(indcount, nsubs_t(1,1) + i) = as_A(as_notB_tmp(as_C(i)))
          end do
!         Store sp3 positions: remaining (in as_A but not sp1 or sp2)
          k = nsubs_t(1,1) + nsubs_t(1,2) + 1
          do i = 1, nsubs
            found = .false.
            do jj = 1, nsubs_t(1,1)
              if (as_B(jj) == i) then
                found = .true.
                exit
              end if
            end do
            if (.not. found) then
              do jj = 1, nsubs_t(1,2)
                if (as_notB_tmp(as_C(jj)) == i) then
                  found = .true.
                  exit
                end if
              end do
            end if
            if (.not. found) then
              indconf(indcount, k) = as_A(i)
              k = k + 1
            end if
          end do

          visited(r) = .true.
          write (47, *) "List of operators for configuration: ", indcount
          opc = 1
          write (47, *) opc
          do i = 1, 3
            write (47, *) (mgroup(1, i, j), j=1, 3), vgroup(1, i)
          end do

!         Apply symmetry operators
          do op = 2, nop
!           Image of all combined positions
            do concurrent (i = 1:nsubs)
              newconf_A(i) = eqmt(1, as_A(i), op)
            end do
            call bubble(newconf_A, nsubs)

!           Image of sp1 positions
            do concurrent (i = 1:nsubs_t(1,1))
              newconf_B(i) = eqmt(1, as_A(as_B(i)), op)
            end do
            call bubble(newconf_B, nsubs_t(1,1))

!           Find indices of sp1 in new combined (two-pointer: both sorted)
            jj = 1
            do i = 1, nsubs_t(1,1)
              do while (newconf_A(jj) /= newconf_B(i))
                jj = jj + 1
              end do
              new_as_B(i) = jj
              jj = jj + 1
            end do
            call bubble(new_as_B, nsubs_t(1,1))

!           Build new_notB_vals: sorted values in newconf_A not in newconf_B
            jj = 0
            j = 1
            do i = 1, nsubs
              if (j <= nsubs_t(1,1) .and. new_as_B(j) == i) then
                j = j + 1
              else
                jj = jj + 1
                new_notB_tmp(jj) = newconf_A(i)
              end if
            end do

!           Image of sp2 positions
            do concurrent (i = 1:nsubs_t(1,2))
              newconf_C(i) = eqmt(1, as_A(as_notB_tmp(as_C(i))), op)
            end do
            call bubble(newconf_C, nsubs_t(1,2))

!           Find indices of sp2 in new_notB_vals (two-pointer: both sorted)
            jj = 1
            do i = 1, nsubs_t(1,2)
              do while (new_notB_tmp(jj) /= newconf_C(i))
                jj = jj + 1
              end do
              new_as_C(i) = jj
              jj = jj + 1
            end do
            call bubble(new_as_C, nsubs_t(1,2))

!           Compute new ranks
            new_r_A = binom_A(nsubs, npos)
            do i = 1, nsubs
              if (npos - newconf_A(i) >= nsubs - i + 1) &
                new_r_A = new_r_A - binom_A(nsubs - i + 1, npos - newconf_A(i))
            end do
            new_r_B = binom_B(nsubs_t(1,1), nsubs)
            do i = 1, nsubs_t(1,1)
              if (nsubs - new_as_B(i) >= nsubs_t(1,1) - i + 1) &
                new_r_B = new_r_B - binom_B(nsubs_t(1,1) - i + 1, nsubs - new_as_B(i))
            end do
            new_r_C = binom_C(nsubs_t(1,2), n2rem)
            do i = 1, nsubs_t(1,2)
              if (n2rem - new_as_C(i) >= nsubs_t(1,2) - i + 1) &
                new_r_C = new_r_C - binom_C(nsubs_t(1,2) - i + 1, n2rem - new_as_C(i))
            end do

            new_r_A = (new_r_A - 1_int64) * ntc_B * ntc_C + (new_r_B - 1_int64) * ntc_C + new_r_C

            if (all(newconf_A(1:nsubs) == as_A(1:nsubs)) .and. &
                all(new_as_B(1:nsubs_t(1,1)) == as_B(1:nsubs_t(1,1))) .and. &
                all(new_as_C(1:nsubs_t(1,2)) == as_C(1:nsubs_t(1,2)))) then
              opc = opc + 1
              write (47, *) opc
              do i = 1, 3
                write (47, *) (mgroup(op, i, j), j=1, 3), vgroup(op, i)
              end do
            else if (.not. visited(new_r_A)) then
              visited(new_r_A) = .true.
              degen(indcount) = degen(indcount) + 1
              equivcount = equivcount + 1
            end if
          end do  ! op loop
          write (47, *) 0

          if ((100.0d0*equivcount/ntc_joint - perc > 5.0d0) .or. (equivcount == ntc_joint)) then
            perc = 100.0d0*equivcount/ntc_joint
            write (*, '(4x,i6,7x,f5.1,a2)') indcount, perc, "% "
          end if
        end if  ! not visited

        if (.not. mores_C) exit
!       Advance as_C to next nsubs_t(1,2)-subset of {1..n2rem}
        do i = nsubs_t(1,2), 1, -1
          if (as_C(i) < n2rem - nsubs_t(1,2) + i) then
            as_C(i) = as_C(i) + 1
            do j_remove = i+1, nsubs_t(1,2)
              as_C(j_remove) = as_C(j_remove-1) + 1
            end do
            mores_C = (as_C(1) /= n2rem - nsubs_t(1,2) + 1)
            exit
          end if
        end do
      end do  ! inner colouring loop (sp2)

      if (.not. mores_B) exit
!     Advance as_B to next nsubs_t(1,1)-subset of {1..nsubs}
      do i = nsubs_t(1,1), 1, -1
        if (as_B(i) < nsubs - nsubs_t(1,1) + i) then
          as_B(i) = as_B(i) + 1
          do j_remove = i+1, nsubs_t(1,1)
            as_B(j_remove) = as_B(j_remove-1) + 1
          end do
          mores_B = (as_B(1) /= nsubs - nsubs_t(1,1) + 1)
          exit
        end if
      end do
    end do  ! middle colouring loop (sp1)

    if (.not. mores_A) exit
!   Advance as_A to next nsubs-subset of {1..npos}
    do i = nsubs, 1, -1
      if (as_A(i) < npos - nsubs + i) then
        as_A(i) = as_A(i) + 1
        do j_remove = i+1, nsubs
          as_A(j_remove) = as_A(j_remove-1) + 1
        end do
        mores_A = (as_A(1) /= npos - nsubs + 1)
        exit
      end if
    end do
  end do  ! outer combined loop

  nic = indcount

! Trim arrays to actual size
  block
    integer, dimension(:,:), allocatable :: tmp2d
    integer, dimension(:), allocatable :: tmp1d
    allocate (tmp2d(1:nic, 1:nsubs))
    tmp2d = indconf(1:nic, 1:nsubs)
    call move_alloc(tmp2d, indconf)
    allocate (tmp1d(1:nic))
    tmp1d = degen(1:nic)
    call move_alloc(tmp1d, degen)
  end block

  write (*, *) " "
  write (*, '(a,i0)') "       Number of inequivalent configurations: ", nic
  write (*, *) " "

! Write OUTSOD header and data (absolute positions)
  write (30, *) nsubs_t(1,1), nsubs_t(1,2), nsubs_t(1,3), " substitutions in ", npos, "sites"
  write (30, *) nic, " configurations"
  do indcount = 1, nic
    write (30, 331) indcount, degen(indcount), indconf(indcount, 1:nsubs) + atini - 1
  end do

  deallocate (as_A, as_B, as_C, newconf_A, newconf_B, newconf_C)
  deallocate (new_as_B, new_as_C, as_notB_tmp, new_notB_tmp)
  deallocate (degen)
  deallocate (indconf)
  deallocate (binom_A, binom_B, binom_C)
  deallocate (visited)

  close (30)
  close (47)

  else if (ntarget >= 2 .and. all(nk(1:ntarget) == 1)) then  ! Stage D: multi-target binary

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Multi-target enumeration (ntarget >= 2): direct only
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  nsubs_tot   = sum(nsubs_t(1:ntarget, 1))
  nsubs_max_t = maxval(nsubs_t(1:ntarget, 1))

! Build directory name nXX[_YY[_ZZ...]] dynamically
  write (nxx_dir, '("n", i2.2)') nsubs_t(1,1)
  do t = 2, ntarget
    write (nxx_dir, '(a, "_", i2.2)') trim(nxx_dir), nsubs_t(t,1)
  end do
  call execute_command_line("mkdir -p " // trim(nxx_dir), wait=.true.)
  open (unit=30, file=trim(nxx_dir) // "/OUTSOD")
  write (30, '(a)') "# SOD OUTSOD format version 2"
  open (unit=47, file=trim(nxx_dir) // "/cSGO")

  write (*, *) "----------------------------------------------------------------------------"
  write (*, '(a)', advance='no') " Multi-target substitution: nsubs = ["
  do t = 1, ntarget
    if (t > 1) write (*, '(a)', advance='no') ", "
    write (*, '(i0)', advance='no') nsubs_t(t,1)
  end do
  write (*, '(a)') "]"
  write (*, *) "----------------------------------------------------------------------------"
  write (*, *) ""
  write (*, *) " > Composition of the substituted supercell:"
  do sp = 1, nsp
    found = .false.
    do t = 1, ntarget
      if (sp == sptarget(t)) then
        write (*, '(a, a, i10)') "    - ", newsymbol(t,1), nsubs_t(t,1)
        write (*, '(a, a, i10)') "    - ", newsymbol(t,2), natsp(sp) - nsubs_t(t,1)
        found = .true.
        exit
      end if
    end do
    if (.not. found) write (*, '(a, a, i10)') "    - ", symbol(sp), natsp(sp)
  end do
  write (*, *) ""

! Precompute 3D binomial table: binom_t(t, k, n) = C(n,k) for each target t
  npos_max = maxval(npos_t(1:ntarget))
  allocate (binom_t(1:ntarget, 0:nsubs_max_t, 0:npos_max))
  binom_t(:, 0, :) = 1_int64
  binom_t(:, 1:, :) = 0_int64
  do t = 1, ntarget
    do k = 1, nsubs_t(t,1)
      do j = k, npos_t(t)
        binom_t(t, k, j) = binom_t(t, k-1, j-1) + binom_t(t, k, j-1)
      end do
    end do
  end do

  do t = 1, ntarget
    ntc_t(t) = binom_t(t, nsubs_t(t,1), npos_t(t))
  end do
  ntc_joint = product(ntc_t(1:ntarget))

  do t = 1, ntarget
    write (*, '(a,i0,a,i0,a,i0,a,i0)') &
      "       Configurations for target ", t, ": C(", npos_t(t), ", ", nsubs_t(t,1), ") = ", ntc_t(t)
  end do
  write (*, '(a,i0)') "       Total joint configurations: ", ntc_joint
  maxentropy = kb * log(real(ntc_joint, real64))
  ientropy = 0.d0
  do t = 1, ntarget
    x = real(nsubs_t(t,1), real64) / real(npos_t(t), real64)
    if (x > 0.d0 .and. x < 1.d0) ientropy = ientropy &
      - kb * (real(npos_t(t), real64) / real(sum(npos_t(1:ntarget)), real64)) &
           * (x * log(x) + (1.d0 - x) * log(1.d0 - x))
  end do
  ientropy = ientropy * real(sum(npos_t(1:ntarget)), real64)
  if (ientropy > 0.d0) then
    perc = maxentropy / ientropy * 100.d0
  else
    perc = 0.d0
  end if
  write (*, '(a,f12.6,a)') "       Maximum entropy (supercell ensemble):    ", maxentropy * 1000.d0, " meV/K"
  write (*, '(a,f12.6,a)') "       Ideal mixing entropy (per cell):          ", ientropy * 1000.d0, " meV/K"
  write (*, '(a,f7.2,a)')  "       Supercell captures                        ", perc, "% of ideal entropy"
  write (*, *) " "

  allocate (visited(1:ntc_joint), stat=ierr)
  if (ierr /= 0) then
    write (*, *) "Error: insufficient memory for visited array (size", ntc_joint, ")"
    stop 1
  end if
  visited(:) = .false.

  allocate (as_t(1:ntarget, 1:nsubs_max_t), stat=ierr)
  allocate (newconf_t(1:ntarget, 1:nsubs_max_t), stat=ierr)
  allocate (tmp_conf(1:nsubs_tot), stat=ierr)
  allocate (indconf(1:ntc_joint, 1:nsubs_tot), stat=ierr)
  if (ierr /= 0) then
    write (*, *) "Error: insufficient memory for indconf"
    stop 1
  end if
  allocate (degen(1:ntc_joint), stat=ierr)
  indconf(:,:) = 0
  degen(:) = 1

  write (*, *) "Finding the inequivalent configurations..."
  write (*, *) " "
  write (*, *) "       Found    Completion "
  write (*, *) "       =====    ========== "

  indcount = 0
  perc = 0.0d0
  equivcount = 0

! Initialise all targets to first subset [1, 2, ..., nsubs_t(t,1)]
  do t = 1, ntarget
    do i = 1, nsubs_t(t,1)
      as_t(t, i) = i
    end do
    mores_t(t) = (nsubs_t(t,1) < npos_t(t))
  end do

  main_loop: do

!   Joint rank of current state (mixed-radix, 1-based)
    r = 0_int64
    do t = 1, ntarget
      r_t(t) = binom_t(t, nsubs_t(t,1), npos_t(t))
      do i = 1, nsubs_t(t,1)
        if (npos_t(t) - as_t(t,i) >= nsubs_t(t,1) - i + 1) &
          r_t(t) = r_t(t) - binom_t(t, nsubs_t(t,1) - i + 1, npos_t(t) - as_t(t,i))
      end do
      r = r * ntc_t(t) + (r_t(t) - 1_int64)
    end do
    r = r + 1_int64

    if (.not. visited(r)) then
      equivcount = equivcount + 1
      indcount = indcount + 1
      offset = 0
      do t = 1, ntarget
        indconf(indcount, offset+1:offset+nsubs_t(t,1)) = as_t(t, 1:nsubs_t(t,1))
        offset = offset + nsubs_t(t,1)
      end do
      visited(r) = .true.
      write (47, *) "List of operators for configuration: ", indcount
      opc = 1
      write (47, *) opc
      do i = 1, 3
        write (47, *) (mgroup(1, i, j), j=1, 3), vgroup(1, i)
      end do

!     Apply symmetry operators
      do op = 2, nop
!       Map each target's positions under op; sort result
        do t = 1, ntarget
          do concurrent (i = 1:nsubs_t(t,1))
            newconf_t(t, i) = eqmt(t, as_t(t, i), op)
          end do
          call bubble(newconf_t(t, 1:nsubs_t(t,1)), nsubs_t(t,1))
        end do
!       Joint rank of image (mixed-radix)
        new_r = 0_int64
        do t = 1, ntarget
          new_r_tt = binom_t(t, nsubs_t(t,1), npos_t(t))
          do i = 1, nsubs_t(t,1)
            if (npos_t(t) - newconf_t(t,i) >= nsubs_t(t,1) - i + 1) &
              new_r_tt = new_r_tt - binom_t(t, nsubs_t(t,1) - i + 1, npos_t(t) - newconf_t(t,i))
          end do
          new_r = new_r * ntc_t(t) + (new_r_tt - 1_int64)
        end do
        new_r = new_r + 1_int64
!       Stabilizer (same rank = same joint config) or new equivalent?
        if (new_r == r) then
          opc = opc + 1
          write (47, *) opc
          do i = 1, 3
            write (47, *) (mgroup(op, i, j), j=1, 3), vgroup(op, i)
          end do
        else if (.not. visited(new_r)) then
          visited(new_r) = .true.
          degen(indcount) = degen(indcount) + 1
          equivcount = equivcount + 1
        end if
      end do
      write (47, *) 0

      if ((100.0d0*equivcount/ntc_joint - perc > 5.0d0) .or. (equivcount == ntc_joint)) then
        perc = 100.0d0*equivcount/ntc_joint
        write (*, '(4x,i6,7x,f5.1,a2)') indcount, perc, "% "
      end if
    end if

!   Advance to next joint configuration (carry propagation, innermost level first)
    do t = ntarget, 1, -1
      if (mores_t(t)) then
!       Advance level t (inline next-k-subset)
        do i = nsubs_t(t,1), 1, -1
          if (as_t(t, i) < npos_t(t) - nsubs_t(t,1) + i) then
            as_t(t, i) = as_t(t, i) + 1
            do j = i+1, nsubs_t(t,1)
              as_t(t, j) = as_t(t, j-1) + 1
            end do
            mores_t(t) = (as_t(t, 1) /= npos_t(t) - nsubs_t(t,1) + 1)
            exit
          end if
        end do
!       Reset all inner levels to initial subsets
        do tt = t+1, ntarget
          do i = 1, nsubs_t(tt,1)
            as_t(tt, i) = i
          end do
          mores_t(tt) = (nsubs_t(tt,1) < npos_t(tt))
        end do
        exit
      else if (t == 1) then
        exit main_loop  ! all joint configurations exhausted
      end if
    end do

  end do main_loop

  nic = indcount

! Trim indconf and degen to actual size
  block
    integer, dimension(:,:), allocatable :: tmp2d
    integer, dimension(:), allocatable :: tmp1d
    allocate (tmp2d(1:nic, 1:nsubs_tot))
    tmp2d = indconf(1:nic, 1:nsubs_tot)
    call move_alloc(tmp2d, indconf)
    allocate (tmp1d(1:nic))
    tmp1d = degen(1:nic)
    call move_alloc(tmp1d, degen)
  end block

  write (*, *) " "
  write (*, '(a,i0)') "       Number of inequivalent configurations: ", nic
  write (*, *) " "

  write (30, *) (nsubs_t(t,1), t=1,ntarget), " substitutions in ", (npos_t(t), t=1,ntarget), "sites"
  write (30, *) nic, " configurations"
  do indcount = 1, nic
    offset = 0
    do t = 1, ntarget
      tmp_conf(offset+1:offset+nsubs_t(t,1)) = &
        indconf(indcount, offset+1:offset+nsubs_t(t,1)) + atini_t(t) - 1
      offset = offset + nsubs_t(t,1)
    end do
    write (30, 331) indcount, degen(indcount), tmp_conf(1:nsubs_tot)
  end do
331 format(i6, 1x, i6, *(1x, i4))

  deallocate (newconf_t, as_t, tmp_conf)
  deallocate (degen)
  deallocate (indconf)
  deallocate (binom_t)
  deallocate (visited)

  close (30)
  close (47)

  else  ! Phase 4: ntarget >= 2, at least one nk(t) >= 2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Multi-target multi-nary enumeration (Phase 4)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! Per-target derived counts
  do t = 1, ntarget
    nsubs_tot_t(t) = sum(nsubs_t(t, 1:nk(t)))
    nrem2_t(t) = nsubs_tot_t(t) - nsubs_t(t,1)
  end do
  nsubs_tot    = sum(nsubs_tot_t(1:ntarget))
  nsubs_pos_max = maxval(nsubs_tot_t(1:ntarget))
  npos_max      = maxval(npos_t(1:ntarget))

! Build directory name: nAA[_BB[_CC]]_DD[_EE[_FF]]...
  write (nxx_dir, '("n", i2.2)') nsubs_t(1,1)
  do j = 2, nk(1)
    write (nxx_dir, '(a,"_",i2.2)') trim(nxx_dir), nsubs_t(1,j)
  end do
  do t = 2, ntarget
    do j = 1, nk(t)
      write (nxx_dir, '(a,"_",i2.2)') trim(nxx_dir), nsubs_t(t,j)
    end do
  end do
  call execute_command_line("mkdir -p " // trim(nxx_dir), wait=.true.)
  open (unit=30, file=trim(nxx_dir) // "/OUTSOD")
  write (30, '(a)') "# SOD OUTSOD format version 2"
  open (unit=47, file=trim(nxx_dir) // "/cSGO")

  write (*, *) "----------------------------------------------------------------------------"
  write (*, '(a)', advance='no') " Multi-target multi-nary substitution: nsubs = ["
  do t = 1, ntarget
    if (t > 1) write (*, '(a)', advance='no') " | "
    do j = 1, nk(t)
      if (j > 1) write (*, '(a)', advance='no') ","
      write (*, '(i0)', advance='no') nsubs_t(t,j)
    end do
  end do
  write (*, '(a)') "]"
  write (*, *) "----------------------------------------------------------------------------"
  write (*, *) ""
  write (*, *) " > Composition of the substituted supercell:"
  do sp = 1, nsp
    found = .false.
    do t = 1, ntarget
      if (sp == sptarget(t)) then
        do j = 1, nk(t)
          write (*, '(a, a, i10)') "    - ", newsymbol(t,j), nsubs_t(t,j)
        end do
        write (*, '(a, a, i10)') "    - ", newsymbol(t, nk(t)+1), natsp(sp) - nsubs_tot_t(t)
        found = .true.
        exit
      end if
    end do
    if (.not. found) write (*, '(a, a, i10)') "    - ", symbol(sp), natsp(sp)
  end do
  write (*, *) ""

! Binomial tables:
!   binom_pos_t(t, k, n) = C(n,k) for position ranks: k=0..nsubs_tot_t(t), n=0..npos_t(t)
!   binom_col_t(t, k, n) = C(n,k) for colouring ranks: k,n in 0..nsubs_pos_max
  allocate (binom_pos_t(1:ntarget, 0:nsubs_pos_max, 0:npos_max))
  allocate (binom_col_t(1:ntarget, 0:nsubs_pos_max, 0:nsubs_pos_max))
  binom_pos_t(:, 0, :) = 1_int64;  binom_pos_t(:, 1:, :) = 0_int64
  binom_col_t(:, 0, :) = 1_int64;  binom_col_t(:, 1:, :) = 0_int64
  do t = 1, ntarget
    do k = 1, nsubs_tot_t(t)
      do j = k, npos_t(t)
        binom_pos_t(t, k, j) = binom_pos_t(t, k-1, j-1) + binom_pos_t(t, k, j-1)
      end do
    end do
    do k = 1, nsubs_pos_max
      do j = k, nsubs_pos_max
        binom_col_t(t, k, j) = binom_col_t(t, k-1, j-1) + binom_col_t(t, k, j-1)
      end do
    end do
  end do

! Per-target space sizes
  do t = 1, ntarget
    ntc_pos_t(t) = binom_pos_t(t, nsubs_tot_t(t), npos_t(t))
    if (nk(t) >= 2) then
      ntc_col1_t(t) = binom_col_t(t, nsubs_t(t,1), nsubs_tot_t(t))
    else
      ntc_col1_t(t) = 1_int64
    end if
    if (nk(t) >= 3) then
      ntc_col2_t(t) = binom_col_t(t, nsubs_t(t,2), nrem2_t(t))
    else
      ntc_col2_t(t) = 1_int64
    end if
    ntc_t(t) = ntc_pos_t(t) * ntc_col1_t(t) * ntc_col2_t(t)
  end do
  ntc_joint = product(ntc_t(1:ntarget))

  do t = 1, ntarget
    write (*, '(a,i0,a,i0,a,i0,a,i0)') &
      "       Configurations for target ", t, &
      ": C(", npos_t(t), ",", nsubs_tot_t(t), ") = ", ntc_pos_t(t)
    if (nk(t) >= 2) write (*, '(a,i0,a,i0,a,i0)') &
      "         col1 space: C(", nsubs_tot_t(t), ",", nsubs_t(t,1), ") = ", ntc_col1_t(t)
    if (nk(t) >= 3) write (*, '(a,i0,a,i0,a,i0)') &
      "         col2 space: C(", nrem2_t(t), ",", nsubs_t(t,2), ") = ", ntc_col2_t(t)
    write (*, '(a,i0,a,i0)') "         ntc(t=", t, ") = ", ntc_t(t)
  end do
  write (*, '(a,i0)') "       Total joint configurations: ", ntc_joint
  write (*, *) ""

  allocate (visited(1:ntc_joint), stat=ierr)
  if (ierr /= 0) then
    write (*, *) "Error: insufficient memory for visited array (size", ntc_joint, ")"
    stop 1
  end if
  visited(:) = .false.

  allocate (as_pos_t(1:ntarget,      1:nsubs_pos_max))
  allocate (as_col1_t(1:ntarget,     1:nsubs_pos_max))
  allocate (as_col2_t(1:ntarget,     1:nsubs_pos_max))
  allocate (notcol1_t(1:ntarget,     1:nsubs_pos_max))
  allocate (newpos_t(1:ntarget,      1:nsubs_pos_max))
  allocate (newval_col1_t(1:ntarget, 1:nsubs_pos_max))
  allocate (new_col1_t(1:ntarget,    1:nsubs_pos_max))
  allocate (new_notcol1_t(1:ntarget, 1:nsubs_pos_max))
  allocate (newval_col2_t(1:ntarget, 1:nsubs_pos_max))
  allocate (new_col2_t(1:ntarget,    1:nsubs_pos_max))
  allocate (tmp_conf(1:nsubs_tot))
  allocate (indconf(1:ntc_joint, 1:nsubs_tot), stat=ierr)
  if (ierr /= 0) then
    write (*, *) "Error: insufficient memory for indconf"
    stop 1
  end if
  allocate (degen(1:ntc_joint))
  indconf(:,:) = 0
  degen(:) = 1

! Initialise all targets to first position subset and colouring
  do t = 1, ntarget
    do i = 1, nsubs_tot_t(t)
      as_pos_t(t, i) = i
    end do
    mores_pos_t(t) = (nsubs_tot_t(t) < npos_t(t))
    if (nk(t) >= 2) then
      do i = 1, nsubs_t(t,1)
        as_col1_t(t, i) = i
      end do
      mores_col1_t(t) = (nsubs_t(t,1) < nsubs_tot_t(t))
      call complement_indices(as_col1_t(t,:), nsubs_t(t,1), nsubs_tot_t(t), notcol1_t(t,:), j)
    end if
    if (nk(t) >= 3) then
      do i = 1, nsubs_t(t,2)
        as_col2_t(t, i) = i
      end do
      mores_col2_t(t) = (nsubs_t(t,2) < nrem2_t(t))
    end if
  end do

  write (*, *) "Finding the inequivalent configurations..."
  write (*, *) " "
  write (*, *) "       Found    Completion "
  write (*, *) "       =====    ========== "
  indcount = 0
  perc     = 0.0d0
  equivcount = 0

  p4_main: do

!   --- Compute joint rank (mixed-radix, 1-based) ---
    r = 0_int64
    do t = 1, ntarget
      r_pos = binom_pos_t(t, nsubs_tot_t(t), npos_t(t))
      do i = 1, nsubs_tot_t(t)
        if (npos_t(t) - as_pos_t(t,i) >= nsubs_tot_t(t) - i + 1) &
          r_pos = r_pos - binom_pos_t(t, nsubs_tot_t(t)-i+1, npos_t(t)-as_pos_t(t,i))
      end do
      if (nk(t) >= 2) then
        r_col1 = binom_col_t(t, nsubs_t(t,1), nsubs_tot_t(t))
        do i = 1, nsubs_t(t,1)
          if (nsubs_tot_t(t) - as_col1_t(t,i) >= nsubs_t(t,1) - i + 1) &
            r_col1 = r_col1 - binom_col_t(t, nsubs_t(t,1)-i+1, nsubs_tot_t(t)-as_col1_t(t,i))
        end do
      else
        r_col1 = 1_int64
      end if
      if (nk(t) >= 3) then
        r_col2 = binom_col_t(t, nsubs_t(t,2), nrem2_t(t))
        do i = 1, nsubs_t(t,2)
          if (nrem2_t(t) - as_col2_t(t,i) >= nsubs_t(t,2) - i + 1) &
            r_col2 = r_col2 - binom_col_t(t, nsubs_t(t,2)-i+1, nrem2_t(t)-as_col2_t(t,i))
        end do
      else
        r_col2 = 1_int64
      end if
      r_t(t) = (r_pos - 1_int64)*ntc_col1_t(t)*ntc_col2_t(t) &
             + (r_col1 - 1_int64)*ntc_col2_t(t) + r_col2
      r = r * ntc_t(t) + (r_t(t) - 1_int64)
    end do
    r = r + 1_int64

    if (.not. visited(r)) then
      equivcount = equivcount + 1
      indcount   = indcount + 1
      visited(r) = .true.

!     Build indconf row
      offset = 0
      do t = 1, ntarget
        if (nk(t) == 1) then
          do i = 1, nsubs_tot_t(t)
            indconf(indcount, offset+i) = as_pos_t(t, i)
          end do
        else if (nk(t) == 2) then
          do i = 1, nsubs_t(t,1)
            indconf(indcount, offset+i) = as_pos_t(t, as_col1_t(t,i))
          end do
          do i = 1, nrem2_t(t)
            indconf(indcount, offset+nsubs_t(t,1)+i) = as_pos_t(t, notcol1_t(t,i))
          end do
        else  ! nk(t) == 3
          do i = 1, nsubs_t(t,1)
            indconf(indcount, offset+i) = as_pos_t(t, as_col1_t(t,i))
          end do
          do i = 1, nsubs_t(t,2)
            indconf(indcount, offset+nsubs_t(t,1)+i) = as_pos_t(t, notcol1_t(t, as_col2_t(t,i)))
          end do
          k = offset + nsubs_t(t,1) + nsubs_t(t,2) + 1
          do i = 1, nrem2_t(t)
            found = .false.
            do j = 1, nsubs_t(t,2)
              if (as_col2_t(t,j) == i) then; found = .true.; exit; end if
            end do
            if (.not. found) then
              indconf(indcount, k) = as_pos_t(t, notcol1_t(t,i))
              k = k + 1
            end if
          end do
        end if
        offset = offset + nsubs_tot_t(t)
      end do

      write (47, *) "List of operators for configuration: ", indcount
      opc = 1
      write (47, *) opc
      do i = 1, 3
        write (47, *) (mgroup(1, i, j), j=1, 3), vgroup(1, i)
      end do

!     Apply symmetry operators
      do op = 2, nop
        do t = 1, ntarget
!         Image of combined position set
          do i = 1, nsubs_tot_t(t)
            newpos_t(t, i) = eqmt(t, as_pos_t(t,i), op)
          end do
          call bubble(newpos_t(t, 1:nsubs_tot_t(t)), nsubs_tot_t(t))

          if (nk(t) >= 2) then
!           Image of sp1 positions (absolute values)
            do i = 1, nsubs_t(t,1)
              newval_col1_t(t, i) = eqmt(t, as_pos_t(t, as_col1_t(t,i)), op)
            end do
            call bubble(newval_col1_t(t, 1:nsubs_t(t,1)), nsubs_t(t,1))
!           Find indices of sp1 values in newpos_t (two-pointer)
            call find_in_sorted(newval_col1_t(t, 1:nsubs_t(t,1)), nsubs_t(t,1), &
                                newpos_t(t, 1:nsubs_tot_t(t)), new_col1_t(t, :))
!           Build new_notcol1_t: values in newpos_t not at new_col1_t indices
            jj = 0; j = 1
            do i = 1, nsubs_tot_t(t)
              if (j <= nsubs_t(t,1) .and. new_col1_t(t,j) == i) then
                j = j + 1
              else
                jj = jj + 1
                new_notcol1_t(t, jj) = newpos_t(t, i)
              end if
            end do
          end if

          if (nk(t) >= 3) then
!           Image of sp2 positions
            do i = 1, nsubs_t(t,2)
              newval_col2_t(t, i) = eqmt(t, as_pos_t(t, notcol1_t(t, as_col2_t(t,i))), op)
            end do
            call bubble(newval_col2_t(t, 1:nsubs_t(t,2)), nsubs_t(t,2))
!           Find indices of sp2 values in new_notcol1_t (two-pointer)
            call find_in_sorted(newval_col2_t(t, 1:nsubs_t(t,2)), nsubs_t(t,2), &
                                new_notcol1_t(t, 1:nrem2_t(t)), new_col2_t(t, :))
          end if
        end do  ! t loop

!       Compute joint rank of image
        new_r = 0_int64
        do t = 1, ntarget
          new_r_pos = binom_pos_t(t, nsubs_tot_t(t), npos_t(t))
          do i = 1, nsubs_tot_t(t)
            if (npos_t(t) - newpos_t(t,i) >= nsubs_tot_t(t) - i + 1) &
              new_r_pos = new_r_pos - binom_pos_t(t, nsubs_tot_t(t)-i+1, npos_t(t)-newpos_t(t,i))
          end do
          if (nk(t) >= 2) then
            new_r_col1 = binom_col_t(t, nsubs_t(t,1), nsubs_tot_t(t))
            do i = 1, nsubs_t(t,1)
              if (nsubs_tot_t(t) - new_col1_t(t,i) >= nsubs_t(t,1) - i + 1) &
                new_r_col1 = new_r_col1 - binom_col_t(t, nsubs_t(t,1)-i+1, nsubs_tot_t(t)-new_col1_t(t,i))
            end do
          else
            new_r_col1 = 1_int64
          end if
          if (nk(t) >= 3) then
            new_r_col2 = binom_col_t(t, nsubs_t(t,2), nrem2_t(t))
            do i = 1, nsubs_t(t,2)
              if (nrem2_t(t) - new_col2_t(t,i) >= nsubs_t(t,2) - i + 1) &
                new_r_col2 = new_r_col2 - binom_col_t(t, nsubs_t(t,2)-i+1, nrem2_t(t)-new_col2_t(t,i))
            end do
          else
            new_r_col2 = 1_int64
          end if
          new_r_tt = (new_r_pos - 1_int64)*ntc_col1_t(t)*ntc_col2_t(t) &
                   + (new_r_col1 - 1_int64)*ntc_col2_t(t) + new_r_col2
          new_r = new_r * ntc_t(t) + (new_r_tt - 1_int64)
        end do
        new_r = new_r + 1_int64

        if (new_r == r) then
          opc = opc + 1
          write (47, *) opc
          do i = 1, 3
            write (47, *) (mgroup(op, i, j), j=1, 3), vgroup(op, i)
          end do
        else if (.not. visited(new_r)) then
          visited(new_r) = .true.
          degen(indcount) = degen(indcount) + 1
          equivcount = equivcount + 1
        end if
      end do  ! op loop
      write (47, *) 0

      if ((100.0d0*equivcount/ntc_joint - perc > 5.0d0) .or. (equivcount == ntc_joint)) then
        perc = 100.0d0*equivcount/ntc_joint
        write (*, '(4x,i6,7x,f5.1,a2)') indcount, perc, "% "
      end if
    end if  ! .not. visited(r)

!   --- Carry propagation (innermost = col2 of tN, outermost = pos of t1) ---
    p4_carry: do t = ntarget, 1, -1
      if (nk(t) >= 3) then
        if (mores_col2_t(t)) then
          call comb_next(as_col2_t(t,:), nsubs_t(t,2), nrem2_t(t), mores_col2_t(t))
          exit p4_carry
        else
          do i = 1, nsubs_t(t,2); as_col2_t(t,i) = i; end do
          mores_col2_t(t) = (nsubs_t(t,2) < nrem2_t(t))
        end if
      end if
      if (nk(t) >= 2) then
        if (mores_col1_t(t)) then
          call comb_next(as_col1_t(t,:), nsubs_t(t,1), nsubs_tot_t(t), mores_col1_t(t))
          call complement_indices(as_col1_t(t,:), nsubs_t(t,1), nsubs_tot_t(t), notcol1_t(t,:), j)
          if (nk(t) >= 3) then
            do i = 1, nsubs_t(t,2); as_col2_t(t,i) = i; end do
            mores_col2_t(t) = (nsubs_t(t,2) < nrem2_t(t))
          end if
          exit p4_carry
        else
          do i = 1, nsubs_t(t,1); as_col1_t(t,i) = i; end do
          mores_col1_t(t) = (nsubs_t(t,1) < nsubs_tot_t(t))
          call complement_indices(as_col1_t(t,:), nsubs_t(t,1), nsubs_tot_t(t), notcol1_t(t,:), j)
          if (nk(t) >= 3) then
            do i = 1, nsubs_t(t,2); as_col2_t(t,i) = i; end do
            mores_col2_t(t) = (nsubs_t(t,2) < nrem2_t(t))
          end if
        end if
      end if
      if (mores_pos_t(t)) then
        call comb_next(as_pos_t(t,:), nsubs_tot_t(t), npos_t(t), mores_pos_t(t))
        if (nk(t) >= 2) then
          do i = 1, nsubs_t(t,1); as_col1_t(t,i) = i; end do
          mores_col1_t(t) = (nsubs_t(t,1) < nsubs_tot_t(t))
          call complement_indices(as_col1_t(t,:), nsubs_t(t,1), nsubs_tot_t(t), notcol1_t(t,:), j)
        end if
        if (nk(t) >= 3) then
          do i = 1, nsubs_t(t,2); as_col2_t(t,i) = i; end do
          mores_col2_t(t) = (nsubs_t(t,2) < nrem2_t(t))
        end if
        exit p4_carry
      else if (t == 1) then
        exit p4_main
      else
        do i = 1, nsubs_tot_t(t); as_pos_t(t,i) = i; end do
        mores_pos_t(t) = (nsubs_tot_t(t) < npos_t(t))
        if (nk(t) >= 2) then
          do i = 1, nsubs_t(t,1); as_col1_t(t,i) = i; end do
          mores_col1_t(t) = (nsubs_t(t,1) < nsubs_tot_t(t))
          call complement_indices(as_col1_t(t,:), nsubs_t(t,1), nsubs_tot_t(t), notcol1_t(t,:), j)
        end if
        if (nk(t) >= 3) then
          do i = 1, nsubs_t(t,2); as_col2_t(t,i) = i; end do
          mores_col2_t(t) = (nsubs_t(t,2) < nrem2_t(t))
        end if
      end if
    end do p4_carry

  end do p4_main

  nic = indcount

  block
    integer, dimension(:,:), allocatable :: tmp2d
    integer, dimension(:), allocatable   :: tmp1d
    allocate (tmp2d(1:nic, 1:nsubs_tot))
    tmp2d = indconf(1:nic, 1:nsubs_tot)
    call move_alloc(tmp2d, indconf)
    allocate (tmp1d(1:nic))
    tmp1d = degen(1:nic)
    call move_alloc(tmp1d, degen)
  end block

  write (*, *) " "
  write (*, '(a,i0)') "       Number of inequivalent configurations: ", nic
  write (*, *) " "

! OUTSOD: per-target nsubs counts, then per-target npos
  write (30, *) ((nsubs_t(t,j), j=1,nk(t)), t=1,ntarget), &
                " substitutions in ", (npos_t(t), t=1,ntarget), "sites"
  write (30, *) nic, " configurations"
  do indcount = 1, nic
    offset = 0
    do t = 1, ntarget
      tmp_conf(offset+1:offset+nsubs_tot_t(t)) = &
        indconf(indcount, offset+1:offset+nsubs_tot_t(t)) + atini_t(t) - 1
      offset = offset + nsubs_tot_t(t)
    end do
    write (30, 331) indcount, degen(indcount), tmp_conf(1:nsubs_tot)
  end do

  deallocate (binom_pos_t, binom_col_t)
  deallocate (as_pos_t, as_col1_t, as_col2_t, notcol1_t)
  deallocate (newpos_t, newval_col1_t, new_col1_t, new_notcol1_t, newval_col2_t, new_col2_t)
  deallocate (tmp_conf, degen, indconf, visited)

  close (30)
  close (47)

  end if  ! ntarget == 1 / ntarget >= 2

!!!!!!! Write FILER to file, to be read by the shell script
  write (43, *) filer
! Calculation input file generation is handled by genersod,
! called automatically by sod_comb.sh when FILER is not -1.

  deallocate (eqmt)

!!!!!!!Reporting the end
  call system_clock(t_now)
  elapsed_total = real(t_now - t_start_total, real64) / real(clock_rate, real64)
  write (*, *) " "
  if (ntarget == 1 .and. nk(1) == 1) then
  write (*, *) " > Timing summary:"
  write (*, *) "   nsubs    Wall time (s)"
  write (*, *) "   -----    -------------"
  do nsubs_loop = 0, nsubs_max - nsubs_min
    write (*, '(3x, i5, 4x, f12.2)') nsubs_level_list(nsubs_loop), elapsed_per_level(nsubs_loop)
  end do
  write (*, *) "   -----    -------------"
  write (*, '(a, f12.2)') "   Total         ", elapsed_total
  write (*, *) ""
  deallocate (elapsed_per_level, nsubs_level_list)
  end if
  write (*, *) " > Done!"
  write (*, *) ""

contains

  ! ---------------------------------------------------------------------------
  ! Helper: combinatorial rank of sorted subset as(1:nchoose) ⊆ {1..nelem}.
  ! Uses precomputed binomial table btab(0:nchoose, 0:nelem).
  ! ---------------------------------------------------------------------------
  pure integer(int64) function comb_rank(as, nchoose, nelem, btab)
    integer,        intent(in) :: as(:), nchoose, nelem
    integer(int64), intent(in) :: btab(0:, 0:)
    integer :: i
    comb_rank = btab(nchoose, nelem)
    do i = 1, nchoose
      if (nelem - as(i) >= nchoose - i + 1) &
        comb_rank = comb_rank - btab(nchoose - i + 1, nelem - as(i))
    end do
  end function comb_rank

  ! ---------------------------------------------------------------------------
  ! Helper: advance sorted subset as(1:nchoose) to next nchoose-subset of {1..nelem}.
  ! Updates mores: .false. when the new value is the last subset.
  ! Must only be called when mores is currently .true.
  ! ---------------------------------------------------------------------------
  subroutine comb_next(as, nchoose, nelem, mores)
    integer, intent(inout) :: as(:)
    integer, intent(in)    :: nchoose, nelem
    logical, intent(inout) :: mores
    integer :: i, j
    do i = nchoose, 1, -1
      if (as(i) < nelem - nchoose + i) then
        as(i) = as(i) + 1
        do j = i + 1, nchoose
          as(j) = as(j-1) + 1
        end do
        mores = (as(1) /= nelem - nchoose + 1)
        return
      end if
    end do
  end subroutine comb_next

  ! ---------------------------------------------------------------------------
  ! Helper: given sorted as(1:k) ⊆ {1..n}, return sorted complement in notAs.
  ! m returns n-k (number of complement elements written).
  ! ---------------------------------------------------------------------------
  subroutine complement_indices(as, k, n, notAs, m)
    integer, intent(in)  :: as(:), k, n
    integer, intent(out) :: notAs(:), m
    integer :: i, j
    m = 0
    j = 1
    do i = 1, n
      if (j <= k .and. as(j) == i) then
        j = j + 1
      else
        m = m + 1
        notAs(m) = i
      end if
    end do
  end subroutine complement_indices

  ! ---------------------------------------------------------------------------
  ! Helper: for each vals(i) (sorted), find its 1-based index in sorted ref,
  ! store in idx_out(i). Output is automatically sorted (two-pointer).
  ! ---------------------------------------------------------------------------
  subroutine find_in_sorted(vals, nvals, ref, idx_out)
    integer, intent(in)  :: vals(:), nvals, ref(:)
    integer, intent(out) :: idx_out(:)
    integer :: i, j
    j = 1
    do i = 1, nvals
      do while (ref(j) /= vals(i))
        j = j + 1
      end do
      idx_out(i) = j
      j = j + 1
    end do
  end subroutine find_in_sorted

end program combsod
