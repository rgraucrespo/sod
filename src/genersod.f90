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

program genersod
  implicit none

  integer, parameter :: nspmax = 10, natmax = 10000, nlineamax = 200
  integer, parameter :: maxgulplines = 1000

  integer :: i, l, m
  integer :: ifound
  integer :: sp, nsp, sptarget, ssp
  integer :: at0, nat0, at, nat, att, filer, ndigits
  character :: outfilename*20, fmtstr*20
  integer :: na, nb, nc, nsubs, nsubs_min, nsubs_max, atini, atfin, cumnatsp, ierr_rd
  character :: insod_line*256
  integer :: npos, nic, indcount
  integer, dimension(:), allocatable :: newconf, degen
  integer, dimension(nspmax) :: natsp0, natsp, snatsp
  integer, dimension(natmax) :: spat
  real, dimension(natmax, 3) :: coords0, coords
  real, dimension(3, 3) :: cellvector
  integer, dimension(:, :), allocatable :: indconf
  real :: a1, b1, c1, alpha, beta, gamma, a, b, c, xc, yc, zc
  character, dimension(nspmax) :: symbol*3, ssymbol*3
  character, dimension(2) :: newsymbol*3
  character :: linea*85, title*15, runtitle*40, trashtext*20, symboltrash*3
  character :: cifline*200, atmlabel*20, atmsymbol*3
  logical :: in_atom_loop

  ! Variables for GULP logic (case 1)
  integer :: ios, ngulplines, struc_line_idx, nmap_gulp, all_nsp_gulp, idx_token
  character :: gulpline_buf*200, sline_trimmed*200, outbuf_gulp*200
  character, dimension(maxgulplines) :: gulptemplate*200
  character, dimension(nspmax) :: map_sod_arr*3, map_gulp_arr*3
  character, dimension(nspmax+1) :: all_sym*3, all_gulptype_arr*3
  integer :: all_ishell_arr(nspmax+1)
  character :: libname_gulp*80, gtype_buf*3, gline_type*3, gline_coreshel*4, pending_section_header*200
  logical :: haslibrary_gulp, hasinline_ff, relevant_line
  character :: numfmt*20, numstr*10, confdir_gulp*30, inpfile_gulp*60, ndir_gulp*10

  ! Variables for LAMMPS logic (case 2)
  character :: lammps_atom_style*20
  integer :: lammps_nmap, lmptypeval, lmpbondval, lmpidx
  character, dimension(nspmax+1) :: lammps_map_sod*3, lammps_map_role*5
  integer :: lammps_map_type(nspmax+1), lammps_map_bond(nspmax+1)
  integer :: lammps_core_type(nspmax+1), lammps_shell_type(nspmax+1), lammps_shell_bond(nspmax+1)
  integer :: natoms_lammps, nbonds_lammps, natom_types_lammps, nbond_types_lammps
  integer :: mol_id_lammps, atom_id_lammps, bond_id_lammps
  real :: xy_tilt, xz_tilt, yz_tilt
  logical :: has_shells_lammps, is_triclinic_lammps
  character :: lmptok1_buf*20, lmptok2_buf*10, lmpslinetail*200, lmptmpline*200
  integer, dimension(natmax) :: lmp_core_id, lmp_shell_id, lmp_mol_id, lmp_sym_idx
  integer :: lammps_sym_idx, lmp_i
  character :: lammps_sym_cur*3

! Input files

  open (unit=9, file="INSOD")
  open (unit=31, file="supercell.cif")

! Output files

  open (unit=43, file="filer")

!
! DEFINITION OF VARIABLES:
!
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

  write (*, *) "Reading INSOD, supercell.cif, and OUTSOD to generate calculation input files"
  write (*, *) " "

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
  read (9, '(A)') insod_line
  read (insod_line, *, IOSTAT=ierr_rd) nsubs_min, nsubs_max
  if (ierr_rd /= 0) then
    read (insod_line, *) nsubs_min
    nsubs_max = nsubs_min
  end if
  if (nsubs_max == 0) then
    write (*, *) "Illegal number of substitutions"
    stop
  end if
  read (9, *)
  read (9, *)
  read (9, *)
  read (9, *) (newsymbol(i), i=1, 2)
  read (9, *)
  read (9, *)
  read (9, *) filer

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      [OUTSOD is now read per level inside the main loop below]
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading supercell.cif; builds natsp and spat from atom loop
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  natsp(1:nsp) = 0
  nat = 0
  in_atom_loop = .false.

  do
    read (31, '(a)', end=801) cifline
    cifline = adjustl(cifline)
    if (cifline(1:14) == '_cell_length_a') then
      read (cifline(15:), *) a
    else if (cifline(1:14) == '_cell_length_b') then
      read (cifline(15:), *) b
    else if (cifline(1:14) == '_cell_length_c') then
      read (cifline(15:), *) c
    else if (cifline(1:17) == '_cell_angle_alpha') then
      read (cifline(18:), *) alpha
    else if (cifline(1:16) == '_cell_angle_beta') then
      read (cifline(17:), *) beta
    else if (cifline(1:17) == '_cell_angle_gamma') then
      read (cifline(18:), *) gamma
    else if (cifline(1:18) == '_atom_site_fract_z') then
      in_atom_loop = .true.
    else if (in_atom_loop .and. len_trim(cifline) > 0 .and. &
             cifline(1:1) /= '_' .and. cifline(1:1) /= '#' .and. &
             cifline(1:4) /= 'loop' .and. cifline(1:4) /= 'data' .and. &
             cifline(1:1) /= "'") then
      nat = nat + 1
      read (cifline, *) atmlabel, atmsymbol, coords(nat, 1), coords(nat, 2), coords(nat, 3)
      do sp = 1, nsp
        if (trim(atmsymbol) == trim(symbol(sp))) then
          natsp(sp) = natsp(sp) + 1
          spat(nat) = sp
          exit
        end if
      end do
    end if
  end do
801 continue

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!    GENERATE INPUT FILES !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!! Write a temporary file with FILER, to be read later by the script
  write (43, *) filer

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Loop over substitution levels
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  select case (filer)
  case (-1)
    write (*, *) "Calculation files not created.&
               & Change filer value in INSOD if you want to create calculation files."
    write (*, *) ""
  case (0)
    write (*, *) " "
    write (*, *) "Creating CIF files for each configuration..."
    write (*, *) " "
  case (1)
    write (*, *) " "
    write (*, *) "Creating input files for GULP in configuration directories..."
    write (*, *) " "
  case (2)
    write (*, *) " "
    write (*, *) "Creating input files for LAMMPS in configuration directories..."
    write (*, *) " "
  case (11)
    write (*, *) " "
    write (*, *) "Creating input files for VASP in configuration directories..."
    write (*, *) " "
  case (12)
    write (*, *) " "
    write (*, *) "Creating input files for CASTEP in configuration directories..."
    write (*, *) " "
  case (13)
    write (*, *) " "
    write (*, *) "Creating input files for Quantum ESPRESSO in configuration directories..."
    write (*, *) " "
  end select

  do nsubs = nsubs_min, nsubs_max

    write (ndir_gulp, '("n", i2.2)') nsubs
    open (unit=30, file=trim(ndir_gulp) // '/OUTSOD', status='old', IOSTAT=ios)
    if (ios /= 0) then
      write (*, *) "Warning: ", trim(ndir_gulp), "/OUTSOD not found, skipping."
      cycle
    end if

    read (30, *) m, trashtext, trashtext, npos
    read (30, *) nic

    ndigits = max(1, int(log10(real(nic))) + 1)
    write (fmtstr, '(a,i0,a,i0,a)') '(a,i', ndigits, '.', ndigits, ')'

    allocate (degen(1:nic))
    allocate (newconf(1:nsubs))
    allocate (indconf(1:nic, 1:nsubs))

    do indcount = 1, nic
      read (30, *) m, degen(indcount), indconf(indcount, 1:nsubs)
      if (m /= indcount) then
        write (*, *) "Error in configuration numbering in OUTSOD. Aborting..."
        stop
      end if
    end do

    close (30)

    select case (filer)

  case (-1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!  CIF (P1, one file per config) !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case (0)

    write (numfmt, '(a,i0,a,i0,a)') '(i', ndigits, '.', ndigits, ')'

    write (ndir_gulp, '(a,i2.2)') 'n', nsubs
    call execute_command_line('mkdir -p ' // trim(ndir_gulp), exitstat=ios)
    if (ios /= 0) then
      write (*, *) "Error: could not create directory ", trim(ndir_gulp)
      stop
    end if

    do indcount = 1, nic
      newconf(1:nsubs) = indconf(indcount, 1:nsubs)

      write (numstr, numfmt) indcount
      confdir_gulp = trim(ndir_gulp) // '/c' // trim(numstr)
      call execute_command_line('mkdir -p ' // trim(confdir_gulp), exitstat=ios)
      if (ios /= 0) then
        write (*, *) "Error: could not create directory ", trim(confdir_gulp)
        stop
      end if

      inpfile_gulp = trim(confdir_gulp) // '/configuration.cif'
      open (unit=72, file=trim(inpfile_gulp), status='replace')

      write (72, '(a, i0)') "data_c", indcount
      write (72, '(a)') " "
      write (72, '(a, f12.6)') "_cell_length_a     ", a
      write (72, '(a, f12.6)') "_cell_length_b     ", b
      write (72, '(a, f12.6)') "_cell_length_c     ", c
      write (72, '(a, f12.6)') "_cell_angle_alpha  ", alpha
      write (72, '(a, f12.6)') "_cell_angle_beta   ", beta
      write (72, '(a, f12.6)') "_cell_angle_gamma  ", gamma
      write (72, '(a)') " "
      write (72, '(a)') "_symmetry_space_group_name_H-M  'P 1'"
      write (72, '(a)') "_symmetry_Int_Tables_number      1"
      write (72, '(a)') " "
      write (72, '(a)') "loop_"
      write (72, '(a)') "_symmetry_equiv_pos_as_xyz"
      write (72, '(a)') "'x, y, z'"
      write (72, '(a)') " "
      write (72, '(a)') "loop_"
      write (72, '(a)') "_atom_site_label"
      write (72, '(a)') "_atom_site_type_symbol"
      write (72, '(a)') "_atom_site_fract_x"
      write (72, '(a)') "_atom_site_fract_y"
      write (72, '(a)') "_atom_site_fract_z"

      do at = 1, nat
        sp = spat(at)
        if (sp /= sptarget) then
          write (72, '(a, i0, 2x, a, 3(2x, f11.7))') &
                trim(symbol(sp)), at, trim(symbol(sp)), coords(at, 1), coords(at, 2), coords(at, 3)
        else
          att = at - atini + 1
          call member(nsubs, newconf, att, ifound)
          if (ifound == 1) then
            write (72, '(a, i0, 2x, a, 3(2x, f11.7))') &
                  trim(newsymbol(1)), at, trim(newsymbol(1)), coords(at, 1), coords(at, 2), coords(at, 3)
          else
            write (72, '(a, i0, 2x, a, 3(2x, f11.7))') &
                  trim(newsymbol(2)), at, trim(newsymbol(2)), coords(at, 1), coords(at, 2), coords(at, 3)
          end if
        end if
      end do

      close (unit=72)

    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!  GULP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case (1)

    ! --- Format statements used by structure block and by CASTEP/QE cases ---
211 format(6(f10.4, 2x))
331 format(a3, 2x, a4, 2x, 3(f11.7, 2x))
332 format(a85)
333 format(a85)

    ! --- Read template_input.gin into memory ---
    open (unit=70, file="template_input.gin", status="old", iostat=ios)
    if (ios /= 0) then
      write (*, *) "Error: template_input.gin not found. Aborting."
      stop
    end if
    ngulplines = 0
    do
      read (70, '(a)', end=1001) gulpline_buf
      ngulplines = ngulplines + 1
      if (ngulplines > maxgulplines) then
        write (*, *) "Error: template_input.gin exceeds ", maxgulplines, " lines. Aborting."
        stop
      end if
      gulptemplate(ngulplines) = gulpline_buf
    end do
1001 continue
    close (70)

    ! --- Validate @configuration_structure@ ---
    struc_line_idx = 0
    do l = 1, ngulplines
      sline_trimmed = adjustl(gulptemplate(l))
      if (trim(sline_trimmed) == "@configuration_structure@") then
        if (struc_line_idx /= 0) then
          write (*, *) "Error: @configuration_structure@ appears more than once in template_input.gin. Aborting."
          stop
        end if
        struc_line_idx = l
      else if (index(gulptemplate(l), "@configuration_structure@") /= 0) then
        write (*, *) "Error: @configuration_structure@ is not alone on its line. Aborting."
        stop
      end if
    end do
    if (struc_line_idx == 0) then
      write (*, *) "Error: @configuration_structure@ not found in template_input.gin. Aborting."
      stop
    end if

    ! --- Build all_sym array ---
    all_nsp_gulp = 0
    do ssp = 1, nsp
      if (ssp < sptarget) then
        all_nsp_gulp = all_nsp_gulp + 1
        all_sym(all_nsp_gulp) = symbol(ssp)
      else if (ssp == sptarget) then
        all_nsp_gulp = all_nsp_gulp + 1
        all_sym(all_nsp_gulp) = newsymbol(1)
        all_nsp_gulp = all_nsp_gulp + 1
        all_sym(all_nsp_gulp) = newsymbol(2)
      else
        all_nsp_gulp = all_nsp_gulp + 1
        all_sym(all_nsp_gulp) = symbol(ssp)
      end if
    end do

    ! --- Parse sod_type_map comment lines ---
    nmap_gulp = 0
    do l = 1, ngulplines
      sline_trimmed = adjustl(gulptemplate(l))
      if (sline_trimmed(1:15) == "# sod_type_map ") then
        nmap_gulp = nmap_gulp + 1
        read (sline_trimmed(16:), *, iostat=ios) map_sod_arr(nmap_gulp), map_gulp_arr(nmap_gulp)
        if (ios /= 0) then
          write (*, *) "Error parsing sod_type_map line: ", trim(sline_trimmed)
          stop
        end if
      end if
    end do

    ! --- Apply mapping to get GULP types (default: SOD name = GULP type) ---
    do i = 1, all_nsp_gulp
      all_gulptype_arr(i) = all_sym(i)
      do m = 1, nmap_gulp
        if (trim(all_sym(i)) == trim(map_sod_arr(m))) then
          all_gulptype_arr(i) = map_gulp_arr(m)
          exit
        end if
      end do
    end do

    ! --- Detect library directive ---
    haslibrary_gulp = .false.
    libname_gulp = ""
    do l = 1, ngulplines
      sline_trimmed = adjustl(gulptemplate(l))
      if (sline_trimmed(1:8) == "library ") then
        haslibrary_gulp = .true.
        libname_gulp = trim(adjustl(sline_trimmed(9:)))
        exit
      end if
    end do

    ! --- Initialize core/shell flags to core-only ---
    all_ishell_arr(1:all_nsp_gulp) = 0

    ! --- Scan settings.gulp for inline species/shel entries ---
    hasinline_ff = .false.
    do l = 1, ngulplines
      sline_trimmed = adjustl(gulptemplate(l))
      if (len_trim(sline_trimmed) == 0) cycle
      if (sline_trimmed(1:1) == "#") cycle
      read (sline_trimmed, *, iostat=ios) gline_type, gline_coreshel
      if (ios == 0) then
        if (trim(gline_coreshel) == "core" .or. trim(gline_coreshel) == "shel") then
          hasinline_ff = .true.
          if (trim(gline_coreshel) == "shel") then
            do i = 1, all_nsp_gulp
              if (trim(all_gulptype_arr(i)) == trim(gline_type)) then
                all_ishell_arr(i) = 1
              end if
            end do
          end if
        end if
      end if
    end do

    ! --- Scan library file for shel entries if present ---
    if (haslibrary_gulp) then
      open (unit=71, file=trim(libname_gulp), status="old", iostat=ios)
      if (ios /= 0) then
        write (*, *) "Error: library file '", trim(libname_gulp), "' not found. Aborting."
        stop
      end if
      do
        read (71, '(a)', end=1002) gulpline_buf
        sline_trimmed = adjustl(gulpline_buf)
        if (len_trim(sline_trimmed) == 0) cycle
        if (sline_trimmed(1:1) == "#") cycle
        read (sline_trimmed, *, iostat=ios) gline_type, gline_coreshel
        if (ios == 0 .and. trim(gline_coreshel) == "shel") then
          do i = 1, all_nsp_gulp
            if (trim(all_gulptype_arr(i)) == trim(gline_type)) then
              all_ishell_arr(i) = 1
            end if
          end do
        end if
      end do
1002  continue
      close (71)
    end if

    ! --- Validate that force-field information is available ---
    if (.not. hasinline_ff .and. .not. haslibrary_gulp) then
      write (*, *) "Error: template_input.gin has no inline force-field information and no library directive. Aborting."
      stop
    end if

    ! --- Build number-only format string ---
    write (numfmt, '(a,i0,a,i0,a)') '(i', ndigits, '.', ndigits, ')'

    ! --- Build nXX parent directory and create it ---
    write (ndir_gulp, '(a,i2.2)') 'n', nsubs
    call execute_command_line('mkdir -p ' // trim(ndir_gulp), exitstat=ios)
    if (ios /= 0) then
      write (*, *) "Error: could not create directory ", trim(ndir_gulp)
      stop
    end if

    ! --- Generate one configuration directory per configuration ---
    do indcount = 1, nic
      newconf(1:nsubs) = indconf(indcount, 1:nsubs)

      write (numstr, numfmt) indcount
      confdir_gulp = trim(ndir_gulp) // '/c' // trim(numstr)

      call execute_command_line('mkdir -p ' // trim(confdir_gulp), exitstat=ios)
      if (ios /= 0) then
        write (*, *) "Error: could not create directory ", trim(confdir_gulp)
        stop
      end if

      inpfile_gulp = trim(confdir_gulp) // '/input.gin'
      open (unit=72, file=trim(inpfile_gulp), status="replace")

      do l = 1, ngulplines
        if (l == struc_line_idx) then

          write (72, '(a)') "cell"
          write (72, 211) a, b, c, alpha, beta, gamma
          write (72, '(a)') "frac"

          do at = 1, nat
            sp = spat(at)
            if (sp /= sptarget) then
              gtype_buf = symbol(sp)
              do m = 1, all_nsp_gulp
                if (trim(all_sym(m)) == trim(symbol(sp))) then
                  gtype_buf = all_gulptype_arr(m)
                  exit
                end if
              end do
              write (72, 331) trim(gtype_buf), "core", coords(at, 1), coords(at, 2), coords(at, 3)
              do m = 1, all_nsp_gulp
                if (trim(all_gulptype_arr(m)) == trim(gtype_buf) .and. all_ishell_arr(m) == 1) then
                  write (72, 331) trim(gtype_buf), "shel", coords(at, 1), coords(at, 2), coords(at, 3)
                  exit
                end if
              end do
            else
              att = at - atini + 1
              call member(nsubs, newconf, att, ifound)
              if (ifound == 1) then
                gtype_buf = newsymbol(1)
                do m = 1, all_nsp_gulp
                  if (trim(all_sym(m)) == trim(newsymbol(1))) then
                    gtype_buf = all_gulptype_arr(m)
                    exit
                  end if
                end do
                write (72, 331) trim(gtype_buf), "core", coords(at, 1), coords(at, 2), coords(at, 3)
                do m = 1, all_nsp_gulp
                  if (trim(all_gulptype_arr(m)) == trim(gtype_buf) .and. all_ishell_arr(m) == 1) then
                    write (72, 331) trim(gtype_buf), "shel", coords(at, 1), coords(at, 2), coords(at, 3)
                    exit
                  end if
                end do
              else
                gtype_buf = newsymbol(2)
                do m = 1, all_nsp_gulp
                  if (trim(all_sym(m)) == trim(newsymbol(2))) then
                    gtype_buf = all_gulptype_arr(m)
                    exit
                  end if
                end do
                write (72, 331) trim(gtype_buf), "core", coords(at, 1), coords(at, 2), coords(at, 3)
                do m = 1, all_nsp_gulp
                  if (trim(all_gulptype_arr(m)) == trim(gtype_buf) .and. all_ishell_arr(m) == 1) then
                    write (72, 331) trim(gtype_buf), "shel", coords(at, 1), coords(at, 2), coords(at, 3)
                    exit
                  end if
                end do
              end if
            end if
          end do

        else

          outbuf_gulp = gulptemplate(l)
          sline_trimmed = adjustl(outbuf_gulp)
          if (sline_trimmed(1:15) == "# sod_type_map ") cycle
          do while (index(outbuf_gulp, "@configuration_number@") /= 0)
            idx_token = index(outbuf_gulp, "@configuration_number@")
            outbuf_gulp = outbuf_gulp(1:idx_token-1) // trim(numstr) // &
                          outbuf_gulp(idx_token+22:)
          end do
          write (72, '(a)') trim(outbuf_gulp)

        end if
      end do

      close (72)

      if (haslibrary_gulp) then
        inpfile_gulp = trim(confdir_gulp) // '/' // trim(libname_gulp)
        open (unit=73, file=trim(inpfile_gulp), status="replace")
        open (unit=71, file=trim(libname_gulp), status="old")
        pending_section_header = ''
        do
          read (71, '(a)', end=1003) gulpline_buf
          sline_trimmed = adjustl(gulpline_buf)
          ! Blank lines and comments: pass through directly
          if (len_trim(sline_trimmed) == 0 .or. sline_trimmed(1:1) == "#") then
            write (73, '(a)') trim(gulpline_buf)
            cycle
          end if
          ! Single-word line: potential section keyword (e.g. morse, spring, three)
          ! Hold it as pending; only write it when a relevant content line follows
          if (index(trim(sline_trimmed), ' ') == 0) then
            pending_section_header = gulpline_buf
            cycle
          end if
          ! Content line: check whether first token is a relevant atom type
          read (sline_trimmed, *, iostat=ios) gline_type
          relevant_line = .false.
          do i = 1, all_nsp_gulp
            if (trim(gline_type) == trim(all_gulptype_arr(i))) then
              relevant_line = .true.
              exit
            end if
          end do
          if (relevant_line) then
            if (len_trim(pending_section_header) > 0) then
              write (73, '(a)') trim(pending_section_header)
              pending_section_header = ''
            end if
            write (73, '(a)') trim(gulpline_buf)
          end if
        end do
1003    continue
        close (71)
        close (73)
      end if

    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!  LAMMPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case (2)

    ! --- Read template_in.lammps into template buffer ---
    open (unit=70, file="template_in.lammps", status="old", iostat=ios)
    if (ios /= 0) then
      write (*, *) "Error: template_in.lammps not found. Aborting."
      stop
    end if
    ngulplines = 0
    do
      read (70, '(a)', end=2001) gulpline_buf
      ngulplines = ngulplines + 1
      if (ngulplines > maxgulplines) then
        write (*, *) "Error: template_in.lammps exceeds ", maxgulplines, " lines. Aborting."
        stop
      end if
      gulptemplate(ngulplines) = gulpline_buf
    end do
2001 continue
    close (70)

    ! --- Build all_sym array ---
    all_nsp_gulp = 0
    do ssp = 1, nsp
      if (ssp < sptarget) then
        all_nsp_gulp = all_nsp_gulp + 1
        all_sym(all_nsp_gulp) = symbol(ssp)
      else if (ssp == sptarget) then
        all_nsp_gulp = all_nsp_gulp + 1
        all_sym(all_nsp_gulp) = newsymbol(1)
        all_nsp_gulp = all_nsp_gulp + 1
        all_sym(all_nsp_gulp) = newsymbol(2)
      else
        all_nsp_gulp = all_nsp_gulp + 1
        all_sym(all_nsp_gulp) = symbol(ssp)
      end if
    end do

    ! --- Parse atom_style ---
    lammps_atom_style = ""
    do l = 1, ngulplines
      sline_trimmed = adjustl(gulptemplate(l))
      if (len_trim(sline_trimmed) == 0 .or. sline_trimmed(1:1) == "#") cycle
      read (sline_trimmed, *, iostat=ios) lmptok1_buf, lmptok2_buf
      if (ios == 0 .and. trim(lmptok1_buf) == "atom_style") then
        lammps_atom_style = trim(lmptok2_buf)
        exit
      end if
    end do
    if (len_trim(lammps_atom_style) == 0) then
      write (*, *) "Error: atom_style not found in template_in.lammps. Aborting."
      stop
    end if

    ! --- Parse # sod_type_map lines ---
    lammps_nmap = 0
    do l = 1, ngulplines
      sline_trimmed = adjustl(gulptemplate(l))
      if (sline_trimmed(1:15) /= "# sod_type_map ") cycle
      lammps_nmap = lammps_nmap + 1
      lmpslinetail = sline_trimmed(16:)
      read (lmpslinetail, *, iostat=ios) lmptok1_buf, lmptok2_buf
      if (ios /= 0) then
        write (*, *) "Error parsing sod_type_map: ", trim(sline_trimmed)
        stop
      end if
      lammps_map_sod(lammps_nmap) = trim(lmptok1_buf)
      lammps_map_role(lammps_nmap) = trim(lmptok2_buf)
      ! Extract bond_type= (if present)
      lmpbondval = 0
      lmpidx = index(lmpslinetail, "bond_type=")
      if (lmpidx > 0) then
        read (lmpslinetail(lmpidx+10:), *, iostat=ios) lmpbondval
        if (ios /= 0) then
          write (*, *) "Error: invalid bond_type= in: ", trim(sline_trimmed)
          stop
        end if
      end if
      lammps_map_bond(lammps_nmap) = lmpbondval
      ! Extract type= (masking bond_type= to avoid substring collision)
      lmptmpline = lmpslinetail
      if (lmpidx > 0) then
        do lmp_i = lmpidx, lmpidx + 9
          lmptmpline(lmp_i:lmp_i) = ' '
        end do
      end if
      lmpidx = index(lmptmpline, "type=")
      if (lmpidx == 0) then
        write (*, *) "Error: no type= in sod_type_map: ", trim(sline_trimmed)
        stop
      end if
      read (lmptmpline(lmpidx+5:), *, iostat=ios) lmptypeval
      if (ios /= 0) then
        write (*, *) "Error: invalid type= in: ", trim(sline_trimmed)
        stop
      end if
      lammps_map_type(lammps_nmap) = lmptypeval
    end do

    ! --- Build per-species core/shell type arrays ---
    lammps_core_type(1:all_nsp_gulp) = 0
    lammps_shell_type(1:all_nsp_gulp) = 0
    lammps_shell_bond(1:all_nsp_gulp) = 0

    do i = 1, all_nsp_gulp
      do lmp_i = 1, lammps_nmap
        if (trim(all_sym(i)) == trim(lammps_map_sod(lmp_i))) then
          if (trim(lammps_map_role(lmp_i)) == "core") then
            lammps_core_type(i) = lammps_map_type(lmp_i)
          else if (trim(lammps_map_role(lmp_i)) == "shell") then
            lammps_shell_type(i) = lammps_map_type(lmp_i)
            lammps_shell_bond(i) = lammps_map_bond(lmp_i)
          end if
        end if
      end do
      if (lammps_core_type(i) == 0) then
        write (*, *) "Error: SOD species '", trim(all_sym(i)), &
                     "' has no core mapping in template_in.lammps. Aborting."
        stop
      end if
      if (lammps_shell_type(i) > 0 .and. lammps_shell_bond(i) == 0) then
        write (*, *) "Error: SOD species '", trim(all_sym(i)), &
                     "' has shell mapping but no bond_type. Aborting."
        stop
      end if
    end do

    ! --- Check atom_style supports shells ---
    has_shells_lammps = .false.
    do i = 1, all_nsp_gulp
      if (lammps_shell_type(i) > 0) then
        has_shells_lammps = .true.
        exit
      end if
    end do
    if (has_shells_lammps .and. trim(lammps_atom_style) /= "full") then
      write (*, *) "Error: shell species require atom_style full. Aborting."
      stop
    end if

    ! --- Compute max atom and bond type counts ---
    natom_types_lammps = 0
    nbond_types_lammps = 0
    do lmp_i = 1, lammps_nmap
      if (lammps_map_type(lmp_i) > natom_types_lammps) natom_types_lammps = lammps_map_type(lmp_i)
      if (lammps_map_bond(lmp_i) > nbond_types_lammps) nbond_types_lammps = lammps_map_bond(lmp_i)
    end do

    ! --- Compute cell vectors for box definition ---
    call cell(cellvector, a, b, c, alpha, beta, gamma)
    xy_tilt = cellvector(1, 2)
    xz_tilt = cellvector(1, 3)
    yz_tilt = cellvector(2, 3)
    is_triclinic_lammps = (abs(xy_tilt) > 1.0e-6 .or. &
                           abs(xz_tilt) > 1.0e-6 .or. &
                           abs(yz_tilt) > 1.0e-6)

    ! --- Build number format ---
    write (numfmt, '(a,i0,a,i0,a)') '(i', ndigits, '.', ndigits, ')'

    ! --- Build nXX parent directory and create it ---
    write (ndir_gulp, '(a,i2.2)') 'n', nsubs
    call execute_command_line('mkdir -p ' // trim(ndir_gulp), exitstat=ios)
    if (ios /= 0) then
      write (*, *) "Error: could not create directory ", trim(ndir_gulp)
      stop
    end if

    ! --- Generate one configuration directory per configuration ---
    do indcount = 1, nic
      newconf(1:nsubs) = indconf(indcount, 1:nsubs)

      write (numstr, numfmt) indcount
      confdir_gulp = trim(ndir_gulp) // '/c' // trim(numstr)

      call execute_command_line('mkdir -p ' // trim(confdir_gulp), exitstat=ios)
      if (ios /= 0) then
        write (*, *) "Error: could not create directory ", trim(confdir_gulp)
        stop
      end if

      ! --- Precompute LAMMPS atom IDs and species indices ---
      atom_id_lammps = 0
      mol_id_lammps = 0
      do at = 1, nat
        sp = spat(at)
        if (sp /= sptarget) then
          lammps_sym_cur = symbol(sp)
        else
          att = at - atini + 1
          call member(nsubs, newconf, att, ifound)
          if (ifound == 1) then
            lammps_sym_cur = newsymbol(1)
          else
            lammps_sym_cur = newsymbol(2)
          end if
        end if
        lammps_sym_idx = 0
        do lmp_i = 1, all_nsp_gulp
          if (trim(all_sym(lmp_i)) == trim(lammps_sym_cur)) then
            lammps_sym_idx = lmp_i
            exit
          end if
        end do
        lmp_sym_idx(at) = lammps_sym_idx
        atom_id_lammps = atom_id_lammps + 1
        lmp_core_id(at) = atom_id_lammps
        if (lammps_shell_type(lammps_sym_idx) > 0) then
          mol_id_lammps = mol_id_lammps + 1
          lmp_mol_id(at) = mol_id_lammps
          atom_id_lammps = atom_id_lammps + 1
          lmp_shell_id(at) = atom_id_lammps
        else
          lmp_mol_id(at) = 0
          lmp_shell_id(at) = 0
        end if
      end do
      natoms_lammps = atom_id_lammps
      nbonds_lammps = mol_id_lammps

      ! --- Write conf.data ---
      inpfile_gulp = trim(confdir_gulp) // '/conf.data'
      open (unit=72, file=trim(inpfile_gulp), status="replace")

      write (72, '(a, i0)') "LAMMPS data file generated by SOD - configuration ", indcount
      write (72, '(a)') " "
      write (72, '(i0, a)') natoms_lammps, " atoms"
      if (nbonds_lammps > 0) write (72, '(i0, a)') nbonds_lammps, " bonds"
      write (72, '(i0, a)') natom_types_lammps, " atom types"
      if (nbond_types_lammps > 0) write (72, '(i0, a)') nbond_types_lammps, " bond types"
      write (72, '(a)') " "
      write (72, '(2(f14.6, 2x), a)') 0.0, cellvector(1, 1), "xlo xhi"
      write (72, '(2(f14.6, 2x), a)') 0.0, cellvector(2, 2), "ylo yhi"
      write (72, '(2(f14.6, 2x), a)') 0.0, cellvector(3, 3), "zlo zhi"
      if (is_triclinic_lammps) then
        write (72, '(3(f14.6, 2x), a)') xy_tilt, xz_tilt, yz_tilt, "xy xz yz"
      end if
      write (72, '(a)') " "

      if (trim(lammps_atom_style) == "atomic") then
        write (72, '(a)') "Atoms  # atomic"
      else if (trim(lammps_atom_style) == "charge") then
        write (72, '(a)') "Atoms  # charge"
      else
        write (72, '(a)') "Atoms  # full"
      end if
      write (72, '(a)') " "

      do at = 1, nat
        lammps_sym_idx = lmp_sym_idx(at)
        xc = cellvector(1,1)*coords(at,1) + cellvector(1,2)*coords(at,2) + cellvector(1,3)*coords(at,3)
        yc = cellvector(2,2)*coords(at,2) + cellvector(2,3)*coords(at,3)
        zc = cellvector(3,3)*coords(at,3)
        if (trim(lammps_atom_style) == "atomic") then
          write (72, '(i0, 2x, i0, 3(2x, f14.6))') &
                lmp_core_id(at), lammps_core_type(lammps_sym_idx), xc, yc, zc
        else if (trim(lammps_atom_style) == "charge") then
          write (72, '(i0, 2x, i0, 2x, f8.4, 3(2x, f14.6))') &
                lmp_core_id(at), lammps_core_type(lammps_sym_idx), 0.0, xc, yc, zc
        else  ! full
          write (72, '(i0, 2x, i0, 2x, i0, 2x, f8.4, 3(2x, f14.6))') &
                lmp_core_id(at), lmp_mol_id(at), lammps_core_type(lammps_sym_idx), 0.0, xc, yc, zc
          if (lmp_shell_id(at) > 0) then
            write (72, '(i0, 2x, i0, 2x, i0, 2x, f8.4, 3(2x, f14.6))') &
                  lmp_shell_id(at), lmp_mol_id(at), lammps_shell_type(lammps_sym_idx), 0.0, xc, yc, zc
          end if
        end if
      end do

      if (nbonds_lammps > 0) then
        write (72, '(a)') " "
        write (72, '(a)') "Bonds"
        write (72, '(a)') " "
        bond_id_lammps = 0
        do at = 1, nat
          if (lmp_shell_id(at) > 0) then
            lammps_sym_idx = lmp_sym_idx(at)
            bond_id_lammps = bond_id_lammps + 1
            write (72, '(i0, 2x, i0, 2x, i0, 2x, i0)') &
                  bond_id_lammps, lammps_shell_bond(lammps_sym_idx), &
                  lmp_core_id(at), lmp_shell_id(at)
          end if
        end do
      end if

      close (72)

      ! --- Write template_in.lammps from template with token replacement ---
      inpfile_gulp = trim(confdir_gulp) // '/in.lammps'
      open (unit=72, file=trim(inpfile_gulp), status="replace")
      do l = 1, ngulplines
        outbuf_gulp = gulptemplate(l)
        sline_trimmed = adjustl(outbuf_gulp)
        if (sline_trimmed(1:15) == "# sod_type_map ") cycle
        do while (index(outbuf_gulp, "@configuration_number@") /= 0)
          idx_token = index(outbuf_gulp, "@configuration_number@")
          outbuf_gulp = outbuf_gulp(1:idx_token-1) // trim(numstr) // &
                        outbuf_gulp(idx_token+22:)
        end do
        write (72, '(a)') trim(outbuf_gulp)
      end do
      close (72)

    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!  VASP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case (11)

    title = 'vasp'
    call cell(cellvector, a, b, c, alpha, beta, gamma)

    write (numfmt, '(a,i0,a,i0,a)') '(i', ndigits, '.', ndigits, ')'

    write (ndir_gulp, '(a,i2.2)') 'n', nsubs
    call execute_command_line('mkdir -p ' // trim(ndir_gulp), exitstat=ios)
    if (ios /= 0) then
      write (*, *) "Error: could not create directory ", trim(ndir_gulp)
      stop
    end if

    do indcount = 1, nic
      newconf(1:nsubs) = indconf(indcount, 1:nsubs)

      write (numstr, numfmt) indcount
      confdir_gulp = trim(ndir_gulp) // '/c' // trim(numstr)
      call execute_command_line('mkdir -p ' // trim(confdir_gulp), exitstat=ios)
      if (ios /= 0) then
        write (*, *) "Error: could not create directory ", trim(confdir_gulp)
        stop
      end if

      inpfile_gulp = trim(confdir_gulp) // '/POSCAR'
      open (unit=72, file=trim(inpfile_gulp), status="replace")

      write (72, *) title
      write (72, *) '1.00000000'
335   format(3(f10.6, 2x))
      write (72, 335) cellvector(1, 1), cellvector(2, 1), cellvector(3, 1)
      write (72, 335) cellvector(1, 2), cellvector(2, 2), cellvector(3, 2)
      write (72, 335) cellvector(1, 3), cellvector(2, 3), cellvector(3, 3)

      do ssp = 1, nsp
        if (ssp < sptarget) then
          snatsp(ssp) = natsp(ssp)
          ssymbol(ssp) = symbol(ssp)
        end if
        if (ssp == sptarget) then
          snatsp(ssp) = nsubs
          snatsp(ssp + 1) = npos - nsubs
          ssymbol(ssp) = newsymbol(1)
          ssymbol(ssp + 1) = newsymbol(2)
        end if
        if (ssp > sptarget) then
          snatsp(ssp + 1) = natsp(ssp)
          ssymbol(ssp + 1) = symbol(ssp)
        end if
      end do

      write (72, *) ssymbol(1:nsp + 1)
336   format(10(i4, 1x))
      write (72, 336) snatsp(1:nsp + 1)
      write (72, *) 'Direct'
      do at = 1, atini - 1
        write (72, 335) coords(at, 1), coords(at, 2), coords(at, 3)
      end do
      do at = atini, atfin
        att = at - atini + 1
        call member(nsubs, newconf, att, ifound)
        if (ifound == 1) write (72, 335) coords(at, 1), coords(at, 2), coords(at, 3)
      end do
      do at = atini, atfin
        att = at - atini + 1
        call member(nsubs, newconf, att, ifound)
        if (ifound == 0) write (72, 335) coords(at, 1), coords(at, 2), coords(at, 3)
      end do
      do at = atfin + 1, nat
        write (72, 335) coords(at, 1), coords(at, 2), coords(at, 3)
      end do

      close (72)

    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!  CASTEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case (12)

    ! --- Read template_castep.cell into memory ---
    open (unit=70, file="template_castep.cell", status="old", iostat=ios)
    if (ios /= 0) then
      write (*, *) "Error: template_castep.cell not found. Aborting."
      stop
    end if
    ngulplines = 0
    do
      read (70, '(a)', end=1201) gulpline_buf
      ngulplines = ngulplines + 1
      if (ngulplines > maxgulplines) then
        write (*, *) "Error: template_castep.cell exceeds ", maxgulplines, " lines. Aborting."
        stop
      end if
      gulptemplate(ngulplines) = gulpline_buf
    end do
1201 continue
    close (70)

    ! --- Validate @configuration_structure@ ---
    struc_line_idx = 0
    do l = 1, ngulplines
      sline_trimmed = adjustl(gulptemplate(l))
      if (trim(sline_trimmed) == "@configuration_structure@") then
        if (struc_line_idx /= 0) then
          write (*, *) "Error: @configuration_structure@ appears more than once in template_castep.cell. Aborting."
          stop
        end if
        struc_line_idx = l
      else if (index(gulptemplate(l), "@configuration_structure@") /= 0) then
        write (*, *) "Error: @configuration_structure@ is not alone on its line. Aborting."
        stop
      end if
    end do
    if (struc_line_idx == 0) then
      write (*, *) "Error: @configuration_structure@ not found in template_castep.cell. Aborting."
      stop
    end if

    call cell(cellvector, a, b, c, alpha, beta, gamma)

    ! --- Build number format ---
    write (numfmt, '(a,i0,a,i0,a)') '(i', ndigits, '.', ndigits, ')'

    ! --- Build nXX parent directory and create it ---
    write (ndir_gulp, '(a,i2.2)') 'n', nsubs
    call execute_command_line('mkdir -p ' // trim(ndir_gulp), exitstat=ios)
    if (ios /= 0) then
      write (*, *) "Error: could not create directory ", trim(ndir_gulp)
      stop
    end if

    do indcount = 1, nic
      newconf(1:nsubs) = indconf(indcount, 1:nsubs)

      write (numstr, numfmt) indcount
      confdir_gulp = trim(ndir_gulp) // '/c' // trim(numstr)

      call execute_command_line('mkdir -p ' // trim(confdir_gulp), exitstat=ios)
      if (ios /= 0) then
        write (*, *) "Error: could not create directory ", trim(confdir_gulp)
        stop
      end if

      inpfile_gulp = trim(confdir_gulp) // '/castep.cell'
      open (unit=72, file=trim(inpfile_gulp), status="replace")

      do l = 1, ngulplines
        if (l == struc_line_idx) then

          write (72, '(a)') "%BLOCK lattice_cart"
          write (72, 335) cellvector(1, 1), cellvector(2, 1), cellvector(3, 1)
          write (72, 335) cellvector(1, 2), cellvector(2, 2), cellvector(3, 2)
          write (72, 335) cellvector(1, 3), cellvector(2, 3), cellvector(3, 3)
          write (72, '(a)') "%ENDBLOCK lattice_cart"
          write (72, '(a)') "%BLOCK positions_frac"

          do at = 1, nat
            sp = spat(at)
            if (sp /= sptarget) then
              write (72, 337) symbol(sp), coords(at, 1), coords(at, 2), coords(at, 3)
            else
              att = at - atini + 1
              call member(nsubs, newconf, att, ifound)
              if (ifound == 1) then
                write (72, 337) newsymbol(1), coords(at, 1), coords(at, 2), coords(at, 3)
              else
                write (72, 337) newsymbol(2), coords(at, 1), coords(at, 2), coords(at, 3)
              end if
            end if
337         format(a3, 3(f11.7, 2x))
          end do

          write (72, '(a)') "%ENDBLOCK positions_frac"

        else

          outbuf_gulp = gulptemplate(l)
          do while (index(outbuf_gulp, "@configuration_number@") /= 0)
            idx_token = index(outbuf_gulp, "@configuration_number@")
            outbuf_gulp = outbuf_gulp(1:idx_token-1) // trim(numstr) // &
                          outbuf_gulp(idx_token+22:)
          end do
          write (72, '(a)') trim(outbuf_gulp)

        end if
      end do

      close (72)

    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!  QUANTUM ESPRESSO !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case (13)

    ! --- Read template_pw.in into memory ---
    open (unit=70, file="template_pw.in", status="old", iostat=ios)
    if (ios /= 0) then
      write (*, *) "Error: template_pw.in not found. Aborting."
      stop
    end if
    ngulplines = 0
    do
      read (70, '(a)', end=1301) gulpline_buf
      ngulplines = ngulplines + 1
      if (ngulplines > maxgulplines) then
        write (*, *) "Error: template_pw.in exceeds ", maxgulplines, " lines. Aborting."
        stop
      end if
      gulptemplate(ngulplines) = gulpline_buf
    end do
1301 continue
    close (70)

    ! --- Validate @configuration_structure@ ---
    struc_line_idx = 0
    do l = 1, ngulplines
      sline_trimmed = adjustl(gulptemplate(l))
      if (trim(sline_trimmed) == "@configuration_structure@") then
        if (struc_line_idx /= 0) then
          write (*, *) "Error: @configuration_structure@ appears more than once in template_pw.in. Aborting."
          stop
        end if
        struc_line_idx = l
      else if (index(gulptemplate(l), "@configuration_structure@") /= 0) then
        write (*, *) "Error: @configuration_structure@ is not alone on its line. Aborting."
        stop
      end if
    end do
    if (struc_line_idx == 0) then
      write (*, *) "Error: @configuration_structure@ not found in template_pw.in. Aborting."
      stop
    end if

    call cell(cellvector, a, b, c, alpha, beta, gamma)

    ! --- Build number format ---
    write (numfmt, '(a,i0,a,i0,a)') '(i', ndigits, '.', ndigits, ')'

    ! --- Build nXX parent directory and create it ---
    write (ndir_gulp, '(a,i2.2)') 'n', nsubs
    call execute_command_line('mkdir -p ' // trim(ndir_gulp), exitstat=ios)
    if (ios /= 0) then
      write (*, *) "Error: could not create directory ", trim(ndir_gulp)
      stop
    end if

    do indcount = 1, nic
      newconf(1:nsubs) = indconf(indcount, 1:nsubs)

      write (numstr, numfmt) indcount
      confdir_gulp = trim(ndir_gulp) // '/c' // trim(numstr)

      call execute_command_line('mkdir -p ' // trim(confdir_gulp), exitstat=ios)
      if (ios /= 0) then
        write (*, *) "Error: could not create directory ", trim(confdir_gulp)
        stop
      end if

      inpfile_gulp = trim(confdir_gulp) // '/pw.in'
      open (unit=72, file=trim(inpfile_gulp), status="replace")

      do l = 1, ngulplines
        if (l == struc_line_idx) then

          write (72, '(a)') "CELL_PARAMETERS {angstrom}"
          write (72, 335) cellvector(1, 1), cellvector(2, 1), cellvector(3, 1)
          write (72, 335) cellvector(1, 2), cellvector(2, 2), cellvector(3, 2)
          write (72, 335) cellvector(1, 3), cellvector(2, 3), cellvector(3, 3)

          write (72, '(a)') "ATOMIC_POSITIONS {crystal}"
          do at = 1, nat
            sp = spat(at)
            if (sp /= sptarget) then
              write (72, 337) symbol(sp), coords(at, 1), coords(at, 2), coords(at, 3)
            else
              att = at - atini + 1
              call member(nsubs, newconf, att, ifound)
              if (ifound == 1) then
                write (72, 337) newsymbol(1), coords(at, 1), coords(at, 2), coords(at, 3)
              else
                write (72, 337) newsymbol(2), coords(at, 1), coords(at, 2), coords(at, 3)
              end if
            end if
          end do

        else

          outbuf_gulp = gulptemplate(l)
          do while (index(outbuf_gulp, "@configuration_number@") /= 0)
            idx_token = index(outbuf_gulp, "@configuration_number@")
            outbuf_gulp = outbuf_gulp(1:idx_token-1) // trim(numstr) // &
                          outbuf_gulp(idx_token+22:)
          end do
          write (72, '(a)') trim(outbuf_gulp)

        end if
      end do

      close (72)

    end do

    end select

    deallocate (newconf)
    deallocate (degen)
    deallocate (indconf)

  end do  ! nsubs

  close (43)

!!!!!!!Reporting the end
  write (*, *) "Done!!!"
  write (*, *) ""
  write (*, *) ""

end program genersod
