!*******************************************************************************
!    Copyright (c) 2018 Ricardo Grau-Crespo and co-authors
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
  use iso_fortran_env, only: real64
  use insod_reader,    only: insod_t, read_insod
  use ensemble_io,     only: read_ensemble
  implicit none

  integer, parameter :: nspmax = 10, natmax = 10000, nlineamax = 200
  integer, parameter :: maxgulplines = 1000
  integer, parameter :: ntargetmax = 5, nkmax = 5

  integer :: i, j, l, m
  integer :: arg_count, max_generate, gen_limit
  integer :: mode, n_chosen, gen_count, ic
  character(len=32) :: choose_label
  integer :: itest, choose_count
  integer :: ifound
  integer :: sp_slot, col_off, col_off_t  ! for multi-nary species-slot determination
  integer :: sp, nsp, ssp
  integer :: ntarget, sptarget(ntargetmax)
  integer :: nk(ntargetmax)
  integer :: nsubs_t(ntargetmax, nkmax), atini_t(ntargetmax), atfin_t(ntargetmax), npos_t(ntargetmax)
  integer :: nsubs_A, nsubs_B, nsubs_tot, nsubs_tot_t(ntargetmax), t_A, t
  integer :: at0, nat0, at, nat, att, filer, ndigits
  character(len=20) :: outfilename, fmtstr
  integer :: na, nb, nc, nsubs, nsubs_min, nsubs_max, atini, atfin, cumnatsp
  integer :: nic, indcount
  integer, dimension(:), allocatable :: newconf, degen
  integer, allocatable :: chosen_indices(:), gen_indices(:)
  integer, dimension(nspmax) :: natsp0, natsp, snatsp
  integer, dimension(natmax) :: spat
  real(real64), dimension(natmax, 3) :: coords0, coords
  real(real64), dimension(3, 3) :: cellvector
  integer, dimension(:, :), allocatable :: indconf
  real(real64) :: a1, b1, c1, alpha, beta, gamma, a, b, c, xc, yc, zc
  character(len=3), dimension(nspmax) :: symbol, ssymbol
  character(len=5), dimension(ntargetmax, nkmax+1) :: newsymbol
  character(len=85) :: linea
  character(len=15) :: title
  character(len=40) :: runtitle
  type(insod_t) :: insod
  character(len=32) :: arg_text
  character(len=3) :: symboltrash
  character(len=200) :: cifline
  character(len=20) :: atmlabel
  character(len=3) :: atmsymbol
  logical :: in_atom_loop
  logical :: found
  logical :: fileexists

! Variables for molecule (@NAME) and vacancy (%NAME) handling
  integer, parameter :: nmolmax_atoms = 500
  integer, parameter :: nmoltypes = 10
  integer :: nmol_types, im, imt
  character(len=4) :: mol_name(nmoltypes)
  integer :: mol_natoms(nmoltypes)
  character(len=3) :: mol_sym_arr(nmoltypes, nmolmax_atoms)
  real(real64) :: mol_xyz_arr(nmoltypes, nmolmax_atoms, 3)
  ! Molecule (@NAME) / vacancy (%NAME) status of every newsymbol slot, indexed
  ! by (target, species slot).  Slot nk(t)+1 is the symbol kept for the sites of
  ! target t that are not substituted.
  logical :: is_mol_t(ntargetmax, nkmax+1), is_vac_t(ntargetmax, nkmax+1)
  integer :: mol_idx_t(ntargetmax, nkmax+1)
  ! Resolution of one target-site atom; see resolve_target_atom.
  integer, parameter :: site_plain = 0, site_molecule = 1, site_vacancy = 2
  integer :: site_action, site_mol
  character(len=5) :: site_sym
  logical :: has_molecules, has_vacancies
! Parent-structure molecules: placeholder species expanded to a rigid molecule
  logical :: sp_is_mol(nspmax)
  integer :: sp_mol_idx(nspmax)
  integer :: ip, sp_loc
! Expanded atom list (used for all filers when molecules/vacancies present)
  integer :: nat_exp, nsp_exp, isp_exp
  character(len=3), dimension(natmax) :: sym_exp
  real(real64), dimension(natmax, 3) :: coords_exp
  character(len=3), dimension(nspmax*2+nmoltypes*5+5) :: uniq_sym_exp
  integer :: uniq_cnt_exp(nspmax*2+nmoltypes*5+5)
  real(real64) :: mol_frac_buf(nmolmax_atoms, 3)
  real(real64) :: Rmat(3,3)

  ! Variables for GULP logic (FILER=1)
  integer :: ios, ngulplines, struc_line_idx, nmap_gulp, all_nsp_gulp, idx_token
  character(len=200) :: gulpline_buf, sline_trimmed, outbuf_gulp
  character(len=200), dimension(maxgulplines) :: gulptemplate
  character(len=3), dimension(nspmax) :: map_sod_arr, map_gulp_arr
  character(len=3), dimension(nspmax+nmoltypes*5+5) :: all_sym, all_gulptype_arr
  integer :: all_ishell_arr(nspmax+nmoltypes*5+5)
  character(len=80) :: libname_gulp
  character(len=3) :: gtype_buf, gline_type
  character(len=4) :: gline_coreshel
  character(len=200) :: pending_section_header
  logical :: haslibrary_gulp, hasinline_ff, haskimmodel_gulp, relevant_line
  character(len=20) :: numfmt
  character(len=10) :: numstr
  character(len=256) :: confdir_gulp
  character(len=256) :: inpfile_gulp
  character(len=256) :: ndir_gulp

  ! Variables for LAMMPS logic (FILER=2)
  character(len=20) :: lammps_atom_style
  integer :: lammps_nmap, lmptypeval, lmpbondval, lmpidx
  character(len=3), dimension(nspmax+nmoltypes*5+5) :: lammps_map_sod
  character(len=5), dimension(nspmax+nmoltypes*5+5) :: lammps_map_role
  integer :: lammps_map_type(nspmax+nmoltypes*5+5), lammps_map_bond(nspmax+nmoltypes*5+5)
  integer :: lammps_core_type(nspmax+nmoltypes*5+5), lammps_shell_type(nspmax+nmoltypes*5+5)
  integer :: lammps_shell_bond(nspmax+nmoltypes*5+5)
  integer :: natoms_lammps, nbonds_lammps, natom_types_lammps, nbond_types_lammps
  integer :: mol_id_lammps, atom_id_lammps, bond_id_lammps
  real(real64) :: xy_tilt, xz_tilt, yz_tilt
  logical :: has_shells_lammps, is_triclinic_lammps
  character(len=20) :: lmptok1_buf
  character(len=10) :: lmptok2_buf
  character(len=200) :: lmpslinetail, lmptmpline
  integer, dimension(natmax) :: lmp_core_id, lmp_shell_id, lmp_mol_id, lmp_sym_idx
  integer :: lammps_sym_idx, lmp_i
  character(len=3) :: lammps_sym_cur
  integer, dimension(natmax) :: exp_sym_idx, exp_core_id, exp_shell_id, exp_mol_id
  logical, dimension(natmax) :: exp_has_shell

  mode = 0
  max_generate = -1
  n_chosen = 0
  choose_label = ''
  choose_count = 1
  arg_count = command_argument_count()
  if (arg_count >= 1) then
    call get_command_argument(1, arg_text)
    if (trim(arg_text) == '-first') then
      if (arg_count /= 2) then
        write (*, *) "Error: -first requires exactly one argument N."
        stop 1
      end if
      call get_command_argument(2, arg_text)
      read (arg_text, *, iostat=ios) max_generate
      if (ios /= 0 .or. max_generate <= 0) then
        write (*, *) "Error: -first N requires N to be a positive integer."
        stop 1
      end if
      mode = 1
    else if (trim(arg_text) == '-choose') then
      if (arg_count < 2) then
        write (*, *) "Error: -choose requires a configuration index or a selection label."
        call write_choose_usage()
        stop 1
      end if
      ! -choose takes either explicit configuration indices or a single label
      ! naming a rule that picks one.  A label cannot be resolved yet: it needs
      ! the level directory and the configuration count, which are only known
      ! once the ENSEMBLE has been located, so it is validated here and applied
      ! further below.
      call get_command_argument(2, arg_text)
      read (arg_text, *, iostat=ios) itest
      if (ios /= 0) then
        if (arg_count > 3) then
          write (*, *) "Error: -choose <label> takes at most one count."
          call write_choose_usage()
          stop 1
        end if
        choose_label = to_lower(trim(arg_text))
        if (choose_label /= 'bestsqs' .and. choose_label /= 'lowestenergy' &
            .and. choose_label /= 'highestenergy' .and. choose_label /= 'random') then
          write (*, '(3A)') " Error: '", trim(arg_text), &
              "' is neither a configuration index nor a known selection label."
          call write_choose_usage()
          stop 1
        end if
        ! Optional count: how many of the best to take.  Default 1.
        if (arg_count == 3) then
          call get_command_argument(3, arg_text)
          read (arg_text, *, iostat=ios) choose_count
          if (ios /= 0 .or. choose_count <= 0) then
            write (*, '(2A)') " Error: the count after a -choose label must be a", &
                " positive integer."
            stop 1
          end if
        end if
      else
        n_chosen = arg_count - 1
        allocate(chosen_indices(n_chosen))
        do i = 1, n_chosen
          call get_command_argument(i + 1, arg_text)
          read (arg_text, *, iostat=ios) chosen_indices(i)
          if (ios /= 0 .or. chosen_indices(i) <= 0) then
            write (*, *) "Error: -choose indices must be positive integers."
            stop 1
          end if
        end do
      end if
      mode = 2
    else
      write (*, '(2A)') "Error: unknown argument: ", trim(arg_text)
      write (*, *) "Usage: genersod [-first N | -choose <index ...|label [N]>]"
      call write_choose_usage()
      stop 1
    end if
  end if

! Input files

  open (unit=31, file="supercell.cif")

! Output files

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

  write (*, '(A)') "SOD (Site-Occupancy Disorder) version 0.92 - genersod"
  write (*, *) "Reading INSOD, supercell.cif, and ENSEMBLE to generate calculation input files"
  write (*, *) " "

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the INSOD file
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  call read_insod('INSOD', insod)
  runtitle = insod%runtitle
  a1 = insod%a1; b1 = insod%b1; c1 = insod%c1
  alpha = insod%alpha; beta = insod%beta; gamma = insod%gamma
  nsp = insod%nsp
  symbol(1:nsp)  = insod%symbol
  natsp0(1:nsp)  = insod%natsp0
  nat0           = insod%nat0
  coords0(1:nat0, 1:3) = insod%coords0
  na = insod%na; nb = insod%nb; nc = insod%nc
  ntarget        = insod%ntarget
  sptarget(1:ntarget) = insod%sptarget(1:ntarget)
  nk(1:ntarget)  = insod%nk(1:ntarget)
  nsubs_t(1:ntarget, 1:nkmax) = insod%nsubs_t(1:ntarget, 1:nkmax)
  nsubs_min      = insod%nsubs_min
  nsubs_max      = insod%nsubs_max
  newsymbol      = insod%newsymbol
  filer          = insod%filer

  ! Parse @ and % prefixes from newsymbol for every target and every species
  ! slot.  Doing it generically is what lets the writers below handle 1, 2 or 5
  ! targets, binary or multi-nary, through a single code path.
  is_mol_t = .false.
  is_vac_t = .false.
  mol_idx_t = 0
  do t = 1, ntarget
    do j = 1, nk(t) + 1
      if (newsymbol(t,j)(1:1) == '@') then
        is_mol_t(t,j) = .true.
      else if (newsymbol(t,j)(1:1) == '%') then
        is_vac_t(t,j) = .true.
      end if
    end do
  end do
  has_molecules = any(is_mol_t(1:ntarget, :))
  has_vacancies = any(is_vac_t(1:ntarget, :))
! Parent-structure molecules (see mapping block at end of INSOD)
  sp_is_mol(:) = .false.
  sp_mol_idx(:) = 0
  has_molecules = has_molecules .or. (insod%n_parentmol > 0)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      [ENSEMBLE is now read per level inside the main loop below]
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

  do t_A = 1, ntarget
    atini_t(t_A) = 1
    do sp = 1, sptarget(t_A) - 1
      atini_t(t_A) = atini_t(t_A) + natsp(sp)
    end do
    atfin_t(t_A) = atini_t(t_A) + natsp(sptarget(t_A)) - 1
    npos_t(t_A) = natsp(sptarget(t_A))
  end do
  atini = atini_t(1)
  atfin = atfin_t(1)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Read molecule .xyz files for any @NAME symbols
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  nmol_types = 0

  if (has_molecules .or. has_vacancies) then
    call cell(cellvector, a, b, c, alpha, beta, gamma)
  end if

  if (has_molecules) call random_seed()

  do t = 1, ntarget
    do j = 1, nk(t) + 1
      if (.not. is_mol_t(t,j)) cycle
      call get_molecule(newsymbol(t,j)(2:), imt)
      mol_idx_t(t,j) = imt
    end do
  end do

! Load molecules for parent-structure placeholder species (@NAME parent mapping)
  do ip = 1, insod%n_parentmol
    sp_loc = 0
    do sp = 1, nsp
      if (trim(symbol(sp)) == trim(insod%parentmol_sym(ip))) then
        sp_loc = sp; exit
      end if
    end do
    if (sp_loc == 0) then
      write (*, *) "Error: parent-molecule placeholder '", trim(insod%parentmol_sym(ip)), &
                   "' is not a species listed in INSOD. Aborting."
      stop 1
    end if
    if (any(sp_loc == sptarget(1:ntarget))) then
      write (*, *) "Error: parent-molecule placeholder '", trim(insod%parentmol_sym(ip)), &
                   "' is also a substitution target; a site cannot be both. Aborting."
      stop 1
    end if
    call get_molecule(insod%parentmol_name(ip), imt)
    sp_is_mol(sp_loc) = .true.
    sp_mol_idx(sp_loc) = imt
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!    GENERATE INPUT FILES !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Loop over substitution levels
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  select case (filer)
  case (-1)
    write (*, *) " > FILER = -1: No calculation files will be created."
  case (0)
    write (*, *) " > Creating CIF files..."
  case (1)
    write (*, *) " > Creating GULP input files..."
  case (2)
    write (*, *) " > Creating LAMMPS input files..."
  case (11)
    write (*, *) " > Creating VASP POSCAR files..."
  case (12)
    write (*, *) " > Creating CASTEP input files..."
  case (13)
    write (*, *) " > Creating Quantum ESPRESSO input files..."
  end select
  write (*, *) ""

  if (ntarget >= 2 .and. any(nk(1:ntarget) >= 2)) then
!   Phase 4: at least one target is multi-nary
    do t = 1, ntarget
      nsubs_tot_t(t) = sum(nsubs_t(t, 1:nk(t)))
    end do
    nsubs_tot = sum(nsubs_tot_t(1:ntarget))
  else if (ntarget >= 2) then
!   Stage D: multi-target binary (all nk==1)
    do t = 1, ntarget
      nsubs_tot_t(t) = nsubs_t(t, 1)
    end do
    nsubs_tot = sum(nsubs_tot_t(1:ntarget))
  else if (ntarget == 1 .and. nk(1) >= 2) then
    nsubs_tot = sum(nsubs_t(1, 1:nk(1)))
  else
    nsubs_tot = nsubs_min  ! will be updated in loop
  end if

  do nsubs = nsubs_min, nsubs_max

    if (ntarget == 1 .and. nk(1) == 1) then
      ndir_gulp = 'n'//trim(npad(nsubs))
      nsubs_tot = nsubs
    else if (ntarget == 1 .and. nk(1) == 2) then
      ndir_gulp = 'n'//trim(npad(nsubs_t(1,1)))//'_'//trim(npad(nsubs_t(1,2)))
    else if (ntarget == 1 .and. nk(1) == 3) then
      ndir_gulp = 'n'//trim(npad(nsubs_t(1,1)))//'_'//trim(npad(nsubs_t(1,2)))//'_'//trim(npad(nsubs_t(1,3)))
    else if (ntarget >= 2 .and. all(nk(1:ntarget) == 1)) then
      ! Stage D: multi-target binary
      ndir_gulp = 'n'//trim(npad(nsubs_t(1,1)))
      do t = 2, ntarget
        ndir_gulp = trim(ndir_gulp)//'_'//trim(npad(nsubs_t(t,1)))
      end do
    else
      ! Phase 4: ntarget >= 2, at least one nk(t) >= 2
      ndir_gulp = 'n'//trim(npad(nsubs_t(1,1)))
      do j = 2, nk(1)
        ndir_gulp = trim(ndir_gulp)//'_'//trim(npad(nsubs_t(1,j)))
      end do
      do t = 2, ntarget
        do j = 1, nk(t)
          ndir_gulp = trim(ndir_gulp)//'_'//trim(npad(nsubs_t(t,j)))
        end do
      end do
    end if

    ! sod_gener.sh sets SOD_ENSEMBLE_DIR when launched from a subdir containing ENSEMBLE
    ! (e.g. nXX/random/), so the binary reads from there instead of the INSOD-derived ndir.
    block
      character(len=256) :: sod_edir
      integer :: edir_len
      integer :: ntarget_rd, nic_rd, iflat
      integer, allocatable :: nsfl_rd(:), npst_rd(:), degen_rd(:), indconf_rd(:,:)
      real(real64) :: tsampling_rd
      logical :: ensemble_ok

      call get_environment_variable('SOD_ENSEMBLE_DIR', sod_edir, length=edir_len)
      if (edir_len > 0) ndir_gulp = trim(adjustl(sod_edir))

      call read_ensemble(trim(ndir_gulp) // '/ENSEMBLE', ntarget_rd, nsfl_rd, npst_rd, &
                         nic_rd, indconf_rd, degen_rd, tsampling_rd, ensemble_ok)
      if (.not. ensemble_ok) then
        ! An explicit SOD_ENSEMBLE_DIR names one directory and nothing else can
        ! stand in for it, so a missing ENSEMBLE there is an error, not a skip.
        if (edir_len > 0) then
          write (*, *) "Error: ", trim(ndir_gulp), "/ENSEMBLE not found or invalid."
          stop 1
        end if
        write (*, *) "Warning: ", trim(ndir_gulp), "/ENSEMBLE not found or invalid, skipping."
        if (ntarget >= 2 .or. (ntarget == 1 .and. nk(1) >= 2)) exit
        cycle
      end if
      ! Unpack flat nsubs back into nsubs_t and npos_t

      iflat = 0
      do t = 1, ntarget
        do j = 1, nk(t)
          nsubs_t(t, j) = nsfl_rd(iflat + j)
        end do
        iflat = iflat + nk(t)
      end do
      npos_t(1:ntarget_rd) = npst_rd(1:ntarget_rd)
      nic = nic_rd
      ! Trust the local ENSEMBLE over INSOD for substitution counts (they may differ,
      ! e.g. when SOD_ENSEMBLE_DIR points to a subdir with a different composition).
      if (sum(nsfl_rd) /= nsubs_tot) then
        write (*, '(A,I0,A,I0,A)') &
          " Warning: INSOD specifies ", nsubs_tot, " substitutions but the local ENSEMBLE has ", &
          sum(nsfl_rd), ". Using the local value — no need to change INSOD."
      end if
      nsubs_tot = sum(nsfl_rd)
      do t = 1, ntarget
        nsubs_tot_t(t) = sum(nsubs_t(t, 1:nk(t)))
      end do
      allocate(degen(nic))
      degen = degen_rd(1:nic)
      allocate(newconf(max(1, nsubs_tot)))
      allocate(indconf(nic, max(1, nsubs_tot)))
      if (nsubs_tot > 0) indconf(1:nic, 1:nsubs_tot) = indconf_rd(1:nic, 1:nsubs_tot)
    end block
    ndigits = max(1, int(log10(real(nic))) + 1)

    ! Resolve a -choose selection label now that the level directory and the
    ! configuration count are known.  Each label names a file in that same
    ! directory, so the answer always belongs to the ENSEMBLE just loaded --
    ! including when SOD_ENSEMBLE_DIR points at nXX/random/.  The optional count
    ! takes the best choose_count of them, in rank order.
    if (len_trim(choose_label) > 0) then
      block
        integer :: iu_sel, rank_read, idx_read, nread, kk, jbest, nfound, ndup
        integer :: itmp
        real(real64) :: rnd
        integer, allocatable :: picked(:)
        real(real64), allocatable :: cand_e(:)
        logical, allocatable :: taken(:), have_e(:)
        real(real64) :: e_read, e_best
        character(len=256) :: line, selfile
        character(len=16) :: label_shown
        ! Echo the label in its canonical spelling, not as the user cased it.
        select case (trim(choose_label))
        case ('bestsqs');       label_shown = 'bestSQS'
        case ('lowestenergy');  label_shown = 'lowestENERGY'
        case ('random');        label_shown = 'random'
        case default;           label_shown = 'highestENERGY'
        end select

        allocate (picked(choose_count))
        picked = 0
        nfound = 0

        if (choose_label == 'random') then
          ! A random sample of the ENSEMBLE, weighted by degeneracy, drawn
          ! without replacement.
          !
          ! Weighted, because the rows of ENSEMBLE are symmetry-inequivalent
          ! configurations standing for degen(i) microstates each.  Sampling the
          ! rows uniformly would over-represent the rare ones, and a property
          ! averaged over such a subset would not estimate the ensemble average.
          ! Drawing with probability proportional to degeneracy makes the subset
          ! a fair sample of configuration space, so a plain mean over it is an
          ! unbiased estimate.
          !
          ! Without replacement, because each selection becomes a directory
          ! cNNN: a repeat would collapse two selections into one and silently
          ! yield fewer structures than asked for.
          !
          ! Efraimidis-Spirakis exponential race: give configuration i the key
          ! -ln(u_i)/degen(i) with u_i uniform on (0,1], and take the N smallest.
          ! That draws exactly the weighted-without-replacement sample, in one
          ! pass, with no rejection loop that could stall as N approaches nic.
          selfile = 'ENSEMBLE'
          if (choose_count > nic) then
            write (*, '(A,I0,A,I0,A)') " Error: ", choose_count, &
                " configurations requested but the ENSEMBLE has only ", nic, "."
            stop 1
          end if
          call random_seed()
          allocate (cand_e(nic))
          do i = 1, nic
            if (degen(i) <= 0) then
              ! Cannot represent any microstate; never selected.
              cand_e(i) = huge(1.0_real64)
              cycle
            end if
            call random_number(rnd)
            ! random_number can return exactly 0; -ln(0) is not finite.
            if (rnd <= 0.0_real64) rnd = tiny(1.0_real64)
            cand_e(i) = -log(rnd) / real(degen(i), real64)
          end do
          allocate (taken(nic))
          taken = .false.
          do kk = 1, choose_count
            jbest = 0
            do i = 1, nic
              if (taken(i)) cycle
              if (jbest == 0) then
                jbest = i
              else if (cand_e(i) < cand_e(jbest)) then
                jbest = i
              end if
            end do
            taken(jbest) = .true.
            picked(kk) = jbest
          end do
          deallocate (cand_e, taken)
          ! No meaningful order among random draws; ascending reads best.
          do kk = 2, choose_count
            itmp = picked(kk)
            jbest = kk - 1
            do while (jbest >= 1)
              if (picked(jbest) <= itmp) exit
              picked(jbest + 1) = picked(jbest)
              jbest = jbest - 1
            end do
            picked(jbest + 1) = itmp
          end do
          nfound = choose_count

        else if (choose_label == 'bestsqs') then
          ! OUTSQS is written sorted by ascending Q with ties broken on
          ! configuration index, so its data rows are already ranks 1, 2, 3...
          ! The SECOND column is the configuration number; the first is the
          ! rank, and confusing the two is what this label exists to prevent.
          selfile = trim(ndir_gulp)//'/OUTSQS'
          open (newunit=iu_sel, file=trim(selfile), status='old', iostat=ios)
          if (ios /= 0) then
            write (*, '(3A)') " Error: -choose bestSQS requires ", trim(selfile), &
                ", which was not found."
            write (*, '(A)') " Run sod_sqs.sh in that directory first."
            stop 1
          end if
          do
            if (nfound >= choose_count) exit
            read (iu_sel, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle
            read (line, *, iostat=ios) rank_read, idx_read
            if (ios /= 0) cycle
            if (idx_read < 1) cycle
            nfound = nfound + 1
            picked(nfound) = idx_read
          end do
          close (iu_sel)
          if (nfound < choose_count) then
            write (*, '(A,I0,3A,I0,A)') " Error: ", choose_count, " configurations requested but ", &
                trim(selfile), " lists only ", nfound, "."
            if (nfound > 0) then
              write (*, '(A)') " sqssod writes the top n_top_sqs rows only (default 10)."
              write (*, '(A)') " Raise n_top_sqs in INSQS, or set it to 0 to rank every configuration."
            end if
            stop 1
          end if

        else
          ! ENERGIES holds one "index energy" pair per configuration, the same
          ! file statsod averages over.
          selfile = trim(ndir_gulp)//'/ENERGIES'
          open (newunit=iu_sel, file=trim(selfile), status='old', iostat=ios)
          if (ios /= 0) then
            write (*, '(4A)') " Error: -choose ", trim(label_shown), " requires ", &
                trim(selfile)//", which was not found."
            write (*, '(A)') " Collect energies first (see the sod_*_ener.sh scripts)."
            stop 1
          end if
          ! Store by configuration index, not by line number.  ENERGIES is
          ! written by the sod_*_ener.sh scripts and edited by hand, so it can
          ! carry repeated or partial entries; counting lines into a
          ! nic-sized buffer would run off the end of it.
          allocate (cand_e(nic), have_e(nic), taken(nic))
          have_e = .false.
          taken = .false.
          cand_e = 0.0_real64
          ndup = 0
          do
            read (iu_sel, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle
            read (line, *, iostat=ios) idx_read, e_read
            if (ios /= 0) cycle
            if (idx_read < 1 .or. idx_read > nic) then
              write (*, '(A,I0,A,I0,A)') " Error: ENERGIES names configuration ", idx_read, &
                  " which is outside 1..", nic, " for this level."
              stop 1
            end if
            if (have_e(idx_read)) ndup = ndup + 1
            have_e(idx_read) = .true.
            cand_e(idx_read) = e_read       ! a repeated index keeps its last value
          end do
          close (iu_sel)
          nread = count(have_e)
          if (ndup > 0) then
            write (*, '(A,I0,3A)') " Warning: ", ndup, " repeated configuration index(es) in ", &
                trim(selfile), "; the last value of each was used."
          end if
          if (nread == 0) then
            write (*, '(3A)') " Error: no usable entries in ", trim(selfile), "."
            stop 1
          end if
          if (nread /= nic) then
            write (*, '(A,I0,A,I0,A)') " Warning: ENERGIES covers ", nread, &
                " of the ", nic, " configurations in the ENSEMBLE."
          end if
          if (choose_count > nread) then
            write (*, '(A,I0,3A,I0,A)') " Error: ", choose_count, " configurations requested but ", &
                trim(selfile), " has energies for only ", nread, "."
            stop 1
          end if
          ! Repeated extremum selection: choose_count is small in practice.
          do kk = 1, choose_count
            jbest = 0
            do i = 1, nic
              if (.not. have_e(i)) cycle
              if (taken(i)) cycle
              if (jbest == 0) then
                jbest = i
              else if (choose_label == 'lowestenergy') then
                if (cand_e(i) < cand_e(jbest)) jbest = i
              else
                if (cand_e(i) > cand_e(jbest)) jbest = i
              end if
            end do
            taken(jbest) = .true.
            picked(kk) = jbest
            if (kk == 1) e_best = cand_e(jbest)
          end do
          nfound = choose_count
          write (*, '(3A,ES16.8)') " > -choose ", trim(label_shown), ": best energy ", e_best
          deallocate (cand_e, have_e, taken)
        end if

        if (choose_label == 'random') then
          write (*, '(3A,I0,A,I0,A)') " > -choose ", trim(label_shown), ": ", choose_count, &
              " of ", nic, " configurations, drawn without replacement"
          write (*, '(A)') "   with probability proportional to degeneracy:"
          write (*, '(A,*(1X,I0))') "  ", (picked(kk), kk = 1, choose_count)
          write (*, '(A)') "   (repeat this exact selection with -choose followed by those indices)"
        else if (choose_count == 1) then
          write (*, '(3A,I0,3A)') " > -choose ", trim(label_shown), " selects configuration ", &
              picked(1), " (from ", trim(selfile), ")"
        else
          write (*, '(3A,I0,3A)') " > -choose ", trim(label_shown), " selects the best ", &
              choose_count, " configurations (from ", trim(selfile), "):"
          write (*, '(A,*(1X,I0))') "  ", (picked(kk), kk = 1, choose_count)
        end if
        if (allocated(chosen_indices)) deallocate (chosen_indices)
        n_chosen = choose_count
        allocate (chosen_indices(n_chosen))
        chosen_indices(1:n_chosen) = picked(1:n_chosen)
        deallocate (picked)
      end block
    end if

    if (mode == 2) then
      do i = 1, n_chosen
        if (chosen_indices(i) > nic) then
          write (*, '(A,I0,A,I0,A)') "Error: -choose index ", chosen_indices(i), &
              " exceeds number of configurations (", nic, ")."
          stop 1
        end if
      end do
      gen_count = n_chosen
      allocate(gen_indices(gen_count))
      gen_indices(1:gen_count) = chosen_indices(1:gen_count)
      write (*, '(A,I0,A,I0,A)') " > Generating ", gen_count, " selected configurations of ", nic, "."
    else
      gen_limit = nic
      if (mode == 1) gen_limit = min(nic, max_generate)
      gen_count = gen_limit
      allocate(gen_indices(gen_count))
      do ic = 1, gen_count
        gen_indices(ic) = ic
      end do
      if (gen_count < nic) then
        write (*, '(A,I0,A,I0,A)') " > Generating first ", gen_count, " of ", nic, " configurations."
      end if
    end if

    block
      integer :: d, probe_status
      character(len=20) :: probe_fmt, probe_str
      do d = 1, 5
        if (d == ndigits) cycle
        write (probe_fmt, '(a,i0,a,i0,a)') '(i', d, '.', d, ')'
        write (probe_str, probe_fmt) 1
        call execute_command_line('test -d "' // trim(ndir_gulp) // '/c' // trim(probe_str) // '"', &
                                  wait=.true., exitstat=probe_status)
        if (probe_status == 0) then
          write (*, '(A)') " Warning: " // trim(ndir_gulp) // &
            " contains cXXX directories with different zero-padding (stale from a previous run)." // &
            " Remove them to avoid confusion."
          exit
        end if
      end do
    end block

    select case (filer)

  case (-1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!  CIF (P1, one file per config) !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case (0)

    write (numfmt, '(a,i0,a,i0,a)') '(i', ndigits, '.', ndigits, ')'

    call execute_command_line('mkdir -p ' // trim(ndir_gulp), exitstat=ios)
    if (ios /= 0) then
      write (*, *) "Error: could not create directory ", trim(ndir_gulp)
      stop 1
    end if

    do ic = 1, gen_count
      indcount = gen_indices(ic)
      newconf(1:nsubs_tot) = indconf(indcount, 1:nsubs_tot)

      write (numstr, numfmt) indcount
      confdir_gulp = trim(ndir_gulp) // '/c' // trim(numstr)
      call execute_command_line('mkdir -p ' // trim(confdir_gulp), exitstat=ios)
      if (ios /= 0) then
        write (*, *) "Error: could not create directory ", trim(confdir_gulp)
        stop 1
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

      nat_exp = 0
      do at = 1, nat
        sp = spat(at)
        if (sp_is_mol(sp)) then
          ! Parent-structure molecule: expand placeholder into rigid molecule
          call mol_rotate_frac(sp_mol_idx(sp), coords(at,1), coords(at,2), coords(at,3), &
                               cellvector, mol_frac_buf, mol_natoms(sp_mol_idx(sp)))
          do im = 1, mol_natoms(sp_mol_idx(sp))
            nat_exp = nat_exp + 1
            write (72, '(a, i0, 2x, a, 3(2x, f11.7))') &
                  trim(mol_sym_arr(sp_mol_idx(sp), im)), nat_exp, &
                  trim(mol_sym_arr(sp_mol_idx(sp), im)), &
                  mol_frac_buf(im, 1), mol_frac_buf(im, 2), mol_frac_buf(im, 3)
          end do
          cycle
        end if
        if (.not. any(sp == sptarget(1:ntarget))) then
          nat_exp = nat_exp + 1
          write (72, '(a, i0, 2x, a, 3(2x, f11.7))') &
                trim(symbol(sp)), nat_exp, trim(symbol(sp)), &
                coords(at, 1), coords(at, 2), coords(at, 3)
        else
          call resolve_target_atom(at, sp, site_action, site_sym, site_mol)
          if (site_action == site_molecule) then
            call mol_rotate_frac(site_mol, coords(at,1), coords(at,2), coords(at,3), &
                                 cellvector, mol_frac_buf, mol_natoms(site_mol))
            do im = 1, mol_natoms(site_mol)
              nat_exp = nat_exp + 1
              write (72, '(a, i0, 2x, a, 3(2x, f11.7))') &
                    trim(mol_sym_arr(site_mol, im)), nat_exp, &
                    trim(mol_sym_arr(site_mol, im)), &
                    mol_frac_buf(im, 1), mol_frac_buf(im, 2), mol_frac_buf(im, 3)
            end do
          else if (site_action == site_plain) then
            nat_exp = nat_exp + 1
            write (72, '(a, i0, 2x, a, 3(2x, f11.7))') &
                  trim(site_sym), nat_exp, trim(site_sym), &
                  coords(at, 1), coords(at, 2), coords(at, 3)
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
      stop 1
    end if
    ngulplines = 0
    do
      read (70, '(a)', end=1001) gulpline_buf
      ngulplines = ngulplines + 1
      if (ngulplines > maxgulplines) then
        write (*, *) "Error: template_input.gin exceeds ", maxgulplines, " lines. Aborting."
        stop 1
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
          stop 1
        end if
        struc_line_idx = l
      else if (index(gulptemplate(l), "@configuration_structure@") /= 0) then
        write (*, *) "Error: @configuration_structure@ is not alone on its line. Aborting."
        stop 1
      end if
    end do
    if (struc_line_idx == 0) then
      write (*, *) "Error: @configuration_structure@ not found in template_input.gin. Aborting."
      stop 1
    end if

    ! --- Build all_sym array ---
    all_nsp_gulp = 0
    do ssp = 1, nsp
      t = 0
      do i = 1, ntarget
        if (ssp == sptarget(i)) then
          t = i
          exit
        end if
      end do
      if (t == 0) then
        ! Not a substitution target: the parent species (or its molecule)
        if (sp_is_mol(ssp)) then
          call add_mol_types(sp_mol_idx(ssp))
        else
          call add_unique_sym(symbol(ssp))
        end if
      else
        ! Target species: every slot this site can carry, including the
        ! remaining-species slot nk(t)+1.  Vacancies contribute no type.
        do j = 1, nk(t) + 1
          if (is_vac_t(t,j)) cycle
          if (is_mol_t(t,j)) then
            call add_mol_types(mol_idx_t(t,j))
          else
            call add_unique_sym(newsymbol(t,j)(1:3))
          end if
        end do
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
          stop 1
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

    ! --- Detect kim_model directive (OpenKIM potential, no library file needed) ---
    haskimmodel_gulp = .false.
    do l = 1, ngulplines
      sline_trimmed = adjustl(gulptemplate(l))
      if (sline_trimmed(1:9) == "kim_model") then
        haskimmodel_gulp = .true.
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
        stop 1
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
    if (.not. hasinline_ff .and. .not. haslibrary_gulp .and. .not. haskimmodel_gulp) then
      write (*, *) "Error: template_input.gin has no inline force-field or library information. Aborting."
      stop 1
    end if

    ! --- Build number-only format string ---
    write (numfmt, '(a,i0,a,i0,a)') '(i', ndigits, '.', ndigits, ')'

    ! --- Build nXX parent directory and create it ---
    call execute_command_line('mkdir -p ' // trim(ndir_gulp), exitstat=ios)
    if (ios /= 0) then
      write (*, *) "Error: could not create directory ", trim(ndir_gulp)
      stop 1
    end if

    ! --- Generate one configuration directory per configuration ---
    do ic = 1, gen_count
      indcount = gen_indices(ic)
      newconf(1:nsubs_tot) = indconf(indcount, 1:nsubs_tot)

      write (numstr, numfmt) indcount
      confdir_gulp = trim(ndir_gulp) // '/c' // trim(numstr)

      call execute_command_line('mkdir -p ' // trim(confdir_gulp), exitstat=ios)
      if (ios /= 0) then
        write (*, *) "Error: could not create directory ", trim(confdir_gulp)
        stop 1
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
            if (sp_is_mol(sp)) then
              ! Parent-structure molecule: expand placeholder into rigid molecule
              call mol_rotate_frac(sp_mol_idx(sp), coords(at,1), coords(at,2), coords(at,3), &
                                   cellvector, mol_frac_buf, mol_natoms(sp_mol_idx(sp)))
              do im = 1, mol_natoms(sp_mol_idx(sp))
                gtype_buf = mol_sym_arr(sp_mol_idx(sp), im)
                do m = 1, all_nsp_gulp
                  if (trim(all_sym(m)) == trim(mol_sym_arr(sp_mol_idx(sp), im))) then
                    gtype_buf = all_gulptype_arr(m)
                    exit
                  end if
                end do
                write (72, 331) trim(gtype_buf), "core", &
                                mol_frac_buf(im,1), mol_frac_buf(im,2), mol_frac_buf(im,3)
                do m = 1, all_nsp_gulp
                  if (trim(all_gulptype_arr(m)) == trim(gtype_buf) .and. all_ishell_arr(m) == 1) then
                    write (72, 331) trim(gtype_buf), "shel", &
                                    mol_frac_buf(im,1), mol_frac_buf(im,2), mol_frac_buf(im,3)
                    exit
                  end if
                end do
              end do
              cycle
            end if
            if (.not. any(sp == sptarget(1:ntarget))) then
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
              call resolve_target_atom(at, sp, site_action, site_sym, site_mol)
              if (site_action == site_molecule) then
                call mol_rotate_frac(site_mol, coords(at,1), coords(at,2), coords(at,3), &
                                     cellvector, mol_frac_buf, mol_natoms(site_mol))
                do im = 1, mol_natoms(site_mol)
                  gtype_buf = mol_sym_arr(site_mol, im)
                  do m = 1, all_nsp_gulp
                    if (trim(all_sym(m)) == trim(mol_sym_arr(site_mol, im))) then
                      gtype_buf = all_gulptype_arr(m)
                      exit
                    end if
                  end do
                  write (72, 331) trim(gtype_buf), "core", &
                                  mol_frac_buf(im,1), mol_frac_buf(im,2), mol_frac_buf(im,3)
                  do m = 1, all_nsp_gulp
                    if (trim(all_gulptype_arr(m)) == trim(gtype_buf) .and. all_ishell_arr(m) == 1) then
                      write (72, 331) trim(gtype_buf), "shel", &
                                      mol_frac_buf(im,1), mol_frac_buf(im,2), mol_frac_buf(im,3)
                      exit
                    end if
                  end do
                end do
              else if (site_action == site_plain) then
                gtype_buf = site_sym(1:3)
                do m = 1, all_nsp_gulp
                  if (trim(all_sym(m)) == trim(site_sym(1:3))) then
                    gtype_buf = all_gulptype_arr(m); exit
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
      stop 1
    end if
    ngulplines = 0
    do
      read (70, '(a)', end=2001) gulpline_buf
      ngulplines = ngulplines + 1
      if (ngulplines > maxgulplines) then
        write (*, *) "Error: template_in.lammps exceeds ", maxgulplines, " lines. Aborting."
        stop 1
      end if
      gulptemplate(ngulplines) = gulpline_buf
    end do
2001 continue
    close (70)

    ! --- Build all_sym array (molecule-aware) ---
    all_nsp_gulp = 0
    do ssp = 1, nsp
      t = 0
      do i = 1, ntarget
        if (ssp == sptarget(i)) then
          t = i
          exit
        end if
      end do
      if (t == 0) then
        ! Not a substitution target: the parent species (or its molecule)
        if (sp_is_mol(ssp)) then
          call add_mol_types(sp_mol_idx(ssp))
        else
          call add_unique_sym(symbol(ssp))
        end if
      else
        ! Target species: every slot this site can carry, including the
        ! remaining-species slot nk(t)+1.  Vacancies contribute no type.
        do j = 1, nk(t) + 1
          if (is_vac_t(t,j)) cycle
          if (is_mol_t(t,j)) then
            call add_mol_types(mol_idx_t(t,j))
          else
            call add_unique_sym(newsymbol(t,j)(1:3))
          end if
        end do
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
      stop 1
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
        stop 1
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
          stop 1
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
        stop 1
      end if
      read (lmptmpline(lmpidx+5:), *, iostat=ios) lmptypeval
      if (ios /= 0) then
        write (*, *) "Error: invalid type= in: ", trim(sline_trimmed)
        stop 1
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
        stop 1
      end if
      if (lammps_shell_type(i) > 0 .and. lammps_shell_bond(i) == 0) then
        write (*, *) "Error: SOD species '", trim(all_sym(i)), &
                     "' has shell mapping but no bond_type. Aborting."
        stop 1
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
      stop 1
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
    call execute_command_line('mkdir -p ' // trim(ndir_gulp), exitstat=ios)
    if (ios /= 0) then
      write (*, *) "Error: could not create directory ", trim(ndir_gulp)
      stop 1
    end if

    ! --- Generate one configuration directory per configuration ---
    do ic = 1, gen_count
      indcount = gen_indices(ic)
      newconf(1:nsubs_tot) = indconf(indcount, 1:nsubs_tot)

      write (numstr, numfmt) indcount
      confdir_gulp = trim(ndir_gulp) // '/c' // trim(numstr)

      call execute_command_line('mkdir -p ' // trim(confdir_gulp), exitstat=ios)
      if (ios /= 0) then
        write (*, *) "Error: could not create directory ", trim(confdir_gulp)
        stop 1
      end if

      ! --- Build expanded atom list (handles molecules, vacancies, shells) ---
      nat_exp = 0
      atom_id_lammps = 0
      mol_id_lammps = 0

      do at = 1, nat
        sp = spat(at)
        if (sp_is_mol(sp)) then
          ! Parent-structure molecule: expand placeholder into rigid molecule
          call mol_rotate_frac(sp_mol_idx(sp), coords(at,1), coords(at,2), coords(at,3), &
                               cellvector, mol_frac_buf, mol_natoms(sp_mol_idx(sp)))
          mol_id_lammps = mol_id_lammps + 1
          do im = 1, mol_natoms(sp_mol_idx(sp))
            nat_exp = nat_exp + 1
            sym_exp(nat_exp) = mol_sym_arr(sp_mol_idx(sp), im)
            coords_exp(nat_exp, 1:3) = mol_frac_buf(im, 1:3)
            lammps_sym_idx = 0
            do lmp_i = 1, all_nsp_gulp
              if (trim(all_sym(lmp_i)) == trim(mol_sym_arr(sp_mol_idx(sp), im))) then
                lammps_sym_idx = lmp_i; exit
              end if
            end do
            exp_sym_idx(nat_exp) = lammps_sym_idx
            exp_mol_id(nat_exp) = mol_id_lammps
            exp_has_shell(nat_exp) = .false.
            exp_shell_id(nat_exp) = 0
            atom_id_lammps = atom_id_lammps + 1
            exp_core_id(nat_exp) = atom_id_lammps
          end do
          cycle
        end if
        if (.not. any(sp == sptarget(1:ntarget))) then
          ! Framework atom
          lammps_sym_cur = symbol(sp)
        else
          call resolve_target_atom(at, sp, site_action, site_sym, site_mol)
          if (site_action == site_vacancy) cycle
          if (site_action == site_molecule) then
            call mol_rotate_frac(site_mol, coords(at,1), coords(at,2), coords(at,3), &
                                 cellvector, mol_frac_buf, mol_natoms(site_mol))
            mol_id_lammps = mol_id_lammps + 1
            do im = 1, mol_natoms(site_mol)
              nat_exp = nat_exp + 1
              sym_exp(nat_exp) = mol_sym_arr(site_mol, im)
              coords_exp(nat_exp, 1:3) = mol_frac_buf(im, 1:3)
              lammps_sym_idx = 0
              do lmp_i = 1, all_nsp_gulp
                if (trim(all_sym(lmp_i)) == trim(mol_sym_arr(site_mol, im))) then
                  lammps_sym_idx = lmp_i; exit
                end if
              end do
              exp_sym_idx(nat_exp) = lammps_sym_idx
              exp_mol_id(nat_exp) = mol_id_lammps
              exp_has_shell(nat_exp) = .false.
              exp_shell_id(nat_exp) = 0
              atom_id_lammps = atom_id_lammps + 1
              exp_core_id(nat_exp) = atom_id_lammps
            end do
            cycle
          end if
          lammps_sym_cur = site_sym(1:3)
        end if
        ! Regular (non-molecule) atom: append to expanded list
        nat_exp = nat_exp + 1
        sym_exp(nat_exp) = lammps_sym_cur
        coords_exp(nat_exp, 1:3) = coords(at, 1:3)
        lammps_sym_idx = 0
        do lmp_i = 1, all_nsp_gulp
          if (trim(all_sym(lmp_i)) == trim(lammps_sym_cur)) then
            lammps_sym_idx = lmp_i; exit
          end if
        end do
        exp_sym_idx(nat_exp) = lammps_sym_idx
        atom_id_lammps = atom_id_lammps + 1
        exp_core_id(nat_exp) = atom_id_lammps
        if (lammps_shell_type(lammps_sym_idx) > 0) then
          mol_id_lammps = mol_id_lammps + 1
          exp_mol_id(nat_exp) = mol_id_lammps
          exp_has_shell(nat_exp) = .true.
          atom_id_lammps = atom_id_lammps + 1
          exp_shell_id(nat_exp) = atom_id_lammps
        else
          exp_mol_id(nat_exp) = 0
          exp_has_shell(nat_exp) = .false.
          exp_shell_id(nat_exp) = 0
        end if
      end do
      natoms_lammps = atom_id_lammps
      nbonds_lammps = 0
      do at = 1, nat_exp
        if (exp_has_shell(at)) nbonds_lammps = nbonds_lammps + 1
      end do

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

      do at = 1, nat_exp
        lammps_sym_idx = exp_sym_idx(at)
        xc = cellvector(1,1)*coords_exp(at,1) + cellvector(1,2)*coords_exp(at,2) + cellvector(1,3)*coords_exp(at,3)
        yc = cellvector(2,2)*coords_exp(at,2) + cellvector(2,3)*coords_exp(at,3)
        zc = cellvector(3,3)*coords_exp(at,3)
        if (trim(lammps_atom_style) == "atomic") then
          write (72, '(i0, 2x, i0, 3(2x, f14.6))') &
                exp_core_id(at), lammps_core_type(lammps_sym_idx), xc, yc, zc
        else if (trim(lammps_atom_style) == "charge") then
          write (72, '(i0, 2x, i0, 2x, f8.4, 3(2x, f14.6))') &
                exp_core_id(at), lammps_core_type(lammps_sym_idx), 0.0, xc, yc, zc
        else  ! full
          write (72, '(i0, 2x, i0, 2x, i0, 2x, f8.4, 3(2x, f14.6))') &
                exp_core_id(at), exp_mol_id(at), lammps_core_type(lammps_sym_idx), 0.0, xc, yc, zc
          if (exp_has_shell(at)) then
            write (72, '(i0, 2x, i0, 2x, i0, 2x, f8.4, 3(2x, f14.6))') &
                  exp_shell_id(at), exp_mol_id(at), lammps_shell_type(lammps_sym_idx), 0.0, xc, yc, zc
          end if
        end if
      end do

      if (nbonds_lammps > 0) then
        write (72, '(a)') " "
        write (72, '(a)') "Bonds"
        write (72, '(a)') " "
        bond_id_lammps = 0
        do at = 1, nat_exp
          if (exp_has_shell(at)) then
            lammps_sym_idx = exp_sym_idx(at)
            bond_id_lammps = bond_id_lammps + 1
            write (72, '(i0, 2x, i0, 2x, i0, 2x, i0)') &
                  bond_id_lammps, lammps_shell_bond(lammps_sym_idx), &
                  exp_core_id(at), exp_shell_id(at)
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

    ! Check for INCAR file (recommended for VASP, but not required)
    inquire(file='INCAR', exist=fileexists)
    if (.not. fileexists) then
      write (*, *) 'Warning: INCAR file not found in current directory.'
      write (*, *) 'VASP input files will be created without an INCAR; place one'
      write (*, *) 'in each c* directory before running VASP.'
    end if

    write (numfmt, '(a,i0,a,i0,a)') '(i', ndigits, '.', ndigits, ')'

    call execute_command_line('mkdir -p ' // trim(ndir_gulp), exitstat=ios)
    if (ios /= 0) then
      write (*, *) "Error: could not create directory ", trim(ndir_gulp)
      stop 1
    end if

    do ic = 1, gen_count
      indcount = gen_indices(ic)
      newconf(1:nsubs_tot) = indconf(indcount, 1:nsubs_tot)

      write (numstr, numfmt) indcount
      confdir_gulp = trim(ndir_gulp) // '/c' // trim(numstr)
      call execute_command_line('mkdir -p ' // trim(confdir_gulp), exitstat=ios)
      if (ios /= 0) then
        write (*, *) "Error: could not create directory ", trim(confdir_gulp)
        stop 1
      end if

      ! Copy INCAR to configuration directory (only if present)
      if (fileexists) then
        call execute_command_line('cp INCAR ' // trim(confdir_gulp) // '/', exitstat=ios)
        if (ios /= 0) then
          write (*, *) "Warning: could not copy INCAR to ", trim(confdir_gulp)
        end if
      end if

      inpfile_gulp = trim(confdir_gulp) // '/POSCAR'
      open (unit=72, file=trim(inpfile_gulp), status="replace")

      write (72, *) title
      write (72, *) '1.00000000'
335   format(3(f10.6, 2x))
      write (72, 335) cellvector(1, 1), cellvector(2, 1), cellvector(3, 1)
      write (72, 335) cellvector(1, 2), cellvector(2, 2), cellvector(3, 2)
      write (72, 335) cellvector(1, 3), cellvector(2, 3), cellvector(3, 3)

      ! --- Two-pass expansion: build sym_exp / coords_exp ---
      nat_exp = 0
      do at = 1, nat
        sp = spat(at)
        if (sp_is_mol(sp)) then
          ! Parent-structure molecule: expand placeholder into rigid molecule
          call mol_rotate_frac(sp_mol_idx(sp), coords(at,1), coords(at,2), coords(at,3), &
                               cellvector, mol_frac_buf, mol_natoms(sp_mol_idx(sp)))
          do im = 1, mol_natoms(sp_mol_idx(sp))
            nat_exp = nat_exp + 1
            sym_exp(nat_exp) = mol_sym_arr(sp_mol_idx(sp), im)
            coords_exp(nat_exp, 1:3) = mol_frac_buf(im, 1:3)
          end do
          cycle
        end if
        if (.not. any(sp == sptarget(1:ntarget))) then
          nat_exp = nat_exp + 1
          sym_exp(nat_exp) = symbol(sp)
          coords_exp(nat_exp, 1:3) = coords(at, 1:3)
        else
          call resolve_target_atom(at, sp, site_action, site_sym, site_mol)
          if (site_action == site_molecule) then
            call mol_rotate_frac(site_mol, coords(at,1), coords(at,2), coords(at,3), &
                                 cellvector, mol_frac_buf, mol_natoms(site_mol))
            do im = 1, mol_natoms(site_mol)
              nat_exp = nat_exp + 1
              sym_exp(nat_exp) = mol_sym_arr(site_mol, im)
              coords_exp(nat_exp, 1:3) = mol_frac_buf(im, 1:3)
            end do
          else if (site_action == site_plain) then
            nat_exp = nat_exp + 1
            sym_exp(nat_exp) = site_sym(1:3)
            coords_exp(nat_exp, 1:3) = coords(at, 1:3)
          end if
        end if
      end do

      ! --- Collect unique species (order of first appearance) ---
      nsp_exp = 0
      uniq_cnt_exp = 0
      do at = 1, nat_exp
        isp_exp = 0
        do i = 1, nsp_exp
          if (trim(sym_exp(at)) == trim(uniq_sym_exp(i))) then
            isp_exp = i
            exit
          end if
        end do
        if (isp_exp == 0) then
          nsp_exp = nsp_exp + 1
          uniq_sym_exp(nsp_exp) = sym_exp(at)
          uniq_cnt_exp(nsp_exp) = 1
        else
          uniq_cnt_exp(isp_exp) = uniq_cnt_exp(isp_exp) + 1
        end if
      end do

      ! --- Write POSCAR ---
      outbuf_gulp = trim(uniq_sym_exp(1))
      do i = 2, nsp_exp
        outbuf_gulp = trim(outbuf_gulp) // ' ' // trim(uniq_sym_exp(i))
      end do
      write (72, '(a)') trim(outbuf_gulp)
336   format(10(i4, 1x))
      write (72, 336) (uniq_cnt_exp(i), i=1, nsp_exp)
      write (72, *) 'Direct'
      do isp_exp = 1, nsp_exp
        do at = 1, nat_exp
          if (trim(sym_exp(at)) == trim(uniq_sym_exp(isp_exp))) then
            write (72, 335) coords_exp(at, 1), coords_exp(at, 2), coords_exp(at, 3)
          end if
        end do
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
      stop 1
    end if
    ngulplines = 0
    do
      read (70, '(a)', end=1201) gulpline_buf
      ngulplines = ngulplines + 1
      if (ngulplines > maxgulplines) then
        write (*, *) "Error: template_castep.cell exceeds ", maxgulplines, " lines. Aborting."
        stop 1
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
          stop 1
        end if
        struc_line_idx = l
      else if (index(gulptemplate(l), "@configuration_structure@") /= 0) then
        write (*, *) "Error: @configuration_structure@ is not alone on its line. Aborting."
        stop 1
      end if
    end do
    if (struc_line_idx == 0) then
      write (*, *) "Error: @configuration_structure@ not found in template_castep.cell. Aborting."
      stop 1
    end if

    call cell(cellvector, a, b, c, alpha, beta, gamma)

    ! --- Build number format ---
    write (numfmt, '(a,i0,a,i0,a)') '(i', ndigits, '.', ndigits, ')'

    ! --- Build nXX parent directory and create it ---
    call execute_command_line('mkdir -p ' // trim(ndir_gulp), exitstat=ios)
    if (ios /= 0) then
      write (*, *) "Error: could not create directory ", trim(ndir_gulp)
      stop 1
    end if

    do ic = 1, gen_count
      indcount = gen_indices(ic)
      newconf(1:nsubs_tot) = indconf(indcount, 1:nsubs_tot)

      write (numstr, numfmt) indcount
      confdir_gulp = trim(ndir_gulp) // '/c' // trim(numstr)

      call execute_command_line('mkdir -p ' // trim(confdir_gulp), exitstat=ios)
      if (ios /= 0) then
        write (*, *) "Error: could not create directory ", trim(confdir_gulp)
        stop 1
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

337       format(a3, 3(f11.7, 2x))
          do at = 1, nat
            sp = spat(at)
            if (sp_is_mol(sp)) then
              ! Parent-structure molecule: expand placeholder into rigid molecule
              call mol_rotate_frac(sp_mol_idx(sp), coords(at,1), coords(at,2), coords(at,3), &
                                   cellvector, mol_frac_buf, mol_natoms(sp_mol_idx(sp)))
              do im = 1, mol_natoms(sp_mol_idx(sp))
                write (72, 337) mol_sym_arr(sp_mol_idx(sp), im), &
                                mol_frac_buf(im,1), mol_frac_buf(im,2), mol_frac_buf(im,3)
              end do
              cycle
            end if
            if (.not. any(sp == sptarget(1:ntarget))) then
              write (72, 337) symbol(sp), coords(at, 1), coords(at, 2), coords(at, 3)
            else
              call resolve_target_atom(at, sp, site_action, site_sym, site_mol)
              if (site_action == site_molecule) then
                call mol_rotate_frac(site_mol, coords(at,1), coords(at,2), coords(at,3), &
                                     cellvector, mol_frac_buf, mol_natoms(site_mol))
                do im = 1, mol_natoms(site_mol)
                  write (72, 337) mol_sym_arr(site_mol, im), &
                                  mol_frac_buf(im,1), mol_frac_buf(im,2), mol_frac_buf(im,3)
                end do
              else if (site_action == site_plain) then
                write (72, 337) site_sym(1:3), coords(at, 1), coords(at, 2), coords(at, 3)
              end if
            end if
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
      stop 1
    end if
    ngulplines = 0
    do
      read (70, '(a)', end=1301) gulpline_buf
      ngulplines = ngulplines + 1
      if (ngulplines > maxgulplines) then
        write (*, *) "Error: template_pw.in exceeds ", maxgulplines, " lines. Aborting."
        stop 1
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
          stop 1
        end if
        struc_line_idx = l
      else if (index(gulptemplate(l), "@configuration_structure@") /= 0) then
        write (*, *) "Error: @configuration_structure@ is not alone on its line. Aborting."
        stop 1
      end if
    end do
    if (struc_line_idx == 0) then
      write (*, *) "Error: @configuration_structure@ not found in template_pw.in. Aborting."
      stop 1
    end if

    call cell(cellvector, a, b, c, alpha, beta, gamma)

    ! --- Build number format ---
    write (numfmt, '(a,i0,a,i0,a)') '(i', ndigits, '.', ndigits, ')'

    ! --- Build nXX parent directory and create it ---
    call execute_command_line('mkdir -p ' // trim(ndir_gulp), exitstat=ios)
    if (ios /= 0) then
      write (*, *) "Error: could not create directory ", trim(ndir_gulp)
      stop 1
    end if

    do ic = 1, gen_count
      indcount = gen_indices(ic)
      newconf(1:nsubs_tot) = indconf(indcount, 1:nsubs_tot)

      write (numstr, numfmt) indcount
      confdir_gulp = trim(ndir_gulp) // '/c' // trim(numstr)

      call execute_command_line('mkdir -p ' // trim(confdir_gulp), exitstat=ios)
      if (ios /= 0) then
        write (*, *) "Error: could not create directory ", trim(confdir_gulp)
        stop 1
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
            if (sp_is_mol(sp)) then
              ! Parent-structure molecule: expand placeholder into rigid molecule
              call mol_rotate_frac(sp_mol_idx(sp), coords(at,1), coords(at,2), coords(at,3), &
                                   cellvector, mol_frac_buf, mol_natoms(sp_mol_idx(sp)))
              do im = 1, mol_natoms(sp_mol_idx(sp))
                write (72, 337) mol_sym_arr(sp_mol_idx(sp), im), &
                                mol_frac_buf(im,1), mol_frac_buf(im,2), mol_frac_buf(im,3)
              end do
              cycle
            end if
            if (.not. any(sp == sptarget(1:ntarget))) then
              write (72, 337) symbol(sp), coords(at, 1), coords(at, 2), coords(at, 3)
            else
              call resolve_target_atom(at, sp, site_action, site_sym, site_mol)
              if (site_action == site_molecule) then
                call mol_rotate_frac(site_mol, coords(at,1), coords(at,2), coords(at,3), &
                                     cellvector, mol_frac_buf, mol_natoms(site_mol))
                do im = 1, mol_natoms(site_mol)
                  write (72, 337) mol_sym_arr(site_mol, im), &
                                  mol_frac_buf(im,1), mol_frac_buf(im,2), mol_frac_buf(im,3)
                end do
              else if (site_action == site_plain) then
                write (72, 337) site_sym(1:3), coords(at, 1), coords(at, 2), coords(at, 3)
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

    if (ntarget == 1 .and. nk(1) >= 2) exit  ! multi-nary: no range loop
    if (ntarget >= 2 .and. any(nk(1:ntarget) >= 2)) exit  ! Phase 4: no range loop

  end do  ! nsubs

  close (43)

!!!!!!!Reporting the end
  write (*, *) "Done!!!"
  write (*, *) ""
  write (*, *) ""

contains

  function npad(n) result(s)
    integer, intent(in) :: n
    character(len=10) :: s
    integer :: nd, tmp
    character(len=20) :: fmt

    tmp = max(1, abs(maxval(npos_t(1:ntarget))))
    nd = 1
    do while (tmp >= 10)
      nd = nd + 1
      tmp = tmp / 10
    end do
    nd = max(2, nd)

    write(fmt, '(A,I0,A,I0,A)') '(i', nd, '.', nd, ')'
    write(s, fmt) n
    s = adjustl(s)
  end function npad

!---------------------------------------------------------------------------
! Return the molecule-type index for NAME, loading NAME.xyz on first request.
  subroutine get_molecule(name, imt_out)
    character(len=*), intent(in) :: name
    integer, intent(out) :: imt_out
    integer :: m_loc, im_loc, ios_loc
    character(len=200) :: comment_loc
    ! Already loaded?
    do m_loc = 1, nmol_types
      if (trim(mol_name(m_loc)) == trim(name)) then
        imt_out = m_loc
        return
      end if
    end do
    ! Load a new molecule type
    nmol_types = nmol_types + 1
    if (nmol_types > nmoltypes) then
      write (*, *) "Error: too many molecule types (max ", nmoltypes, "). Aborting."
      stop 1
    end if
    imt_out = nmol_types
    mol_name(imt_out) = trim(name)
    open (unit=32, file=trim(mol_name(imt_out)) // '.xyz', status='old', iostat=ios_loc)
    if (ios_loc /= 0) then
      write (*, *) "Error: ", trim(mol_name(imt_out)), ".xyz not found. Aborting."
      stop 1
    end if
    read (32, *) mol_natoms(imt_out)
    if (mol_natoms(imt_out) > nmolmax_atoms) then
      write (*, *) "Error: ", trim(mol_name(imt_out)), ".xyz has more atoms than nmolmax_atoms =", &
                   nmolmax_atoms, ". Aborting."
      stop 1
    end if
    read (32, '(a)') comment_loc   ! comment line
    do im_loc = 1, mol_natoms(imt_out)
      read (32, *) mol_sym_arr(imt_out, im_loc), &
                   mol_xyz_arr(imt_out, im_loc, 1), mol_xyz_arr(imt_out, im_loc, 2), &
                   mol_xyz_arr(imt_out, im_loc, 3)
    end do
    close (32)
    ! Shift to CoM = origin
    call shift_to_com(imt_out)
    write (*, '(a,a,a,i0,a)') " Molecule @", trim(mol_name(imt_out)), &
                                " read from ", mol_natoms(imt_out), &
                                " atoms in " // trim(mol_name(imt_out)) // ".xyz"
  end subroutine get_molecule

!---------------------------------------------------------------------------
! Append the (unique) atom types of molecule MIDX to the all_sym registry.
  ! Append a symbol to the all_sym table unless it is already there.
  ! Lowercase a string, so -choose labels can be typed in any case.
  function to_lower(str) result(out)
    character(len=*), intent(in) :: str
    character(len=32) :: out
    integer :: k, code
    out = ''
    do k = 1, min(len_trim(str), len(out))
      code = iachar(str(k:k))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        out(k:k) = achar(code + 32)
      else
        out(k:k) = str(k:k)
      end if
    end do
  end function to_lower

  ! The -choose labels, listed in one place so the usage text cannot drift from
  ! the set the parser accepts.
  subroutine write_choose_usage()
    write (*, '(A)') " -choose takes either configuration indices or a selection label:"
    write (*, '(A)') "   -choose i1 i2 ...        those configurations"
    write (*, '(A)') "   -choose <label> [N]      the best N by that label (default 1)"
    write (*, '(A)') " Labels, case-insensitive:"
    write (*, '(A)') "   bestSQS         best-ranked in OUTSQS   (written by sod_sqs.sh)"
    write (*, '(A)') "   lowestENERGY    lowest in ENERGIES"
    write (*, '(A)') "   highestENERGY   highest in ENERGIES"
    write (*, '(A)') "   random          a degeneracy-weighted sample, no repeats"
    write (*, '(A)') " OUTSQS and ENERGIES are read from the directory holding the ENSEMBLE."
  end subroutine write_choose_usage

  subroutine add_unique_sym(symstr)
    character(len=*), intent(in) :: symstr
    integer :: m_loc
    do m_loc = 1, all_nsp_gulp
      if (trim(all_sym(m_loc)) == trim(symstr)) return
    end do
    all_nsp_gulp = all_nsp_gulp + 1
    all_sym(all_nsp_gulp) = symstr
  end subroutine add_unique_sym

  subroutine add_mol_types(midx)
    integer, intent(in) :: midx
    integer :: im_loc, m_loc
    logical :: found_loc
    do im_loc = 1, mol_natoms(midx)
      found_loc = .false.
      do m_loc = 1, all_nsp_gulp
        if (trim(all_sym(m_loc)) == trim(mol_sym_arr(midx, im_loc))) then
          found_loc = .true.; exit
        end if
      end do
      if (.not. found_loc) then
        all_nsp_gulp = all_nsp_gulp + 1
        all_sym(all_nsp_gulp) = mol_sym_arr(midx, im_loc)
      end if
    end do
  end subroutine add_mol_types

!---------------------------------------------------------------------------
  subroutine shift_to_com(imt)
    integer, intent(in) :: imt
    integer :: im_loc
    real(real64) :: com(3), totmass, amass
    com = 0.0
    totmass = 0.0
    do im_loc = 1, mol_natoms(imt)
      amass = get_mass(trim(mol_sym_arr(imt, im_loc)))
      com(1) = com(1) + amass * mol_xyz_arr(imt, im_loc, 1)
      com(2) = com(2) + amass * mol_xyz_arr(imt, im_loc, 2)
      com(3) = com(3) + amass * mol_xyz_arr(imt, im_loc, 3)
      totmass = totmass + amass
    end do
    com = com / totmass
    do im_loc = 1, mol_natoms(imt)
      mol_xyz_arr(imt, im_loc, 1) = mol_xyz_arr(imt, im_loc, 1) - com(1)
      mol_xyz_arr(imt, im_loc, 2) = mol_xyz_arr(imt, im_loc, 2) - com(2)
      mol_xyz_arr(imt, im_loc, 3) = mol_xyz_arr(imt, im_loc, 3) - com(3)
    end do
  end subroutine shift_to_com

!---------------------------------------------------------------------------
  function get_mass(sym) result(mass)
    character(len=*), intent(in) :: sym
    real(real64) :: mass
    select case (trim(sym))
    case ('H');  mass = 1.008
    case ('He'); mass = 4.003
    case ('Li'); mass = 6.941
    case ('Be'); mass = 9.012
    case ('B');  mass = 10.811
    case ('C');  mass = 12.011
    case ('N');  mass = 14.007
    case ('O');  mass = 15.999
    case ('F');  mass = 18.998
    case ('Ne'); mass = 20.180
    case ('Na'); mass = 22.990
    case ('Mg'); mass = 24.305
    case ('Al'); mass = 26.982
    case ('Si'); mass = 28.086
    case ('P');  mass = 30.974
    case ('S');  mass = 32.065
    case ('Cl'); mass = 35.453
    case ('Ar'); mass = 39.948
    case ('K');  mass = 39.098
    case ('Ca'); mass = 40.078
    case ('Sc'); mass = 44.956
    case ('Ti'); mass = 47.867
    case ('V');  mass = 50.942
    case ('Cr'); mass = 51.996
    case ('Mn'); mass = 54.938
    case ('Fe'); mass = 55.845
    case ('Co'); mass = 58.933
    case ('Ni'); mass = 58.693
    case ('Cu'); mass = 63.546
    case ('Zn'); mass = 65.38
    case ('Ga'); mass = 69.723
    case ('Ge'); mass = 72.630
    case ('As'); mass = 74.922
    case ('Se'); mass = 78.971
    case ('Br'); mass = 79.904
    case ('Kr'); mass = 83.798
    case ('Rb'); mass = 85.468
    case ('Sr'); mass = 87.620
    case ('Y');  mass = 88.906
    case ('Zr'); mass = 91.224
    case ('Nb'); mass = 92.906
    case ('Mo'); mass = 95.960
    case ('Tc'); mass = 98.0
    case ('Ru'); mass = 101.07
    case ('Rh'); mass = 102.91
    case ('Pd'); mass = 106.42
    case ('Ag'); mass = 107.87
    case ('Cd'); mass = 112.41
    case ('In'); mass = 114.82
    case ('Sn'); mass = 118.71
    case ('Sb'); mass = 121.76
    case ('Te'); mass = 127.60
    case ('I');  mass = 126.90
    case ('Xe'); mass = 131.29
    case ('Cs'); mass = 132.91
    case ('Ba'); mass = 137.33
    case ('La'); mass = 138.91
    case ('Ce'); mass = 140.12
    case ('Pr'); mass = 140.91
    case ('Nd'); mass = 144.24
    case ('Pm'); mass = 145.0
    case ('Sm'); mass = 150.36
    case ('Eu'); mass = 151.96
    case ('Gd'); mass = 157.25
    case ('Tb'); mass = 158.93
    case ('Dy'); mass = 162.50
    case ('Ho'); mass = 164.93
    case ('Er'); mass = 167.26
    case ('Tm'); mass = 168.93
    case ('Yb'); mass = 173.04
    case ('Lu'); mass = 174.97
    case ('Hf'); mass = 178.49
    case ('Ta'); mass = 180.95
    case ('W');  mass = 183.84
    case ('Re'); mass = 186.21
    case ('Os'); mass = 190.23
    case ('Ir'); mass = 192.22
    case ('Pt'); mass = 195.08
    case ('Au'); mass = 196.97
    case ('Hg'); mass = 200.59
    case ('Tl'); mass = 204.38
    case ('Pb'); mass = 207.2
    case ('Bi'); mass = 208.98
    case ('Po'); mass = 209.0
    case ('At'); mass = 210.0
    case ('Rn'); mass = 222.0
    case ('Fr'); mass = 223.0
    case ('Ra'); mass = 226.0
    case ('Ac'); mass = 227.0
    case ('Th'); mass = 232.04
    case ('Pa'); mass = 231.04
    case ('U');  mass = 238.03
    case ('Np'); mass = 237.0
    case ('Pu'); mass = 244.0
    case ('Am'); mass = 243.0
    case ('Cm'); mass = 247.0
    case ('Bk'); mass = 247.0
    case ('Cf'); mass = 251.0
    case ('Es'); mass = 252.0
    case ('Fm'); mass = 257.0
    case ('Md'); mass = 258.0
    case ('No'); mass = 259.0
    case ('Lr'); mass = 262.0
    case ('Rf'); mass = 265.0
    case ('Db'); mass = 268.0
    case ('Sg'); mass = 271.0
    case ('Bh'); mass = 270.0
    case ('Hs'); mass = 277.0
    case ('Mt'); mass = 276.0
    case ('Ds'); mass = 281.0
    case ('Rg'); mass = 280.0
    case ('Cn'); mass = 285.0
    case ('Nh'); mass = 284.0
    case ('Fl'); mass = 289.0
    case ('Mc'); mass = 288.0
    case ('Lv'); mass = 293.0
    case ('Ts'); mass = 294.0
    case ('Og'); mass = 294.0
    case default
      write (*, *) "Warning: unknown element '", trim(sym), "', using mass = 1.0"
      mass = 1.0
    end select
  end function get_mass

!---------------------------------------------------------------------------
  subroutine random_rotation_matrix(R)
    real(real64), intent(out) :: R(3,3)
    real(real64) :: u1, u2, u3, q0, q1, q2, q3, twopi
    twopi = 8.0_real64 * atan(1.0_real64)
    call random_number(u1)
    call random_number(u2)
    call random_number(u3)
    q0 = sqrt(1.0-u1) * sin(twopi*u2)
    q1 = sqrt(1.0-u1) * cos(twopi*u2)
    q2 = sqrt(u1)     * sin(twopi*u3)
    q3 = sqrt(u1)     * cos(twopi*u3)
    ! Quaternion to rotation matrix
    R(1,1) = 1.0 - 2.0*(q2**2 + q3**2)
    R(1,2) = 2.0*(q1*q2 - q0*q3)
    R(1,3) = 2.0*(q1*q3 + q0*q2)
    R(2,1) = 2.0*(q1*q2 + q0*q3)
    R(2,2) = 1.0 - 2.0*(q1**2 + q3**2)
    R(2,3) = 2.0*(q2*q3 - q0*q1)
    R(3,1) = 2.0*(q1*q3 - q0*q2)
    R(3,2) = 2.0*(q2*q3 + q0*q1)
    R(3,3) = 1.0 - 2.0*(q1**2 + q2**2)
  end subroutine random_rotation_matrix

!---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------
  !  Resolve one atom sitting on a substitution target.
  !
  !  The configuration row lists substituted sites target by target and, within
  !  a target, species by species, so walking that layout gives the (target,
  !  slot) this atom belongs to; slot nk(t)+1 means "not substituted".  Doing it
  !  here rather than in each writer is what makes every writer correct for any
  !  number of targets and any number of species per target.
  ! ---------------------------------------------------------------------------
  subroutine resolve_target_atom(at_in, sp_in, action, sym_out, mol_out)
    integer, intent(in)  :: at_in, sp_in
    integer, intent(out) :: action, mol_out
    character(len=*), intent(out) :: sym_out
    integer :: tt, jj, ii, slot, off, off_t

    action  = site_plain
    mol_out = 0
    sym_out = ''

    off = 0
    do tt = 1, ntarget
      if (sp_in == sptarget(tt)) then
        slot  = nk(tt) + 1
        off_t = off
        do jj = 1, nk(tt)
          do ii = 1, nsubs_t(tt, jj)
            if (newconf(off_t + ii) == at_in) then
              slot = jj
              exit
            end if
          end do
          if (slot <= nk(tt)) exit
          off_t = off_t + nsubs_t(tt, jj)
        end do
        if (is_vac_t(tt, slot)) then
          action = site_vacancy
        else if (is_mol_t(tt, slot)) then
          action  = site_molecule
          mol_out = mol_idx_t(tt, slot)
        else
          sym_out = newsymbol(tt, slot)
        end if
        return
      end if
      off = off + nsubs_tot_t(tt)
    end do
  end subroutine resolve_target_atom

  subroutine mol_rotate_frac(mol_idx, fx, fy, fz, cv, frac_out, nm)
    integer, intent(in) :: mol_idx, nm
    real(real64), intent(in)    :: fx, fy, fz, cv(3,3)
    real(real64), intent(out)   :: frac_out(nmolmax_atoms, 3)
    integer :: im_loc, k, j
    real(real64) :: rc(3), d_rot(3), r_abs(3), f(3)

    ! Convert site fractional to Cartesian: r = cv * f
    do k = 1, 3
      rc(k) = cv(k,1)*fx + cv(k,2)*fy + cv(k,3)*fz
    end do

    ! Generate random rotation matrix
    call random_rotation_matrix(Rmat)

    do im_loc = 1, nm
      ! Rotate molecule atom displacement
      do k = 1, 3
        d_rot(k) = 0.0
        do j = 1, 3
          d_rot(k) = d_rot(k) + Rmat(k,j) * mol_xyz_arr(mol_idx, im_loc, j)
        end do
      end do
      ! Absolute Cartesian position
      do k = 1, 3
        r_abs(k) = rc(k) + d_rot(k)
      end do
      ! Back-substitute to get fractional (cv is upper triangular)
      f(3) = r_abs(3) / cv(3,3)
      f(2) = (r_abs(2) - cv(2,3)*f(3)) / cv(2,2)
      f(1) = (r_abs(1) - cv(1,2)*f(2) - cv(1,3)*f(3)) / cv(1,1)
      ! Wrap into [0,1)
      do k = 1, 3
        f(k) = f(k) - floor(f(k))
      end do
      frac_out(im_loc, 1:3) = f(1:3)
    end do
  end subroutine mol_rotate_frac

end program genersod
