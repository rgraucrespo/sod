!*******************************************************************************
!    Copyright (c) 2022 Ricardo Grau-Crespo and co-authors
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

program sqssod
  use iso_fortran_env, only: real64
  use insod_reader,    only: insod_t, parse_insod => read_insod
  use ensemble_io,     only: parse_ensemble => read_ensemble
  use config_sampling, only: load_eqmatrix_all
  implicit none

  ! ================================================================
  ! Parameters
  ! ================================================================
  integer, parameter :: max_k_param = 6
  integer, parameter :: natmax = 10000
  integer, parameter :: max_pair_orbits_param = 20000
  integer, parameter :: n_top_default = 10
  real(real64), parameter :: dist_tol = 1.0e-6_real64

  ! ================================================================
  ! Data structures
  ! ================================================================
  type :: target_t
    integer :: nsites = 0
    integer :: nspecies = 0
    integer, allocatable :: global_sites(:)
    character(len=10), allocatable :: species(:)
    integer, allocatable :: counts(:)
    real(real64), allocatable :: concentration(:)
  end type target_t

  type :: pair_orbit_t
    integer :: n_inst = 0
    integer, allocatable :: site_i(:)
    integer, allocatable :: site_j(:)
    integer :: target_i = 0
    integer :: target_j = 0
    integer :: family = 0
    real(real64) :: distance = 0.0_real64
    real(real64) :: weight = 1.0_real64
  end type pair_orbit_t

  type :: pair_channel_t
    integer :: orbit = 0
    integer :: species_i = 0
    integer :: species_j = 0
    real(real64) :: target_corr = 0.0_real64
    real(real64) :: eps_grid = 0.0_real64
  end type pair_channel_t

  type :: pair_family_t
    integer :: target_i = 0
    integer :: target_j = 0
  end type pair_family_t

  ! ================================================================
  ! INSQS input
  ! ================================================================
  character(len=256) :: ensemble_dir
  integer :: max_order
  real(real64) :: cutoff_r(max_k_param)
  real(real64) :: weight_k(max_k_param)
  real(real64) :: omega
  real(real64) :: eps_tol
  integer :: n_top_sqs

  ! ================================================================
  ! Disorder model
  ! ================================================================
  type(insod_t) :: insod
  integer :: ntarget, ndis
  integer, allocatable :: nk_t(:)
  integer, allocatable :: npos_t(:), atini_t(:), atfin_t(:), site_offset(:)
  integer, allocatable :: target_of_site(:), local_site_index(:), global_atom_index(:)
  type(target_t), allocatable :: targets(:)

  ! ================================================================
  ! Symmetry data (from EQMATRIX)
  ! ================================================================
  integer :: nop
  integer, allocatable :: eqmt_blocks(:,:,:)  ! (ntarget, nop, max(npos_t))
  integer, allocatable :: eqm_all(:,:)        ! (ndis, nop), concatenated sites

  ! ================================================================
  ! Supercell geometry (from supercell.cif)
  ! ================================================================
  real(real64) :: sc_a, sc_b, sc_c, sc_alpha, sc_beta, sc_gamma
  real(real64) :: cellvec(3, 3)
  integer :: nat_total
  real(real64), allocatable :: tcoords(:,:)   ! (ndis, 3) fractional
  real(real64), allocatable :: distmat(:,:)   ! (ndis, ndis)

  ! ================================================================
  ! Pair orbits and species-resolved channels
  ! ================================================================
  integer :: n_pair_orbits, n_channels, n_families
  type(pair_orbit_t), allocatable :: pair_orbits(:)
  type(pair_channel_t), allocatable :: channels(:)
  type(pair_family_t), allocatable :: families(:)
  integer, allocatable :: family_index(:,:)
  integer, allocatable :: channel_start(:), channel_end(:)
  integer, allocatable :: pair_order(:)

  ! ================================================================
  ! Configuration data (from ENSEMBLE)
  ! ================================================================
  integer :: nic, nsubs_tot
  integer, allocatable :: nsubs_flat(:)
  integer, allocatable :: indconf(:,:)        ! raw global site indices from ENSEMBLE
  integer, allocatable :: degen(:)
  integer, allocatable :: occupation(:,:)     ! (ndis, nic), target-local species code

  ! ================================================================
  ! Scores and correlations
  ! ================================================================
  real(real64), allocatable :: scores(:)
  real(real64), allocatable :: corr_all(:,:)       ! (nic, n_channels)
  real(real64), allocatable :: L_val(:)
  real(real64), allocatable :: abs_residual(:)
  real(real64), allocatable :: eps_min_cfg(:)      ! (n_channels)
  real(real64), allocatable :: orbit_error(:,:)    ! (nic, n_pair_orbits)
  real(real64), allocatable :: family_error(:,:)   ! (nic, n_families)

  ! ================================================================
  ! Main program
  ! ================================================================
  write (*, '(A)') "SOD (Site-Occupancy Disorder) version 0.90 - sqssod"

  call read_insqs()
  call read_insod_model()
  call read_supercell_cif()
  call read_eqmatrix()
  call compute_distances()
  call generate_pair_orbits()
  call read_ensemble_data()
  call generate_pair_channels()
  call score_configurations()
  call report_results()

  write (*, *)
  write (*, *) " > Done."

contains

  ! ================================================================
  ! Read INSQS input file
  ! ================================================================
  subroutine read_insqs()
    implicit none
    integer :: iu, ios, k
    character(len=256) :: insqs_path
    character(len=256) :: line
    logical :: insqs_exists
    iu = 20

    ensemble_dir = "."
    if (command_argument_count() >= 1) then
      call get_command_argument(1, ensemble_dir)
      ensemble_dir = adjustl(ensemble_dir)
    end if

    insqs_path = trim(ensemble_dir) // "/INSQS"
    inquire (file=trim(insqs_path), exist=insqs_exists)
    if (.not. insqs_exists) insqs_path = "INSQS"

    open (unit=iu, file=trim(insqs_path), status='old', IOSTAT=ios)
    if (ios /= 0) then
      write (*, '(a)') "Error: INSQS not found in nXX/ or SODPROJECT/."
      stop 1
    end if

    read (iu, *)
    read (iu, *) max_order
    if (max_order < 2 .or. max_order > max_k_param) then
      write (*, '(a,i0,a,i0)') " Error: MaxOrder must be 2..", max_k_param, &
        ", got ", max_order
      stop 1
    end if

    read (iu, *)
    read (iu, *)
    cutoff_r = 0.0_real64
    read (iu, *) (cutoff_r(k), k = 2, max_order)

    read (iu, *)
    read (iu, *)
    weight_k = 0.0_real64
    read (iu, *) (weight_k(k), k = 2, max_order)

    read (iu, *)
    read (iu, *)
    read (iu, *) omega, eps_tol

    ! Optional trailing line: n_top_sqs (0 = rank and list all configurations).
    ! Older INSQS files without this line default to n_top_default. Blank and
    ! '#' comment lines before the value are skipped.
    n_top_sqs = n_top_default
    do
      read (iu, '(a)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      read (line, *, iostat=ios) k
      if (ios == 0) n_top_sqs = k
      exit
    end do
    if (n_top_sqs < 0) then
      write (*, '(a,i0)') " Error: n_top_sqs must be >= 0, got ", n_top_sqs
      stop 1
    end if
    close (iu)

    write (*, *) " > Input parameters (INSQS):"
    write (*, '(a,a)') "    ENSEMBLE directory:   ", trim(ensemble_dir)
    write (*, '(a,i0)') "    Max cluster order:  ", max_order
    write (*, '(a)', advance='no') "    Cutoff radii (A):  "
    do k = 2, max_order
      write (*, '(f8.3)', advance='no') cutoff_r(k)
    end do
    write (*, *)
    write (*, '(a)', advance='no') "    Weights:           "
    do k = 2, max_order
      write (*, '(f8.4)', advance='no') weight_k(k)
    end do
    write (*, *)
    write (*, '(a,es10.3,a,es10.3)') "    Scoring:            van de Walle  omega=", &
      omega, "  eps=", eps_tol
    if (n_top_sqs == 0) then
      write (*, '(a)') "    Ranking:            all configurations (n_top_sqs=0)"
    else
      write (*, '(a,i0,a)') "    Ranking:            top ", n_top_sqs, " configurations (n_top_sqs)"
    end if
    if (max_order > 2) then
      write (*, '(a)') "    Note: generalized SQS currently uses species-resolved pair correlations only."
      write (*, '(a)') "          INSQS orders above 2 are read for compatibility and ignored here."
    end if
    write (*, *)
  end subroutine read_insqs

  ! ================================================================
  ! Read INSOD with the central parser and initialise target species.
  ! Counts/concentrations are set later from the ENSEMBLE target line.
  ! ================================================================
  subroutine read_insod_model()
    implicit none
    integer :: t, j

    call parse_insod('INSOD', insod)

    ntarget = insod%ntarget
    if (ntarget < 1) then
      write (*, '(a)') " Error: INSOD contains no substitution targets."
      stop 1
    end if

    allocate (targets(ntarget), nk_t(ntarget))
    do t = 1, ntarget
      if (insod%sptarget(t) < 1 .or. insod%sptarget(t) > insod%nsp) then
        write (*, '(a,i0,a,i0)') " Error: invalid sptarget for target ", t, &
          ": ", insod%sptarget(t)
        stop 1
      end if
      nk_t(t) = insod%nk(t)
      targets(t)%nspecies = nk_t(t) + 1
      if (targets(t)%nspecies < 2) then
        write (*, '(a,i0,a)') " Error: target ", t, " has fewer than two species."
        stop 1
      end if

      allocate (targets(t)%species(targets(t)%nspecies))
      allocate (targets(t)%counts(targets(t)%nspecies))
      allocate (targets(t)%concentration(targets(t)%nspecies))
      targets(t)%species(1) = adjustl(insod%newsymbol(t, nk_t(t) + 1))
      do j = 1, nk_t(t)
        targets(t)%species(j + 1) = adjustl(insod%newsymbol(t, j))
      end do
      targets(t)%counts = 0
      targets(t)%concentration = 0.0_real64
    end do

    write (*, '(a,i0)') " > Substitution targets (from INSOD): ", ntarget
    do t = 1, ntarget
      write (*, '(a,i0,a,a,a)', advance='no') "    Target ", t, ": parent ", &
        trim(insod%symbol(insod%sptarget(t))), " ->"
      do j = 1, targets(t)%nspecies
        write (*, '(1x,a)', advance='no') trim(targets(t)%species(j))
      end do
      write (*, *)
    end do
    write (*, *)
  end subroutine read_insod_model

  ! ================================================================
  ! Read supercell.cif and extract coordinates of every target block.
  ! ================================================================
  subroutine read_supercell_cif()
    implicit none
    integer :: iu, ios, iat, t, loc, sp, n_this
    character(len=512) :: line
    character(len=20) :: label, typesym
    character(len=10) :: target_sym
    real(real64) :: fx, fy, fz
    character(len=10) :: all_sym(natmax)
    real(real64) :: all_frac(natmax, 3)
    logical :: in_atoms

    iu = 22
    open (unit=iu, file="supercell.cif", status='old', IOSTAT=ios)
    if (ios /= 0) then
      write (*, *) "Error: supercell.cif not found. Run combsod first."
      stop 1
    end if

    sc_a = 0.0_real64; sc_b = 0.0_real64; sc_c = 0.0_real64
    sc_alpha = 90.0_real64; sc_beta = 90.0_real64; sc_gamma = 90.0_real64
    nat_total = 0
    in_atoms = .false.

    do
      read (iu, '(A)', IOSTAT=ios) line
      if (ios /= 0) exit
      line = adjustl(line)

      if (.not. in_atoms) then
        if (index(line, '_cell_length_a') == 1) then
          read (line(15:), *) sc_a
        else if (index(line, '_cell_length_b') == 1) then
          read (line(15:), *) sc_b
        else if (index(line, '_cell_length_c') == 1) then
          read (line(15:), *) sc_c
        else if (index(line, '_cell_angle_alpha') == 1) then
          read (line(18:), *) sc_alpha
        else if (index(line, '_cell_angle_beta') == 1) then
          read (line(17:), *) sc_beta
        else if (index(line, '_cell_angle_gamma') == 1) then
          read (line(18:), *) sc_gamma
        else if (index(line, '_atom_site_fract_z') == 1) then
          in_atoms = .true.
        end if
      else
        if (len_trim(line) == 0) exit
        read (line, *, IOSTAT=ios) label, typesym, fx, fy, fz
        if (ios /= 0) exit
        nat_total = nat_total + 1
        if (nat_total > natmax) then
          write (*, '(a,i0)') " Error: too many atoms in supercell.cif; natmax=", natmax
          stop 1
        end if
        all_sym(nat_total) = adjustl(typesym)
        all_frac(nat_total, 1) = fx
        all_frac(nat_total, 2) = fy
        all_frac(nat_total, 3) = fz
      end if
    end do
    close (iu)

    allocate (npos_t(ntarget), atini_t(ntarget), atfin_t(ntarget), site_offset(ntarget))

    ndis = 0
    do t = 1, ntarget
      target_sym = adjustl(insod%symbol(insod%sptarget(t)))
      n_this = 0
      atini_t(t) = 0
      atfin_t(t) = 0
      do iat = 1, nat_total
        if (trim(all_sym(iat)) == trim(target_sym)) then
          n_this = n_this + 1
          if (n_this == 1) atini_t(t) = iat
          atfin_t(t) = iat
        end if
      end do
      if (n_this <= 0) then
        write (*, '(a,i0,a,a,a)') " Error: found no atoms for target ", t, &
          " parent species ", trim(target_sym), " in supercell.cif."
        stop 1
      end if
      if (atfin_t(t) - atini_t(t) + 1 /= n_this) then
        write (*, '(a,i0,a)') " Error: target ", t, &
          " atoms are not contiguous in supercell.cif; ENSEMBLE offsets require contiguous blocks."
        stop 1
      end if
      do sp = t + 1, ntarget
        if (insod%sptarget(sp) == insod%sptarget(t)) then
          write (*, '(a)') " Error: duplicate target parent species are not supported by sqssod."
          stop 1
        end if
      end do

      npos_t(t) = n_this
      targets(t)%nsites = n_this
      allocate (targets(t)%global_sites(n_this))
      do loc = 1, n_this
        targets(t)%global_sites(loc) = atini_t(t) + loc - 1
      end do
      site_offset(t) = ndis
      ndis = ndis + n_this
    end do

    allocate (target_of_site(ndis), local_site_index(ndis), global_atom_index(ndis))
    allocate (tcoords(ndis, 3))
    do t = 1, ntarget
      do loc = 1, npos_t(t)
        iat = site_offset(t) + loc
        target_of_site(iat) = t
        local_site_index(iat) = loc
        global_atom_index(iat) = atini_t(t) + loc - 1
        tcoords(iat, :) = all_frac(global_atom_index(iat), :)
      end do
    end do

    call cell(cellvec, sc_a, sc_b, sc_c, sc_alpha, sc_beta, sc_gamma)

    write (*, '(a,3f10.4,3f9.2)') " > Cell: ", sc_a, sc_b, sc_c, &
      sc_alpha, sc_beta, sc_gamma
    write (*, '(a,i0,a,i0,a)') " > Atoms: ", nat_total, " total, ", ndis, " disordered sites"
    do t = 1, ntarget
      write (*, '(a,i0,a,i0,a,i0,a,i0)') "    Target ", t, ": ", npos_t(t), &
        " sites, global atom range ", atini_t(t), "-", atfin_t(t)
    end do
    write (*, *)
  end subroutine read_supercell_cif

  ! ================================================================
  ! Read every target block from EQMATRIX and build concatenated maps.
  ! ================================================================
  subroutine read_eqmatrix()
    implicit none
    integer :: iop, t, loc, isite, image_site

    call load_eqmatrix_all('EQMATRIX', ntarget, npos_t, eqmt_blocks, nop)

    allocate (eqm_all(ndis, nop))
    do iop = 1, nop
      do t = 1, ntarget
        do loc = 1, npos_t(t)
          isite = site_offset(t) + loc
          image_site = site_offset(t) + eqmt_blocks(t, iop, loc)
          eqm_all(isite, iop) = image_site
        end do
      end do
    end do

    call validate_eqmatrix()

    write (*, '(a,i0,a,i0,a,i0,a)') " > EQMATRIX: ", nop, " operators, ", &
      ntarget, " target block(s), ", ndis, " concatenated sites"
    write (*, *)
  end subroutine read_eqmatrix

  subroutine validate_eqmatrix()
    implicit none
    integer :: iop, isite, image_site
    logical, allocatable :: seen(:)

    allocate (seen(ndis))
    do iop = 1, nop
      seen = .false.
      do isite = 1, ndis
        image_site = eqm_all(isite, iop)
        if (image_site < 1 .or. image_site > ndis) then
          write (*, '(a,i0,a,i0)') " Error: EQMATRIX operator ", iop, &
            " maps a site out of range: ", image_site
          stop 1
        end if
        if (seen(image_site)) then
          write (*, '(a,i0,a)') " Error: EQMATRIX operator ", iop, &
            " is not a permutation."
          stop 1
        end if
        if (target_of_site(image_site) /= target_of_site(isite)) then
          write (*, '(a,i0,a)') " Error: EQMATRIX operator ", iop, &
            " maps a site outside its target block."
          stop 1
        end if
        seen(image_site) = .true.
      end do
    end do
    deallocate (seen)
  end subroutine validate_eqmatrix

  ! ================================================================
  ! Compute minimum-image distance matrix between disordered sites.
  ! ================================================================
  subroutine compute_distances()
    implicit none
    integer :: ia, ja
    real(real64) :: df(3), dr(3), dd

    allocate (distmat(ndis, ndis))
    distmat = 0.0_real64

    do ia = 1, ndis
      do ja = ia + 1, ndis
        df(:) = tcoords(ia, :) - tcoords(ja, :)
        df(:) = df(:) - nint(df(:))
        dr = matmul(cellvec, df)
        dd = sqrt(sum(dr**2))
        distmat(ia, ja) = dd
        distmat(ja, ia) = dd
      end do
    end do

    if (ndis > 1) then
      write (*, '(a,f8.4,a)') " > Distances: min = ", &
        minval(distmat, mask = (distmat > 0.0_real64)), " A"
      write (*, '(a,f8.4,a)') "              max = ", maxval(distmat), " A"
    else
      write (*, '(a)') " > Distances: only one disordered site"
    end if
    write (*, *)
  end subroutine compute_distances

  ! ================================================================
  ! Generate symmetry-distinct pair orbits over all disordered targets.
  ! Same-target orbits store both orientations. Cross-target orbits are
  ! oriented from the lower target index to the higher target index.
  ! ================================================================
  subroutine generate_pair_orbits()
    implicit none
    integer :: ia, ja, iop, a, b, tmp, m, n_temp
    integer, allocatable :: temp_i(:), temp_j(:)
    logical :: is_canonical, is_dup

    allocate (pair_orbits(max_pair_orbits_param))
    allocate (temp_i(nop), temp_j(nop))
    n_pair_orbits = 0

    do ia = 1, ndis - 1
      do ja = ia + 1, ndis
        if (distmat(ia, ja) > cutoff_r(2) + dist_tol) cycle

        is_canonical = .true.
        do iop = 1, nop
          a = eqm_all(ia, iop)
          b = eqm_all(ja, iop)
          if (a > b) then
            tmp = a; a = b; b = tmp
          end if
          if (pair_lex_lt(a, b, ia, ja)) then
            is_canonical = .false.
            exit
          end if
        end do
        if (.not. is_canonical) cycle

        n_temp = 0
        do iop = 1, nop
          a = eqm_all(ia, iop)
          b = eqm_all(ja, iop)
          if (a > b) then
            tmp = a; a = b; b = tmp
          end if

          if (abs(distmat(a, b) - distmat(ia, ja)) > 1.0e-5_real64) then
            write (*, '(a)') " Error: EQMATRIX does not preserve pair distances."
            stop 1
          end if

          is_dup = .false.
          do m = 1, n_temp
            if (temp_i(m) == a .and. temp_j(m) == b) then
              is_dup = .true.
              exit
            end if
          end do
          if (.not. is_dup) then
            n_temp = n_temp + 1
            temp_i(n_temp) = a
            temp_j(n_temp) = b
          end if
        end do

        n_pair_orbits = n_pair_orbits + 1
        if (n_pair_orbits > max_pair_orbits_param) then
          write (*, '(a,i0,a)') " Error: more than ", max_pair_orbits_param, &
            " pair orbits. Reduce the pair cutoff or increase max_pair_orbits_param."
          stop 1
        end if

        if (target_of_site(ia) == target_of_site(ja)) then
          pair_orbits(n_pair_orbits)%n_inst = 2 * n_temp
          allocate (pair_orbits(n_pair_orbits)%site_i(2 * n_temp))
          allocate (pair_orbits(n_pair_orbits)%site_j(2 * n_temp))
          do m = 1, n_temp
            pair_orbits(n_pair_orbits)%site_i(m) = temp_i(m)
            pair_orbits(n_pair_orbits)%site_j(m) = temp_j(m)
            pair_orbits(n_pair_orbits)%site_i(n_temp + m) = temp_j(m)
            pair_orbits(n_pair_orbits)%site_j(n_temp + m) = temp_i(m)
          end do
        else
          pair_orbits(n_pair_orbits)%n_inst = n_temp
          allocate (pair_orbits(n_pair_orbits)%site_i(n_temp))
          allocate (pair_orbits(n_pair_orbits)%site_j(n_temp))
          pair_orbits(n_pair_orbits)%site_i(1:n_temp) = temp_i(1:n_temp)
          pair_orbits(n_pair_orbits)%site_j(1:n_temp) = temp_j(1:n_temp)
        end if

        pair_orbits(n_pair_orbits)%target_i = target_of_site(pair_orbits(n_pair_orbits)%site_i(1))
        pair_orbits(n_pair_orbits)%target_j = target_of_site(pair_orbits(n_pair_orbits)%site_j(1))
        pair_orbits(n_pair_orbits)%distance = distmat(ia, ja)
        pair_orbits(n_pair_orbits)%weight = weight_k(2)
      end do
    end do

    deallocate (temp_i, temp_j)

    call assign_pair_families()
    call sort_pair_order()

    write (*, '(a)') " > Pair orbits generated:"
    if (n_pair_orbits == 0) then
      write (*, '(a)') "    No pair orbits survived the distance filter."
    else
      write (*, '(a)') " Orbit  T_i  T_j  Instances  Dist(A)  Weight"
      write (*, '(a)') "  ---- ---- ---- ---------- -------- ------"
      do ia = 1, n_pair_orbits
        m = pair_order(ia)
        write (*, '(i6, i5, i5, i10, f9.4, f8.4)') m, pair_orbits(m)%target_i, &
          pair_orbits(m)%target_j, pair_orbits(m)%n_inst, pair_orbits(m)%distance, &
          pair_orbits(m)%weight
      end do
    end if
    write (*, '(a,i0)') "    Total pair orbits: ", n_pair_orbits
    write (*, '(a,i0)') "    Pair families:     ", n_families
    write (*, *)
  end subroutine generate_pair_orbits

  subroutine assign_pair_families()
    implicit none
    integer :: p, t, u, f

    allocate (family_index(ntarget, ntarget))
    family_index = 0
    allocate (families(max(1, n_pair_orbits)))
    n_families = 0

    do p = 1, n_pair_orbits
      t = pair_orbits(p)%target_i
      u = pair_orbits(p)%target_j
      if (t > u) then
        write (*, '(a)') " Error: internal cross-target pair orientation is inconsistent."
        stop 1
      end if
      if (family_index(t, u) == 0) then
        n_families = n_families + 1
        family_index(t, u) = n_families
        families(n_families)%target_i = t
        families(n_families)%target_j = u
      end if
      f = family_index(t, u)
      pair_orbits(p)%family = f
    end do
  end subroutine assign_pair_families

  subroutine sort_pair_order()
    implicit none
    integer :: i, j, imin, tmp

    allocate (pair_order(max(1, n_pair_orbits)))
    do i = 1, n_pair_orbits
      pair_order(i) = i
    end do

    do i = 1, n_pair_orbits - 1
      imin = i
      do j = i + 1, n_pair_orbits
        if (pair_before(pair_order(j), pair_order(imin))) imin = j
      end do
      if (imin /= i) then
        tmp = pair_order(i)
        pair_order(i) = pair_order(imin)
        pair_order(imin) = tmp
      end if
    end do
  end subroutine sort_pair_order

  logical function pair_before(pa, pb) result(res)
    implicit none
    integer, intent(in) :: pa, pb
    res = .false.
    if (pair_orbits(pa)%distance < pair_orbits(pb)%distance - dist_tol) then
      res = .true.
    else if (abs(pair_orbits(pa)%distance - pair_orbits(pb)%distance) <= dist_tol) then
      if (pair_orbits(pa)%target_i < pair_orbits(pb)%target_i) then
        res = .true.
      else if (pair_orbits(pa)%target_i == pair_orbits(pb)%target_i) then
        if (pair_orbits(pa)%target_j < pair_orbits(pb)%target_j) res = .true.
      end if
    end if
  end function pair_before

  ! ================================================================
  ! Read ENSEMBLE and construct dense target-local occupations.
  ! ================================================================
  subroutine read_ensemble_data()
    implicit none
    integer :: ntarget_ens, ic, t, j, k, flat_off, conf_off
    integer :: gsite, local_idx, isite, expected_nsubs
    integer, allocatable :: npos_t_ens(:), counts_seen(:)
    real(real64) :: tsampling
    logical :: ok
    logical, allocatable :: assigned(:)
    character(len=256) :: path

    path = trim(ensemble_dir) // "/ENSEMBLE"
    call parse_ensemble(trim(path), ntarget_ens, nsubs_flat, npos_t_ens, &
                        nic, indconf, degen, tsampling, ok)
    if (.not. ok) then
      write (*, '(a,a)') " Error: cannot read ENSEMBLE: ", trim(path)
      stop 1
    end if

    if (ntarget_ens /= ntarget) then
      write (*, '(a,i0,a,i0)') " Error: ENSEMBLE target count ", ntarget_ens, &
        " does not match INSOD target count ", ntarget
      stop 1
    end if
    if (size(nsubs_flat) /= sum(nk_t(1:ntarget))) then
      write (*, '(a,i0,a,i0)') " Error: ENSEMBLE has ", size(nsubs_flat), &
        " substituent count(s), but INSOD target species require ", sum(nk_t(1:ntarget))
      stop 1
    end if
    do t = 1, ntarget
      if (npos_t_ens(t) /= npos_t(t)) then
        write (*, '(a,i0,a,i0,a,i0)') " Error: ENSEMBLE target ", t, &
          " has ", npos_t_ens(t), " sites but supercell/EQMATRIX has ", npos_t(t)
        stop 1
      end if
    end do

    flat_off = 0
    do t = 1, ntarget
      targets(t)%counts(1) = npos_t(t) - sum(nsubs_flat(flat_off + 1:flat_off + nk_t(t)))
      do j = 1, nk_t(t)
        targets(t)%counts(j + 1) = nsubs_flat(flat_off + j)
      end do
      if (any(targets(t)%counts < 0) .or. sum(targets(t)%counts) /= npos_t(t)) then
        write (*, '(a,i0,a)') " Error: invalid species counts for target ", t, "."
        stop 1
      end if
      targets(t)%concentration = real(targets(t)%counts, real64) / real(npos_t(t), real64)
      if (abs(sum(targets(t)%concentration) - 1.0_real64) > 1.0e-10_real64) then
        write (*, '(a,i0,a)') " Error: concentrations do not sum to one for target ", t, "."
        stop 1
      end if
      flat_off = flat_off + nk_t(t)
    end do

    nsubs_tot = sum(nsubs_flat)
    expected_nsubs = sum(nsubs_flat)
    if (size(indconf, 2) < max(1, expected_nsubs)) then
      write (*, '(a)') " Error: ENSEMBLE configuration array is shorter than expected."
      stop 1
    end if

    allocate (occupation(ndis, nic))
    allocate (assigned(ndis), counts_seen(maxval(nk_t(1:ntarget)) + 1))

    do ic = 1, nic
      occupation(:, ic) = 1
      assigned = .false.
      conf_off = 0
      flat_off = 0

      do t = 1, ntarget
        do j = 1, nk_t(t)
          do k = 1, targets(t)%counts(j + 1)
            gsite = indconf(ic, conf_off + k)
            local_idx = gsite - atini_t(t) + 1
            if (local_idx < 1 .or. local_idx > npos_t(t)) then
              write (*, '(a,i0,a,i0,a,i0,a)') " Error: config ", ic, &
                " contains site ", gsite, " outside target ", t, "."
              stop 1
            end if
            isite = site_offset(t) + local_idx
            if (assigned(isite)) then
              write (*, '(a,i0,a,i0,a)') " Error: config ", ic, &
                " assigns site ", gsite, " more than once."
              stop 1
            end if
            occupation(isite, ic) = j + 1
            assigned(isite) = .true.
          end do
          conf_off = conf_off + targets(t)%counts(j + 1)
        end do
        flat_off = flat_off + nk_t(t)
      end do

      if (conf_off /= expected_nsubs) then
        write (*, '(a)') " Error: internal ENSEMBLE offset mismatch."
        stop 1
      end if

      do t = 1, ntarget
        counts_seen = 0
        do k = 1, npos_t(t)
          isite = site_offset(t) + k
          if (occupation(isite, ic) < 1 .or. occupation(isite, ic) > targets(t)%nspecies) then
            write (*, '(a,i0,a)') " Error: invalid occupation code in config ", ic, "."
            stop 1
          end if
          counts_seen(occupation(isite, ic)) = counts_seen(occupation(isite, ic)) + 1
        end do
        do j = 1, targets(t)%nspecies
          if (counts_seen(j) /= targets(t)%counts(j)) then
            write (*, '(a,i0,a,i0,a)') " Error: species counts do not match for config ", &
              ic, ", target ", t, "."
            stop 1
          end if
        end do
      end do
    end do

    deallocate (assigned, counts_seen, npos_t_ens)

    write (*, '(a,i0,a,i0,a)') " > ENSEMBLE: ", nic, " configurations, ", &
      nsubs_tot, " substituted sites"
    write (*, '(a)') "    Target compositions:"
    do t = 1, ntarget
      write (*, '(a,i0,a)', advance='no') "      Target ", t, ":"
      do j = 1, targets(t)%nspecies
        write (*, '(1x,a,a,i0,a,f8.5,a)', advance='no') trim(targets(t)%species(j)), &
          "=", targets(t)%counts(j), " (", targets(t)%concentration(j), ")"
      end do
      write (*, *)
    end do
    write (*, *)
  end subroutine read_ensemble_data

  ! ================================================================
  ! Create one chemical channel for every species pair on every pair orbit.
  ! ================================================================
  subroutine generate_pair_channels()
    implicit none
    integer :: p, ch, a, b, t, u, m, n_near
    real(real64) :: nearest

    n_channels = 0
    do p = 1, n_pair_orbits
      t = pair_orbits(p)%target_i
      u = pair_orbits(p)%target_j
      n_channels = n_channels + targets(t)%nspecies * targets(u)%nspecies
    end do

    allocate (channels(max(1, n_channels)))
    allocate (channel_start(max(1, n_pair_orbits)), channel_end(max(1, n_pair_orbits)))

    ch = 0
    do p = 1, n_pair_orbits
      t = pair_orbits(p)%target_i
      u = pair_orbits(p)%target_j
      channel_start(p) = ch + 1
      do a = 1, targets(t)%nspecies
        do b = 1, targets(u)%nspecies
          ch = ch + 1
          channels(ch)%orbit = p
          channels(ch)%species_i = a
          channels(ch)%species_j = b
          channels(ch)%target_corr = targets(t)%concentration(a) * targets(u)%concentration(b)
          m = pair_orbits(p)%n_inst
          n_near = nint(real(m, real64) * channels(ch)%target_corr)
          n_near = max(0, min(m, n_near))
          nearest = real(n_near, real64) / real(m, real64)
          channels(ch)%eps_grid = abs(nearest - channels(ch)%target_corr)
        end do
      end do
      channel_end(p) = ch
    end do

    write (*, '(a,i0)') " > Chemical pair channels: ", n_channels
    write (*, *)
  end subroutine generate_pair_channels

  ! ================================================================
  ! Score all configurations using species-resolved pair probabilities.
  ! ================================================================
  subroutine score_configurations()
    implicit none
    integer :: ic, ch, p, m, i, j, a, b, f, t, u
    integer :: count_ab, n_ch_orbit
    real(real64) :: sum_prob, fam_norm(max(1, n_pair_orbits))

    allocate (scores(nic), L_val(nic), abs_residual(nic))

    if (n_channels == 0) then
      scores = 0.0_real64
      L_val = 0.0_real64
      abs_residual = 0.0_real64
      allocate (corr_all(nic, 0), eps_min_cfg(0), orbit_error(nic, 0), family_error(nic, 0))
      write (*, '(a)') " > Scoring complete: no pair channels were generated."
      write (*, *)
      return
    end if

    allocate (corr_all(nic, n_channels))
    allocate (eps_min_cfg(n_channels))
    allocate (orbit_error(nic, n_pair_orbits))
    allocate (family_error(nic, n_families))

    do ic = 1, nic
      do ch = 1, n_channels
        p = channels(ch)%orbit
        a = channels(ch)%species_i
        b = channels(ch)%species_j
        count_ab = 0
        do m = 1, pair_orbits(p)%n_inst
          i = pair_orbits(p)%site_i(m)
          j = pair_orbits(p)%site_j(m)
          if (occupation(i, ic) == a .and. occupation(j, ic) == b) count_ab = count_ab + 1
        end do
        corr_all(ic, ch) = real(count_ab, real64) / real(pair_orbits(p)%n_inst, real64)
      end do

      do p = 1, n_pair_orbits
        sum_prob = sum(corr_all(ic, channel_start(p):channel_end(p)))
        if (abs(sum_prob - 1.0_real64) > 1.0e-10_real64) then
          write (*, '(a,i0,a,i0,a,es12.4)') " Error: probabilities for config ", ic, &
            ", pair orbit ", p, " sum to ", sum_prob
          stop 1
        end if
      end do
    end do

    do ch = 1, n_channels
      eps_min_cfg(ch) = minval(abs(corr_all(:, ch) - channels(ch)%target_corr))
    end do

    orbit_error = 0.0_real64
    family_error = 0.0_real64
    fam_norm = 0.0_real64
    do p = 1, n_pair_orbits
      f = pair_orbits(p)%family
      fam_norm(f) = fam_norm(f) + pair_orbits(p)%weight
      n_ch_orbit = channel_end(p) - channel_start(p) + 1
      do ic = 1, nic
        orbit_error(ic, p) = sum(abs(corr_all(ic, channel_start(p):channel_end(p)) - &
          channels(channel_start(p):channel_end(p))%target_corr)) / real(n_ch_orbit, real64)
        family_error(ic, f) = family_error(ic, f) + pair_orbits(p)%weight * orbit_error(ic, p)
      end do
    end do

    do f = 1, n_families
      if (fam_norm(f) > 0.0_real64) then
        family_error(:, f) = family_error(:, f) / fam_norm(f)
      else
        family_error(:, f) = 0.0_real64
      end if
    end do

    call compute_vdw_scores()

    write (*, '(a)') " > Target pair-channel probabilities and empirical errors:"
    write (*, '(a)') "   Ch Orbit  T_i  T_j  Sp_i      Sp_j        Target    eps_grid     eps_cfg"
    write (*, '(a)') "  ---- ---- ---- ---- ---------- ---------- ---------- ---------- ----------"
    do ch = 1, n_channels
      p = channels(ch)%orbit
      t = pair_orbits(p)%target_i
      u = pair_orbits(p)%target_j
      write (*, '(i6,i5,i5,i5,1x,a10,1x,a10,f11.6,es11.3,es11.3)') ch, p, t, u, &
        trim(targets(t)%species(channels(ch)%species_i)), &
        trim(targets(u)%species(channels(ch)%species_j)), &
        channels(ch)%target_corr, channels(ch)%eps_grid, eps_min_cfg(ch)
    end do
    write (*, *)

    write (*, '(a)') " > Scoring complete."
    write (*, *)
  end subroutine score_configurations

  ! ================================================================
  ! Van de Walle-style score: Q = -omega*L + normalized pair error.
  ! L advances by complete distance shells only when every channel in
  ! that shell is within its empirical minimum plus eps_tol.
  ! ================================================================
  subroutine compute_vdw_scores()
    implicit none
    integer :: ic, jj, kk, p, ch, f
    real(real64) :: cur_diam
    logical :: shell_ok

    do ic = 1, nic
      L_val(ic) = 0.0_real64
      jj = 1
      do while (jj <= n_pair_orbits)
        cur_diam = pair_orbits(pair_order(jj))%distance
        shell_ok = .true.
        kk = jj
        do while (kk <= n_pair_orbits)
          p = pair_order(kk)
          if (abs(pair_orbits(p)%distance - cur_diam) > dist_tol) exit
          do ch = channel_start(p), channel_end(p)
            if (abs(corr_all(ic, ch) - channels(ch)%target_corr) > eps_min_cfg(ch) + eps_tol) then
              shell_ok = .false.
              exit
            end if
          end do
          if (.not. shell_ok) exit
          kk = kk + 1
        end do
        if (.not. shell_ok) exit
        L_val(ic) = cur_diam
        jj = kk
      end do

      abs_residual(ic) = 0.0_real64
      do f = 1, n_families
        abs_residual(ic) = abs_residual(ic) + family_error(ic, f)
      end do
      if (n_families > 0) abs_residual(ic) = abs_residual(ic) / real(n_families, real64)
      scores(ic) = -omega * L_val(ic) + abs_residual(ic)
    end do
  end subroutine compute_vdw_scores

  ! ================================================================
  ! Report results to stdout and write OUTSQS/SQS_CORRELATIONS.
  ! With n_top_sqs > 0 only the best configurations are selected
  ! (single O(nic*n_top) pass; the full ensemble is never sorted);
  ! n_top_sqs = 0 sorts and lists the whole ensemble.
  ! ================================================================
  subroutine report_results()
    implicit none
    integer :: rank, n_show, n_stdout, f, p, ch, best_cfg
    integer, allocatable :: ridx(:)
    integer :: iu_out, iu_corr
    character(len=256) :: path_out, path_corr

    call rank_configs(ridx, n_show)
    best_cfg = ridx(1)

    n_stdout = min(n_show, 20)
    write (*, '(a,i0,a)') " > Top ", n_stdout, " configurations:"
    write (*, *)
    write (*, '(a)') "    Rank  Config  Degen      L(A)       Error                Q"
    write (*, '(a)') "    ----  ------  -----  --------  ----------  ---------------"
    do rank = 1, n_stdout
      write (*, '(i7, i8, i7, f10.4, es12.4, es16.6)') rank, ridx(rank), degen(ridx(rank)), &
        L_val(ridx(rank)), abs_residual(ridx(rank)), scores(ridx(rank))
    end do
    write (*, *)

    write (*, '(a,i0)') " > Best SQS: configuration ", best_cfg
    write (*, '(a,i0)') "    Degeneracy:  ", degen(best_cfg)
    write (*, '(a,f10.4,a)') "    Matched L:   ", L_val(best_cfg), " Angstroms"
    write (*, '(a,es14.6)') "    Error:       ", abs_residual(best_cfg)
    write (*, '(a,es14.6)') "    Q (score):   ", scores(best_cfg)
    if (n_families > 0) then
      write (*, '(a)') "    Family errors:"
      do f = 1, n_families
        write (*, '(a,i0,i0,a,es14.6)') "      E", families(f)%target_i, families(f)%target_j, &
          ": ", family_error(best_cfg, f)
      end do
    end if
    write (*, *)

    path_out = trim(ensemble_dir) // "/OUTSQS"
    iu_out = 30
    open (unit=iu_out, file=trim(path_out))
    write (iu_out, '(a)') "# SQS scoring results from sqssod"
    write (iu_out, '(a,i0,a,i0,a,i0,a,i0)') "# ntarget=", ntarget, &
      " ndis=", ndis, " n_pair_orbits=", n_pair_orbits, " n_channels=", n_channels
    write (iu_out, '(a,es10.3,a,es10.3)') "# Mode: species_resolved_pairs  omega=", &
      omega, "  eps=", eps_tol
    tblock: block
      integer :: t, j
      do t = 1, ntarget
        write (iu_out, '(a,i0,a,i0,a)', advance='no') "# Target ", t, ": sites=", targets(t)%nsites, " species="
        do j = 1, targets(t)%nspecies
          write (iu_out, '(1x,a,a,i0,a,f8.5,a)', advance='no') trim(targets(t)%species(j)), &
            "=", targets(t)%counts(j), "(", targets(t)%concentration(j), ")"
        end do
        write (iu_out, *)
      end do
    end block tblock
    if (n_show == nic) then
      write (iu_out, '(a,i0,a)') "# All ", nic, &
        " configurations (ascending Q; ties rank by configuration index)"
    else
      write (iu_out, '(a,i0,a,i0,a)') "# Top ", n_show, " of ", nic, &
        " configurations (ascending Q; ties rank by configuration index)"
    end if
    write (iu_out, '(a)', advance='no') "# Rank  Config   Degen      L(A)          Error              Q"
    do f = 1, n_families
      write (iu_out, '(1x,a,i0,i0)', advance='no') "E", families(f)%target_i, families(f)%target_j
    end do
    write (iu_out, *)
    do rank = 1, n_show
      write (iu_out, '(i6,1x,i8,1x,i7,1x,f12.4,1x,es15.7,1x,es15.7)', advance='no') &
        rank, ridx(rank), degen(ridx(rank)), L_val(ridx(rank)), abs_residual(ridx(rank)), scores(ridx(rank))
      do f = 1, n_families
        write (iu_out, '(1x,es15.7)', advance='no') family_error(ridx(rank), f)
      end do
      write (iu_out, *)
    end do
    close (iu_out)
    write (*, '(a,a)') " > Results written to ", trim(path_out)

    path_corr = trim(ensemble_dir) // "/SQS_CORRELATIONS"
    iu_corr = 31
    open (unit=iu_corr, file=trim(path_corr))
    write (iu_corr, '(a)') "# Detailed species-resolved pair correlations for best SQS"
    write (iu_corr, '(a,i0)') "# Config ", best_cfg
    write (iu_corr, '(a)') "# Config Orbit Distance(A) Target_i Target_j Species_i Species_j P_cell P_random Delta"
    do p = 1, n_pair_orbits
      do ch = channel_start(p), channel_end(p)
        write (iu_corr, '(i8,1x,i6,1x,f12.5,1x,i8,1x,i8,1x,a10,1x,a10,1x,f12.8,1x,f12.8,1x,es13.5)') &
          best_cfg, p, pair_orbits(p)%distance, pair_orbits(p)%target_i, pair_orbits(p)%target_j, &
          trim(targets(pair_orbits(p)%target_i)%species(channels(ch)%species_i)), &
          trim(targets(pair_orbits(p)%target_j)%species(channels(ch)%species_j)), &
          corr_all(best_cfg, ch), channels(ch)%target_corr, corr_all(best_cfg, ch) - channels(ch)%target_corr
      end do
    end do
    close (iu_corr)
    write (*, '(a,a)') " > Correlations written to ", trim(path_corr)
  end subroutine report_results

  ! ================================================================
  ! Rank configurations by ascending score, breaking ties by
  ! configuration index. n_top_sqs = 0 sorts the whole ensemble;
  ! otherwise only the n_top_sqs best are selected in one pass.
  ! ridx(1:nsel) holds the result in rank order.
  ! ================================================================
  subroutine rank_configs(ridx, nsel)
    implicit none
    integer, allocatable, intent(out) :: ridx(:)
    integer, intent(out) :: nsel
    integer :: ic

    if (n_top_sqs == 0) then
      nsel = nic
      allocate (ridx(nic))
      do ic = 1, nic
        ridx(ic) = ic
      end do
      call qsort_idx(ridx, 1, nic, scores)
    else
      nsel = min(nic, n_top_sqs)
      allocate (ridx(nsel))
      call select_top_configs(ridx, nsel)
    end if
  end subroutine rank_configs

  ! ================================================================
  ! Select the ncap best configurations by ascending score without
  ! sorting the full ensemble (insertion into a small sorted list).
  ! ================================================================
  subroutine select_top_configs(ridx, ncap)
    implicit none
    integer, intent(in) :: ncap
    integer, intent(out) :: ridx(ncap)
    integer :: ic, pos, m, nsel

    ridx = 0
    nsel = 0
    do ic = 1, nic
      if (nsel == ncap) then
        if (.not. idx_before(ic, ridx(nsel), scores)) cycle
        nsel = nsel - 1
      end if
      pos = nsel
      do while (pos >= 1)
        if (.not. idx_before(ic, ridx(pos), scores)) exit
        pos = pos - 1
      end do
      do m = nsel, pos + 1, -1
        ridx(m + 1) = ridx(m)
      end do
      ridx(pos + 1) = ic
      nsel = nsel + 1
    end do
  end subroutine select_top_configs

  ! ================================================================
  ! Quicksort index array by (score, index) ascending; used only for
  ! the n_top_sqs = 0 full-ranking mode.
  ! ================================================================
  recursive subroutine qsort_idx(idx, lo, hi, key)
    implicit none
    integer, intent(inout) :: idx(:)
    integer, intent(in) :: lo, hi
    real(real64), intent(in) :: key(:)
    integer :: i, j, tmp, pidx

    if (lo >= hi) return

    pidx = idx(lo + (hi - lo) / 2)
    i = lo - 1
    j = hi + 1
    do
      do
        i = i + 1
        if (.not. idx_before(idx(i), pidx, key)) exit
      end do
      do
        j = j - 1
        if (.not. idx_before(pidx, idx(j), key)) exit
      end do
      if (i >= j) exit
      tmp = idx(i); idx(i) = idx(j); idx(j) = tmp
    end do
    call qsort_idx(idx, lo, j, key)
    call qsort_idx(idx, j + 1, hi, key)
  end subroutine qsort_idx

  logical function idx_before(ia, ib, key) result(res)
    implicit none
    integer, intent(in) :: ia, ib
    real(real64), intent(in) :: key(:)
    res = key(ia) < key(ib) .or. (key(ia) == key(ib) .and. ia < ib)
  end function idx_before

  logical function pair_lex_lt(a1, a2, b1, b2) result(res)
    implicit none
    integer, intent(in) :: a1, a2, b1, b2
    res = .false.
    if (a1 < b1) then
      res = .true.
    else if (a1 == b1 .and. a2 < b2) then
      res = .true.
    end if
  end function pair_lex_lt

end program sqssod
