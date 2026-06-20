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
  implicit none

  ! ================================================================
  ! Parameters
  ! ================================================================
  integer, parameter :: max_k_param = 6       ! max supported cluster order
  integer, parameter :: natmax = 10000        ! max atoms in supercell
  integer, parameter :: max_clusters_param = 20000

  ! ================================================================
  ! INSQS input
  ! ================================================================
  character(len=256) :: ensemble_dir
  character(len=10)  :: target_species
  integer :: max_order
  real(real64) :: cutoff_r(max_k_param)
  real(real64) :: weight_k(max_k_param)
  real(real64) :: omega                        ! weight for matched-diameter term
  real(real64) :: eps_tol                      ! tolerance for "exact match"

  ! ================================================================
  ! Symmetry data (from EQMATRIX)
  ! ================================================================
  integer :: nop, npos
  integer, allocatable :: eqmt(:,:)           ! (npos, nop)

  ! ================================================================
  ! Supercell geometry (from supercell.cif)
  ! ================================================================
  real(real64) :: sc_a, sc_b, sc_c, sc_alpha, sc_beta, sc_gamma
  real(real64) :: cellvec(3, 3)
  integer :: nat_total, atini
  real(real64), allocatable :: tcoords(:,:)   ! (npos, 3) fractional

  ! ================================================================
  ! Pairwise distance matrix for target sites
  ! ================================================================
  real(real64), allocatable :: distmat(:,:)   ! (npos, npos)

  ! ================================================================
  ! Cluster types
  ! ================================================================
  type :: cluster_t
    integer :: order
    integer :: n_inst
    integer, allocatable :: inst(:,:)         ! (order, n_inst)
    real(real64) :: target_corr
    real(real64) :: w                         ! weight in score
    real(real64) :: char_dist                 ! max pairwise distance
    real(real64) :: eps_min                   ! min achievable |Δρ| in this supercell
  end type cluster_t

  integer :: n_clusters
  type(cluster_t), allocatable :: clusters(:)
  integer, allocatable :: pi_order(:)           ! pi_order(j) = cluster index for Pi_j output

  ! ================================================================
  ! Configuration data (from ENSEMBLE)
  ! ================================================================
  integer :: nic, nsubs
  integer, allocatable :: indconf(:,:)        ! (nic, nsubs) global site indices
  integer, allocatable :: degen(:)            ! (nic)
  real(real64) :: x_comp                      ! composition of substituted species

  ! ================================================================
  ! Scores and correlations
  ! ================================================================
  real(real64), allocatable :: scores(:)        ! (nic) — sort key for ranking
  real(real64), allocatable :: corr_all(:,:)    ! (nic, n_clusters)
  real(real64), allocatable :: L_val(:)         ! (nic) matched diameter L, Angstroms
  real(real64), allocatable :: abs_residual(:)  ! (nic) normalized weighted mean of |Δρ|
  real(real64), allocatable :: eps_min_cfg(:)   ! (n_clusters) min |Δρ| achieved over all configs

  ! ================================================================
  ! Main program
  ! ================================================================
  write (*, '(A)') "SOD (Site-Occupancy Disorder) version 0.82 - sqssod"

  call read_insqs()
  call read_eqmatrix()
  call read_insod()
  call read_supercell_cif()
  call compute_distances()
  call generate_clusters()
  call read_ensemble()
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
    logical :: insqs_exists
    iu = 20

    ! Get ENSEMBLE directory from command-line argument (default ".")
    ensemble_dir = "."
    target_species = ""
    if (command_argument_count() >= 1) then
      call get_command_argument(1, ensemble_dir)
      ensemble_dir = adjustl(ensemble_dir)
    end if

    ! INSQS: nXX/ takes priority over SODPROJECT/
    insqs_path = trim(ensemble_dir) // "/INSQS"
    inquire (file=trim(insqs_path), exist=insqs_exists)
    if (.not. insqs_exists) insqs_path = "INSQS"

    open (unit=iu, file=trim(insqs_path), status='old', IOSTAT=ios)
    if (ios /= 0) then
      write (*, '(a)') "Error: INSQS not found in nXX/ or SODPROJECT/."
      stop 1
    end if

    ! Block 1: Maximum cluster order
    read (iu, *)                                       ! comment
    read (iu, *) max_order                             ! data
    if (max_order < 2 .or. max_order > max_k_param) then
      write (*, '(a,i0,a,i0)') " Error: MaxOrder must be 2..", max_k_param, &
        ", got ", max_order
      stop 1
    end if

    ! Block 2: Cutoff radii for orders 2..MaxOrder (Angstroms)
    read (iu, *)                                       ! blank
    read (iu, *)                                       ! comment
    read (iu, *) (cutoff_r(k), k = 2, max_order)      ! data

    ! Block 3: Weights for orders 2..MaxOrder (order 1 is always composition, fixed)
    read (iu, *)                                       ! blank
    read (iu, *)                                       ! comment
    weight_k(1) = 0.0_real64
    read (iu, *) (weight_k(k), k = 2, max_order)      ! data

    ! Block 4: omega and eps_tol
    read (iu, *)                                       ! blank
    read (iu, *)                                       ! comment
    read (iu, *) omega, eps_tol                        ! data

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
    write (*, *)
  end subroutine read_insqs

  ! ================================================================
  ! Read INSOD to extract target species
  ! ================================================================
  subroutine read_insod()
    implicit none
    integer :: iu, ios, nsp, sptarget, isp
    character(len=256) :: line
    character(len=10) :: symbols(10)

    iu = 19
    open (unit=iu, file="INSOD", status='old', IOSTAT=ios)
    if (ios /= 0) then
      write (*, *) "Error: INSOD file not found."
      stop 1
    end if

    ! Skip to "nsp: Number of species in the parent structure"
    do
      read (iu, '(A)') line
      if (index(line, '# nsp:') > 0) exit
    end do
    read (iu, *) nsp

    ! Skip to "symbol(1:nsp): Atom symbols"
    do
      read (iu, '(A)') line
      if (index(line, '# symbol(1:nsp):') > 0) exit
    end do
    read (iu, *) (symbols(isp), isp = 1, nsp)

    ! Skip to "sptarget: Species to be substituted"
    do
      read (iu, '(A)') line
      if (index(line, '# sptarget:') > 0) exit
    end do
    read (iu, *) sptarget

    close (iu)

    if (sptarget < 1 .or. sptarget > nsp) then
      write (*, '(a,i0,a,i0)') " Error: sptarget=", sptarget, &
        " out of range [1..", nsp, "]"
      stop 1
    end if

    target_species = symbols(sptarget)

    write (*, '(a,a)') " > Target species (from INSOD): ", trim(target_species)
    write (*, *)
  end subroutine read_insod

  ! ================================================================
  ! Read EQMATRIX (symmetry permutation table for target-species sites)
  ! ================================================================
  subroutine read_eqmatrix()
    implicit none
    integer :: iu, ios, iop, j
    iu = 21
    open (unit=iu, file="EQMATRIX", status='old', IOSTAT=ios)
    if (ios /= 0) then
      write (*, *) "Error: EQMATRIX not found. Run combsod first."
      stop 1
    end if

    ! For single-target binary: one block (nop, npos) then nop lines.
    ! Multi-target support: would need to read and skip to the correct block.
    read (iu, *) nop, npos

    allocate (eqmt(npos, nop))
    do iop = 1, nop
      read (iu, *) (eqmt(j, iop), j = 1, npos)
    end do
    close (iu)

    write (*, '(a,i0,a,i0,a)') " > EQMATRIX: ", nop, " operators, ", npos, " sites"
    write (*, *)
  end subroutine read_eqmatrix

  ! ================================================================
  ! Read supercell.cif — cell parameters and target-species coordinates
  ! ================================================================
  subroutine read_supercell_cif()
    implicit none
    integer :: iu, ios, iat, itarg
    character(len=256) :: line
    character(len=20) :: label, typesym
    real(real64) :: fx, fy, fz
    integer :: all_nat
    character(len=10) :: all_sym(natmax)
    real(real64) :: all_frac(natmax, 3)
    logical :: in_atoms

    iu = 22
    open (unit=iu, file="supercell.cif", status='old', IOSTAT=ios)
    if (ios /= 0) then
      write (*, *) "Error: supercell.cif not found. Run combsod first."
      stop 1
    end if

    sc_a = 0; sc_b = 0; sc_c = 0
    sc_alpha = 90; sc_beta = 90; sc_gamma = 90
    all_nat = 0
    in_atoms = .false.

    do
      read (iu, '(A)', IOSTAT=ios) line
      if (ios /= 0) exit
      line = adjustl(line)

      if (.not. in_atoms) then
        if (line(1:14) == '_cell_length_a') then
          read (line(15:), *) sc_a
        else if (line(1:14) == '_cell_length_b') then
          read (line(15:), *) sc_b
        else if (line(1:14) == '_cell_length_c') then
          read (line(15:), *) sc_c
        else if (line(1:17) == '_cell_angle_alpha') then
          read (line(18:), *) sc_alpha
        else if (line(1:16) == '_cell_angle_beta') then
          read (line(17:), *) sc_beta
        else if (line(1:17) == '_cell_angle_gamma') then
          read (line(18:), *) sc_gamma
        else if (line(1:18) == '_atom_site_fract_z') then
          in_atoms = .true.
        end if
      else
        if (len_trim(line) == 0) exit
        read (line, *, IOSTAT=ios) label, typesym, fx, fy, fz
        if (ios /= 0) exit
        all_nat = all_nat + 1
        all_sym(all_nat) = adjustl(typesym)
        all_frac(all_nat, 1) = fx
        all_frac(all_nat, 2) = fy
        all_frac(all_nat, 3) = fz
      end if
    end do
    close (iu)

    nat_total = all_nat

    ! Identify target-species atoms and determine atini
    atini = 0
    itarg = 0
    do iat = 1, nat_total
      if (trim(all_sym(iat)) == trim(target_species)) then
        itarg = itarg + 1
        if (itarg == 1) atini = iat
      end if
    end do

    if (itarg /= npos) then
      write (*, '(a,i0,a,a,a,i0)') " Error: found ", itarg, " atoms of ", &
        trim(target_species), " but EQMATRIX expects ", npos
      stop 1
    end if

    ! Extract fractional coordinates of target-species sites
    allocate (tcoords(npos, 3))
    itarg = 0
    do iat = 1, nat_total
      if (trim(all_sym(iat)) == trim(target_species)) then
        itarg = itarg + 1
        tcoords(itarg, :) = all_frac(iat, :)
      end if
    end do

    call cell(cellvec, sc_a, sc_b, sc_c, sc_alpha, sc_beta, sc_gamma)

    write (*, '(a,3f10.4,3f9.2)') " > Cell: ", sc_a, sc_b, sc_c, &
      sc_alpha, sc_beta, sc_gamma
    write (*, '(a,i0,a,i0,a,a,a,i0)') " > Atoms: ", nat_total, " total, ", &
      npos, " ", trim(target_species), " starting at index ", atini
    write (*, *)
  end subroutine read_supercell_cif

  ! ================================================================
  ! Compute minimum-image distance matrix between target-species sites
  ! ================================================================
  subroutine compute_distances()
    implicit none
    integer :: ia, ja
    real(real64) :: df(3), dr(3), dd

    allocate (distmat(npos, npos))
    distmat = 0.0_real64

    do ia = 1, npos
      do ja = ia + 1, npos
        df(:) = tcoords(ia, :) - tcoords(ja, :)
        df(:) = df(:) - nint(df(:))
        dr = matmul(cellvec, df)
        dd = sqrt(sum(dr**2))
        distmat(ia, ja) = dd
        distmat(ja, ia) = dd
      end do
    end do

    write (*, '(a,f8.4,a)') " > Distances: min = ", &
      minval(distmat, mask = (distmat > 0.0_real64)), " A"
    write (*, '(a,f8.4,a)') "              max = ", maxval(distmat), " A"
    write (*, *)
  end subroutine compute_distances

  ! ================================================================
  ! Generate cluster types by internal enumeration of k-tuples,
  ! geometric filtering, and symmetry reduction via EQMATRIX
  ! ================================================================
  subroutine generate_clusters()
    implicit none
    integer :: kk, iop, ia, ja, jc
    integer :: tuple_k(max_k_param), trans_k(max_k_param)
    real(real64) :: max_pd
    logical :: more_sub, is_canonical, is_dup
    integer, allocatable :: temp_inst(:,:)
    integer :: n_temp
    integer :: n_per_order(max_k_param)

    allocate (clusters(max_clusters_param))
    allocate (temp_inst(max_k_param, nop))
    n_clusters = 0
    n_per_order = 0

    ! ----------------------------------------------------------
    ! Order 1: one type — all single sites form one orbit.
    ! Singlet correlation equals (1-2x) for every configuration
    ! at the same composition, so it cannot distinguish configs.
    ! Included only when the user assigns nonzero weight.
    ! ----------------------------------------------------------
    if (weight_k(1) > 0.0_real64) then
      n_clusters = 1
      clusters(1)%order = 1
      clusters(1)%n_inst = npos
      allocate (clusters(1)%inst(1, npos))
      do ia = 1, npos
        clusters(1)%inst(1, ia) = ia
      end do
      clusters(1)%char_dist = 0.0_real64
      clusters(1)%w = weight_k(1)
      n_per_order(1) = 1
    end if

    ! ----------------------------------------------------------
    ! Orders 2..max_order: enumerate all k-subsets of sites,
    ! apply distance cutoff, keep only canonical representatives
    ! (lex-smallest member of each symmetry orbit), then generate
    ! the full orbit for each retained type.
    ! ----------------------------------------------------------
    do kk = 2, max_order
      more_sub = .false.
      do
        call ksubset(npos, kk, tuple_k, more_sub)

        ! Distance filter: max pairwise distance within the tuple
        max_pd = 0.0_real64
        do ia = 1, kk - 1
          do ja = ia + 1, kk
            if (distmat(tuple_k(ia), tuple_k(ja)) > max_pd) &
              max_pd = distmat(tuple_k(ia), tuple_k(ja))
          end do
        end do
        if (max_pd > cutoff_r(kk)) then
          if (.not. more_sub) exit
          cycle
        end if

        ! Canonical check: is this the lex-smallest tuple in its orbit?
        ! Since ksubset generates tuples in lex order, the first tuple
        ! encountered from each orbit IS the canonical representative.
        ! Non-canonical tuples have at least one symmetry transform
        ! that is lex-smaller, so we detect and skip them early.
        is_canonical = .true.
        do iop = 1, nop
          do ia = 1, kk
            trans_k(ia) = eqmt(tuple_k(ia), iop)
          end do
          call isort(trans_k, kk)
          if (ilex_lt(trans_k, tuple_k, kk)) then
            is_canonical = .false.
            exit
          end if
        end do

        if (is_canonical) then
          ! Generate the full orbit: all distinct transforms under symmetry
          n_temp = 0
          do iop = 1, nop
            do ia = 1, kk
              trans_k(ia) = eqmt(tuple_k(ia), iop)
            end do
            call isort(trans_k, kk)

            is_dup = .false.
            do jc = 1, n_temp
              if (all(temp_inst(1:kk, jc) == trans_k(1:kk))) then
                is_dup = .true.
                exit
              end if
            end do
            if (.not. is_dup) then
              n_temp = n_temp + 1
              temp_inst(1:kk, n_temp) = trans_k(1:kk)
            end if
          end do

          ! Store cluster type
          n_clusters = n_clusters + 1
          if (n_clusters > max_clusters_param) then
            write (*, *) "Error: too many cluster types. Increase max_clusters_param."
            stop 1
          end if
          clusters(n_clusters)%order = kk
          clusters(n_clusters)%n_inst = n_temp
          allocate (clusters(n_clusters)%inst(kk, n_temp))
          clusters(n_clusters)%inst(1:kk, 1:n_temp) = temp_inst(1:kk, 1:n_temp)
          clusters(n_clusters)%char_dist = max_pd
          clusters(n_clusters)%w = weight_k(kk)
          n_per_order(kk) = n_per_order(kk) + 1
        end if

        if (.not. more_sub) exit
      end do
    end do

    deallocate (temp_inst)

    write (*, *) " > Cluster types generated:"
    do kk = 1, max_order
      if (n_per_order(kk) > 0) &
        write (*, '(a,i0,a,i0,a)') "    Order ", kk, ": ", n_per_order(kk), " types"
    end do
    write (*, '(a,i0)') "    Total: ", n_clusters
    write (*, *)

    if (n_clusters == 0) then
      write (*, *) "Warning: no cluster types survived the distance filter."
      write (*, *) "         All configurations will have zero score."
      write (*, *)
      return
    end if

    ! Build pi_order: sort by (char_dist ascending, order ascending)
    allocate (pi_order(n_clusters))
    do jc = 1, n_clusters
      pi_order(jc) = jc
    end do
    call sort_pi_order(pi_order, n_clusters)

    write (*, '(a)') "   Pi   Order  Instances  MaxDist(A)  Weight"
    write (*, '(a)') "  ----  -----  ---------  ----------  ------"
    do jc = 1, n_clusters
      write (*, '(i6, i7, i11, f12.4, f9.4)') jc, clusters(pi_order(jc))%order, &
        clusters(pi_order(jc))%n_inst, clusters(pi_order(jc))%char_dist, &
        clusters(pi_order(jc))%w
    end do
    write (*, *)
  end subroutine generate_clusters

  ! ================================================================
  ! Sort pi_order by (char_dist ascending, order ascending)
  ! ================================================================
  subroutine sort_pi_order(idx, nn)
    implicit none
    integer, intent(in) :: nn
    integer, intent(inout) :: idx(nn)
    integer :: ia, ja, imin, tmp
    real(real64) :: dmin
    integer :: omin
    do ia = 1, nn - 1
      imin = ia
      dmin = clusters(idx(ia))%char_dist
      omin = clusters(idx(ia))%order
      do ja = ia + 1, nn
        if (clusters(idx(ja))%char_dist < dmin - 1.0e-8_real64 .or. &
            (abs(clusters(idx(ja))%char_dist - dmin) < 1.0e-8_real64 .and. &
             clusters(idx(ja))%order < omin)) then
          imin = ja
          dmin = clusters(idx(ja))%char_dist
          omin = clusters(idx(ja))%order
        end if
      end do
      if (imin /= ia) then
        tmp = idx(ia); idx(ia) = idx(imin); idx(imin) = tmp
      end if
    end do
  end subroutine sort_pi_order

  ! ================================================================
  ! Read ENSEMBLE — configurations to score (supports v2 and v3 formats)
  ! ================================================================
  subroutine read_ensemble()
    implicit none
    integer :: iu, ios, ic, m_idx, npos_check
    character(len=256) :: line, path
    character(len=20) :: word1

    path = trim(ensemble_dir) // "/ENSEMBLE"
    iu = 23
    open (unit=iu, file=trim(path), status='old', IOSTAT=ios)
    if (ios /= 0) then
      write (*, '(a,a)') " Error: cannot open ", trim(path)
      stop 1
    end if

    ! Read first non-blank line to detect v2 vs v3
    do
      read (iu, '(A)', iostat=ios) line
      if (ios /= 0) then
        write (*, *) " Error: cannot read ENSEMBLE."
        stop 1
      end if
      if (len_trim(line) > 0) exit
    end do

    if (line(1:1) == '#') then
      ! v2: skip comment lines
      do while (line(1:1) == '#')
        read (iu, '(A)', iostat=ios) line
        if (ios /= 0) then
          write (*, *) " Error: cannot read ENSEMBLE."
          stop 1
        end if
      end do
      ! line is now "nsubs substitutions in npos sites"
      read (line, *) nsubs, word1
      if (trim(word1) /= 'substitutions') then
        write (*, *) "Error: only binary single-target ENSEMBLE is currently supported."
        stop 1
      end if
      m_idx = index(line, ' in ')
      read (line(m_idx + 4:), *) npos_check
      if (npos_check /= npos) then
        write (*, '(a,i0,a,i0)') " Error: ENSEMBLE npos=", npos_check, &
          " does not match EQMATRIX npos=", npos
        stop 1
      end if
      read (iu, *) nic
    else
      ! v3: first line is "... ensemble ...: nic configurations"
      block
        integer :: cp, kp
        cp = index(line, ':', back=.true.)
        kp = index(line, 'configurations')
        if (cp > 0 .and. kp > cp) then
          read (line(cp+1:kp-1), *, iostat=ios) nic
          if (ios /= 0) then
            write (*, *) " Error: cannot parse ENSEMBLE configuration count."
            stop 1
          end if
        else
          write (*, *) " Error: cannot parse ENSEMBLE header."
          stop 1
        end if
      end block
      ! Read target line(s) and column-header comment; verify binary single-target
      nsubs = -1
      npos_check = -1
      do
        read (iu, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) == 0) cycle
        if (line(1:1) == '#') exit   ! column-header comment — file now at first data row
        if (index(line, 'sites') > 0 .and. index(line, '->') > 0) then
          if (nsubs >= 0) then
            write (*, *) "Error: only binary single-target ENSEMBLE is currently supported."
            stop 1
          end if
          read (line, *, iostat=ios) npos_check
          m_idx = index(line, '->')
          read (line(m_idx+2:), *, iostat=ios) nsubs
        end if
      end do
      if (nsubs < 0) then
        write (*, *) " Error: no target line found in ENSEMBLE."
        stop 1
      end if
      if (npos_check /= npos) then
        write (*, '(a,i0,a,i0)') " Error: ENSEMBLE npos=", npos_check, &
          " does not match EQMATRIX npos=", npos
        stop 1
      end if
    end if

    allocate (degen(nic), indconf(nic, nsubs))
    do ic = 1, nic
      read (iu, *) m_idx, degen(ic), indconf(ic, 1:nsubs)
    end do
    close (iu)

    x_comp = real(nsubs, real64) / real(npos, real64)

    ! Set target correlations and minimum achievable errors
    do ic = 1, n_clusters
      clusters(ic)%target_corr = (1.0_real64 - 2.0_real64 * x_comp) ** clusters(ic)%order
      ! eps_min: smallest |Δρ| achievable in this supercell.
      ! Valid correlations: ρ = (2n₊ - m)/m, n₊ ∈ {0,...,m}, spacing = 2/m.
      ! Nearest n₊ to target: round(m*(1+target)/2), clamped to [0,m].
      block
        integer :: m, n_plus
        real(real64) :: rho_nearest
        m = clusters(ic)%n_inst
        n_plus = nint(real(m, real64) * (1.0_real64 + clusters(ic)%target_corr) / 2.0_real64)
        n_plus = max(0, min(m, n_plus))
        rho_nearest = real(2*n_plus - m, real64) / real(m, real64)
        clusters(ic)%eps_min = abs(rho_nearest - clusters(ic)%target_corr)
      end block
    end do

    write (*, '(a,i0,a,i0,a)') " > ENSEMBLE: ", nic, " configurations, ", &
      nsubs, " substitutions"
    write (*, '(a,f8.5)') "    Composition x = ", x_comp
    write (*, *)
  end subroutine read_ensemble

  ! ================================================================
  ! Score all configurations
  ! ================================================================
  subroutine score_configurations()
    implicit none
    integer :: ic, jc, js, ia2, local_idx
    integer :: sigma(npos)
    real(real64) :: prod_val, corr_val

    allocate (scores(nic), corr_all(nic, n_clusters))

    do ic = 1, nic
      ! Build sigma: +1 for original species, -1 for substituted
      sigma = 1
      do js = 1, nsubs
        local_idx = indconf(ic, js) - atini + 1
        if (local_idx < 1 .or. local_idx > npos) then
          write (*, '(a,i0,a,i0,a,i0)') " Error: config ", ic, ", site ", &
            indconf(ic, js), " maps to local index ", local_idx
          stop 1
        end if
        sigma(local_idx) = -1
      end do

      ! Evaluate cluster correlations
      do jc = 1, n_clusters
        corr_val = 0.0_real64
        do js = 1, clusters(jc)%n_inst
          prod_val = 1.0_real64
          do ia2 = 1, clusters(jc)%order
            prod_val = prod_val * sigma(clusters(jc)%inst(ia2, js))
          end do
          corr_val = corr_val + prod_val
        end do
        corr_val = corr_val / real(clusters(jc)%n_inst, real64)
        corr_all(ic, jc) = corr_val
      end do
    end do

    ! Compute eps_min_cfg: minimum |Δρ| actually achieved over all configurations.
    ! Used by van de Walle scoring as the match tolerance (so the best-achieving
    ! configuration for each cluster always counts as "matched").
    allocate (eps_min_cfg(n_clusters))
    do jc = 1, n_clusters
      eps_min_cfg(jc) = minval(abs(corr_all(:, jc) - clusters(jc)%target_corr))
    end do

    ! Print target correlations table (now that both eps values are available)
    write (*, '(a)') " > Target correlations and achievable errors:"
    write (*, '(a)') "   Pi   Order  MaxDist(A)      Target    eps_grid     eps_cfg"
    write (*, '(a)') "  ----  -----  ----------  ----------  ----------  ----------"
    do jc = 1, n_clusters
      associate (cl => clusters(pi_order(jc)))
        write (*, '(i6, i7, f12.4, f12.8, es12.3, es12.3)') jc, cl%order, &
          cl%char_dist, cl%target_corr, cl%eps_min, eps_min_cfg(pi_order(jc))
      end associate
    end do
    write (*, *)

    call compute_vdw_scores()

    write (*, '(a)') " > Scoring complete."
    write (*, *)
  end subroutine score_configurations

  ! ================================================================
  ! Van de Walle 2013 scoring: Q = -omega*L + sum(w_k*|Δρ|)/sum(w_k)
  ! L = largest diameter (Å) such that ALL clusters with
  ! char_dist <= L satisfy |Δρ| <= eps_min + eps_tol.
  ! Overwrites scores() with Q for ranking.
  ! ================================================================
  subroutine compute_vdw_scores()
    implicit none
    integer :: ic, jj, jc
    integer :: sorted_c(n_clusters)
    real(real64) :: dev_val, cur_diam, L_this, weight_norm

    allocate (L_val(nic), abs_residual(nic))

    ! Sort cluster indices by char_dist ascending
    do jc = 1, n_clusters
      sorted_c(jc) = jc
    end do
    call sort_clusters_by_dist(sorted_c, n_clusters)

    weight_norm = 0.0_real64
    do jc = 1, n_clusters
      weight_norm = weight_norm + clusters(jc)%w
    end do

    do ic = 1, nic
      abs_residual(ic) = 0.0_real64

      ! Determine L: sweep sorted clusters; L advances only while every
      ! cluster at the current diameter is matched within eps_min + eps_tol.
      L_this = 0.0_real64
      jj = 1
      do while (jj <= n_clusters)
        ! Collect all clusters sharing the same char_dist (same shell)
        cur_diam = clusters(sorted_c(jj))%char_dist
        ! Check all clusters in this shell
        do while (jj <= n_clusters)
          if (abs(clusters(sorted_c(jj))%char_dist - cur_diam) >= 1.0e-8_real64) exit
          jc = sorted_c(jj)
          dev_val = corr_all(ic, jc) - clusters(jc)%target_corr
          if (abs(dev_val) > eps_min_cfg(jc) + eps_tol) goto 10
          jj = jj + 1
        end do
        L_this = cur_diam   ! whole shell matched
      end do
10    continue

      L_val(ic) = L_this

      ! Normalized weighted absolute residual over all clusters
      do jj = 1, n_clusters
        jc = sorted_c(jj)
        abs_residual(ic) = abs_residual(ic) + &
          clusters(jc)%w * abs(corr_all(ic, jc) - clusters(jc)%target_corr)
      end do
      if (weight_norm > 0.0_real64) abs_residual(ic) = abs_residual(ic) / weight_norm

      scores(ic) = -omega * L_val(ic) + abs_residual(ic)
    end do
  end subroutine compute_vdw_scores

  ! ================================================================
  ! Sort cluster indices by characteristic distance (selection sort;
  ! n_clusters is small so this is fine)
  ! ================================================================
  subroutine sort_clusters_by_dist(idx, nn)
    implicit none
    integer, intent(in) :: nn
    integer, intent(inout) :: idx(nn)
    integer :: ia, ja, imin, tmp
    real(real64) :: dmin
    do ia = 1, nn - 1
      imin = ia
      dmin = clusters(idx(ia))%char_dist
      do ja = ia + 1, nn
        if (clusters(idx(ja))%char_dist < dmin) then
          imin = ja
          dmin = clusters(idx(ja))%char_dist
        end if
      end do
      if (imin /= ia) then
        tmp = idx(ia); idx(ia) = idx(imin); idx(imin) = tmp
      end if
    end do
  end subroutine sort_clusters_by_dist

  ! ================================================================
  ! Report results to stdout and write OUTSQS
  ! ================================================================
  subroutine report_results()
    implicit none
    integer :: ic, jc, n_show
    integer, allocatable :: ridx(:)
    integer :: iu_out
    character(len=256) :: path_out

    ! Build ranking (sort by ascending scores)
    allocate (ridx(nic))
    do ic = 1, nic
      ridx(ic) = ic
    end do
    call qsort_idx(ridx, 1, nic, scores)

    ! ----------------------------------------------------------
    ! Top configurations
    ! ----------------------------------------------------------
    n_show = min(nic, 20)
    write (*, '(a,i0,a)') " > Top ", n_show, " configurations:"
    write (*, *)
    write (*, '(a)') "    Rank  Config  Degen      L(A)       AbsErr                Q"
    write (*, '(a)') "    ----  ------  -----  --------  ----------  ---------------"
    do ic = 1, n_show
      write (*, '(i7, i8, i7, f10.4, es12.4, es16.6)') ic, ridx(ic), degen(ridx(ic)), &
        L_val(ridx(ic)), abs_residual(ridx(ic)), scores(ridx(ic))
    end do
    write (*, *)

    ! ----------------------------------------------------------
    ! Best SQS details
    ! ----------------------------------------------------------
    write (*, '(a,i0)') " > Best SQS: configuration ", ridx(1)
    write (*, '(a,i0)') "    Degeneracy:  ", degen(ridx(1))
    write (*, '(a,f10.4,a)') "    Matched L:   ", L_val(ridx(1)), " Angstroms"
    write (*, '(a,es14.6)') "    AbsErr:      ", abs_residual(ridx(1))
    write (*, '(a,es14.6)') "    Q (score):   ", scores(ridx(1))
    write (*, *)

    write (*, '(a)') "    Cluster correlations for best SQS:"
    write (*, '(a)') "   Pi   Order  Dist(A)       Target       Actual    Deviation"
    write (*, '(a)') "  ----  -----  --------  ----------  ----------  -----------"
    do jc = 1, n_clusters
      associate (cl => clusters(pi_order(jc)))
        write (*, '(i6, i7, f10.4, f12.6, f12.6, es13.4)') jc, &
          cl%order, cl%char_dist, cl%target_corr, &
          corr_all(ridx(1), pi_order(jc)), corr_all(ridx(1), pi_order(jc)) - cl%target_corr
      end associate
    end do
    write (*, *)

    write (*, '(a)', advance='no') "    Substituted sites: "
    do ic = 1, nsubs
      write (*, '(i0,1x)', advance='no') indconf(ridx(1), ic)
    end do
    write (*, *)
    write (*, *)

    ! ----------------------------------------------------------
    ! Write OUTSQS file
    ! ----------------------------------------------------------
    path_out = trim(ensemble_dir) // "/OUTSQS"
    iu_out = 30
    open (unit=iu_out, file=trim(path_out))
    write (iu_out, '(a)') "# SQS scoring results from sqssod"
    write (iu_out, '(a,i0,a,i0,a,f10.6,a,i0)') "# nsubs=", nsubs, &
      " npos=", npos, " x=", x_comp, " n_clusters=", n_clusters
    write (iu_out, '(a,es10.3,a,es10.3)') "# Mode: van_de_Walle  omega=", &
      omega, "  eps=", eps_tol
    ! Header: fixed columns + one 10-char label per Pi (sorted by dist, order)
    write (iu_out, '(a)', advance='no') "# Rank  Config   Degen      L(A)        AbsErr             Q"
    do jc = 1, n_clusters
      write (iu_out, '("    Pi_",i3)', advance='no') jc
    end do
    write (iu_out, *)
    ! Ideal disorder row (target correlations)
    write (iu_out, '("# Ideal disorder:",45x,*(f10.5))') &
      (clusters(pi_order(jc))%target_corr, jc=1, n_clusters)
    ! Ranked configuration rows
    do ic = 1, nic
      write (iu_out, '(i6, i8, i7, f12.4, es15.7, es15.7, *(f10.5))') ic, ridx(ic), degen(ridx(ic)), &
        L_val(ridx(ic)), abs_residual(ridx(ic)), scores(ridx(ic)), &
        (corr_all(ridx(ic), pi_order(jc)), jc=1, n_clusters)
    end do
    close (iu_out)
    write (*, '(a,a)') " > Results written to ", trim(path_out)

    deallocate (ridx)
  end subroutine report_results

  ! ================================================================
  ! Quicksort index array by key (ascending) — standard Hoare partition
  ! ================================================================
  recursive subroutine qsort_idx(idx, lo, hi, key)
    implicit none
    integer, intent(inout) :: idx(:)
    integer, intent(in) :: lo, hi
    real(real64), intent(in) :: key(:)
    integer :: i, j, tmp
    real(real64) :: pivot

    if (lo >= hi) return

    pivot = key(idx(lo))
    i = lo - 1
    j = hi + 1
    do
      do
        i = i + 1
        if (key(idx(i)) >= pivot) exit
      end do
      do
        j = j - 1
        if (key(idx(j)) <= pivot) exit
      end do
      if (i >= j) exit
      tmp = idx(i); idx(i) = idx(j); idx(j) = tmp
    end do
    call qsort_idx(idx, lo, j, key)
    call qsort_idx(idx, j + 1, hi, key)
  end subroutine qsort_idx

  ! ================================================================
  ! Sort small integer array in place (bubble sort — used for
  ! k-tuples with at most max_k_param ~ 6 elements)
  ! ================================================================
  subroutine isort(arr, nn)
    implicit none
    integer, intent(inout) :: arr(*)
    integer, intent(in) :: nn
    integer :: ia, ja, tmp
    do ia = 1, nn - 1
      do ja = ia + 1, nn
        if (arr(ja) < arr(ia)) then
          tmp = arr(ia); arr(ia) = arr(ja); arr(ja) = tmp
        end if
      end do
    end do
  end subroutine isort

  ! ================================================================
  ! Lexicographic less-than for integer arrays of length nn
  ! ================================================================
  function ilex_lt(a, b, nn) result(res)
    implicit none
    integer, intent(in) :: a(*), b(*), nn
    logical :: res
    integer :: ia
    res = .false.
    do ia = 1, nn
      if (a(ia) < b(ia)) then
        res = .true.; return
      else if (a(ia) > b(ia)) then
        return
      end if
    end do
  end function ilex_lt

end program sqssod
