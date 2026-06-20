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

program gqssod
  use iso_fortran_env, only: real64
  use ensemble_io,     only: read_energies_file
  implicit none

  integer, parameter :: max_k_param = 6, natmax = 10000, max_clusters_param = 20000
  integer, parameter :: nclusmax = 1000, ntemp_max = 1000

  ! Input parameters
  character(len=256) :: ensemble_dir
  character(len=10) :: target_species
  integer :: max_order
  real(real64) :: cutoff_r(max_k_param)
  real(real64) :: weight_k(max_k_param)
  real(real64) :: omega
  real(real64) :: eps_tol

  ! Symmetry data
  integer :: nop, npos
  integer, allocatable :: eqmt(:,:)

  ! Supercell geometry
  real(real64) :: sc_a, sc_b, sc_c, sc_alpha, sc_beta, sc_gamma
  real(real64) :: cellvec(3, 3)
  integer :: nat_total, atini
  real(real64), allocatable :: tcoords(:,:)

  ! Distance matrix
  real(real64), allocatable :: distmat(:,:)

  ! Clusters
  type :: cluster_t
    integer :: order
    integer :: n_inst
    integer, allocatable :: inst(:,:)
    real(real64) :: target_corr
    real(real64) :: w
    real(real64) :: char_dist
    real(real64) :: eps_min
  end type cluster_t

  integer :: n_clusters
  type(cluster_t), allocatable :: clusters(:)
  integer, allocatable :: pi_order(:)
  integer, allocatable :: cluster_active(:)  ! maps cluster index to OUTSQS column (0 if weight=0)

  ! Configuration data
  integer :: nic, nsubs
  integer, allocatable :: indconf(:,:)
  integer, allocatable :: degen(:)
  real(real64) :: x_comp

  ! Scores and correlations
  real(real64), allocatable :: scores(:)
  real(real64), allocatable :: corr_all(:,:)
  real(real64), allocatable :: L_val(:)
  real(real64), allocatable :: abs_residual(:)
  real(real64), allocatable :: eps_min_cfg(:)
  real(real64), allocatable :: thermal_avg(:,:)

  ! Thermal data
  integer :: ntemp, tt
  real(real64) :: temperature
  real(real64), allocatable :: temps(:)
  real(real64), allocatable :: ene(:)
  real(real64), allocatable :: avg_target(:)
  real(real64), parameter :: kb = 8.61734e-5_real64

  write (*, '(A)') "SOD (Site-Occupancy Disorder) version 0.82 - gqssod"

  call read_insqs()
  call read_eqmatrix()
  call read_insod()
  call read_supercell_cif()
  call compute_distances()
  call generate_clusters()
  call read_ensemble()
  call read_temperatures()
  call read_energies_and_correlations()

  ! Compute and store thermal averages for all temperatures
  allocate(thermal_avg(ntemp, n_clusters))
  do tt = 1, ntemp
    temperature = temps(tt)
    call set_thermal_targets(temperature)
    thermal_avg(tt, 1:n_clusters) = avg_target(1:n_clusters)
    call score_configurations()
    call report_results(temperature, tt)
  end do

  ! Write thermal averages to file
  call write_thermal_averages()

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

    read (iu, *)
    read (iu, *) max_order
    if (max_order < 2 .or. max_order > max_k_param) then
      write (*, '(a,i0,a,i0)') " Error: MaxOrder must be 2..", max_k_param, &
        ", got ", max_order
      stop 1
    end if

    read (iu, *)
    read (iu, *)
    read (iu, *) (cutoff_r(k), k = 2, max_order)

    read (iu, *)
    read (iu, *)
    weight_k(1) = 0.0_real64
    read (iu, *) (weight_k(k), k = 2, max_order)

    read (iu, *)
    read (iu, *)
    read (iu, *) omega, eps_tol

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

    do
      read (iu, '(A)') line
      if (index(line, '# nsp:') > 0) exit
    end do
    read (iu, *) nsp

    do
      read (iu, '(A)') line
      if (index(line, '# symbol(1:nsp):') > 0) exit
    end do
    read (iu, *) (symbols(isp), isp = 1, nsp)

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
  ! Read EQMATRIX
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
  ! Read supercell.cif
  ! ================================================================
  subroutine read_supercell_cif()
    implicit none
    integer :: iu, ios, iat, itarg, all_nat
    character(len=256) :: line
    character(len=20) :: label, typesym
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
  ! Compute distance matrix
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
  ! Generate clusters (copy from sqssod)
  ! ================================================================
  subroutine generate_clusters()
    implicit none
    integer :: kk, iop, ia, ja, jc, k
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

    ! Order 1: all single sites
    n_clusters = 1
    allocate (clusters(1)%inst(1, npos))
    clusters(1)%order = 1
    clusters(1)%n_inst = npos
    do ia = 1, npos
      clusters(1)%inst(1, ia) = ia
    end do
    clusters(1)%char_dist = 0.0_real64
    clusters(1)%w = weight_k(1)
    n_per_order(1) = 1

    ! Orders 2 to max_order
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

    if (n_clusters == 0) then
      write (*, *) "Warning: no cluster types survived the distance filter."
      write (*, *)
      return
    end if

    allocate (pi_order(n_clusters))
    do jc = 1, n_clusters
      pi_order(jc) = jc
    end do
    call sort_pi_order(pi_order, n_clusters)

    write (*, '(a)') " > Cluster types generated:"
    do kk = 1, max_order
      if (n_per_order(kk) > 0) &
        write (*, '(a,i0,a,i0,a)') "    Order ", kk, ": ", n_per_order(kk), " types"
    end do
    write (*, '(a,i0)') "    Total: ", n_clusters
    write (*, *)

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
  ! Read ENSEMBLE (supports v2 and v3 formats)
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

    write (*, '(a,i0,a,i0,a)') " > ENSEMBLE: ", nic, " configurations, ", &
      nsubs, " substitutions"
    write (*, '(a,f8.5)') "    Composition x = ", x_comp
    write (*, *)
  end subroutine read_ensemble

  ! ================================================================
  ! Read TEMPERATURES (from nXX/ if present, otherwise from SODPROJECT/)
  ! ================================================================
  subroutine read_temperatures()
    implicit none
    integer :: iu, ios
    logical :: exists
    character(len=256) :: temps_path

    ! Try nXX/ first, then SODPROJECT/
    inquire (file=trim(ensemble_dir) // "/TEMPERATURES", exist=exists)
    if (exists) then
      temps_path = trim(ensemble_dir) // "/TEMPERATURES"
    else
      inquire (file="TEMPERATURES", exist=exists)
      if (.not. exists) then
        write (*, *) "Error: TEMPERATURES file not found in ", trim(ensemble_dir), &
          " or in SODPROJECT/."
        stop 1
      end if
      temps_path = "TEMPERATURES"
    end if

    allocate (temps(ntemp_max))
    iu = 24
    open (unit=iu, file=trim(temps_path))
    ntemp = 0
    do
      read (iu, *, iostat=ios) temps(ntemp + 1)
      if (ios /= 0) exit
      ntemp = ntemp + 1
      if (ntemp >= ntemp_max) exit
    end do
    close (iu)

    write (*, '(a,i0)') " > TEMPERATURES: ", ntemp, " points"
    write (*, *)
  end subroutine read_temperatures

  ! ================================================================
  ! Read ENERGIES and OUTSQS correlations
  ! ================================================================
  subroutine read_energies_and_correlations()
    implicit none
    integer :: ic, jc, ios, idum, idum2, idum3, ic_config
    integer :: outsqs_n_clusters, jc_outsqs, jc_gen
    integer :: n_missing
    character(len=256) :: line
    character(len=2048) :: data_line
    logical :: exists, ene_ok
    real(real64) :: rdum, rdum2, rdum3
    real(real64), allocatable :: temp_corr(:)

    allocate (ene(nic))
    allocate (corr_all(nic, n_clusters))
    allocate (avg_target(n_clusters))

    ! Read ENERGIES from nXX/ folder
    call read_energies_file(trim(ensemble_dir)//'/ENERGIES', nic, ene, ene_ok, n_missing)
    if (.not. ene_ok) then
      write (*, '(a,a)') " Error: could not read ENERGIES from ", trim(ensemble_dir)
      stop 1
    end if
    if (n_missing > 0) then
      write (*, '(a,I0,a,a)') " Error: missing energies for ", n_missing, &
        " configuration(s) in ", trim(ensemble_dir)//'/ENERGIES'
      stop 1
    end if

    ! Map cluster indices to OUTSQS columns (only non-zero weight clusters are in OUTSQS)
    allocate (cluster_active(n_clusters))
    cluster_active = 0
    outsqs_n_clusters = 0
    do jc_gen = 1, n_clusters
      if (clusters(jc_gen)%w > 0.0_real64) then
        outsqs_n_clusters = outsqs_n_clusters + 1
        cluster_active(jc_gen) = outsqs_n_clusters
      end if
    end do
    ! Initialize corr_all to zero for zero-weight clusters
    corr_all = 0.0_real64

    ! Read OUTSQS correlations (skip all comment lines, read nic config rows)
    open (unit=26, file=trim(ensemble_dir) // "/OUTSQS", status='old', IOSTAT=ios)
    if (ios /= 0) then
      write (*, '(a,a)') " Error: cannot open ", trim(ensemble_dir) // "/OUTSQS"
      stop 1
    end if
    ! Skip comment lines; 'line' holds the first config row after the loop
    do
      read (26, '(a)', IOSTAT=ios) line
      if (ios /= 0) then
        write (*, '(a)') " Error: no data rows found in OUTSQS"
        stop 1
      end if
      if (line(1:1) /= '#' .and. len_trim(line) > 0) exit
    end do

    ! Read configuration correlations and map to full cluster array
    allocate (temp_corr(outsqs_n_clusters))
    do ic = 1, nic
      if (ic == 1) then
        data_line = line    ! first config row already read by the comment-skip loop
      else
        read (26, '(a)', iostat=ios) data_line
        if (ios /= 0) then
          write (*, '(a,i0)') " Error reading line for config ", ic
          stop 1
        end if
      end if
      ! rank config degen L AbsErr Q Pi...
      read (data_line, *, iostat=ios) idum, ic_config, idum3, rdum, rdum2, rdum3, temp_corr(1:outsqs_n_clusters)
      if (ios /= 0) then
        write (*, '(a,i0)') " Error parsing config ", ic
        write (*, '(a,i0)') " iostat = ", ios
        stop 1
      end if
      ! ic_config is the actual configuration ID (1-71), use it as the index into corr_all
      ! Map back to full cluster array
      do jc_gen = 1, n_clusters
        if (cluster_active(jc_gen) > 0) then
          jc_outsqs = cluster_active(jc_gen)
          corr_all(ic_config, jc_gen) = temp_corr(jc_outsqs)
        end if
      end do
    end do
    deallocate (temp_corr)
    close (26)
  end subroutine read_energies_and_correlations

  ! ================================================================
  ! Compute thermal target correlations at temperature T using Boltzmann averaging
  ! ================================================================
  subroutine set_thermal_targets(T)
    implicit none
    real(real64), intent(in) :: T
    integer :: ic, jc, kc, ic_min
    real(real64) :: beta, z, weight, emin, erel_ic

    ! Compute thermal averages: ⟨ρα(T)⟩ = Σ_i ωi exp(-Ei/kT) ρα(i) / Z(T)

    if (abs(T) < 1.0e-8_real64) then
      ! T=0 limit: use lowest energy configuration
      emin = minval(ene(1:nic))
      ic_min = 1
      do ic = 1, nic
        if (abs(ene(ic) - emin) < 1.0e-10_real64) then
          ic_min = ic
          exit
        end if
      end do
      write (*, '(a,i0,a,f15.8)') " > Ground state: configuration ", ic_min, " with energy ", emin
      do jc = 1, n_clusters
        avg_target(jc) = corr_all(ic_min, jc)
      end do
    else
      ! Finite T: Boltzmann weighted average
      ! Subtract minimum energy to avoid overflow/underflow in exponentials
      emin = minval(ene(1:nic))
      beta = 1.0_real64 / (kb * T)
      z = 0.0_real64
      do ic = 1, nic
        erel_ic = ene(ic) - emin
        z = z + real(degen(ic), real64) * exp(-beta * erel_ic)
      end do

      do jc = 1, n_clusters
        avg_target(jc) = 0.0_real64
        do ic = 1, nic
          erel_ic = ene(ic) - emin
          avg_target(jc) = avg_target(jc) + real(degen(ic), real64) * exp(-beta * erel_ic) * corr_all(ic, jc)
        end do
        avg_target(jc) = avg_target(jc) / z
      end do
    end if

    do jc = 1, n_clusters
      clusters(jc)%target_corr = avg_target(jc)
    end do

    ! Precompute eps_min_cfg at this temperature
    if (.not. allocated(eps_min_cfg)) then
      allocate (eps_min_cfg(n_clusters))
    end if
    do kc = 1, n_clusters
      eps_min_cfg(kc) = minval(abs(corr_all(:, kc) - clusters(kc)%target_corr))
    end do
  end subroutine set_thermal_targets

  ! ================================================================
  ! Score configurations (copy from sqssod, modified for thermal targets)
  ! ================================================================
  subroutine score_configurations()
    implicit none

    if (allocated(scores)) deallocate(scores)
    allocate (scores(nic))
    if (allocated(L_val)) deallocate(L_val, abs_residual)
    allocate (L_val(nic), abs_residual(nic))

    call compute_vdw_scores()

    write (*, '(a)') " > Scoring complete."
    write (*, *)
  end subroutine score_configurations

  ! ================================================================
  ! Van de Walle scoring (copy from sqssod)
  ! ================================================================
  subroutine compute_vdw_scores()
    implicit none
    integer :: ic, jj, jc
    integer :: sorted_c(n_clusters)
    real(real64) :: dev_val, cur_diam, L_this, weight_norm

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

      L_this = 0.0_real64
      jj = 1
      do while (jj <= n_clusters)
        cur_diam = clusters(sorted_c(jj))%char_dist
        do while (jj <= n_clusters)
          if (abs(clusters(sorted_c(jj))%char_dist - cur_diam) >= 1.0e-8_real64) exit
          jc = sorted_c(jj)
          dev_val = corr_all(ic, jc) - clusters(jc)%target_corr
          if (abs(dev_val) > eps_min_cfg(jc) + eps_tol) goto 10
          jj = jj + 1
        end do
        L_this = cur_diam
      end do
10    continue

      L_val(ic) = L_this

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
  ! Sort clusters by char_dist, order
  ! ================================================================
  subroutine sort_clusters_by_dist(idx, nn)
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
  end subroutine sort_clusters_by_dist

  ! ================================================================
  ! Report results for one temperature
  ! ================================================================
  subroutine report_results(T, temp_index)
    implicit none
    real(real64), intent(in) :: T
    integer, intent(in) :: temp_index
    integer :: ic, n_show
    integer, allocatable :: ridx(:)

    allocate (ridx(nic))
    do ic = 1, nic
      ridx(ic) = ic
    end do
    call qsort_idx(ridx, 1, nic, scores)

    n_show = min(nic, 10)
    write (*, '(a,f10.2,a)') " > Temperature: ", T, " K"
    write (*, '(a,i0,a)') " > Top ", n_show, " GQS candidates:"
    write (*, *)
    write (*, '(a)') "    Rank  Config  Degen      L(A)       AbsErr                Q"
    write (*, '(a)') "    ----  ------  -----  --------  ----------  ---------------"
    do ic = 1, n_show
      write (*, '(i7, i8, i7, f10.4, es12.4, es16.6)') ic, ridx(ic), degen(ridx(ic)), &
        L_val(ridx(ic)), abs_residual(ridx(ic)), scores(ridx(ic))
    end do
    write (*, *)

    deallocate (ridx)
  end subroutine report_results

  ! ================================================================
  ! Write thermal averages to OUTGQS
  ! ================================================================
  subroutine write_thermal_averages()
    implicit none
    integer :: iu, tt, jc, ic
    real(real64) :: beta, z, weight, T_inf_avg
    character(len=256) :: path_out

    path_out = trim(ensemble_dir) // "/OUTGQS"
    iu = 31
    open (unit=iu, file=trim(path_out))
    write (iu, '(a)') "# Thermal averages of cluster correlations (GQS targets)"
    write (iu, '(a,i0)') "# Number of configurations: ", nic
    write (iu, '(a,i0)') "# Number of clusters: ", n_clusters
    write (iu, '(a)') "#"

    ! Write thermal averages for each temperature
    do tt = 1, ntemp
      write (iu, '(a,f10.2,a)') "# T = ", temps(tt), " K"
      do jc = 1, n_clusters
        write (iu, '(a,i0,a,f12.8)') "Pi_", jc, " : ", thermal_avg(tt, jc)
      end do
      write (iu, '(a)') ""
    end do

    ! T→∞ limit: equiprobable average
    write (iu, '(a)') "# T -> infinity (equiprobable configurations)"
    do jc = 1, n_clusters
      T_inf_avg = 0.0_real64
      do ic = 1, nic
        T_inf_avg = T_inf_avg + real(degen(ic), real64) * corr_all(ic, jc)
      end do
      T_inf_avg = T_inf_avg / sum(real(degen(1:nic), real64))
      write (iu, '(a,i0,a,f12.8)') "Pi_", jc, " : ", T_inf_avg
    end do

    close (iu)
    write (*, '(a,a)') " > Thermal averages written to ", trim(path_out)
    write (*, *)
  end subroutine write_thermal_averages

  ! ================================================================
  ! Quicksort (copy from sqssod)
  ! ================================================================
  recursive subroutine qsort_idx(idx, lo, hi, key)
    implicit none
    integer, intent(inout) :: idx(:)
    integer, intent(in) :: lo, hi
    real(real64), intent(in) :: key(:)
    integer :: i, j, pi
    integer :: tmp

    if (lo >= hi) return
    if (hi - lo < 2) then
      ! For small arrays, use insertion sort
      if (key(idx(lo)) > key(idx(hi))) then
        tmp = idx(lo); idx(lo) = idx(hi); idx(hi) = tmp
      end if
      return
    end if

    ! Partition using Lomuto scheme (simpler, safer)
    pi = partition_lomuto(idx, lo, hi, key)
    call qsort_idx(idx, lo, pi - 1, key)
    call qsort_idx(idx, pi + 1, hi, key)
  end subroutine qsort_idx

  function partition_lomuto(idx, lo, hi, key) result(pi)
    implicit none
    integer, intent(inout) :: idx(:)
    integer, intent(in) :: lo, hi
    real(real64), intent(in) :: key(:)
    integer :: pi, i, tmp
    real(real64) :: pivot

    pivot = key(idx(hi))
    pi = lo
    do i = lo, hi - 1
      if (key(idx(i)) <= pivot) then
        tmp = idx(pi); idx(pi) = idx(i); idx(i) = tmp
        pi = pi + 1
      end if
    end do
    tmp = idx(pi); idx(pi) = idx(hi); idx(hi) = tmp
  end function partition_lomuto

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

end program gqssod
