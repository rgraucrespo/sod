!*******************************************************************************
!    mcsod — Monte Carlo sampling with SOD effective Hamiltonian (PME).
!
!    Runs one Metropolis MC chain at the temperature supplied as a positional
!    argument and writes output to nXX/MCT_TTTK/PMEx/.
!    For uniform random sampling (energy-free, all moves accepted) use randomsod.
!
!    Required input files: INMC, INSOD, EQMATRIX, SGO, the PME training data
!    (n00/ENERGIES, …) and a temperature: mcsod <temperature/K>.
!
!    Part of the SOD package (v0.83) — GNU GPL v3+.
!*******************************************************************************

program mcsod
  use iso_fortran_env, only: real64, int64, error_unit
  use pmemod
  use ensemble_io,     only: write_ensemble
  use config_sampling, only: mc_random_subset, seed_rng_clock, seed_rng_fixed, &
                             read_next_data_line, to_lower
  implicit none

  integer :: target_level, ios, iseed, isym_reduction, write_trace
  integer :: arg_count, mkdir_status, copy_status, iarg
  real(real64) :: restart_prob_in
  real(real64) :: requested_temperature
  integer :: n_equil, n_prod
  logical :: use_symmetry_reduction, recal_ok
  character(len=32) :: seed_label, arg_text, level_dir, pme_variant, temp_dir
  integer :: unit_in
  character(len=256) :: line, start_config_line, model_filename_arg
  character(len=256) :: output_dir, energies_dir

  write (*, '(A)') "SOD (Site-Occupancy Disorder) version 0.83 - mcsod"

  ! --- Read INMC ---
  open(newunit=unit_in, file='INMC', status='old', action='read', iostat=ios)
  if (ios /= 0) then
    write(error_unit, '(A)') ' Error: INMC file not found.'
    stop 1
  end if

  ! symmetry_reduction (0 = off, 1 = on)
  call read_next_data_line(unit_in, line, ios)
  if (ios /= 0) then
    write(error_unit, '(A)') ' Error: could not read symmetry_reduction from INMC'; stop 1
  end if
  read(line, *, iostat=ios) isym_reduction
  if (ios /= 0 .or. (isym_reduction /= 0 .and. isym_reduction /= 1)) then
    write(error_unit, '(A,A)') ' Error: invalid symmetry_reduction in INMC (0=off, 1=on): ', trim(line); stop 1
  end if
  use_symmetry_reduction = (isym_reduction == 1)

  ! n_prod (production steps)
  call read_next_data_line(unit_in, line, ios)
  if (ios /= 0) then
    write(error_unit, '(A)') ' Error: could not read n_prod from INMC'; stop 1
  end if
  read(line, *, iostat=ios) n_prod
  if (ios /= 0 .or. n_prod <= 0) then
    write(error_unit, '(A,A)') ' Error: invalid n_prod in INMC (must be > 0): ', trim(line); stop 1
  end if

  ! starting configuration ('random' or space-separated site indices)
  call read_next_data_line(unit_in, line, ios)
  if (ios /= 0) then
    write(error_unit, '(A)') ' Error: could not read starting configuration from INMC'; stop 1
  end if
  start_config_line = trim(adjustl(line))

  ! write_trace (0 = off, 1 = write output_dir/MCTRACE)
  call read_next_data_line(unit_in, line, ios)
  if (ios /= 0) then
    write(error_unit, '(A)') ' Error: could not read write_trace from INMC'; stop 1
  end if
  read(line, *, iostat=ios) write_trace
  if (ios /= 0 .or. (write_trace /= 0 .and. write_trace /= 1)) then
    write(error_unit, '(A,A)') ' Error: invalid write_trace in INMC (0 or 1): ', trim(line); stop 1
  end if

  ! n_equil (equilibration steps)
  call read_next_data_line(unit_in, line, ios)
  if (ios /= 0) then
    write(error_unit, '(A)') ' Error: could not read n_equil from INMC'; stop 1
  end if
  read(line, *, iostat=ios) n_equil
  if (ios /= 0 .or. n_equil < 0) then
    write(error_unit, '(A,A)') ' Error: invalid n_equil in INMC (must be >= 0): ', trim(line); stop 1
  end if

  ! restart_prob
  call read_next_data_line(unit_in, line, ios)
  if (ios /= 0) then
    write(error_unit, '(A)') ' Error: could not read restart_prob from INMC'; stop 1
  end if
  read(line, *, iostat=ios) restart_prob_in
  if (ios /= 0 .or. restart_prob_in < 0.0_real64 .or. restart_prob_in >= 1.0_real64) then
    write(error_unit, '(A,A)') ' Error: invalid restart_prob in INMC (0 <= p < 1): ', trim(line); stop 1
  end if

  ! random_seed (-1 = system clock, >0 = fixed)
  call read_next_data_line(unit_in, line, ios)
  if (ios /= 0) then
    write(error_unit, '(A)') ' Error: could not read random_seed from INMC'; stop 1
  end if
  read(line, *, iostat=ios) iseed
  if (ios /= 0 .or. iseed == 0 .or. iseed < -1) then
    write(error_unit, '(A,A)') &
      ' Error: invalid random_seed in INMC (-1 = system clock, >0 = fixed integer): ', trim(line)
    stop 1
  end if

  close(unit_in)

  ! Parse command-line arguments: collect positional args; handle -model flag
  requested_temperature = 0.0_real64
  arg_text = ''
  arg_count = command_argument_count()
  iarg = 1
  do while (iarg <= arg_count)
    call get_command_argument(iarg, model_filename_arg)
    if (trim(model_filename_arg) == '-model') then
      if (iarg + 1 > arg_count) then
        write(error_unit,'(A)') ' Error: -model requires a filename argument.'
        stop 1
      end if
      iarg = iarg + 1
      call get_command_argument(iarg, model_filename_arg)
      call pme_set_model_filename(trim(model_filename_arg))
      write(*,'(A,A)') ' > Using model file: ', trim(model_filename_arg)
    else
      ! Treat as positional; only one positional arg expected (temperature)
      if (arg_text == '') then
        arg_text = trim(model_filename_arg)
      else
        write(error_unit,'(A,A)') ' Error: unexpected argument: ', trim(model_filename_arg)
        stop 1
      end if
    end if
    iarg = iarg + 1
  end do

  if (arg_text == '') then
    write(error_unit, '(A)') ' Error: mcsod requires a temperature argument, e.g. mcsod 600.0'
    stop 1
  end if
  read(arg_text, *, iostat=ios) requested_temperature
  if (ios /= 0 .or. requested_temperature <= 0.0_real64) then
    write(error_unit, '(A,A)') ' Error: invalid mcsod temperature argument: ', trim(arg_text)
    stop 1
  end if

  ! --- Print parameters ---
  write(*, '(A)') ' --- MC parameters ----------------------------------------------------------'
  write(*, '(A)')          '  Sampler           : Metropolis'
  write(*, '(A,I0)')       '  Equil. steps      : ', n_equil
  write(*, '(A,F6.4)')     '  Restart prob.     : ', restart_prob_in
  if (iseed < 0) then
    seed_label = 'system clock'
  else
    write(seed_label, '(I0)') iseed
  end if
  write(*, '(A,A)')        '  Random seed       : ', trim(seed_label)
  write(*, '(A,F10.2,A)')  '  Temperature       : ', requested_temperature, ' K'
  if (use_symmetry_reduction) then
    write(*, '(A)')          '  Symmetry          : on'
  else
    write(*, '(A)')          '  Symmetry          : off'
  end if
  write(*, '(A,I0)')         '  Production steps  : ', n_prod
  write(*, '(A,A)')          '  Start config      : ', trim(start_config_line)
  if (write_trace == 1) then
    write(*, '(A)')          '  Trace             : on (MCTRACE)'
  else
    write(*, '(A)')          '  Trace             : off'
  end if
  write(*, '(A)') ' ---------------------------------------------------------------------------'
  write(*, *)

  ! --- Seed the RNG (once before the temperature loop) ---
  if (iseed > 0) then
    call seed_rng_fixed(iseed)
  else
    call seed_rng_clock()
  end if

  ! --- Initialize PME model ---
  call pme_get_target_level_from_insod(target_level)
  call pme_initialize_model(target_level)
  recal_ok = .true.
  if (pme_energies_available()) call pme_preload_recalibration(target_level, recal_ok)
  call pme_print_model_summary()

  call format_level_directory(target_level, level_dir)
  call format_pme_variant_directory(pme_get_choice(), pme_variant)

  ! Output layout: nXX/MCT_TK/PMEx/{ENSEMBLE,ENERGIES} (the walk is Hamiltonian-driven).
  call format_metropolis_directory(requested_temperature, temp_dir)
  output_dir   = trim(level_dir)//'/'//trim(temp_dir)//'/'//trim(pme_variant)
  energies_dir = output_dir

  call execute_command_line('mkdir -p '//trim(output_dir), exitstat=mkdir_status)
  if (mkdir_status /= 0) then
    write(error_unit,'(A,A)') ' Error: could not create output directory ', trim(output_dir)
    stop 1
  end if
  call execute_command_line('cp INMC '//trim(output_dir)//'/INMC', exitstat=copy_status)
  if (copy_status /= 0) then
    write(*,'(A,A)') ' Warning: could not copy INMC to ', trim(output_dir)
  end if

  ! --- Run MC ---
  write(*, '(A,F10.2,A)') '  Running Metropolis MC at ', requested_temperature, ' K'
  write(*, *)
  call run_mc(target_level, requested_temperature, n_equil, n_prod, restart_prob_in, &
              use_symmetry_reduction, write_trace, &
              start_config_line, output_dir, energies_dir)

  ! --- Cleanup ---
  call pme_finalize_model()

  write(*, '(A)') ' ============================================================================'
  write(*, '(A)') '  mcsod completed.'
  write(*, '(A)') ' ============================================================================'
  write(*, *)

contains

  subroutine format_pme_variant_directory(choice, dirname)
    integer, intent(in) :: choice
    character(len=*), intent(out) :: dirname

    select case (choice)
    case (0)
      dirname = 'PME0'
    case (1)
      dirname = 'PME1'
    case default
      dirname = 'PMEh'
    end select
  end subroutine format_pme_variant_directory

  subroutine format_metropolis_directory(temperature, dirname)
    real(real64), intent(in) :: temperature
    character(len=*), intent(out) :: dirname
    integer :: temp_integer

    temp_integer = int(temperature)
    write(dirname, '(A,I0,A)') 'MCT_', temp_integer, 'K'
  end subroutine format_metropolis_directory

  ! =========================================================================
  !  Main MC procedure
  ! =========================================================================

  subroutine run_mc(target_level, temperature, n_equil, n_prod, restart_prob, &
                    use_symmetry_reduction, write_trace, &
                    start_config_line, output_dir, energies_dir)
    !  Runs n_equil equilibration steps (chain advances, no storage) followed by
    !  n_prod production steps (configurations stored, contribute to averages).
    !  Writes output_dir/{ENSEMBLE,OUTMC,ENERGIES} (energies_dir == output_dir).
    !  Optionally writes output_dir/MCTRACE (equil + prod steps, phase-labelled).
    integer,          intent(in) :: target_level, n_equil, n_prod, write_trace
    real(real64),     intent(in) :: temperature, restart_prob
    logical,          intent(in) :: use_symmetry_reduction
    character(len=*), intent(in) :: start_config_line, output_dir, energies_dir

    real(real64), parameter :: kB = 8.617333262e-5_real64  ! eV/K

    integer :: k, lev, npos_val, atini_val, max_lo, max_hi
    logical :: has_high
    real(real64) :: v0_lo, v0_hi, mu_lo(0:4), mu_hi(0:4), alpha_hybrid_val, eta_hybrid_val
    integer, allocatable :: subset(:), trial_subset(:), best_subset(:)
    real(real64) :: energy, trial_energy, best_energy
    real(real64) :: energy_low, energy_high, trial_low, trial_high
    real(real64) :: low_terms(4), high_terms(4)             ! running terms of current config
    real(real64) :: trial_low_terms(4), trial_high_terms(4) ! terms of the trial config
    real(real64) :: dlow_terms(4), dhigh_terms(4)           ! incremental change for a swap
    real(real64) :: delta_e, beta, rand_num
    integer :: n_accepted_equil, n_accepted_prod
    integer :: i_step, i, j
    integer :: removed_site, added_site, hole_pick, hole_count
    integer, allocatable :: cur_holes(:)   ! complement (hole) list of current config
    integer :: steps_since_resync
    ! Periodic full recompute to bound floating-point drift in the running
    ! term vectors (equivalence test B showed zero drift, but resync keeps it
    ! robust for larger models / longer chains). One full eval per n_resync steps.
    integer, parameter :: n_resync = 5000
    logical :: accepted, is_restart
    character(len=256) :: outmc_path, trace_path
    integer :: unit_out, unit_trace

    character(len=32) :: sampler_label
    character(len=64) :: symmetry_label
    character(len=128) :: ensemble_rows_label, omega_label, stats_weighting_label

    integer, allocatable :: current_canonical(:), trial_canonical(:)
    integer :: current_idx

    integer :: n_unique, max_unique
    integer, allocatable :: unique_subsets(:,:)
    integer, allocatable :: degeneracy(:)
    integer, allocatable :: visit_count(:)
    real(real64), allocatable :: unique_low(:,:), unique_high(:,:)
    real(real64), allocatable :: unique_e(:)
    logical :: found_match

    real(real64) :: sum_e, sum_ge, mean_e, std_e
    integer :: best_unique_idx, total_deg

    real(real64), allocatable :: energy_trace(:)
    integer, parameter :: n_block_sizes = 3
    integer, parameter :: block_sizes(n_block_sizes) = [8, 16, 32]
    real(real64) :: sem_vals(n_block_sizes)
    integer :: i_b, n_bs, n_block_steps, i_block
    real(real64) :: block_mean_e, block_mean_bar, block_ss
    logical :: do_block_sem
    integer(int64) :: total_omega

    lev       = target_level
    npos_val  = pme_get_npos()
    atini_val = pme_get_target_atini()
    max_lo    = pme_get_max_low_order()
    max_hi    = pme_get_max_high_order()
    has_high  = pme_has_high_side()
    call pme_get_v0(v0_lo, v0_hi)

    outmc_path = trim(output_dir)//'/OUTMC'

    if (lev == 0) then
      write(*, '(A)') ' > Level 0: single configuration — delegating to regular outputs.'
      call pme_write_level_outputs(target_level)
      return
    end if

    if (lev >= npos_val) then
      write(error_unit,'(A,I0,A,I0,A)') &
        ' Error: substitution level (', lev, ') >= total sites (', npos_val, ').'
      write(error_unit,'(A)') ' MC requires at least one unoccupied site for swap moves.'
      stop 1
    end if

    block
      integer :: mkdir_status
      call execute_command_line('mkdir -p '//trim(output_dir), exitstat=mkdir_status)
      if (mkdir_status /= 0) then
        write(error_unit,'(A,A)') ' Error: could not create directory ', trim(output_dir)
        stop 1
      end if
    end block

    beta = 1.0_real64 / (kB * temperature)

    allocate(subset(lev), trial_subset(lev), best_subset(lev))
    hole_count = npos_val - lev          ! constant: substitution count is fixed
    allocate(cur_holes(max(1, hole_count)))

    ! --- Set mode labels ---
    sampler_label = 'Metropolis'
    stats_weighting_label = 'Omega/sum(Omega) (Metropolis sample already energy-biased)'

    if (use_symmetry_reduction) then
      symmetry_label = 'reduced (symmetry representatives)'
      ensemble_rows_label = 'symmetry representatives'
      omega_label = 'visits to representative, including rejected-step residence'
    else
      symmetry_label = 'full (no symmetry reduction)'
      ensemble_rows_label = 'explicit trajectory states, no global deduplication'
      omega_label = 'residence count for consecutive stays on same explicit configuration'
    end if

    write(*, '(A)') ' --- Sampling ---------------------------------------------------------------'
    write(*, '(A,I0,A,I0,A)') '  Level             : ', lev, ' substitutions in ', npos_val, ' positions'

    max_unique = n_prod
    allocate(current_canonical(lev), trial_canonical(lev))
    allocate(unique_subsets(lev, max_unique), degeneracy(max_unique))
    allocate(unique_low(4, max_unique), unique_high(4, max_unique))
    allocate(unique_e(max_unique), visit_count(max_unique))
    allocate(energy_trace(n_prod))
    visit_count = 0
    n_unique    = 0

    ! --- Parse starting configuration ---
    block
      integer :: buf(lev), ios_cfg
      character(len=256) :: cfg_lower

      cfg_lower = start_config_line
      call to_lower(cfg_lower)
      if (trim(cfg_lower) == 'random') then
        call mc_random_subset(npos_val, lev, subset)
        write(*, '(A)') '  Initial config    : random'
      else
        read(start_config_line, *, iostat=ios_cfg) buf(1:lev)
        if (ios_cfg == 0) then
          subset(1:lev) = buf(1:lev)
          call mc_insertion_sort(subset, lev)
          write(*,'(A,*(1X,I0))') '  Initial config    : from INMC:', subset(1:lev)
        else
          write(*,'(A)') '  Initial config    : INMC start config unreadable — using random'
          call mc_random_subset(npos_val, lev, subset)
        end if
      end if
    end block

    ! Hole list for the starting configuration (maintained incrementally below).
    call rebuild_hole_set(subset, lev, npos_val, cur_holes)

    call pme_evaluate_configuration(subset(1:lev), lev, &
      energy, energy_low, energy_high, low_terms, high_terms)
    call apply_epsilon_energy(lev, low_terms, high_terms, energy)

    n_accepted_equil = 0
    n_accepted_prod  = 0
    steps_since_resync = 0

    if (write_trace == 1) then
      trace_path = trim(output_dir)//'/MCTRACE'
      open(newunit=unit_trace, file=trim(trace_path), status='replace', action='write')
      write(unit_trace,'(A)') '# phase       step        energy(eV)              stored_idx (0 during equil)'
    end if

    write(*, '(A,I0,A,I0,A)') '  Running           : ', n_equil, ' equil. + ', n_prod, ' production steps...'
    write(*, *)

    ! =================================================================
    ! Phase 1a: Equilibration
    ! =================================================================
    do i_step = 1, n_equil
      trial_subset = subset
      call random_number(rand_num)
      if (rand_num < restart_prob) then
        do
          call mc_random_subset(npos_val, lev, trial_subset)
          if (.not. all(trial_subset(1:lev) == subset(1:lev))) exit
        end do
        is_restart = .true.
        ! Restart replaces many sites at once — recompute the trial in full.
        call pme_evaluate_configuration(trial_subset(1:lev), lev, &
          trial_energy, trial_low, trial_high, trial_low_terms, trial_high_terms)
      else
        is_restart = .false.
        call mc_swap_move(trial_subset, lev, cur_holes, hole_count, &
          removed_site, added_site, hole_pick)
        ! Single swap — update the term vectors incrementally from the current
        ! config, reusing the maintained hole list (no O(npos) rebuild).
        call pme_evaluate_swap_delta(subset(1:lev), lev, removed_site, added_site, &
          dlow_terms, dhigh_terms, cur_holes, hole_count)
        trial_low_terms  = low_terms  + dlow_terms
        trial_high_terms = high_terms + dhigh_terms
      end if
      call apply_epsilon_energy(lev, trial_low_terms, trial_high_terms, trial_energy)

      delta_e = trial_energy - energy
      if (delta_e <= 0.0_real64) then
        accepted = .true.
      else
        call random_number(rand_num)
        accepted = rand_num < exp(-beta * delta_e)
      end if

      if (accepted) then
        n_accepted_equil = n_accepted_equil + 1
        subset     = trial_subset
        energy     = trial_energy
        low_terms  = trial_low_terms
        high_terms = trial_high_terms
        if (is_restart) then
          call rebuild_hole_set(subset, lev, npos_val, cur_holes)
        else
          ! O(1) hole-list update: added_site leaves the holes, removed_site joins.
          cur_holes(hole_pick) = removed_site
        end if
      end if

      ! Periodic full recompute to bound floating-point drift.
      steps_since_resync = steps_since_resync + 1
      if (steps_since_resync >= n_resync) then
        call pme_evaluate_configuration(subset(1:lev), lev, &
          energy, energy_low, energy_high, low_terms, high_terms)
        call apply_epsilon_energy(lev, low_terms, high_terms, energy)
        steps_since_resync = 0
      end if

      if (write_trace == 1) then
        write(unit_trace,'(A5,2X,I10,2X,F22.10,2X,I6)') 'equil', i_step, energy, 0
      end if
    end do

    ! =================================================================
    ! Transition: full re-evaluation gives a clean resync of the running
    ! term vectors before production storage begins.
    ! =================================================================
    if (n_equil > 0) then
      call pme_evaluate_configuration(subset(1:lev), lev, &
        energy, energy_low, energy_high, low_terms, high_terms)
      call apply_epsilon_energy(lev, low_terms, high_terms, energy)
      steps_since_resync = 0
    end if

    best_energy      = energy
    best_subset      = subset
    n_unique         = 1
    if (use_symmetry_reduction) then
      call canonicalize_config(subset, lev, current_canonical)
      unique_subsets(:, 1) = current_canonical
    else
      unique_subsets(:, 1) = subset
    end if
    unique_low(:, 1)  = low_terms
    unique_high(:, 1) = high_terms
    visit_count(1)    = 1
    current_idx       = 1

    energy_trace(1) = energy

    if (write_trace == 1) then
      write(unit_trace,'(A4,3X,I10,2X,F22.10,2X,I6)') 'prod', 1, energy, current_idx
    end if

    ! =================================================================
    ! Phase 1b: Production
    ! =================================================================
    do i_step = 2, n_prod
      trial_subset = subset
      call random_number(rand_num)
      if (rand_num < restart_prob) then
        do
          call mc_random_subset(npos_val, lev, trial_subset)
          if (.not. all(trial_subset(1:lev) == subset(1:lev))) exit
        end do
        is_restart = .true.
        ! Restart replaces many sites at once — recompute the trial in full.
        call pme_evaluate_configuration(trial_subset(1:lev), lev, &
          trial_energy, trial_low, trial_high, trial_low_terms, trial_high_terms)
      else
        is_restart = .false.
        call mc_swap_move(trial_subset, lev, cur_holes, hole_count, &
          removed_site, added_site, hole_pick)
        ! Single swap — incremental update from the current config, reusing the
        ! maintained hole list (no O(npos) rebuild).
        call pme_evaluate_swap_delta(subset(1:lev), lev, removed_site, added_site, &
          dlow_terms, dhigh_terms, cur_holes, hole_count)
        trial_low_terms  = low_terms  + dlow_terms
        trial_high_terms = high_terms + dhigh_terms
      end if
      call apply_epsilon_energy(lev, trial_low_terms, trial_high_terms, trial_energy)

      delta_e = trial_energy - energy
      if (delta_e <= 0.0_real64) then
        accepted = .true.
      else
        call random_number(rand_num)
        accepted = rand_num < exp(-beta * delta_e)
      end if

      if (accepted) then
        n_accepted_prod = n_accepted_prod + 1
        subset     = trial_subset
        energy     = trial_energy
        low_terms  = trial_low_terms
        high_terms = trial_high_terms

        ! Keep the hole list in step with the accepted config.
        if (is_restart) then
          call rebuild_hole_set(subset, lev, npos_val, cur_holes)
        else
          cur_holes(hole_pick) = removed_site
        end if

        if (energy < best_energy) then
          best_energy = energy
          best_subset = subset
        end if

        if (use_symmetry_reduction) then
          call canonicalize_config(subset, lev, trial_canonical)
          if (.not. all(trial_canonical(1:lev) == current_canonical(1:lev))) then
            found_match = .false.
            do j = 1, n_unique
              if (all(trial_canonical(1:lev) == unique_subsets(1:lev, j))) then
                found_match = .true.
                current_idx = j
                exit
              end if
            end do
            if (.not. found_match) then
              if (n_unique < max_unique) then
                n_unique = n_unique + 1
                unique_subsets(:, n_unique) = trial_canonical
                unique_low(:, n_unique)     = low_terms
                unique_high(:, n_unique)    = high_terms
                current_idx                 = n_unique
              else
                write(error_unit,'(A)') ' Error: mcsod storage exhausted in symmetry-reduction path.'
                stop 1
              end if
            end if
            current_canonical = trial_canonical
          end if
        else
          if (n_unique >= max_unique) then
            write(error_unit,'(A)') ' Error: mcsod storage exhausted while appending full-space sample.'
            stop 1
          end if
          n_unique = n_unique + 1
          unique_subsets(:, n_unique) = subset
          unique_low(:, n_unique)     = low_terms
          unique_high(:, n_unique)    = high_terms
          visit_count(n_unique)       = 0
          current_idx                 = n_unique
        end if
      end if

      visit_count(current_idx) = visit_count(current_idx) + 1

      energy_trace(i_step) = energy

      if (write_trace == 1) then
        write(unit_trace,'(A4,3X,I10,2X,F22.10,2X,I6)') 'prod', i_step, energy, current_idx
      end if

      ! Periodic full recompute to bound floating-point drift in the running term vectors.
      steps_since_resync = steps_since_resync + 1
      if (steps_since_resync >= n_resync) then
        call pme_evaluate_configuration(subset(1:lev), lev, &
          energy, energy_low, energy_high, low_terms, high_terms)
        call apply_epsilon_energy(lev, low_terms, high_terms, energy)
        steps_since_resync = 0
      end if
    end do

    if (write_trace == 1) close(unit_trace)

    deallocate(current_canonical, trial_canonical)

    ! === Block-average SEM estimation from production energy trajectory ===
    ! SEM = std(block_means) / sqrt(n_blocks); blocks must be >> autocorrelation time.
    ! If SEM converges as n_blocks decreases (block size grows), the estimate is reliable.
    do_block_sem = n_prod >= block_sizes(n_block_sizes)
    if (do_block_sem) then
      do i_b = 1, n_block_sizes
        n_bs = block_sizes(i_b)
        n_block_steps = n_prod / n_bs
        block_mean_bar = 0.0_real64
        do i_block = 1, n_bs
          block_mean_e = sum(energy_trace((i_block-1)*n_block_steps+1 : i_block*n_block_steps)) &
                         / real(n_block_steps, real64)
          block_mean_bar = block_mean_bar + block_mean_e
        end do
        block_mean_bar = block_mean_bar / real(n_bs, real64)
        block_ss = 0.0_real64
        do i_block = 1, n_bs
          block_mean_e = sum(energy_trace((i_block-1)*n_block_steps+1 : i_block*n_block_steps)) &
                         / real(n_block_steps, real64)
          block_ss = block_ss + (block_mean_e - block_mean_bar)**2
        end do
        sem_vals(i_b) = sqrt(block_ss / real(n_bs - 1, real64)) / sqrt(real(n_bs, real64))
      end do
    end if
    deallocate(energy_trace)

    ! === Visit counts as weights ===
    total_omega = sum(int(visit_count(1:n_unique), int64))
    do i = 1, n_unique
      degeneracy(i) = visit_count(i)
    end do

    write(*, '(A)') ' --- Sampling statistics ----------------------------------------------------'
    if (n_equil > 0) then
      write(*, '(A,I0)')       '  Equil. steps      : ', n_equil
      write(*, '(A,F5.1,A)')   '  Equil. acceptance : ', &
        100.0_real64 * real(n_accepted_equil, real64) / real(max(1, n_equil), real64), ' %'
    end if
    write(*, '(A,I0)')         '  Production steps  : ', n_prod
    write(*, '(A,F10.2)')      '  Steps / site      : ', real(n_prod, real64) / real(npos_val, real64)
    write(*, '(A,I0)')         '  Stored configs    : ', n_unique
    write(*, '(A,I0,A,F5.1,A)') '  Accepted (prod.)  : ', n_accepted_prod, &
      '  (', 100.0_real64 * real(n_accepted_prod, real64) / real(max(1, n_prod), real64), ' %)'
    write(*, '(A)') ' ---------------------------------------------------------------------------'
    write(*, *)

    ! === Compute final calibrated energies and report ===
    do i = 1, n_unique
      call apply_epsilon_energy(lev, unique_low(:, i), unique_high(:, i), unique_e(i))
    end do

    ! === Canonical ensemble average (degeneracy-weighted mean) ===
    total_deg = sum(degeneracy(1:n_unique))
    sum_e  = 0.0_real64
    sum_ge = 0.0_real64
    do i = 1, n_unique
      sum_e  = sum_e  + real(degeneracy(i), real64) * unique_e(i)
      sum_ge = sum_ge + real(degeneracy(i), real64) * unique_e(i)**2
    end do
    mean_e = sum_e / real(total_deg, real64)
    std_e  = sqrt(max(0.0_real64, sum_ge / real(total_deg, real64) - mean_e**2))

    best_unique_idx = 1
    do i = 2, n_unique
      if (unique_e(i) < unique_e(best_unique_idx)) best_unique_idx = i
    end do
    best_energy = unique_e(best_unique_idx)
    best_subset = unique_subsets(:, best_unique_idx)

    write(*, '(A)') ' --- Results ----------------------------------------------------------------'
    write(*, '(A,F22.10,A)') '  E_min             = ', best_energy, ' eV'
    write(*, '(A,F22.10,A)') '  E_ave (sample)    = ', mean_e,      ' eV'
    write(*, '(A,F22.10,A)') '  E_std             = ', std_e,       ' eV'
    if (do_block_sem) then
      write(*, '(A)') ' --- Block-average SEM ------------------------------------------------------'
      do i_b = 1, n_block_sizes
        n_bs         = block_sizes(i_b)
        n_block_steps = n_prod / n_bs
        write(*, '(A,I2,A,F22.10,A,I0,A)') &
          '  SEM (', n_bs, ' blocks)    = ', sem_vals(i_b), ' eV  (', n_block_steps, ' steps/block)'
      end do
      write(*, '(A)') '  Note: if SEM(8) ≈ SEM(16), blocks are independent; trust that estimate.'
      write(*, '(A)') '        If SEM grows as n_blocks decreases, increase n_prod or n_equil.'
    end if
    write(*, '(A)') ' ---------------------------------------------------------------------------'
    write(*, *)

    ! === Write ENSEMBLE ===
    block
      character(len=256) :: ensemble_path, energies_path
      integer :: unit_ensemble, unit_energies

      ensemble_path = trim(output_dir)//'/ENSEMBLE'
      energies_path = trim(energies_dir)//'/ENERGIES'

      open(newunit=unit_ensemble, file=trim(ensemble_path), status='replace', action='write')
      block
        integer, allocatable :: mc_ic(:,:), mc_dg(:)
        integer :: ic_i
        character(len=10) :: orig_sym_mc, sym_new_mc, sym_rem_mc
        character(len=10) :: syms_mc(1,2)
        character(len=10) :: orig_arr(1)
        allocate(mc_ic(n_unique, max(1, lev)), mc_dg(n_unique))
        call pme_get_newsymbol(orig_sym_mc, sym_new_mc, sym_rem_mc)
        syms_mc(1,1) = sym_new_mc; syms_mc(1,2) = sym_rem_mc
        orig_arr(1)  = orig_sym_mc
        mc_dg = degeneracy(1:n_unique)
        mc_ic = 0
        do ic_i = 1, n_unique
          if (lev > 0) mc_ic(ic_i, 1:lev) = unique_subsets(1:lev, ic_i)
        end do
        call write_ensemble(unit_ensemble, 1, [1], [lev], [npos_val], [atini_val], &
                            n_unique, mc_ic, mc_dg, 'metropolis', orig_arr, syms_mc, &
                            tsampling=temperature)
      end block
      close(unit_ensemble)

      ! === Write ENERGIES (energies_dir == output_dir, already created) ===
      open(newunit=unit_energies, file=trim(energies_path), status='replace', action='write')
      do i = 1, n_unique
        write(unit_energies,'(I0,2X,F22.10)') i, unique_e(i)
      end do
      close(unit_energies)

    end block

    ! === Write OUTMC (summary only; config table is in ENSEMBLE + ENERGIES) ===
    call pme_get_epsilon(mu_lo, mu_hi)
    call pme_get_alpha_eta(alpha_hybrid_val, eta_hybrid_val)

    open(newunit=unit_out, file=trim(outmc_path), status='replace', action='write')
    write(unit_out,'(A)') ' ================================================================='
    write(unit_out,'(A)') '  SOD PME — Monte Carlo sampling summary'
    write(unit_out,'(A)') ' ================================================================='
    write(unit_out,'(A,A)')       '  Sampler           : ', trim(sampler_label)
    write(unit_out,'(A,A)')       '  Symmetry mode     : ', trim(symmetry_label)
    write(unit_out,'(A,A)')       '  ENSEMBLE rows     : ', trim(ensemble_rows_label)
    write(unit_out,'(A,A)')       '  ENSEMBLE omega    : ', trim(omega_label)
    write(unit_out,'(A,A)')       '  Stats weighting   : ', trim(stats_weighting_label)
    write(unit_out,'(A,I0)')      '  Level (nsubs)     : ', lev
    write(unit_out,'(A,I0)')      '  Total sites       : ', npos_val
    write(unit_out,'(A,F14.4,A)') '  Temperature       : ', temperature, ' K'
    write(unit_out,'(A,I0)')      '  Equil. steps      : ', n_equil
    write(unit_out,'(A,I0)')      '  Production steps  : ', n_prod
    write(unit_out,'(A,I0)')      '  Stored configs    : ', n_unique
    write(unit_out,'(A,I0,A,F5.1,A)') '  Accepted (equil.) : ', n_accepted_equil, &
      '  (', 100.0_real64 * real(n_accepted_equil, real64) / real(max(1, n_equil), real64), ' %)'
    write(unit_out,'(A,I0,A,F5.1,A)') '  Accepted (prod.)  : ', n_accepted_prod, &
      '  (', 100.0_real64 * real(n_accepted_prod, real64) / real(max(1, n_prod), real64), ' %)'
    write(unit_out,'(A)') ' -----------------------------------------------------------------'
    write(unit_out,'(A)',advance='no') '  eps_low           :'
    do k = 0, max(1, max_lo)
      write(unit_out,'(1X,F10.6)',advance='no') mu_lo(k)
    end do
    write(unit_out,*)
    if (has_high) then
      write(unit_out,'(A)',advance='no') '  eps_high          :'
      do k = 0, max(1, max_hi)
        write(unit_out,'(1X,F10.6)',advance='no') mu_hi(k)
      end do
      write(unit_out,*)
      write(unit_out,'(A,F10.6)')       '  alpha (hybrid)    : ', alpha_hybrid_val
      write(unit_out,'(A,F10.6)')       '  eta   (hybrid)    : ', eta_hybrid_val
    end if
    write(unit_out,'(A)') ' -----------------------------------------------------------------'
    write(unit_out,'(A,F22.10,A)') '  E_min             =', best_energy, ' eV'
    write(unit_out,'(A,F22.10,A)') '  E_ave (sample)    =', mean_e,      ' eV'
    write(unit_out,'(A,F22.10,A)') '  E_std             =', std_e,       ' eV'
    if (do_block_sem) then
      write(unit_out,'(A)') ' -----------------------------------------------------------------'
      write(unit_out,'(A)') '  Block-average SEM (standard error of the mean):'
      do i_b = 1, n_block_sizes
        n_bs          = block_sizes(i_b)
        n_block_steps = n_prod / n_bs
        write(unit_out,'(A,I2,A,F22.10,A,I0,A)') &
          '  SEM (', n_bs, ' blocks)    =', sem_vals(i_b), ' eV  (', n_block_steps, ' steps/block)'
      end do
      write(unit_out,'(A)') &
        '  Note: if SEM(8) ~ SEM(16), blocks are independent; trust that estimate.'
      write(unit_out,'(A)') &
        '        If SEM grows as n_blocks decreases, increase n_prod or n_equil.'
    end if
    write(unit_out,'(A)') ' ================================================================='
    close(unit_out)

    ! === Output files summary ===
    block
      character(len=256) :: ensemble_path2, energies_path2

      ensemble_path2    = trim(output_dir)//'/ENSEMBLE'
      energies_path2  = trim(energies_dir)//'/ENERGIES'

      write(*, '(A)') ' --- Output files -----------------------------------------------------------'
      write(*, '(A,A)') '  ENSEMBLE          : ', trim(ensemble_path2)
      write(*, '(A,A)') '  ENERGIES          : ', trim(energies_path2)
      write(*, '(A,A)') '  OUTMC             : ', trim(outmc_path)
      if (write_trace == 1) then
        write(*, '(A,A)') '  MCTRACE           : ', trim(trace_path)
      end if
      write(*, '(A)') ' ---------------------------------------------------------------------------'
      write(*, *)
    end block

    deallocate(subset, trial_subset, best_subset)
    deallocate(cur_holes)
    deallocate(unique_subsets, degeneracy, unique_low, unique_high, unique_e)
    deallocate(visit_count)
  end subroutine run_mc

  ! =========================================================================
  !  MC helpers
  ! =========================================================================

  subroutine mc_swap_move(subset, level, holes, hole_count, removed_site, added_site, hole_pick)
    !  Propose a swap: remove a random occupied site and add a random hole drawn
    !  directly from the maintained complement (hole) list — both O(1), with no
    !  rejection loop and no membership scan. Reports the removed/added sites (for
    !  the incremental energy update) and hole_pick (the index of added_site in
    !  holes, so the caller can update the hole list in O(1) on accept).
    integer, intent(inout) :: subset(:)
    integer, intent(in) :: level, hole_count
    integer, intent(in) :: holes(:)
    integer, intent(out) :: removed_site, added_site, hole_pick
    integer :: remove_idx
    real(real64) :: r

    call random_number(r)
    remove_idx = min(int(r * real(level, real64)) + 1, level)
    removed_site = subset(remove_idx)

    call random_number(r)
    hole_pick = min(int(r * real(hole_count, real64)) + 1, hole_count)
    added_site = holes(hole_pick)

    subset(remove_idx) = added_site
    call mc_insertion_sort(subset, level)
  end subroutine mc_swap_move

  subroutine rebuild_hole_set(subset, level, npos_val, holes)
    !  Rebuild the complement (hole) list = the npos_val-level sites not in
    !  subset(1:level). O(npos). Used for the initial config and after a restart
    !  move (which replaces many sites at once); single swaps update holes in O(1).
    integer, intent(in)  :: subset(:), level, npos_val
    integer, intent(out) :: holes(:)
    logical :: occ(npos_val)
    integer :: i, idx

    occ = .false.
    do i = 1, level
      occ(subset(i)) = .true.
    end do
    idx = 0
    do i = 1, npos_val
      if (.not. occ(i)) then
        idx = idx + 1
        holes(idx) = i
      end if
    end do
  end subroutine rebuild_hole_set

end program mcsod
