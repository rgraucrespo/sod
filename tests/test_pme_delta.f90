!*******************************************************************************
!    test_pme_delta — correctness gate for the incremental swap evaluator.
!
!    Run inside a working directory that already holds a built PME model
!    (INSOD, SGO, EQMATRIX + PME training data — exactly the inputs mcsod
!    needs). Initialises the model the same way mcsod does, then checks the
!    incremental update pme_evaluate_swap_delta against the full
!    pme_evaluate_configuration in two ways:
!
!      A. Per-step equivalence: for random configurations and random single
!         swaps, over a sweep of substitution levels (1 .. npos-1, which
!         exercises the 1-/2-/3-/4-body loops and the high side), assert
!         old_terms + delta == full(new_config) per order and that the
!         apply_epsilon_energy energies agree.
!
!      B. Drift bound: accumulate the term vectors through a long chain of
!         swaps using only the delta update and compare against an independent
!         full recompute at every step; assert the running error stays bounded
!         (this justifies the resync interval used by mcsod).
!
!    Exit status 0 = all checks pass, 1 = failure.
!
!    Part of the SOD package (v0.81) — GNU GPL v3+.
!*******************************************************************************

program test_pme_delta
  use iso_fortran_env, only: real64, error_unit, output_unit
  use pmemod
  implicit none

  real(real64), parameter :: tol_term   = 1.0e-9_real64   ! per-order term agreement
  real(real64), parameter :: tol_energy = 1.0e-9_real64   ! calibrated-energy agreement
  real(real64), parameter :: tol_drift  = 1.0e-8_real64   ! accumulated drift over a chain

  integer :: target_level, npos_val, lev
  logical :: recal_ok, have_energy
  integer :: n_levels_tested, n_swaps_per_level
  integer :: n_fail
  real(real64) :: max_term_err, max_energy_err, max_drift_err
  character(len=32) :: arg1

  call seed_fixed(20240607)

  ! Optional 'bench' mode: report per-call timing of the full recompute vs the
  ! incremental delta across a level sweep (per-step cost; not a pass/fail gate).
  arg1 = ''
  if (command_argument_count() >= 1) call get_command_argument(1, arg1)

  call pme_get_target_level_from_insod(target_level)
  call pme_initialize_model(target_level)
  recal_ok = .true.
  have_energy = pme_energies_available()
  if (have_energy) call pme_preload_recalibration(target_level, recal_ok)

  npos_val = pme_get_npos()
  if (.not. have_energy) then
    write(error_unit,'(A)') ' test_pme_delta: model has no energies — cannot test delta evaluator.'
    stop 1
  end if
  if (npos_val < 2) then
    write(error_unit,'(A)') ' test_pme_delta: need at least 2 sites.'
    stop 1
  end if

  write(output_unit,'(A,I0,A,I0)') ' test_pme_delta: npos = ', npos_val, &
    ', model target level = ', target_level

  if (trim(arg1) == 'bench') then
    call bench_scaling(npos_val)
    call pme_finalize_model()
    stop 0
  end if

  n_fail            = 0
  max_term_err      = 0.0_real64
  max_energy_err    = 0.0_real64
  max_drift_err     = 0.0_real64
  n_levels_tested   = 0
  n_swaps_per_level = 200

  ! ---- A. per-step equivalence across a level sweep ----
  do lev = 1, npos_val - 1
    call check_level_equivalence(lev, n_swaps_per_level, n_fail, max_term_err, max_energy_err)
    n_levels_tested = n_levels_tested + 1
  end do

  ! ---- B. drift bound over a long chain at the model's target level ----
  lev = target_level
  if (lev < 1) lev = max(1, npos_val / 4)
  if (lev > npos_val - 1) lev = npos_val - 1
  call check_drift(lev, 50000, n_fail, max_drift_err)

  write(output_unit,'(A,I0)')        '   levels tested      : ', n_levels_tested
  write(output_unit,'(A,I0)')        '   swaps per level    : ', n_swaps_per_level
  write(output_unit,'(A,ES12.4)')    '   max term error     : ', max_term_err
  write(output_unit,'(A,ES12.4)')    '   max energy error   : ', max_energy_err
  write(output_unit,'(A,ES12.4)')    '   max chain drift    : ', max_drift_err

  if (n_fail == 0) then
    write(output_unit,'(A)') ' test_pme_delta: PASS'
    call pme_finalize_model()
    stop 0
  else
    write(error_unit,'(A,I0,A)') ' test_pme_delta: FAIL (', n_fail, ' mismatches)'
    call pme_finalize_model()
    stop 1
  end if

contains

  subroutine bench_scaling(npos_val)
    !  Per-call timing: full pme_evaluate_configuration vs the incremental
    !  old_terms + pme_evaluate_swap_delta, at several substitution levels.
    !  Both paths do the same per-step work as the MC loop (generate a swap,
    !  evaluate, apply epsilon), so the ratio is the per-step speedup. The
    !  ratio grows with cluster size (L on the low side, H = npos-L on the
    !  high side), illustrating the O(N^k) -> O(N^(k-1)) reduction.
    integer, intent(in) :: npos_val
    integer, parameter :: m_calls = 40000
    integer :: levels(9), nlev, il, level, i, a, b, remove_idx
    integer, allocatable :: subset(:), newset(:)
    real(real64) :: e, elo, ehi, lo(4), hi(4), dlo(4), dhi(4), en
    real(real64) :: t0, t1, t_full, t_delta

    levels = [2, 4, 6, 8, 10, 12, 16, 20, npos_val - 1]
    nlev = size(levels)

    write(output_unit,'(A)') ' --- per-call scaling: full recompute vs incremental delta ---'
    write(output_unit,'(A,I0,A)') '   (', m_calls, ' evaluations per cell, model order as loaded)'
    write(output_unit,'(A)') '   level   H=npos-L     full(s)    delta(s)   speedup'

    do il = 1, nlev
      level = levels(il)
      if (level < 1 .or. level > npos_val - 1) cycle
      if (il > 1) then
        if (level <= levels(il-1)) cycle   ! skip duplicates (e.g. npos-1 already covered)
      end if
      allocate(subset(level), newset(level))
      call random_subset(npos_val, level, subset)

      ! full recompute path
      call cpu_time(t0)
      do i = 1, m_calls
        call pick_swap(subset, level, npos_val, remove_idx, a, b)
        if (a < 0) cycle
        newset = subset; newset(remove_idx) = b
        call mc_insertion_sort(newset, level)
        call pme_evaluate_configuration(newset, level, e, elo, ehi, lo, hi)
        call apply_epsilon_energy(level, lo, hi, en)
      end do
      call cpu_time(t1); t_full = t1 - t0

      ! incremental delta path
      call random_subset(npos_val, level, subset)
      call pme_evaluate_configuration(subset, level, e, elo, ehi, lo, hi)
      call cpu_time(t0)
      do i = 1, m_calls
        call pick_swap(subset, level, npos_val, remove_idx, a, b)
        if (a < 0) cycle
        call pme_evaluate_swap_delta(subset, level, a, b, dlo, dhi)
        call apply_epsilon_energy(level, lo + dlo, hi + dhi, en)
      end do
      call cpu_time(t1); t_delta = t1 - t0

      write(output_unit,'(I8,I11,2F12.4,F10.2,A)') level, npos_val - level, &
        t_full, t_delta, t_full / max(t_delta, 1.0e-12_real64), 'x'
      deallocate(subset, newset)
    end do
  end subroutine bench_scaling

  subroutine check_level_equivalence(level, n_swaps, n_fail, max_term_err, max_energy_err)
    integer, intent(in)    :: level, n_swaps
    integer, intent(inout) :: n_fail
    real(real64), intent(inout) :: max_term_err, max_energy_err

    integer :: subset(level), newset(level)
    integer :: holes(max(1, npos_val - level))
    integer :: a, b, remove_idx, i, s, hc
    real(real64) :: e_old, elo_old, ehi_old, low_old(4), high_old(4)
    real(real64) :: e_new, elo_new, ehi_new, low_new(4), high_new(4)
    real(real64) :: dlow(4), dhigh(4), dlow2(4), dhigh2(4)
    real(real64) :: low_pred(4), high_pred(4)
    real(real64) :: en_full, en_pred, terr, eerr

    call random_subset(npos_val, level, subset)

    do i = 1, n_swaps
      call pick_swap(subset, level, npos_val, remove_idx, a, b)
      if (a < 0) cycle   ! no valid swap (level == npos); skipped by caller range anyway

      ! full evaluation of the current config (term vectors only matter)
      call pme_evaluate_configuration(subset, level, e_old, elo_old, ehi_old, low_old, high_old)

      ! incremental delta for the swap a -> b (build-path: routine rebuilds holes)
      call pme_evaluate_swap_delta(subset, level, a, b, dlow, dhigh)

      ! same delta via the optional holes_in path, with holes in DESCENDING order
      ! to confirm the result is independent of the supplied hole-list ordering
      ! (this is the path the MC loop uses with its maintained hole list).
      hc = 0
      do s = npos_val, 1, -1
        if (.not. any(subset(1:level) == s)) then
          hc = hc + 1
          holes(hc) = s
        end if
      end do
      call pme_evaluate_swap_delta(subset, level, a, b, dlow2, dhigh2, holes, hc)

      ! apply swap and full-evaluate the new config
      newset = subset
      newset(remove_idx) = b
      call mc_insertion_sort(newset, level)
      call pme_evaluate_configuration(newset, level, e_new, elo_new, ehi_new, low_new, high_new)

      low_pred  = low_old  + dlow
      high_pred = high_old + dhigh

      terr = max(maxval(abs(low_pred - low_new)), maxval(abs(high_pred - high_new)))
      terr = max(terr, maxval(abs(low_old + dlow2 - low_new)), maxval(abs(high_old + dhigh2 - high_new)))
      if (terr > max_term_err) max_term_err = terr
      if (terr > tol_term) then
        n_fail = n_fail + 1
        if (n_fail <= 3) then
          write(error_unit,'(A,I0,A,ES12.4)') '   term mismatch at level ', level, ' : ', terr
          write(error_unit,'(A,I0,A,I0,A,I0)') '     a=', a, ' b=', b, ' remove_idx=', remove_idx
          write(error_unit,'(A,4ES14.6)') '     low_old  =', low_old
          write(error_unit,'(A,4ES14.6)') '     dlow     =', dlow
          write(error_unit,'(A,4ES14.6)') '     low_new  =', low_new
          write(error_unit,'(A,4ES14.6)') '     high_old =', high_old
          write(error_unit,'(A,4ES14.6)') '     dhigh    =', dhigh
          write(error_unit,'(A,4ES14.6)') '     high_new =', high_new
        end if
      end if

      ! calibrated-energy agreement (the quantity the Metropolis test uses)
      call apply_epsilon_energy(level, low_new,  high_new,  en_full)
      call apply_epsilon_energy(level, low_pred, high_pred, en_pred)
      eerr = abs(en_full - en_pred)
      if (eerr > max_energy_err) max_energy_err = eerr
      if (eerr > tol_energy) then
        n_fail = n_fail + 1
        if (n_fail <= 10) write(error_unit,'(A,I0,A,ES12.4)') &
          '   energy mismatch at level ', level, ' : ', eerr
      end if

      subset = newset   ! continue the walk from the accepted config
    end do
  end subroutine check_level_equivalence

  subroutine check_drift(level, n_steps, n_fail, max_drift_err)
    !  Maintain the term vectors using ONLY deltas through a long swap chain and
    !  compare to a full recompute every step; this measures FP accumulation.
    integer, intent(in)    :: level, n_steps
    integer, intent(inout) :: n_fail
    real(real64), intent(inout) :: max_drift_err

    integer :: subset(level)
    integer :: a, b, remove_idx, i
    real(real64) :: e0, elo0, ehi0
    real(real64) :: low_run(4), high_run(4)
    real(real64) :: low_full(4), high_full(4)
    real(real64) :: dlow(4), dhigh(4)
    real(real64) :: e_full, elo_f, ehi_f, derr

    call random_subset(npos_val, level, subset)
    call pme_evaluate_configuration(subset, level, e0, elo0, ehi0, low_run, high_run)

    do i = 1, n_steps
      call pick_swap(subset, level, npos_val, remove_idx, a, b)
      if (a < 0) cycle
      call pme_evaluate_swap_delta(subset, level, a, b, dlow, dhigh)
      low_run  = low_run  + dlow
      high_run = high_run + dhigh
      subset(remove_idx) = b
      call mc_insertion_sort(subset, level)

      call pme_evaluate_configuration(subset, level, e_full, elo_f, ehi_f, low_full, high_full)
      derr = max(maxval(abs(low_run - low_full)), maxval(abs(high_run - high_full)))
      if (derr > max_drift_err) max_drift_err = derr
    end do

    if (max_drift_err > tol_drift) then
      n_fail = n_fail + 1
      write(error_unit,'(A,ES12.4)') '   chain drift exceeds tolerance: ', max_drift_err
    end if
  end subroutine check_drift

  subroutine pick_swap(subset, level, npos_val, remove_idx, a, b)
    !  Choose a random occupied site to remove (a) and a random hole to add (b).
    integer, intent(in)  :: subset(:), level, npos_val
    integer, intent(out) :: remove_idx, a, b
    real(real64) :: r
    integer :: tries

    remove_idx = -1; a = -1; b = -1
    if (level <= 0 .or. npos_val <= level) return

    call random_number(r)
    remove_idx = min(int(r * real(level, real64)) + 1, level)
    a = subset(remove_idx)

    tries = 0
    do while (tries < npos_val * 8)
      tries = tries + 1
      call random_number(r)
      b = min(int(r * real(npos_val, real64)) + 1, npos_val)
      if (.not. any(subset(1:level) == b)) return
    end do
    a = -1; b = -1   ! failed to find a hole
  end subroutine pick_swap

  subroutine random_subset(n, k, subset)
    integer, intent(in)  :: n, k
    integer, intent(out) :: subset(:)
    integer :: pool(n), tmp, i, j
    real(real64) :: r
    do i = 1, n
      pool(i) = i
    end do
    do i = 1, k
      call random_number(r)
      j = min(i + int(r * real(n - i + 1, real64)), n)
      tmp = pool(i); pool(i) = pool(j); pool(j) = tmp
    end do
    subset(1:k) = pool(1:k)
    call mc_insertion_sort(subset, k)
  end subroutine random_subset

  subroutine seed_fixed(seed_val)
    integer, intent(in) :: seed_val
    integer :: n
    integer, allocatable :: seed_arr(:)
    call random_seed(size=n)
    allocate(seed_arr(n))
    seed_arr = seed_val
    call random_seed(put=seed_arr)
  end subroutine seed_fixed

end program test_pme_delta
