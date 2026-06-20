!*******************************************************************************
!    randomsod — uniform random sampling of the configuration space.
!
!    Draws nconfigs independent uniform configurations at the INSOD target
!    substitution level and writes them as an ENSEMBLE. No energies are involved:
!    the sample is Hamiltonian-independent, so energies (if wanted) are computed
!    a posteriori by the normal structure-writer -> DFT -> statsod path, exactly
!    as for an enumerated combsod ensemble. randomsod is the sampling counterpart
!    of combsod, for substitution levels too large to enumerate.
!
!    With -symmetry on (default) draws are folded to symmetry representatives and
!    the degeneracy column holds visit counts: the uniform draw already visits
!    each orbit in proportion to its size, so visit counts are the correct
!    importance weights for canonical averages in statsod. -nconfigs is the
!    number of draws (not the number of distinct configurations, which is not
!    known a priori).
!
!    USAGE
!      randomsod -nconfigs <N> [-symmetry on|off] [-seed clock|<int>]
!
!    Required input files: INSOD, SGO (always); EQMATRIX (when -symmetry on).
!    Output: nXX/random/ENSEMBLE (XX = target level).
!
!    Part of the SOD package (v0.82) — GNU GPL v3+.
!*******************************************************************************

program randomsod
  use iso_fortran_env, only: real64, error_unit
  use insod_reader,    only: insod_t, read_insod
  use config_sampling, only: derive_target_geometry, load_eqmatrix, &
                             canonicalize_config, mc_random_subset, &
                             seed_rng_clock, seed_rng_fixed
  use ensemble_io,     only: write_ensemble
  implicit none

  type(insod_t) :: d
  integer :: nconfigs, level, npos, atini, atfin, species_index
  integer :: nop, npos_eqm
  integer, allocatable :: eqmatrix(:, :)
  logical :: use_symmetry
  integer :: iseed                       ! -1 = clock, >0 = fixed
  character(len=64) :: seed_label

  integer :: arg_count, iarg, ios, mkdir_status
  character(len=256) :: arg, val
  logical :: have_nconfigs

  integer, allocatable :: subset(:), canonical(:)
  integer, allocatable :: unique_subsets(:, :), visit_count(:)
  integer :: n_unique, max_unique, i_draw, j
  logical :: found

  character(len=64) :: nxx_dir
  character(len=256) :: ensemble_path
  integer :: unit_ens

  write(*, '(A)') "SOD (Site-Occupancy Disorder) version 0.82 - randomsod"

  ! --- Defaults ---
  nconfigs      = 0
  have_nconfigs = .false.
  use_symmetry  = .true.        ! -symmetry on by default
  iseed         = -1            ! -seed clock by default

  ! --- Parse command line ---
  arg_count = command_argument_count()
  iarg = 1
  do while (iarg <= arg_count)
    call get_command_argument(iarg, arg)
    select case (trim(arg))
    case ('-nconfigs')
      call require_value(iarg, arg_count, '-nconfigs', val)
      read(val, *, iostat=ios) nconfigs
      if (ios /= 0 .or. nconfigs <= 0) then
        write(error_unit,'(A,A)') ' Error: -nconfigs requires a positive integer, got: ', trim(val)
        stop 1
      end if
      have_nconfigs = .true.
    case ('-symmetry')
      call require_value(iarg, arg_count, '-symmetry', val)
      select case (trim(val))
      case ('on');  use_symmetry = .true.
      case ('off'); use_symmetry = .false.
      case default
        write(error_unit,'(A,A)') ' Error: -symmetry must be on or off, got: ', trim(val)
        stop 1
      end select
    case ('-seed')
      call require_value(iarg, arg_count, '-seed', val)
      if (trim(val) == 'clock') then
        iseed = -1
      else
        read(val, *, iostat=ios) iseed
        if (ios /= 0 .or. iseed <= 0) then
          write(error_unit,'(A,A)') ' Error: -seed must be "clock" or a positive integer, got: ', trim(val)
          stop 1
        end if
      end if
    case default
      write(error_unit,'(A,A)') ' Error: unrecognised argument: ', trim(arg)
      write(error_unit,'(A)')   ' Usage: randomsod -nconfigs <N> [-symmetry on|off] [-seed clock|<int>]'
      stop 1
    end select
    iarg = iarg + 1
  end do

  if (.not. have_nconfigs) then
    write(error_unit,'(A)') ' Error: -nconfigs is required.'
    write(error_unit,'(A)') ' Usage: randomsod -nconfigs <N> [-symmetry on|off] [-seed clock|<int>]'
    stop 1
  end if

  ! --- Read INSOD and derive the target geometry (INSOD + SGO) ---
  call read_insod('INSOD', d)
  level = d%nsubs_max
  if (d%nsubs_min /= d%nsubs_max) &
    write(*,'(A,I0)') ' > INSOD range detected; sampling at the upper endpoint: nsubs = ', level

  call derive_target_geometry(d, 'SGO', species_index, npos, atini, atfin)

  if (level < 1 .or. level > npos) then
    write(error_unit,'(A,I0,A,I0,A)') ' Error: target level (', level, &
      ') must be between 1 and the number of target sites (', npos, ').'
    stop 1
  end if

  ! --- Symmetry operators (only when reducing) ---
  if (use_symmetry) then
    call load_eqmatrix('EQMATRIX', eqmatrix, nop, npos_eqm)
    if (npos_eqm /= npos) then
      write(error_unit,'(A,I0,A,I0,A)') ' Error: EQMATRIX target count (', npos_eqm, &
        ') does not match INSOD/SGO target count (', npos, ').'
      stop 1
    end if
  end if

  ! --- Seed the RNG ---
  if (iseed > 0) then
    call seed_rng_fixed(iseed)
    write(seed_label, '(I0)') iseed
  else
    call seed_rng_clock()
    seed_label = 'system clock'
  end if

  ! --- Report parameters ---
  write(*, '(A)') ' --- Random sampling parameters ---------------------------------------------'
  write(*, '(A,I0,A,I0,A)') '  Level             : ', level, ' substitutions in ', npos, ' positions'
  write(*, '(A,I0)')        '  Draws (nconfigs)  : ', nconfigs
  if (use_symmetry) then
    write(*, '(A,I0,A)')    '  Symmetry          : on (', nop, ' operators; visit-count degeneracies)'
  else
    write(*, '(A)')         '  Symmetry          : off (one row per draw)'
  end if
  write(*, '(A,A)')         '  Random seed       : ', trim(seed_label)
  write(*, '(A)') ' ---------------------------------------------------------------------------'
  write(*, *)

  ! --- Draw loop ---
  max_unique = nconfigs
  allocate(subset(level), canonical(level))
  allocate(unique_subsets(level, max_unique), visit_count(max_unique))
  n_unique = 0

  if (use_symmetry) then
    do i_draw = 1, nconfigs
      call mc_random_subset(npos, level, subset)
      call canonicalize_config(subset, level, eqmatrix, nop, canonical)
      found = .false.
      do j = 1, n_unique
        if (all(unique_subsets(1:level, j) == canonical(1:level))) then
          visit_count(j) = visit_count(j) + 1
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        n_unique = n_unique + 1
        unique_subsets(:, n_unique) = canonical
        visit_count(n_unique) = 1
      end if
    end do
  else
    do i_draw = 1, nconfigs
      call mc_random_subset(npos, level, subset)
      unique_subsets(:, i_draw) = subset
      visit_count(i_draw) = 1
    end do
    n_unique = nconfigs
  end if

  ! --- Write nXX/random/ENSEMBLE ---
  nxx_dir = 'n'//trim(npad(level))//'/random'
  call execute_command_line('mkdir -p '//trim(nxx_dir), exitstat=mkdir_status)
  if (mkdir_status /= 0) then
    write(error_unit,'(A,A)') ' Error: could not create output directory ', trim(nxx_dir)
    stop 1
  end if
  ensemble_path = trim(nxx_dir)//'/ENSEMBLE'

  open(newunit=unit_ens, file=trim(ensemble_path), status='replace', action='write')
  block
    integer, allocatable :: mc_ic(:,:)
    integer :: m
    allocate(mc_ic(n_unique, level))
    do m = 1, n_unique
      mc_ic(m, 1:level) = unique_subsets(1:level, m)
    end do
    call write_ensemble(unit_ens, 1, [1], [level], [npos], [atini], &
                        n_unique, mc_ic, visit_count(1:n_unique), 'uniform', &
                        [d%symbol(species_index)], d%newsymbol(1:1, 1:2))
  end block
  close(unit_ens)

  ! --- Summary ---
  write(*, '(A)') ' --- Output -----------------------------------------------------------------'
  write(*, '(A,I0)')  '  Draws             : ', nconfigs
  write(*, '(A,I0)')  '  Distinct configs  : ', n_unique
  write(*, '(A,A)')   '  ENSEMBLE          : ', trim(ensemble_path)
  write(*, '(A)') ' ---------------------------------------------------------------------------'
  write(*, *)
  write(*, '(A)') ' ============================================================================'
  write(*, '(A)') '  randomsod completed.'
  write(*, '(A)') ' ============================================================================'

  deallocate(subset, canonical, unique_subsets, visit_count)
  if (allocated(eqmatrix)) deallocate(eqmatrix)

contains

  function npad(n) result(s)
    integer, intent(in) :: n
    character(len=10) :: s
    integer :: nd, tmp
    character(len=20) :: fmt

    tmp = max(1, abs(npos))
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

  subroutine require_value(iarg, arg_count, flag, val)
    !  Advance to and return the value following a flag; error if missing.
    integer, intent(inout)        :: iarg
    integer, intent(in)           :: arg_count
    character(len=*), intent(in)  :: flag
    character(len=*), intent(out) :: val
    if (iarg + 1 > arg_count) then
      write(error_unit,'(A,A,A)') ' Error: ', trim(flag), ' requires a value.'
      stop 1
    end if
    iarg = iarg + 1
    call get_command_argument(iarg, val)
  end subroutine require_value

end program randomsod
