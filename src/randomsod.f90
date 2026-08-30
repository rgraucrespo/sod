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
!    Supports the general substitution model: multiple target sublattices
!    (ntarget) and multiple substituting species per site (multinary, nk>1).
!    Each draw independently places, on every target sublattice, a uniform
!    colored assignment of the requested per-species counts.
!
!    With -sym on draws are folded to symmetry representatives and
!    the degeneracy column holds visit counts: the uniform draw already visits
!    each orbit in proportion to its size, so visit counts are the correct
!    importance weights for canonical averages in statsod. -nconf is the
!    number of draws (not the number of distinct configurations, which is not
!    known a priori). Symmetry folding uses a colored orbit representative, so it
!    works for the general multinary/multisite case.
!
!    USAGE
!      randomsod -nconf <N> [-sym on|off] [-seed clock|<int>]
!
!    Required input files: INSOD, SGO (always); EQMATRIX (when -symmetry on).
!    Output: nXX/random/ENSEMBLE (XX = target level).
!
!    Part of the SOD package (v0.91) — GNU GPL v3+.
!*******************************************************************************

program randomsod
  use iso_fortran_env, only: real64, error_unit
  use insod_reader,    only: insod_t, read_insod
  use config_sampling, only: derive_target_geometry_all, load_eqmatrix_all, &
                             canonicalize_colored, mc_random_colored, &
                             seed_rng_clock, seed_rng_fixed
  use ensemble_io,     only: write_ensemble
  implicit none

  type(insod_t) :: d
  integer :: nconfigs, ntarget, nsubs_tot, max_nk, npos_max
  integer :: nop, level_t
  integer, allocatable :: npos_t(:), atini_t(:), atfin_t(:)
  integer, allocatable :: nk_t(:), nsubs2(:, :)
  integer, allocatable :: eqmt(:, :, :)
  logical :: use_symmetry, found
  integer :: iseed                       ! -1 = clock, >0 = fixed
  character(len=64) :: seed_label

  integer :: arg_count, iarg, ios, mkdir_status
  character(len=256) :: arg, val
  logical :: have_nconfigs

  integer, allocatable :: conf(:), canonical(:)
  integer, allocatable :: unique_confs(:, :), visit_count(:)
  integer :: n_unique, max_unique, i_draw, j, t, off

  character(len=128) :: nxx_dir
  character(len=256) :: ensemble_path
  integer :: unit_ens

  write(*, '(A)') "SOD (Site-Occupancy Disorder) version 0.92 - randomsod"

  ! --- Defaults ---
  nconfigs      = 0
  have_nconfigs = .false.
  use_symmetry  = .false.       ! -sym off by default
  iseed         = -1            ! -seed clock by default

  ! --- Parse command line ---
  arg_count = command_argument_count()
  iarg = 1
  do while (iarg <= arg_count)
    call get_command_argument(iarg, arg)
    select case (trim(arg))
    case ('-nconf')
      call require_value(iarg, arg_count, '-nconf', val)
      read(val, *, iostat=ios) nconfigs
      if (ios /= 0 .or. nconfigs <= 0) then
        write(error_unit,'(A,A)') ' Error: -nconf requires a positive integer, got: ', trim(val)
        stop 1
      end if
      have_nconfigs = .true.
    case ('-sym')
      call require_value(iarg, arg_count, '-sym', val)
      select case (trim(val))
      case ('on');  use_symmetry = .true.
      case ('off'); use_symmetry = .false.
      case default
        write(error_unit,'(A,A)') ' Error: -sym must be on or off, got: ', trim(val)
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
      write(error_unit,'(A)')   ' Usage: sod_random.sh -nconf <N> [-sym on|off] [-seed clock|<int>]'
      stop 1
    end select
    iarg = iarg + 1
  end do

  if (.not. have_nconfigs) then
    write(error_unit,'(A)') ' Error: -nconf is required.'
    write(error_unit,'(A)') ' Usage: sod_random.sh -nconf <N> [-sym on|off] [-seed clock|<int>]'
    stop 1
  end if

  ! --- Read INSOD ---
  call read_insod('INSOD', d)
  ntarget = d%ntarget

  ! Per-target species-per-site counts (nk_t) and per-species substitution
  ! counts (nsubs2). For the single-target binary case nsubs_max is authoritative
  ! (it also carries the colon-range upper endpoint, where nsubs_t is left 0).
  max_nk = maxval(d%nk(1:ntarget))
  allocate(nk_t(ntarget), nsubs2(ntarget, max_nk))
  nsubs2 = 0
  do t = 1, ntarget
    nk_t(t) = d%nk(t)
    nsubs2(t, 1:nk_t(t)) = d%nsubs_t(t, 1:nk_t(t))
  end do
  if (ntarget == 1 .and. nk_t(1) == 1) nsubs2(1, 1) = d%nsubs_max
  if (d%nsubs_min /= d%nsubs_max) &
    write(*,'(A,I0)') ' > INSOD range detected; sampling at the upper endpoint: nsubs = ', d%nsubs_max

  ! --- Derive per-target geometry (INSOD + SGO) ---
  call derive_target_geometry_all(d, 'SGO', npos_t, atini_t, atfin_t)
  npos_max = maxval(npos_t(1:ntarget))

  ! --- Validate substitution counts against available sites ---
  nsubs_tot = 0
  do t = 1, ntarget
    level_t = sum(nsubs2(t, 1:nk_t(t)))
    if (level_t < 0 .or. level_t > npos_t(t)) then
      write(error_unit,'(A,I0,A,I0,A,I0,A)') ' Error: target ', t, ' requests ', level_t, &
        ' substitutions but has only ', npos_t(t), ' sites.'
      stop 1
    end if
    nsubs_tot = nsubs_tot + level_t
  end do
  if (nsubs_tot < 1) then
    write(error_unit,'(A)') ' Error: no substitutions requested (all counts zero).'
    stop 1
  end if

  ! --- Symmetry operators (only when reducing) ---
  if (use_symmetry) &
    call load_eqmatrix_all('EQMATRIX', ntarget, npos_t, eqmt, nop)

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
  write(*, '(A,I0)')        '  Target sublattices: ', ntarget
  do t = 1, ntarget
    level_t = sum(nsubs2(t, 1:nk_t(t)))
    write(*, '(A,I0,A,I0,A,I0,A)') '   - target ', t, ': ', level_t, &
      ' substitutions in ', npos_t(t), ' sites'
  end do
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
  allocate(conf(nsubs_tot), canonical(nsubs_tot))
  allocate(unique_confs(nsubs_tot, max_unique), visit_count(max_unique))
  n_unique = 0

  do i_draw = 1, nconfigs
    ! Uniform colored draw on every target sublattice.
    off = 0
    do t = 1, ntarget
      level_t = sum(nsubs2(t, 1:nk_t(t)))
      if (level_t > 0) &
        call mc_random_colored(npos_t(t), nk_t(t), nsubs2(t, 1:nk_t(t)), conf(off+1:off+level_t))
      off = off + level_t
    end do

    if (use_symmetry) then
      ! Fold to the colored orbit representative and accumulate visit counts.
      call canonicalize_colored(conf, ntarget, nk_t, nsubs2, eqmt, nop, canonical)
      found = .false.
      do j = 1, n_unique
        if (all(unique_confs(1:nsubs_tot, j) == canonical(1:nsubs_tot))) then
          visit_count(j) = visit_count(j) + 1
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        n_unique = n_unique + 1
        unique_confs(:, n_unique) = canonical
        visit_count(n_unique) = 1
      end if
    else
      unique_confs(:, i_draw) = conf
      visit_count(i_draw) = 1
    end if
  end do

  if (.not. use_symmetry) n_unique = nconfigs

  ! --- Output directory name: 'n' + all per-species counts joined by '_' ---
  nxx_dir = 'n'
  do t = 1, ntarget
    do j = 1, nk_t(t)
      if (t == 1 .and. j == 1) then
        nxx_dir = 'n'//trim(npad(nsubs2(t, j)))
      else
        nxx_dir = trim(nxx_dir)//'_'//trim(npad(nsubs2(t, j)))
      end if
    end do
  end do
  nxx_dir = trim(nxx_dir)//'/random'
  call execute_command_line('mkdir -p '//trim(nxx_dir), exitstat=mkdir_status)
  if (mkdir_status /= 0) then
    write(error_unit,'(A,A)') ' Error: could not create output directory ', trim(nxx_dir)
    stop 1
  end if
  ensemble_path = trim(nxx_dir)//'/ENSEMBLE'

  ! --- Write ENSEMBLE ---
  open(newunit=unit_ens, file=trim(ensemble_path), status='replace', action='write')
  block
    integer, allocatable :: mc_ic(:, :), nsubs_flat(:)
    character(len=3), allocatable :: orig_sym(:)
    integer :: m, nkflat, foff
    allocate(mc_ic(n_unique, nsubs_tot))
    do m = 1, n_unique
      mc_ic(m, 1:nsubs_tot) = unique_confs(1:nsubs_tot, m)
    end do
    nkflat = sum(nk_t(1:ntarget))
    allocate(nsubs_flat(nkflat), orig_sym(ntarget))
    foff = 0
    do t = 1, ntarget
      nsubs_flat(foff+1:foff+nk_t(t)) = nsubs2(t, 1:nk_t(t))
      foff = foff + nk_t(t)
      orig_sym(t) = d%symbol(d%sptarget(t))
    end do
    call write_ensemble(unit_ens, ntarget, nk_t, nsubs_flat, npos_t, atini_t, &
                        n_unique, mc_ic, visit_count(1:n_unique), 'uniform', &
                        orig_sym, d%newsymbol(1:ntarget, 1:max_nk+1))
  end block
  close(unit_ens)

  ! --- Summary ---
  write(*, '(A)') ' --- Output -----------------------------------------------------------------'
  write(*, '(A,I0)')  '  Draws             : ', nconfigs
  if (use_symmetry) &
  write(*, '(A,I0)')  '  Distinct configs  : ', n_unique
  write(*, '(A,A)')   '  ENSEMBLE          : ', trim(ensemble_path)
  write(*, '(A)') ' ---------------------------------------------------------------------------'
  write(*, *)
  write(*, '(A)') ' ============================================================================'
  write(*, '(A)') '  randomsod completed.'
  write(*, '(A)') ' ============================================================================'

  deallocate(conf, canonical, unique_confs, visit_count)
  deallocate(npos_t, atini_t, atfin_t, nk_t, nsubs2)
  if (allocated(eqmt)) deallocate(eqmt)

contains

  function npad(n) result(s)
    integer, intent(in) :: n
    character(len=10) :: s
    integer :: nd, tmp
    character(len=20) :: fmt

    tmp = max(1, abs(npos_max))
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
