!*******************************************************************************
!    config_sampling — shared generator-side utilities for SOD configuration
!    sampling (randomsod and the Metropolis mcsod chain).
!
!    Two groups of routines:
!      (1) Target geometry: derive_target_geometry expands the INSOD asymmetric
!          unit through the SGO space-group operators to the conventional cell,
!          multiplies by the supercell, and returns the number of target sites
!          (npos) and the global atom offset of the target species (atini/atfin).
!          This is the single source of truth shared with pmemod; combsod still
!          carries an independent copy pending a later dedup.
!      (2) Sampling helpers: uniform random subsets, symmetry canonicalization,
!          EQMATRIX loading, RNG seeding and small input/string helpers.
!
!    Part of the SOD package (v0.83) — GNU GPL v3+.
!*******************************************************************************

module config_sampling
  use iso_fortran_env, only: real64, error_unit
  use insod_reader,    only: insod_t
  implicit none
  private

  public :: derive_target_geometry, wrap_fractional
  public :: load_eqmatrix, canonicalize_config
  public :: mc_random_subset, mc_insertion_sort
  public :: seed_rng_clock, seed_rng_fixed
  public :: read_next_data_line, to_lower

contains

  ! =========================================================================
  !  Target geometry
  ! =========================================================================

  subroutine derive_target_geometry(d, sgo_filename, species_index, &
                                     npos_out, atini_out, atfin_out)
    !  Expand the INSOD asymmetric unit (d%coords0, d%natsp0) through the SGO
    !  operators to the conventional cell, deduplicate, multiply by the
    !  d%na*d%nb*d%nc supercell, and report:
    !    species_index = d%sptarget(1)
    !    npos_out      = number of target-species sites in the supercell
    !    atini_out     = global index of the first target-species atom
    !    atfin_out     = global index of the last target-species atom
    !  The caller owns INSOD (passed in as d) and any EQMATRIX cross-check.
    type(insod_t),    intent(in)  :: d
    character(len=*), intent(in)  :: sgo_filename
    integer,          intent(out) :: species_index, npos_out, atini_out, atfin_out

    integer :: unit_sgo, ios, nsp, nat0, at0, op1, nop1
    integer :: na, nb, nc, sp, cumnatsp, at1r, at1, at1i, nat1r, nat1
    integer :: nat_super
    integer, allocatable :: natsp0(:), natsp1(:), natsp(:), spat0(:), spat1(:), spat1r(:)
    real(real64), allocatable :: coords0(:, :), coords1(:, :), coords1r(:, :)
    real(real64), allocatable :: mgroup1(:, :, :), vgroup1(:, :)
    real(real64) :: prod
    logical :: found
    real(real64), parameter :: tol0 = 1.0e-3_real64

    nsp           = d%nsp
    nat0          = d%nat0
    na = d%na; nb = d%nb; nc = d%nc
    species_index = d%sptarget(1)
    if (species_index <= 0) stop 'Error: invalid sptarget in INSOD.'

    allocate(natsp0(nsp), natsp1(nsp), natsp(nsp))
    natsp0 = d%natsp0
    allocate(coords0(nat0, 3), spat0(nat0))
    coords0(1:nat0, 1:3) = d%coords0

    open(newunit=unit_sgo, file=trim(sgo_filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(error_unit,'(A)') 'Error: could not open SGO while deriving target geometry.'
      stop 1
    end if

    cumnatsp = 0
    do sp = 1, nsp
      do at0 = cumnatsp + 1, cumnatsp + natsp0(sp)
        spat0(at0) = sp
      end do
      cumnatsp = cumnatsp + natsp0(sp)
    end do

    read(unit_sgo,*,iostat=ios)
    if (ios /= 0) stop 'Error: malformed SGO header.'
    read(unit_sgo,*,iostat=ios) op1
    if (ios /= 0) stop 'Error: malformed SGO operator header.'
    nop1 = 0
    do while (op1 > 0)
      nop1 = max(nop1, op1)
      do at0 = 1, 3
        read(unit_sgo,*,iostat=ios)
        if (ios /= 0) stop 'Error: malformed SGO operator.'
      end do
      read(unit_sgo,*,iostat=ios) op1
      if (ios /= 0) stop 'Error: malformed SGO operator terminator.'
    end do
    rewind(unit_sgo)
    read(unit_sgo,*)
    allocate(mgroup1(nop1, 3, 3), vgroup1(nop1, 3))
    mgroup1 = 0.0_real64
    vgroup1 = 0.0_real64
    read(unit_sgo,*) op1
    do while (op1 > 0)
      do at0 = 1, 3
        read(unit_sgo,*,iostat=ios) mgroup1(op1, at0, 1:3), vgroup1(op1, at0)
        if (ios /= 0) stop 'Error: malformed SGO operator body.'
      end do
      read(unit_sgo,*,iostat=ios) op1
      if (ios /= 0) stop 'Error: malformed SGO operator terminator.'
    end do
    close(unit_sgo)

    allocate(coords1r(nat0 * nop1, 3), spat1r(nat0 * nop1))
    at1r = 0
    do at0 = 1, nat0
      do op1 = 1, nop1
        at1r = at1r + 1
        coords1r(at1r, 1:3) = matmul(mgroup1(op1, 1:3, 1:3), coords0(at0, 1:3)) + vgroup1(op1, 1:3)
        coords1r(at1r, 1) = wrap_fractional(coords1r(at1r, 1))
        coords1r(at1r, 2) = wrap_fractional(coords1r(at1r, 2))
        coords1r(at1r, 3) = wrap_fractional(coords1r(at1r, 3))
        spat1r(at1r) = spat0(at0)
      end do
    end do
    nat1r = at1r

    allocate(coords1(nat1r, 3), spat1(nat1r))
    nat1 = 0
    do at1r = 1, nat1r
      found = .false.
      do at1i = 1, nat1
        prod = sum((coords1r(at1r, 1:3) - coords1(at1i, 1:3))**2)
        if (prod <= tol0) then
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        nat1 = nat1 + 1
        coords1(nat1, 1:3) = coords1r(at1r, 1:3)
        spat1(nat1) = spat1r(at1r)
      end if
    end do

    natsp1 = 0
    do at1 = 1, nat1
      natsp1(spat1(at1)) = natsp1(spat1(at1)) + 1
    end do
    nat_super = na * nb * nc
    natsp = nat_super * natsp1

    atini_out = 1
    do sp = 1, species_index - 1
      atini_out = atini_out + natsp(sp)
    end do
    atfin_out = atini_out + natsp(species_index) - 1
    npos_out  = natsp(species_index)
    if (npos_out <= 0) then
      write(error_unit,'(A)') 'Error: target species has zero multiplicity in the supercell.'
      stop 1
    end if

    deallocate(natsp0, natsp1, natsp, spat0, spat1, spat1r)
    deallocate(coords0, coords1, coords1r, mgroup1, vgroup1)
  end subroutine derive_target_geometry

  pure real(real64) function wrap_fractional(x) result(corx)
    real(real64), intent(in) :: x
    real(real64), parameter :: tol1 = 1.0e-4_real64

    corx = x - floor(x)
    if (corx > (1.0_real64 - tol1)) corx = 0.0_real64
  end function wrap_fractional

  ! =========================================================================
  !  EQMATRIX loading and symmetry canonicalization
  ! =========================================================================

  subroutine load_eqmatrix(filename, eqmatrix, nop, npos)
    !  Read EQMATRIX: header line "nop npos" followed by nop rows of npos ints.
    character(len=*),     intent(in)  :: filename
    integer, allocatable, intent(out) :: eqmatrix(:, :)
    integer,              intent(out) :: nop, npos
    integer :: unit_in, ios, op

    open(newunit=unit_in, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(error_unit,'(A,1X,A)') 'Error: could not open EQMATRIX file', trim(filename)
      stop 1
    end if

    read(unit_in,*,iostat=ios) nop, npos
    if (ios /= 0 .or. nop <= 0 .or. npos <= 0) then
      close(unit_in)
      write(error_unit,'(A)') 'Error: invalid EQMATRIX header.'
      stop 1
    end if
    allocate(eqmatrix(nop, npos))
    do op = 1, nop
      read(unit_in,*,iostat=ios) eqmatrix(op, 1:npos)
      if (ios /= 0) then
        close(unit_in)
        write(error_unit,'(A)') 'Error: invalid EQMATRIX row.'
        stop 1
      end if
    end do
    close(unit_in)
  end subroutine load_eqmatrix

  subroutine canonicalize_config(conf, level, eqmatrix, nop, canonical)
    !  Map a configuration to its lexicographically minimal form under all nop
    !  symmetry operations in eqmatrix(nop, npos).
    integer, intent(in)  :: conf(:), level, nop
    integer, intent(in)  :: eqmatrix(:, :)
    integer, intent(out) :: canonical(level)
    integer :: op, i, mapped(level)
    logical :: is_less

    canonical = conf(1:level)  ! start with identity
    do op = 1, nop
      do i = 1, level
        mapped(i) = eqmatrix(op, conf(i))
      end do
      call mc_insertion_sort(mapped, level)
      is_less = .false.
      do i = 1, level
        if (mapped(i) < canonical(i)) then
          is_less = .true.; exit
        else if (mapped(i) > canonical(i)) then
          exit
        end if
      end do
      if (is_less) canonical = mapped
    end do
  end subroutine canonicalize_config

  ! =========================================================================
  !  Sampling and sorting
  ! =========================================================================

  subroutine mc_random_subset(n, k, subset)
    !  Uniformly random sorted subset of size k from {1..n} via partial Fisher-Yates.
    integer, intent(in) :: n, k
    integer, intent(out) :: subset(:)
    integer :: pool(n), tmp, i, j
    real(real64) :: r

    do i = 1, n
      pool(i) = i
    end do
    do i = 1, k
      call random_number(r)
      j = i + int(r * real(n - i + 1, real64))
      j = min(j, n)
      tmp     = pool(i)
      pool(i) = pool(j)
      pool(j) = tmp
    end do
    subset(1:k) = pool(1:k)
    call mc_insertion_sort(subset, k)
  end subroutine mc_random_subset

  subroutine mc_insertion_sort(arr, n)
    !  In-place insertion sort for a small integer array.
    integer, intent(inout) :: arr(:)
    integer, intent(in) :: n
    integer :: i, j, key

    do i = 2, n
      key = arr(i)
      j   = i - 1
      do while (j >= 1)
        if (arr(j) <= key) exit
        arr(j+1) = arr(j)
        j = j - 1
      end do
      arr(j+1) = key
    end do
  end subroutine mc_insertion_sort

  ! =========================================================================
  !  RNG seeding
  ! =========================================================================

  subroutine seed_rng_fixed(seed_val)
    integer, intent(in) :: seed_val
    integer :: n
    integer, allocatable :: seed_arr(:)
    call random_seed(size=n)
    allocate(seed_arr(n))
    seed_arr = seed_val
    call random_seed(put=seed_arr)
  end subroutine seed_rng_fixed

  subroutine seed_rng_clock()
    integer :: i, n, clock_val
    integer, allocatable :: seed_arr(:)
    call random_seed(size=n)
    allocate(seed_arr(n))
    call system_clock(count=clock_val)
    do i = 1, n
      seed_arr(i) = clock_val + 37 * (i - 1)
    end do
    call random_seed(put=seed_arr)
  end subroutine seed_rng_clock

  ! =========================================================================
  !  Input / string helpers
  ! =========================================================================

  subroutine read_next_data_line(unit_num, out_line, io_status)
    integer, intent(in) :: unit_num
    character(len=*), intent(out) :: out_line
    integer, intent(out) :: io_status
    do
      read(unit_num, '(A)', iostat=io_status) out_line
      if (io_status /= 0) return
      out_line = adjustl(out_line)
      if (len_trim(out_line) == 0) cycle
      if (out_line(1:1) == '#') cycle
      return
    end do
  end subroutine read_next_data_line

  subroutine to_lower(str)
    character(len=*), intent(inout) :: str
    integer :: i, c
    do i = 1, len(str)
      c = iachar(str(i:i))
      if (c >= 65 .and. c <= 90) str(i:i) = achar(c + 32)
    end do
  end subroutine to_lower

end module config_sampling
