!*******************************************************************************
!    mcstatsod — Thermodynamic integration over MC temperatures.
!
!    Reads E_ave(T) from MCT_TTTK/PMEx/ENSEMBLE + MCT_TTTK/PMEx/ENERGIES for each
!    temperature in ../TEMPERATURES (PMEx is read from pme.model).
!    Uses the Gibbs-Helmholtz relation d(βF)/dβ = U(β) with
!    the infinite-temperature reference F(T→∞) = −kBT ln Ω_total, where
!    Ω_total = C(npos, lev) is the total number of configurations (binomial).
!
!    βF(β) = −ln(Ω_total) + ∫₀^β U(β') dβ'
!
!    The integral is evaluated by trapezoidal rule over the sampled temperatures.
!    The tail from β=0 to the highest sampled β is handled by linear extrapolation
!    of U(β) to β=0.
!
!    Run from nXX/ (where the MCT_*K/ directories and pme.model reside).
!    Output: thermodynamics.dat  (same format as statsod)
!
!    Part of the SOD package (v0.83) — GNU GPL v3+.
!*******************************************************************************

program mcstatsod
  use iso_fortran_env, only: real64, error_unit
  use ensemble_io,     only: read_energies_file
  use pmemod,          only: pme_variant_dir_from_model
  implicit none

  integer,      parameter :: ntempmax  = 1000
  integer,      parameter :: nconfmax  = 2000000
  real(real64), parameter :: kB        = 8.617333262e-5_real64  ! eV/K

  ! --- Temperatures ---
  real(real64) :: temp_raw(ntempmax)
  integer      :: n_temp, i_temp, ios

  ! --- Per-temperature E_ave (same order as TEMPERATURES) ---
  real(real64) :: e_ave_raw(ntempmax)

  ! --- Geometry from first ENSEMBLE ---
  integer :: npos, lev, mm, mm_check, auxm, omega_m
  real(real64) :: ln_omega    ! ln(C(npos, lev)) — the infinite-T entropy / kB

  ! --- Working arrays ---
  real(real64), allocatable :: ene(:)
  integer,      allocatable :: omega(:)
  real(real64) :: sum_omega_r, sum_oe
  logical :: got_geometry

  ! --- TI (sorted arrays, decreasing T / increasing beta) ---
  integer      :: n_ti                   ! = n_temp, after filtering
  real(real64) :: t_sort(ntempmax)       ! sorted temperatures, decreasing
  real(real64) :: u_sort(ntempmax)       ! U(T) correspondingly sorted
  real(real64) :: beta_sort(ntempmax)    ! 1/(kB*T)
  real(real64) :: integral(ntempmax)     ! cumulative ∫₀^β U dβ'
  real(real64) :: f_out(ntempmax)        ! F(T)
  real(real64) :: s_out(ntempmax)        ! S(T)
  real(real64) :: u_inf                  ! extrapolated U at β=0
  real(real64) :: s_inf                  ! kB * ln_omega
  real(real64) :: slope, tail

  ! --- I/O ---
  integer :: unit_temps, unit_ensemble, unit_out
  integer :: n_missing
  logical :: ene_ok
  character(len=32)  :: txx, pme_variant
  logical            :: variant_ok
  character(len=256) :: ensemble_path, energies_path
  character(len=256) :: ensemble_line
  integer :: k, m, j, kpos_mc, kpos2_mc

  write (*, '(A)') "SOD (Site-Occupancy Disorder) version 0.83 - mcstatsod"

  ! -----------------------------------------------------------------------
  ! 0. Resolve the PME variant subdirectory (sampling method first, Hamiltonian
  !    variant second: MCT_*K/PMEx/).  pme.model lives in the main problem
  !    directory (../ from here) and is optional: when absent, mcsod builds the
  !    Hamiltonian directly from reference ENERGIES and defaults to PMEh, so we
  !    default to PMEh here too for consistency.
  ! -----------------------------------------------------------------------
  call pme_variant_dir_from_model('../pme.model', pme_variant, variant_ok)
  if (variant_ok) then
    write(*, '(A,A)') '  PME variant (from ../pme.model): ', trim(pme_variant)
  else
    pme_variant = 'PMEh'
    write(*, '(A)') '  No readable ../pme.model; defaulting to PME variant PMEh.'
  end if

  ! -----------------------------------------------------------------------
  ! 1. Read TEMPERATURES
  ! -----------------------------------------------------------------------
  open(newunit=unit_temps, file='../TEMPERATURES', status='old', action='read', iostat=ios)
  if (ios /= 0) then
    write(error_unit,'(A)') ' Error: ../TEMPERATURES file not found.'
    stop 1
  end if
  n_temp = 0
  do
    read(unit_temps, *, iostat=ios) temp_raw(n_temp + 1)
    if (ios /= 0) exit
    if (temp_raw(n_temp + 1) <= 0.0_real64) then
      write(error_unit,'(A,F12.4)') &
        ' Error: non-positive temperature in TEMPERATURES: ', temp_raw(n_temp + 1)
      stop 1
    end if
    n_temp = n_temp + 1
    if (n_temp >= ntempmax) then
      write(error_unit,'(A)') ' Error: too many temperatures (increase ntempmax).'
      stop 1
    end if
  end do
  close(unit_temps)

  if (n_temp < 2) then
    write(error_unit,'(A)') ' Error: at least 2 temperatures are required for thermodynamic integration.'
    stop 1
  end if

  write(*, '(A,I0,A)') '  Found ', n_temp, ' temperatures in ../TEMPERATURES.'

  ! -----------------------------------------------------------------------
  ! 2. Read TXX/PMEx/ENSEMBLE + TXX/PMEx/ENERGIES for each temperature → E_ave(T)
  ! -----------------------------------------------------------------------
  got_geometry = .false.
  npos = 0; lev = 0

  allocate(omega(nconfmax), ene(nconfmax))

  do i_temp = 1, n_temp
    call format_metropolis_directory(temp_raw(i_temp), txx)
    ensemble_path = trim(txx)//'/'//trim(pme_variant)//'/ENSEMBLE'
    energies_path = trim(txx)//'/'//trim(pme_variant)//'/ENERGIES'

    open(newunit=unit_ensemble, file=trim(ensemble_path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(error_unit,'(A,A)') ' Error: cannot open ', trim(ensemble_path)
      stop 1
    end if

    ! Read first non-blank line; detect v2 ('#') vs v3 ('ensemble:')
    do
      read(unit_ensemble, '(A)', iostat=ios) ensemble_line
      if (ios /= 0) then
        write(error_unit,'(A,A)') ' Error: unexpected end of file in ', trim(ensemble_path)
        stop 1
      end if
      if (len_trim(ensemble_line) > 0) exit
    end do

    if (ensemble_line(1:1) == '#') then
      ! v2: skip comment lines; first non-comment line: "lev substitutions in npos sites"
      do while (ensemble_line(1:1) == '#')
        read(unit_ensemble, '(A)', iostat=ios) ensemble_line
        if (ios /= 0) then
          write(error_unit,'(A,A)') ' Error: unexpected end of file in ', trim(ensemble_path)
          stop 1
        end if
      end do
      if (.not. got_geometry) then
        block
          character(len=64) :: w1, w2
          read(ensemble_line, *, iostat=ios) lev, w1, w2, npos
        end block
        if (ios /= 0) then
          write(error_unit,'(A)') ' Error: could not parse geometry from ENSEMBLE.'
          stop 1
        end if
        got_geometry = .true.
        ln_omega = log_gamma(real(npos + 1, real64)) &
                 - log_gamma(real(lev + 1, real64)) &
                 - log_gamma(real(npos - lev + 1, real64))
        write(*, '(A,I0,A,I0,A)') '  System: ', lev, ' substitutions in ', npos, ' sites'
        write(*, '(A,ES14.6)') '  ln(Omega_total) = ln(C(npos,lev)) = ', ln_omega
      else
        read(ensemble_line, *, iostat=ios) mm_check
        if (ios /= 0 .or. mm_check /= lev) then
          write(error_unit,'(A,A)') &
            ' Warning: substitution level mismatch in ', trim(ensemble_path)
        end if
      end if
      read(unit_ensemble, *, iostat=ios) mm
      if (ios /= 0 .or. mm <= 0 .or. mm > nconfmax) then
        write(error_unit,'(A,A,I0)') ' Error: invalid mm in ', trim(ensemble_path), mm
        stop 1
      end if
    else
      ! v3: first line is "Metropolis ensemble (T K): nic configurations"
      kpos_mc  = index(ensemble_line, ':', back=.true.)
      kpos2_mc = index(ensemble_line, 'configurations')
      if (kpos_mc > 0 .and. kpos2_mc > kpos_mc) then
        read(ensemble_line(kpos_mc+1:kpos2_mc-1), *, iostat=ios) mm
        if (ios /= 0 .or. mm <= 0 .or. mm > nconfmax) then
          write(error_unit,'(A,A)') ' Error: invalid configuration count in ', trim(ensemble_path)
          stop 1
        end if
      end if
      ! Read target line(s) and column-header comment; extract lev and npos
      do
        read(unit_ensemble, '(A)', iostat=ios) ensemble_line
        if (ios /= 0) exit
        if (len_trim(ensemble_line) == 0) cycle
        if (ensemble_line(1:1) == '#') exit  ! column-header comment — at first data row
        if (index(ensemble_line, 'sites') > 0 .and. index(ensemble_line, '->') > 0) then
          if (.not. got_geometry) then
            read(ensemble_line, *, iostat=ios) npos
            kpos_mc = index(ensemble_line, '->')
            read(ensemble_line(kpos_mc+2:), *, iostat=ios) lev
            got_geometry = .true.
            ln_omega = log_gamma(real(npos + 1, real64)) &
                     - log_gamma(real(lev + 1, real64)) &
                     - log_gamma(real(npos - lev + 1, real64))
            write(*, '(A,I0,A,I0,A)') '  System: ', lev, ' substitutions in ', npos, ' sites'
            write(*, '(A,ES14.6)') '  ln(Omega_total) = ln(C(npos,lev)) = ', ln_omega
          else
            read(ensemble_line, *, iostat=ios) npos; read(ensemble_line(index(ensemble_line,'->')+2:), *, iostat=ios) mm_check
            if (mm_check /= lev) write(error_unit,'(A,A)') &
              ' Warning: substitution level mismatch in ', trim(ensemble_path)
          end if
        end if
      end do
    end if

    ! Read omega values (index, omega, site_indices — skip site indices via list-directed I/O)
    do m = 1, mm
      read(unit_ensemble, *, iostat=ios) auxm, omega_m
      if (ios /= 0) then
        write(error_unit,'(A,I0,A,A)') ' Error: could not read config ', m, ' from ', trim(ensemble_path)
        stop 1
      end if
      omega(m) = omega_m
    end do
    close(unit_ensemble)

    ! Read energies
    call read_energies_file(trim(energies_path), mm, ene, ene_ok, n_missing)
    if (.not. ene_ok) then
      write(error_unit,'(A,A)') ' Error: could not open or read ', trim(energies_path)
      stop 1
    end if
    if (n_missing > 0) then
      write(error_unit,'(A,I0,A,A)') ' Error: missing energies for ', n_missing, &
        ' configuration(s) in ', trim(energies_path)
      stop 1
    end if

    ! Compute E_ave = sum(omega * E) / sum(omega)
    sum_omega_r = 0.0_real64
    sum_oe      = 0.0_real64
    do m = 1, mm
      sum_omega_r = sum_omega_r + real(omega(m), real64)
      sum_oe      = sum_oe      + real(omega(m), real64) * ene(m)
    end do
    if (sum_omega_r <= 0.0_real64) then
      write(error_unit,'(A,A)') ' Error: zero total omega in ', trim(ensemble_path)
      stop 1
    end if
    e_ave_raw(i_temp) = sum_oe / sum_omega_r

    write(*, '(A,A,A,F10.2,A,F16.8,A)') &
      '  ', trim(txx), ':  T = ', temp_raw(i_temp), ' K,  E_ave = ', e_ave_raw(i_temp), ' eV'
  end do

  deallocate(omega, ene)
  write(*, *)

  ! -----------------------------------------------------------------------
  ! 3. Sort temperatures descending (ascending beta) for integration
  ! -----------------------------------------------------------------------
  ! Copy to working arrays
  n_ti = n_temp
  do k = 1, n_ti
    t_sort(k) = temp_raw(k)
    u_sort(k) = e_ave_raw(k)
  end do

  ! Insertion sort: decreasing T
  do k = 2, n_ti
    j = k
    do while (j > 1 .and. t_sort(j) > t_sort(j-1))
      call swap_r(t_sort(j), t_sort(j-1))
      call swap_r(u_sort(j), u_sort(j-1))
      j = j - 1
    end do
  end do

  ! Compute beta values
  do k = 1, n_ti
    beta_sort(k) = 1.0_real64 / (kB * t_sort(k))
  end do

  ! Sanity: beta must be strictly increasing
  do k = 2, n_ti
    if (beta_sort(k) <= beta_sort(k-1)) then
      write(error_unit,'(A,I0,A,F10.2,A,F10.2)') &
        ' Warning: non-monotonic temperatures after sorting at k=', k, &
        ': T=', t_sort(k), ' and T=', t_sort(k-1)
    end if
  end do

  ! -----------------------------------------------------------------------
  ! 4. Thermodynamic integration
  !    d(βF)/dβ = U(β)  →  βF(β) = -ln(Ω) + ∫₀^β U(β') dβ'
  ! -----------------------------------------------------------------------

  ! High-T extrapolation of U to β=0 using linear fit from two highest-T points
  if (n_ti >= 2) then
    slope = (u_sort(2) - u_sort(1)) / (beta_sort(2) - beta_sort(1))
    u_inf = u_sort(1) - slope * beta_sort(1)
  else
    u_inf = u_sort(1)
  end if

  s_inf = kB * ln_omega

  write(*, '(A)') ' --- Thermodynamic integration parameters -----------------------------------'
  write(*, '(A,ES14.6)') '  ln(Omega_total)   = ', ln_omega
  write(*, '(A,ES14.6)') '  S(T=inf)/kB       = ', ln_omega
  write(*, '(A,F16.8,A)') '  U(T→∞) extrap.   = ', u_inf, ' eV  (linear extrapolation)'
  write(*, '(A,F10.2,A)') '  Highest sampled T = ', t_sort(1), ' K'
  write(*, '(A,F10.2,A)') '  Lowest  sampled T = ', t_sort(n_ti), ' K'
  if (n_ti < 4) then
    write(*, '(A)') '  Warning: fewer than 4 temperatures; TI accuracy may be limited.'
  end if
  write(*, '(A)') ' ---------------------------------------------------------------------------'
  write(*, *)

  ! Tail integral: trapezoidal from (beta=0, U_inf) to (beta_sort(1), u_sort(1))
  tail = 0.5_real64 * (u_inf + u_sort(1)) * beta_sort(1)
  integral(1) = tail

  ! Cumulative trapezoidal sum from k=1 outward
  do k = 2, n_ti
    integral(k) = integral(k-1) + &
      0.5_real64 * (u_sort(k-1) + u_sort(k)) * (beta_sort(k) - beta_sort(k-1))
  end do

  ! F and S at each sampled temperature
  do k = 1, n_ti
    f_out(k) = (-ln_omega + integral(k)) / beta_sort(k)   ! kB*T * (...)
    s_out(k) = (u_sort(k) - f_out(k)) / t_sort(k)
  end do

  ! -----------------------------------------------------------------------
  ! 5. Write thermodynamics.dat
  ! -----------------------------------------------------------------------
  open(newunit=unit_out, file='thermodynamics.dat', status='replace', action='write')

  write(unit_out,'(A)') &
    '# mcstatsod: thermodynamic integration over MC temperatures'
  write(unit_out,'(A,I0,A,I0,A,ES12.4)') &
    '# System: ', lev, ' subs in ', npos, ' sites;  ln(Omega) = ', ln_omega
  write(unit_out,'(A)') &
    '# F and S from Gibbs-Helmholtz TI with T=inf reference; U from MC trajectory averages'
  write(unit_out,'(A)') &
    '      T/K            E/eV            F/eV         S/(eV/K)'

  ! Write rows from lowest T to highest T (reverse of sort order) to match statsod convention
  do k = n_ti, 1, -1
    write(unit_out, '(F10.1,2X,2(F14.6,2X),ES12.6)') &
      t_sort(k), u_sort(k), f_out(k), s_out(k)
  end do

  ! Infinite-temperature limit
  write(unit_out, '(A10,2X,F14.6,8X,A3,7X,ES12.6)') &
    'Infinite  ', u_inf, ' - ', s_inf

  close(unit_out)

  write(*, '(A)') '  Written: thermodynamics.dat'
  write(*, *)
  write(*, '(A)') ' ============================================================================'
  write(*, '(A)') '  mcstatsod completed.'
  write(*, '(A)') ' ============================================================================'
  write(*, *)

contains

  subroutine format_metropolis_directory(temperature, dirname)
    real(real64), intent(in) :: temperature
    character(len=*), intent(out) :: dirname
    integer :: temp_integer

    temp_integer = int(temperature)
    write(dirname, '(A,I0,A)') 'MCT_', temp_integer, 'K'
  end subroutine format_metropolis_directory

  subroutine swap_r(a, b)
    real(real64), intent(inout) :: a, b
    real(real64) :: tmp
    tmp = a; a = b; b = tmp
  end subroutine swap_r

end program mcstatsod
