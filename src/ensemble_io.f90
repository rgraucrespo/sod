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

module ensemble_io
  use iso_fortran_env, only: real64
  implicit none
  private
  public :: write_ensemble, read_ensemble, read_energies_file
  integer, parameter :: ensemble_line_len = 65536

contains

  ! Write ENSEMBLE content to an already-opened unit (format version 3).
  !
  ! source    : 'enumerated', 'uniform', or 'metropolis'
  ! tsampling : temperature in K (metropolis only)
  ! orig_sym  : original parent species per target, shape (ntarget)
  ! symbols   : new + remaining species per target, shape (ntarget, max_nk+1)
  !             symbols(t, 1:nk(t)) = new species; symbols(t, nk(t)+1) = remaining
  ! indconf   : LOCAL site indices (1..npos_t(t)); converted to global on write.
  ! nsubs_flat: flat list of all nsubs, length = sum(nk(1:ntarget)).
  ! atini_t(t): global index offset; global = local + atini_t(t) - 1.
  subroutine write_ensemble(unit_out, ntarget, nk, nsubs_flat, npos_t, atini_t, &
                            nic, indconf, degen, source, orig_sym, symbols, tsampling)
    integer,          intent(in)           :: unit_out, ntarget, nic
    integer,          intent(in)           :: nk(:)           ! (ntarget)
    integer,          intent(in)           :: nsubs_flat(:)   ! (sum(nk(1:ntarget)))
    integer,          intent(in)           :: npos_t(:)       ! (ntarget)
    integer,          intent(in)           :: atini_t(:)      ! (ntarget)
    integer,          intent(in)           :: indconf(:,:)    ! (nic, nsubs_tot)
    integer,          intent(in)           :: degen(:)        ! (nic)
    character(len=*), intent(in)           :: source          ! 'enumerated','uniform','metropolis'
    character(len=*), intent(in)           :: orig_sym(:)     ! (ntarget)
    character(len=*), intent(in)           :: symbols(:,:)    ! (ntarget, max_nk+1)
    real(real64),     intent(in), optional :: tsampling       ! K, metropolis only

    integer :: m, t, j, flat_off, conf_off, nsubs_t_tot, nsubs_tot, nrem, sum_nk, sum_degen
    integer, allocatable :: gconf(:)
    character(len=20) :: tstr

    sum_degen = sum(degen(1:nic))

    ! Line 1: ensemble type + nic + sum of degeneracies
    select case (trim(source))
    case ('enumerated')
      write(unit_out, '(A,I0,A,I0)') &
        'Enumerated ensemble: ', nic, ' configurations; sum_degeneracies = ', sum_degen
    case ('uniform')
      write(unit_out, '(A,I0,A,I0)') &
        'Uniform random ensemble: ', nic, ' configurations; sum_degeneracies = ', sum_degen
    case ('metropolis')
      if (present(tsampling)) then
        write(tstr, '(F0.1)') tsampling
        write(unit_out, '(A,A,A,I0,A,I0)') &
          'Metropolis ensemble (', trim(adjustl(tstr)), ' K): ', nic, &
          ' configurations; sum_degeneracies = ', sum_degen
      else
        write(unit_out, '(A,I0,A,I0)') &
          'Metropolis ensemble: ', nic, ' configurations; sum_degeneracies = ', sum_degen
      end if
    end select

    ! Target lines: one per target
    flat_off = 0
    do t = 1, ntarget
      nrem = npos_t(t) - sum(nsubs_flat(flat_off+1:flat_off+nk(t)))
      write(unit_out, '(I0,1X,A,A)', advance='no') &
        npos_t(t), trim(orig_sym(t)), ' sites ->'
      do j = 1, nk(t)
        write(unit_out, '(1X,I0,1X,A)', advance='no') &
          nsubs_flat(flat_off+j), trim(symbols(t,j))
      end do
      write(unit_out, '(1X,I0,1X,A)') nrem, trim(symbols(t, nk(t)+1))
      flat_off = flat_off + nk(t)
    end do

    ! Column-header comment (lists substituting new species only, not remaining)
    write(unit_out, '(A)', advance='no') '# Configuration  Degeneracy'
    flat_off = 0
    do t = 1, ntarget
      do j = 1, nk(t)
        write(unit_out, '(2A)', advance='no') '  ', trim(symbols(t,j)) // '_positions'
      end do
      flat_off = flat_off + nk(t)
    end do
    write(unit_out, *)

    ! Data rows
    sum_nk    = sum(nk(1:ntarget))
    nsubs_tot = sum(nsubs_flat(1:sum_nk))
    if (nsubs_tot > 0) allocate(gconf(nsubs_tot))

    do m = 1, nic
      if (nsubs_tot > 0) then
        flat_off = 0
        conf_off = 0
        do t = 1, ntarget
          nsubs_t_tot = sum(nsubs_flat(flat_off+1:flat_off+nk(t)))
          gconf(conf_off+1:conf_off+nsubs_t_tot) = &
            indconf(m, conf_off+1:conf_off+nsubs_t_tot) + atini_t(t) - 1
          flat_off = flat_off + nk(t)
          conf_off = conf_off + nsubs_t_tot
        end do
      end if
      write(unit_out, '(I0,1X,I0)', advance='no') m, degen(m)
      if (nsubs_tot > 0) write(unit_out, '(*(1X,I0))', advance='no') gconf(1:nsubs_tot)
      write(unit_out, *)
    end do

    if (allocated(gconf)) deallocate(gconf)
  end subroutine write_ensemble


  ! Read an ENSEMBLE file (supports format version 3 and backward-compat version 2).
  !
  ! Returns:
  !   ntarget_out    : number of targets
  !   nsubs_flat_out : flat list of substitution counts, length = sum(nk over targets)
  !   npos_t_out     : per-target site count, length = ntarget_out
  !   nic            : number of configurations
  !   indconf        : (nic, nsubs_tot) global indices as stored in file
  !   degen          : (nic) degeneracies
  !   tsampling      : temperature (-1 if not Metropolis)
  !   ok             : .true. on success
  subroutine read_ensemble(filename, ntarget_out, nsubs_flat_out, npos_t_out, &
                           nic, indconf, degen, tsampling, ok)
    character(len=*),              intent(in)  :: filename
    integer,                       intent(out) :: ntarget_out, nic
    integer,          allocatable, intent(out) :: nsubs_flat_out(:)
    integer,          allocatable, intent(out) :: npos_t_out(:)
    integer,          allocatable, intent(out) :: indconf(:,:)
    integer,          allocatable, intent(out) :: degen(:)
    real(real64),                  intent(out) :: tsampling
    logical,                       intent(out) :: ok

    integer :: unit_in, ios, i, aux, nsubs_tot
    character(len=ensemble_line_len) :: line

    ok        = .false.
    ntarget_out = 0
    nic         = 0
    tsampling   = -1.0_real64

    open(newunit=unit_in, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) return

    ! Read first non-blank line
    do
      read(unit_in, '(A)', iostat=ios) line
      if (ios /= 0) then; close(unit_in); return; end if
      if (len_trim(line) > 0) exit
    end do

    if (index(line, 'configurations') > 0 .and. line(1:1) /= '#') then
      ! ── Version 3 (new format: "... ensemble ...: N configurations") ────────
      call read_v3(unit_in, line)
    else
      ! ── Version 2 (old format: # comments then "N substitutions in M sites")
      ! Works for both v2 with and without leading # header.
      call read_v2(unit_in, line)
    end if

    close(unit_in)

  contains

    ! ── Version 2 parser (kept for backward compatibility) ────────────────────
    subroutine read_v2(u, first_line)
      integer,          intent(in) :: u
      character(len=*), intent(in) :: first_line

      integer :: kpos, kpos2, kpos3, n_nsubs, n_npos
      integer :: nsubs_buf(25), npos_buf(5)
      logical :: rows_ok
      character(len=ensemble_line_len) :: ln, rest

      ln = first_line

      ! Skip remaining comment/blank lines; detect optional Metropolis temperature
      do while (ln(1:1) == '#')
        kpos = index(ln, 'Sampling temperature')
        if (kpos > 0) then
          kpos2 = index(ln, ':')
          if (kpos2 > 0) read(ln(kpos2+1:), *, iostat=ios) tsampling
        end if
        read(u, '(A)', iostat=ios) ln
        if (ios /= 0) return
      end do

      ! Parse "N1 [N2 ...] substitutions in M1 [M2 ...] sites"
      kpos = index(ln, 'substitutions')
      if (kpos <= 1) return
      call count_and_read_ints(ln(1:kpos-1), nsubs_buf, n_nsubs)
      if (n_nsubs == 0) return

      rest = ln(kpos+13:)
      kpos2 = index(rest, ' in ')
      if (kpos2 <= 0) return
      rest = rest(kpos2+4:)
      kpos3 = index(rest, 'sites')
      if (kpos3 <= 0) return
      call count_and_read_ints(rest(1:kpos3-1), npos_buf, n_npos)
      if (n_npos == 0) return
      ntarget_out = n_npos

      allocate(nsubs_flat_out(n_nsubs))
      nsubs_flat_out = nsubs_buf(1:n_nsubs)
      allocate(npos_t_out(n_npos))
      npos_t_out = npos_buf(1:n_npos)

      read(u, *, iostat=ios) nic
      if (ios /= 0 .or. nic <= 0) then
        call cleanup_read_arrays()
        return
      end if

      call read_data_rows(u, rows_ok)
      if (.not. rows_ok) return
      ok = .true.
    end subroutine read_v2

    ! ── Version 3 parser ──────────────────────────────────────────────────────
    subroutine read_v3(u, first_line)
      integer,          intent(in) :: u
      character(len=*), intent(in) :: first_line

      integer :: kpos, kpos2, n_vals, n_nsubs_total, t_count
      integer :: vals_buf(30), npos_buf(5)
      integer :: nsubs_acc(25)
      integer :: npos_acc(5)
      real(real64) :: tval
      logical :: rows_ok
      character(len=ensemble_line_len) :: ln, after_arrow

      ! Parse line 1: ensemble type + nic
      ln = first_line
      if (index(ln, 'Metropolis') > 0) then
        ! Extract temperature from "Metropolis ensemble (T K):"
        kpos  = index(ln, '(')
        kpos2 = index(ln, 'K)')
        if (kpos > 0 .and. kpos2 > kpos) then
          read(ln(kpos+1:kpos2-1), *, iostat=ios) tval
          if (ios == 0) tsampling = tval
        end if
      end if
      ! Extract nic from "<nic> configurations" — find last integer before 'configurations'
      kpos = index(ln, 'configurations')
      if (kpos <= 1) return
      block
        integer :: tmp_buf(20), tmp_n
        call extract_ints_from_string(ln(1:kpos-1), tmp_buf, tmp_n)
        if (tmp_n == 0) return
        nic = tmp_buf(tmp_n)
      end block
      if (nic <= 0) return

      ! Read target lines and column-header comment; collect npos_t and nsubs per target
      t_count     = 0
      n_nsubs_total = 0

      do
        read(u, '(A)', iostat=ios) ln
        if (ios /= 0) return
        if (len_trim(ln) == 0) cycle
        if (ln(1:1) == '#') exit          ! column-header comment → data rows next
        if (index(ln, 'sites') == 0) exit ! not a target line → data rows

        ! Parse "<npos_t> <orig_sym> sites -> <n1> <sym1> ... <nrem> <rem_sym>"
        kpos = index(ln, 'sites')
        if (kpos == 0) exit

        ! npos_t is the first integer on the line
        call count_and_read_ints(ln(1:kpos-1), npos_buf, n_vals)
        if (n_vals == 0) return
        t_count = t_count + 1
        npos_acc(t_count) = npos_buf(1)

        ! Extract integers after "->"
        kpos2 = index(ln, '->')
        if (kpos2 == 0) return
        after_arrow = ln(kpos2+2:)
        call extract_ints_from_string(after_arrow, vals_buf, n_vals)
        ! n_vals integers: first (n_vals-1) are nsubs values, last is nrem (not stored)
        if (n_vals < 1) return
        do i = 1, n_vals - 1
          n_nsubs_total = n_nsubs_total + 1
          nsubs_acc(n_nsubs_total) = vals_buf(i)
        end do
      end do
      ! ln now holds the first data row (or was a '#' comment line, in which case
      ! we need to read the first data row separately)

      if (t_count == 0) return
      ntarget_out = t_count

      allocate(npos_t_out(ntarget_out))
      npos_t_out = npos_acc(1:ntarget_out)
      allocate(nsubs_flat_out(n_nsubs_total))
      nsubs_flat_out = nsubs_acc(1:n_nsubs_total)

      ! If we exited the target loop on a '#' comment, read one more line for first data row
      if (ln(1:1) == '#') then
        read(u, '(A)', iostat=ios) ln
        if (ios /= 0) then
          call cleanup_read_arrays()
          return
        end if
      end if

      call read_data_rows_from_line(u, ln, rows_ok)
      if (.not. rows_ok) return
      ok = .true.
    end subroutine read_v3

    ! Read nic data rows; first row is already in first_line
    subroutine read_data_rows_from_line(u, first_line, rows_ok)
      integer,          intent(in) :: u
      character(len=*), intent(in) :: first_line
      logical,          intent(out) :: rows_ok
      integer :: row
      character(len=ensemble_line_len) :: ln

      rows_ok = .false.
      nsubs_tot = sum(nsubs_flat_out)
      allocate(indconf(nic, max(1, nsubs_tot)))
      allocate(degen(nic))
      indconf = 0; degen = 0

      ln = first_line
      do row = 1, nic
        if (nsubs_tot > 0) then
          read(ln, *, iostat=ios) aux, degen(row), indconf(row, 1:nsubs_tot)
        else
          read(ln, *, iostat=ios) aux, degen(row)
        end if
        if (ios /= 0) then
          call cleanup_read_arrays()
          return
        end if
        if (row < nic) then
          read(u, '(A)', iostat=ios) ln
          if (ios /= 0) then
            call cleanup_read_arrays()
            return
          end if
        end if
      end do
      rows_ok = .true.
    end subroutine read_data_rows_from_line

    ! Read nic data rows from unit u (used by v2 path)
    subroutine read_data_rows(u, rows_ok)
      integer, intent(in) :: u
      logical, intent(out) :: rows_ok
      integer :: row

      rows_ok = .false.
      nsubs_tot = sum(nsubs_flat_out)
      allocate(indconf(nic, max(1, nsubs_tot)))
      allocate(degen(nic))
      indconf = 0; degen = 0

      do row = 1, nic
        if (nsubs_tot > 0) then
          read(u, *, iostat=ios) aux, degen(row), indconf(row, 1:nsubs_tot)
        else
          read(u, *, iostat=ios) aux, degen(row)
        end if
        if (ios /= 0) then
          call cleanup_read_arrays()
          return
        end if
      end do
      rows_ok = .true.
    end subroutine read_data_rows

    subroutine cleanup_read_arrays()
      if (allocated(indconf))        deallocate(indconf)
      if (allocated(degen))          deallocate(degen)
      if (allocated(nsubs_flat_out)) deallocate(nsubs_flat_out)
      if (allocated(npos_t_out))     deallocate(npos_t_out)
    end subroutine cleanup_read_arrays

  end subroutine read_ensemble


  ! Extract all integers from a mixed int/string token stream (e.g. "2 Mg  6 La").
  ! Non-integer tokens (species symbols) are silently skipped.
  subroutine extract_ints_from_string(str, vals, n)
    character(len=*), intent(in)  :: str
    integer,          intent(out) :: vals(:)
    integer,          intent(out) :: n

    integer :: pos, len_s, token_start, iv, ios
    character(len=64) :: token

    n     = 0
    len_s = len_trim(str)
    pos   = 1

    do while (pos <= len_s .and. n < size(vals))
      ! skip spaces
      do while (pos <= len_s .and. str(pos:pos) == ' ')
        pos = pos + 1
      end do
      if (pos > len_s) exit
      ! collect token (up to next space)
      token_start = pos
      do while (pos <= len_s .and. str(pos:pos) /= ' ')
        pos = pos + 1
      end do
      token = str(token_start:pos-1)
      ! try to parse as integer; skip if not
      read(token, *, iostat=ios) iv
      if (ios == 0) then
        n = n + 1
        vals(n) = iv
      end if
    end do
  end subroutine extract_ints_from_string


  ! Count how many integers can be read consecutively from str and store in vals(1:n).
  subroutine count_and_read_ints(str, vals, n)
    character(len=*), intent(in)  :: str
    integer,          intent(out) :: vals(:)
    integer,          intent(out) :: n
    integer :: ios
    n = 0
    do
      if (n >= size(vals)) exit
      read(str, *, iostat=ios) vals(1:n+1)
      if (ios /= 0) exit
      n = n + 1
    end do
  end subroutine count_and_read_ints

  ! Read an ENERGIES file in two-column format: "m  E_nm" per line.
  ! Comments (#) and blank lines are skipped.  Missing configs leave
  ! energies(m)=0 and are counted in n_missing.
  ! ok=.false. means the file could not be opened or a line could not be parsed.
  subroutine read_energies_file(filename, mm, energies, ok, n_missing, e_mask)
    character(len=*), intent(in)  :: filename
    integer,          intent(in)  :: mm
    real(real64),     intent(out) :: energies(mm)
    logical,          intent(out) :: ok
    integer,          intent(out) :: n_missing
    logical, optional, intent(out) :: e_mask(mm)

    integer      :: unit_in, ios, m_read
    real(real64) :: e_read
    character(len=256) :: line
    logical      :: e_present(mm)

    ok        = .false.
    n_missing = 0
    energies  = 0.0_real64
    e_present = .false.

    open(newunit=unit_in, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) return

    do
      read(unit_in, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      read(line, *, iostat=ios) m_read, e_read
      if (ios /= 0) then
        close(unit_in)
        return
      end if
      if (m_read >= 1 .and. m_read <= mm) then
        energies(m_read)  = e_read
        e_present(m_read) = .true.
      end if
    end do
    close(unit_in)

    n_missing = count(.not. e_present)
    if (present(e_mask)) e_mask = e_present
    ok        = .true.
  end subroutine read_energies_file

end module ensemble_io
