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

module insod_reader
  use iso_fortran_env, only: real64
  implicit none
  private

  integer, parameter, public :: insod_ntargetmax = 5
  integer, parameter, public :: insod_nkmax      = 5

  type, public :: insod_t
    character(len=40)             :: runtitle
    real(real64)                  :: a1, b1, c1, alpha, beta, gamma
    integer                       :: nsp
    character(len=3), allocatable :: symbol(:)    ! (nsp)
    integer,          allocatable :: natsp0(:)    ! (nsp)
    integer                       :: nat0
    real(real64),     allocatable :: coords0(:,:) ! (nat0, 3)
    integer                       :: na, nb, nc
    integer                       :: ntarget
    integer                       :: sptarget(insod_ntargetmax)
    integer                       :: nk(insod_ntargetmax)
    integer                       :: nsubs_t(insod_ntargetmax, insod_nkmax)
    integer                       :: nsubs_min, nsubs_max
    character(len=5)              :: newsymbol(insod_ntargetmax, insod_nkmax+1)
    integer                       :: filer
    ! Parent molecules (optional): a parent species written as @NAME in the
    ! symbol list is a placeholder expanded to a rigid molecule (NAME.xyz) at
    ! write time. Empty (n_parentmol == 0) when no symbol carries an @ prefix.
    integer                       :: n_parentmol
    character(len=3), allocatable :: parentmol_sym(:)   ! (n_parentmol) placeholder species
    character(len=4), allocatable :: parentmol_name(:)  ! (n_parentmol) molecule name (NAME.xyz)
  end type insod_t

  public :: read_insod

contains

  subroutine read_insod(filename, d)
    character(len=*), intent(in)  :: filename
    type(insod_t),    intent(out) :: d
    integer :: unit_in, ios, at0, t, j, k, pos_colon, ip
    character(len=512) :: line
    character(len=8), allocatable :: rawsym(:)
    logical :: ok

    open(newunit=unit_in, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*, '(a,a)') 'Error: could not open INSOD file: ', trim(filename)
      stop 1
    end if

    ! title
    call next_data_line(unit_in, line, ok)
    if (.not. ok) then; write(*, *) 'Error: malformed INSOD (missing title).'; stop 1; end if
    read(line, *, iostat=ios) d%runtitle
    if (ios /= 0) d%runtitle = ''

    ! a b c alpha beta gamma
    call next_data_line(unit_in, line, ok)
    if (.not. ok) then; write(*, *) 'Error: malformed INSOD (missing cell line).'; stop 1; end if
    read(line, *, iostat=ios) d%a1, d%b1, d%c1, d%alpha, d%beta, d%gamma
    if (ios /= 0) then
      write(*, *) 'Error: could not parse cell parameters from INSOD: ', trim(line)
      stop 1
    end if

    ! nsp
    call next_data_line(unit_in, line, ok)
    if (.not. ok) then; write(*, *) 'Error: malformed INSOD (missing nsp).'; stop 1; end if
    read(line, *, iostat=ios) d%nsp
    if (ios /= 0 .or. d%nsp <= 0) then
      write(*, *) 'Error: invalid nsp entry in INSOD: ', trim(line)
      stop 1
    end if

    ! symbol(1:nsp). A symbol written as @NAME is a parent-molecule placeholder:
    ! the @ is stripped so downstream code sees the ordinary species label NAME,
    ! and NAME is recorded as a molecule to be materialised (from NAME.xyz) at
    ! write time. Read into a wide buffer first so a multi-character NAME after
    ! the @ is not truncated by the len=3 symbol field.
    allocate(d%symbol(d%nsp))
    allocate(rawsym(d%nsp))
    call next_data_line(unit_in, line, ok)
    if (.not. ok) then; write(*, *) 'Error: malformed INSOD (missing symbol list).'; stop 1; end if
    ! Guard against a stale nsp: a list-directed read of rawsym(1:nsp) would
    ! silently ignore extra symbols (the common 'nsp=2 but 3 symbols listed' bug),
    ! so cross-check the token count on the symbol line against nsp explicitly.
    if (count_tokens(line) /= d%nsp) then
      write(*, *) 'Error: INSOD nsp =', d%nsp, 'but the symbol list has', &
                  count_tokens(line), 'entries: ', trim(line)
      stop 1
    end if
    read(line, *, iostat=ios) rawsym(1:d%nsp)
    if (ios /= 0) then
      write(*, *) 'Error: could not parse species symbols from INSOD: ', trim(line)
      stop 1
    end if
    d%n_parentmol = count(rawsym(1:d%nsp)(1:1) == '@')
    if (d%n_parentmol > 0) then
      allocate(d%parentmol_sym(d%n_parentmol))
      allocate(d%parentmol_name(d%n_parentmol))
    end if
    ip = 0
    do j = 1, d%nsp
      if (rawsym(j)(1:1) == '@') then
        ip = ip + 1
        d%symbol(j)            = rawsym(j)(2:)   ! species label, @ stripped
        d%parentmol_sym(ip)    = d%symbol(j)
        d%parentmol_name(ip)   = rawsym(j)(2:)   ! molecule file NAME (NAME.xyz)
      else
        d%symbol(j) = rawsym(j)
      end if
    end do
    deallocate(rawsym)

    ! natsp0(1:nsp)
    allocate(d%natsp0(d%nsp))
    call next_data_line(unit_in, line, ok)
    if (.not. ok) then; write(*, *) 'Error: malformed INSOD (missing natsp0).'; stop 1; end if
    if (count_tokens(line) /= d%nsp) then
      write(*, *) 'Error: INSOD nsp =', d%nsp, 'but the natsp0 list has', &
                  count_tokens(line), 'entries: ', trim(line)
      stop 1
    end if
    read(line, *, iostat=ios) d%natsp0(1:d%nsp)
    if (ios /= 0) then
      write(*, *) 'Error: could not parse natsp0 from INSOD: ', trim(line)
      stop 1
    end if
    d%nat0 = sum(d%natsp0)

    ! coords0(1:nat0, 1:3)
    allocate(d%coords0(d%nat0, 3))
    do at0 = 1, d%nat0
      call next_data_line(unit_in, line, ok)
      if (.not. ok) then
        write(*, *) 'Error: malformed INSOD (missing coordinate lines).'
        stop 1
      end if
      read(line, *, iostat=ios) d%coords0(at0, 1), d%coords0(at0, 2), d%coords0(at0, 3)
      if (ios /= 0) then
        write(*, '(a,i0,a,a)') 'Error: could not parse coords0 line ', at0, ' from INSOD: ', trim(line)
        stop 1
      end if
    end do

    ! na nb nc
    call next_data_line(unit_in, line, ok)
    if (.not. ok) then; write(*, *) 'Error: malformed INSOD (missing supercell multipliers).'; stop 1; end if
    read(line, *, iostat=ios) d%na, d%nb, d%nc
    if (ios /= 0) then
      write(*, *) 'Error: could not parse supercell multipliers from INSOD: ', trim(line)
      stop 1
    end if

    ! sptarget: up to insod_ntargetmax integers on one line; detect count by trying
    call next_data_line(unit_in, line, ok)
    if (.not. ok) then; write(*, *) 'Error: malformed INSOD (missing sptarget line).'; stop 1; end if
    d%ntarget = 0
    do t = 1, insod_ntargetmax
      read(line, *, iostat=ios) d%sptarget(1:t)
      if (ios /= 0) exit
      d%ntarget = t
    end do
    if (d%ntarget == 0) then
      write(*, *) 'Error: could not parse sptarget from INSOD: ', trim(line)
      stop 1
    end if

    ! nsubs: one line per target.
    !   ntarget==1: single integer, space-separated multi-nary ints, or X:Y colon range.
    !   ntarget>=2: first line for target 1, one additional line per remaining target.
    call next_data_line(unit_in, line, ok)
    if (.not. ok) then; write(*, *) 'Error: malformed INSOD (missing nsubs line).'; stop 1; end if
    d%nk(:) = 1
    d%nsubs_t(:, :) = 0

    if (index(trim(line), '/') > 0) then
      write(*, *) "Error: '/' separator in nsubs is no longer supported."
      write(*, *) "  Use one line per target site instead (e.g. first line: 2, second line: 1)."
      stop 1
    else if (d%ntarget == 1 .and. index(trim(line), ':') > 0) then
      pos_colon = index(trim(line), ':')
      read(line(1:pos_colon-1), *, iostat=ios) d%nsubs_min
      if (ios /= 0) then
        write(*, *) 'Error: could not parse nsubs_min from colon notation: ', trim(line)
        stop 1
      end if
      read(line(pos_colon+1:), *, iostat=ios) d%nsubs_max
      if (ios /= 0) then
        write(*, *) 'Error: could not parse nsubs_max from colon notation: ', trim(line)
        stop 1
      end if
      if (d%nsubs_min < 0 .or. d%nsubs_max < 0) then
        write(*, *) 'Error: number of substitutions must be >= 0 (colon range): ', trim(line)
        stop 1
      end if
      if (d%nsubs_min > d%nsubs_max) then
        write(*, *) 'Error: nsubs_min > nsubs_max in colon range: ', trim(line)
        stop 1
      end if
    else if (d%ntarget == 1) then
      d%nk(1) = 0
      do j = 1, insod_nkmax
        read(line, *, iostat=ios) d%nsubs_t(1, 1:j)
        if (ios /= 0) exit
        d%nk(1) = j
      end do
      if (d%nk(1) == 0) then
        write(*, *) 'Error: could not parse nsubs from INSOD: ', trim(line)
        stop 1
      end if
      if (d%nk(1) > 3) then
        write(*, *) 'Error: more than 3 species per site (k=', d%nk(1), ') not yet supported.'
        stop 1
      end if
      if (d%nk(1) == 1 .and. d%nsubs_t(1,1) < 0) then
        write(*, *) 'Error: number of substitutions must be >= 0: ', trim(line)
        stop 1
      end if
      d%nsubs_min = d%nsubs_t(1,1)
      d%nsubs_max = d%nsubs_t(1,1)
    else
      ! ntarget >= 2
      d%nk(1) = 0
      do j = 1, insod_nkmax
        read(line, *, iostat=ios) d%nsubs_t(1, 1:j)
        if (ios /= 0) exit
        d%nk(1) = j
      end do
      if (d%nk(1) == 0) then
        write(*, *) 'Error: could not parse nsubs for target 1 from INSOD: ', trim(line)
        stop 1
      end if
      if (d%nk(1) > 3) then
        write(*, *) 'Error: more than 3 species per site (k=', d%nk(1), ') not yet supported.'
        stop 1
      end if
      do k = 1, d%nk(1)
        if (d%nsubs_t(1,k) < 0) then
          write(*, *) 'Error: number of substitutions must be >= 0.'
          stop 1
        end if
      end do
      do t = 2, d%ntarget
        call next_data_line(unit_in, line, ok)
        if (.not. ok) then
          write(*, '(a,i0,a)') 'Error: malformed INSOD (missing nsubs line for target ', t, ').'
          stop 1
        end if
        d%nk(t) = 0
        do j = 1, insod_nkmax
          read(line, *, iostat=ios) d%nsubs_t(t, 1:j)
          if (ios /= 0) exit
          d%nk(t) = j
        end do
        if (d%nk(t) == 0) then
          write(*, '(a,i0,a,a)') 'Error: could not parse nsubs for target ', t, &
            ' from INSOD: ', trim(line)
          stop 1
        end if
        if (d%nk(t) > 3) then
          write(*, *) 'Error: more than 3 species per site (k=', d%nk(t), ') not yet supported.'
          stop 1
        end if
        do k = 1, d%nk(t)
          if (d%nsubs_t(t,k) < 0) then
            write(*, *) 'Error: number of substitutions must be >= 0.'
            stop 1
          end if
        end do
      end do
      d%nsubs_min = d%nsubs_t(1,1)
      d%nsubs_max = d%nsubs_t(1,1)
    end if

    ! newsymbol: one line per target, nk(t)+1 tokens each
    call next_data_line(unit_in, line, ok)
    if (.not. ok) then; write(*, *) 'Error: malformed INSOD (missing newsymbol line).'; stop 1; end if
    read(line, *, iostat=ios) (d%newsymbol(1, j), j=1, d%nk(1)+1)
    if (ios /= 0) then
      write(*, *) 'Error: could not parse newsymbol for target 1 from INSOD: ', trim(line)
      stop 1
    end if
    do t = 2, d%ntarget
      call next_data_line(unit_in, line, ok)
      if (.not. ok) then
        write(*, '(a,i0,a)') 'Error: malformed INSOD (missing newsymbol line for target ', t, ').'
        stop 1
      end if
      read(line, *, iostat=ios) (d%newsymbol(t, j), j=1, d%nk(t)+1)
      if (ios /= 0) then
        write(*, '(a,i0,a,a)') 'Error: could not parse newsymbol for target ', t, &
          ' from INSOD: ', trim(line)
        stop 1
      end if
    end do

    ! filer
    call next_data_line(unit_in, line, ok)
    if (.not. ok) then; write(*, *) 'Error: malformed INSOD (missing filer line).'; stop 1; end if
    read(line, *, iostat=ios) d%filer
    if (ios /= 0) then
      write(*, *) 'Error: could not parse filer from INSOD: ', trim(line)
      stop 1
    end if

    close(unit_in)
  end subroutine read_insod

  subroutine next_data_line(unit_in, line, ok)
    integer,          intent(in)  :: unit_in
    character(len=*), intent(out) :: line
    logical,          intent(out) :: ok
    integer :: ios
    ok = .false.
    line = ''
    do
      read(unit_in, '(A)', iostat=ios) line
      if (ios /= 0) return
      line = adjustl(trim(line))
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      ok = .true.
      return
    end do
  end subroutine next_data_line

  ! Count whitespace-separated tokens in a line (used to validate list lengths).
  integer function count_tokens(line) result(n)
    character(len=*), intent(in) :: line
    integer :: i
    logical :: in_token
    n = 0
    in_token = .false.
    do i = 1, len_trim(line)
      if (line(i:i) == ' ' .or. line(i:i) == char(9)) then
        in_token = .false.
      else if (.not. in_token) then
        in_token = .true.
        n = n + 1
      end if
    end do
  end function count_tokens

end module insod_reader
