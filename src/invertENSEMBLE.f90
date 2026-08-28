!*******************************************************************************
!    Copyright (c) 2014 Ricardo Grau-Crespo and co-authors
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

program invertensemble
  implicit none

  integer :: m, mm, npos, pos, nsubs, nrem, trash, ios, sum_degen
  integer :: posinverted, kpos, kpos2
  character(len=20) :: trashstr, orig_sym, sub_sym, rem_sym
  character(len=200) :: ensemble_line
  integer, dimension(:), allocatable:: omega
  integer, dimension(:, :), allocatable:: conf, confinverted

  write (*, '(A)') "SOD (Site-Occupancy Disorder) version 0.90 - invertENSEMBLE"
  write (*, *) " > Inverting ENSEMBLE configurations (n -> npos-n)..."
  write (*, *) ""

  open (unit=10, file="ENSEMBLE_original")
  open (unit=11, file="ENSEMBLE_inverted")

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the original ENSEMBLE file (format version 3)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! Line 1 is "<type> ensemble[ (T K)]: mm configurations; ..."
  do
    read (10, '(a)') ensemble_line
    if (len_trim(ensemble_line) > 0) exit
  end do

  kpos  = index(ensemble_line, ':', back=.true.)
  kpos2 = index(ensemble_line, 'configurations')
  if (ensemble_line(1:1) == '#' .or. kpos <= 0 .or. kpos2 <= kpos) then
    write (*, *) 'ERROR: ENSEMBLE_original is not a version 3 ENSEMBLE.'
    write (*, *) '       Version 2 ENSEMBLE files are no longer supported;'
    write (*, *) '       regenerate the level with sod_comb.sh or sod_random.sh.'
    stop 1
  end if
  read (ensemble_line(kpos+1:kpos2-1), *, iostat=ios) mm
  if (ios /= 0 .or. mm < 1) then
    write (*, *) 'ERROR: could not parse the configuration count of ENSEMBLE_original.'
    stop 1
  end if

  ! Target line: "<npos> <orig_sym> sites -> <nsubs> <sub_sym> <nrem> <rem_sym>"
  nsubs = -1
  npos  = -1
  do
    read (10, '(a)') ensemble_line
    if (len_trim(ensemble_line) == 0) cycle
    if (ensemble_line(1:1) == '#') exit   ! column-header comment - at data rows
    if (index(ensemble_line, 'sites') > 0 .and. index(ensemble_line, '->') > 0) then
      if (nsubs >= 0) then
        write (*, *) 'ERROR: invertENSEMBLE only supports single-site binary substitutions.'
        write (*, *) '       Multi-target ENSEMBLE detected. Aborting.'
        stop 1
      end if
      kpos = index(ensemble_line, 'sites')
      read (ensemble_line(1:kpos-1), *, iostat=ios) npos, orig_sym
      if (ios /= 0) then
        write (*, *) 'ERROR: could not parse the target line of ENSEMBLE_original.'
        stop 1
      end if
      kpos = index(ensemble_line, '->')
      read (ensemble_line(kpos+2:), *, iostat=ios) nsubs, sub_sym, nrem, rem_sym
      if (ios /= 0) then
        write (*, *) 'ERROR: could not parse the target line of ENSEMBLE_original.'
        stop 1
      end if
      ! A fifth item means more than one substituting species on this target.
      read (ensemble_line(kpos+2:), *, iostat=ios) nsubs, sub_sym, nrem, rem_sym, trashstr
      if (ios == 0 .or. nrem /= npos - nsubs) then
        write (*, *) 'ERROR: invertENSEMBLE only supports single-site binary substitutions.'
        write (*, *) '       Multi-nary ENSEMBLE detected. Aborting.'
        stop 1
      end if
    end if
  end do
  if (nsubs < 0) then
    write (*, *) 'ERROR: no target line found in ENSEMBLE_original.'
    stop 1
  end if

  allocate (conf(1:mm, 1:nsubs))
  allocate (confinverted(1:mm, 1:npos - nsubs))
  allocate (omega(1:mm))

  do m = 1, mm
    read (10, *) trash, omega(m), conf(m, 1:nsubs)
  end do
  sum_degen = sum(omega(1:mm))

  do m = 1, mm
    posinverted = 0
    do pos = 1, npos
      if (.not. any(conf(m, :) == pos)) then
        posinverted = posinverted + 1
        confinverted(m, posinverted) = pos
      end if
    end do
  end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Writing the inverted ENSEMBLE file (format version 3)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! The substituted set becomes its complement, so the roles of the parent and
  ! the substituting species swap: nrem sites of orig_sym are now listed, and
  ! nsubs sites of sub_sym remain.
  write (11, '(A,I0,A,I0)') &
    'Enumerated ensemble: ', mm, ' configurations; sum_degeneracies = ', sum_degen
  write (11, '(I0,1X,A,A,1X,I0,1X,A,1X,I0,1X,A)') &
    npos, trim(sub_sym), ' sites ->', npos - nsubs, trim(orig_sym), nsubs, trim(sub_sym)
  write (11, '(3A)') '# Configuration  Degeneracy  ', trim(orig_sym), '_positions'

  do m = 1, mm
    write (11, 10) m, omega(m), confinverted(m, 1:npos - nsubs)
10  format(i0, 1x, i0, *(1x, i0))
  end do

  deallocate (conf)
  deallocate (confinverted)
  deallocate (omega)

  write (*, *) " > Inversion completed."
  write (*, *) ""

end program invertensemble
