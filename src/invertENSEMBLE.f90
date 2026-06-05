!*******************************************************************************
!    Copyright (c) 2014 Ricardo Grau-Crespo, Said Hamad
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

  integer :: m, mm, npos, pos, nsubs, trash
  integer :: posinverted, kpos, kpos2
  character(len=20) :: trashstr
  character(len=200) :: ensemble_line
  integer, dimension(:), allocatable:: omega
  integer, dimension(:, :), allocatable:: conf, confinverted

  write (*, '(A)') "SOD (Site-Occupancy Disorder) version 0.80 - invertENSEMBLE"
  write (*, *) " > Inverting ENSEMBLE configurations (n -> npos-n)..."
  write (*, *) ""

  open (unit=10, file="ENSEMBLE_original")
  open (unit=11, file="ENSEMBLE_inverted")

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the original ENSEMBLE file (v2 and v3 formats)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! Read first non-blank line to detect v2 vs v3
  do
    read (10, '(a)') ensemble_line
    if (len_trim(ensemble_line) > 0) exit
  end do

  if (ensemble_line(1:1) == '#') then
    ! v2: skip comment lines; read "nsubs substitutions in npos sites"
    do while (ensemble_line(1:1) == '#')
      read (10, '(a)') ensemble_line
    end do
    read (ensemble_line, *) nsubs, trashstr
    if (trim(trashstr) /= 'substitutions') then
      write (*, *) 'ERROR: invertENSEMBLE only supports single-site binary substitutions.'
      write (*, *) '       Multi-target or multi-nary ENSEMBLE detected. Aborting.'
      stop 1
    end if
    read (ensemble_line, *) nsubs, trashstr, trashstr, npos
    read (10, *) mm
  else
    ! v3: first line is "... ensemble ...: mm configurations"
    kpos  = index(ensemble_line, ':', back=.true.)
    kpos2 = index(ensemble_line, 'configurations')
    if (kpos > 0 .and. kpos2 > kpos) then
      read (ensemble_line(kpos+1:kpos2-1), *) mm
    else
      write (*, *) 'ERROR: could not parse ENSEMBLE_original header.'
      stop 1
    end if
    ! Read target line (must be single binary target) and column-header comment
    nsubs = -1
    npos  = -1
    do
      read (10, '(a)') ensemble_line
      if (len_trim(ensemble_line) == 0) cycle
      if (ensemble_line(1:1) == '#') exit   ! column-header comment — at data rows
      if (index(ensemble_line, 'sites') > 0 .and. index(ensemble_line, '->') > 0) then
        if (nsubs >= 0) then
          write (*, *) 'ERROR: invertENSEMBLE only supports single-site binary substitutions.'
          write (*, *) '       Multi-target ENSEMBLE detected. Aborting.'
          stop 1
        end if
        read (ensemble_line, *) npos
        kpos = index(ensemble_line, '->')
        read (ensemble_line(kpos+2:), *) nsubs
      end if
    end do
    if (nsubs < 0) then
      write (*, *) 'ERROR: no target line found in ENSEMBLE_original.'
      stop 1
    end if
  end if

  ! Write inverted header (v2 format for broad compatibility)
  write (11, '(a)') "# SOD ENSEMBLE format version 2 (inverted)"
  write (11, *) npos - nsubs, "substitutions in", npos, "sites"
  write (11, *) mm, "configurations"

  allocate (conf(1:mm, 1:nsubs))
  allocate (confinverted(1:mm, 1:npos - nsubs))
  allocate (omega(1:mm))

  do m = 1, mm
    read (10, *) trash, omega(m), conf(m, 1:nsubs)
  end do

  do m = 1, mm
    posinverted = 0
    do pos = 1, npos
      if (.not. any(conf(m, :) == pos)) then
        posinverted = posinverted + 1
        confinverted(m, posinverted) = pos
      end if
    end do
  end do

  do m = 1, mm
    write (11, 10) m, omega(m), confinverted(m, 1:npos - nsubs)
10  format(i6, 1x, i6, *(1x, i4))
  end do

  deallocate (conf)
  deallocate (confinverted)
  deallocate (omega)

  write (*, *) " > Inversion completed."
  write (*, *) ""

end program invertensemble
