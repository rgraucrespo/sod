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

program invertoutsod
  implicit none

  integer :: m, mm, npos, pos, nsubs, trash
  integer :: posinverted
  character(len=20) :: trashstr
  character(len=200) :: outsod_line
  integer, dimension(:), allocatable:: omega
  integer, dimension(:, :), allocatable:: conf, confinverted

  write (*, *) "============================================================================"
  write (*, *) "         SOD (Site Occupancy Disorder) version 0.71"
  write (*, *) ""
  write (*, *) "         Authors: R. Grau-Crespo and S. Hamad"
  write (*, *) "         Contact: <r.grau-crespo@qmul.ac.uk>"
  write (*, *) "============================================================================"
  write (*, *) ""
  write (*, *) " > Inverting OUTSOD configurations (n -> npos-n)..."
  write (*, *) ""

!!!!!! Input files

  open (unit=10, file="OUTSOD_original")

!!!!!! Output files

  open (unit=11, file="OUTSOD_inverted")

!
! DEFINITION OF VARIABLES:
!
! nsubs               numbers of substitutions
! npos                Number of atoms of the target species
! pos                 Index for the positions
! conf                List of all independent configurations for the structure with nsubs substitutions
! conf1               List of all independent configurations for the structure with 1 substitutions
! conf2               List of all independent configurations for the structure with 2 substitutions
!                     All configurations are given as a list of the substituted positions
! m                   Index for the configurations (conf)
! Mn                  Total number of configurations in conf (m=1,Mm)
! m1                  Index for the configurations (conf1) for the structure with 1 defect
! Mm1                 Total number of configurations in conf1 (m1=1,Mm1) for the structure with 1 defect
! m2                  Index for the configurations (conf2) for the structure with 2 defects
! Mm2                 Total number of configurations in conf2 (m2=1,Mm2) for the structure with 2 defects

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the original OUTSOD file
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!          10  substitutions in 12 sites
!           5  configurations
!     1      6    1    2    3    4    5    6    7    8    9   10
!     2     24    1    2    3    4    5    6    7    8    9   11
!     3      6    1    2    3    4    5    7    8    9   10   11
!     4      6    1    2    3    4    5    7    8    9   10   12

  do
    read (10, '(a)') outsod_line
    if (outsod_line(1:1) /= '#') exit
  end do
  read (outsod_line, *) nsubs, trashstr
  if (trim(trashstr) /= 'substitutions') then
    write (*, *) 'ERROR: invertOUTSOD only supports single-site binary substitutions.'
    write (*, *) '       Multi-species or multi-nary OUTSOD detected. Aborting.'
    stop 1
  end if
  read (outsod_line, *) nsubs, trashstr, trashstr, npos
  write (11, '(a)') "# SOD OUTSOD format version 2"
  write (11, *) npos - nsubs, "substitutions in", npos, "sites"
  read (10, *) mm
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

!!!!!!Deallocating arrays
  deallocate (conf)
  deallocate (confinverted)
  deallocate (omega)

  write (*, *) " > Inversion completed."
  write (*, *) ""

end program invertoutsod

