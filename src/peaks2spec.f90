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

program peaks2spec

  implicit none

  real(kind=8), parameter :: pi = acos(-1.0d0)
  integer :: iconf, ipeak, nconfs, npeaks, npoints
  real(kind=8)  :: x, xmin, xmax, sigma, range, step, gaussian, spectrum
  integer :: i
  real(kind=8), dimension(:, :), allocatable :: peaks, spectra
  open (10, file='INP2S', status='OLD')
  open (11, file='PEAKS', status='OLD')
  open (12, file='SPECTRA')
  open (13, file='XSPEC')

!      Reading INP2S file
  read (10, *)
  read (10, *) nconfs
  read (10, *)
  read (10, *) npeaks
  read (10, *)
  read (10, *) xmin
  read (10, *)
  read (10, *) xmax
  read (10, *)
  read (10, *) npoints
  read (10, *)
  read (10, *) sigma
  close (10)

  allocate (peaks(npeaks, nconfs))
  allocate (spectra(npoints, nconfs))

!      Reading PEAKS file
  read (11, *) ((peaks(ipeak, iconf), ipeak=1, npeaks), iconf=1, nconfs)

!      Creating spectra
  range = xmax - xmin
  step = range/(npoints - 1)
  x = xmin

  do i = 1, npoints
    do iconf = 1, nconfs
      spectrum = 0.0
      do ipeak = 1, npeaks
        if (abs(x - peaks(ipeak, iconf)) < 5*sigma) then
          gaussian = (1/(sigma*sqrt(2*pi))*exp(-(x - peaks(ipeak, iconf))**2/(2*sigma**2)))
        else
          gaussian = 0.0
        end if
        spectrum = spectrum + gaussian
      end do
      spectra(i, iconf) = spectrum
    end do
    x = x + step
  end do

  write (12, *) npoints
  do iconf = 1, nconfs
    do i = 1, npoints
      write (12, '(F10.4)', advance="no") spectra(i, iconf)
    end do
    write (12, *)
  end do

  do i = 1, npoints
    write (13, '(F8.3)') xmin + step*(i - 1)
  end do

  deallocate (peaks)
  deallocate (spectra)

end program peaks2spec
