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

function cc(x) &
  result(corx)
  use iso_fortran_env, only: real64
  implicit none
  real(real64), intent(in) :: x
  real(real64) :: corx
  real(real64), parameter :: tol1 = 0.0001_real64
  corx = x - floor(x)
  if (corx > (1.0_real64 - tol1)) corx = 0.0_real64
end function cc

