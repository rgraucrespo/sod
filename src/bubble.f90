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

subroutine bubble(a, n)

  implicit none
  integer, parameter :: natmax = 1000
  integer :: i, j, n
  integer, dimension(1:natmax)   :: a

!
!       Make N-1 passes through the array.
!       On pass i, "bubble" the next smallest element
!       up from the end of the array to position i.
!
  do i = 1, n - 1
  do j = n, i + 1, -1
    if (a(j) < a(j - 1)) then
      call swap(a(j), a(j - 1))
    end if
  end do
  end do
  return
end

subroutine swap(i, j)
!       Exchanges the integers I and J.
  implicit none
  integer :: i, j, itemp
  itemp = i
  i = j
  j = itemp
  return
end

