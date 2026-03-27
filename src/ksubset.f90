subroutine ksubset(n, k, a, more)

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
!
!
!  This subroutine generates the subsets of size K from a set of size N.
!
!  Adapted from reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the size of the set from which subsets are drawn.
!
!    Input, integer K, the desired size of the subsets.  K must
!    be between 0 and N.
!
!    Input/output, integer A(K).  A(I) is the I-th element of the
!    subset.  Thus A(I) will be an integer between 1 and N.
!    Note that the routine will return the values in A
!    in sorted order: 1 <= A(1) < A(2) < ... < A(K) <= N
!
!    Input/output, logical MORE.  Set MORE = FALSE before first call
!    for a new sequence of subsets.  It then is set and remains
!    TRUE as long as the subset computed on this call is not the
!    final one.  When the final subset is computed, MORE is set to
!    FALSE as a signal that the computation is done.
!
!*********************************************************************************

  implicit none

  integer k
  integer a(k)
  integer j
  integer, save :: m = 0
  integer, save :: m2 = 0
  logical more
  integer n

  if (k < 1 .or. n < k) then
    write (*, '(a)') ' '
    write (*, '(a)') 'KSUB_NEXT - Fatal error!'
    write (*, '(a,i6)') 'N = ', n
    write (*, '(a,i6)') 'K = ', k
    write (*, '(a)') 'but 1 <= K <= N is required!'
    stop
  end if

  if (.not. more) then
    m2 = 0
    m = k
  else
    if (m2 < n - m) then
      m = 0
    end if
    m = m + 1
    m2 = a(k + 1 - m)
  end if

  do j = 1, m
    a(k + j - m) = m2 + j
  end do

  more = a(1) /= (n - k + 1)

  return
end

