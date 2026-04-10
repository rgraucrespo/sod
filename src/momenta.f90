function momentaa0(nsubsmax, npos, x)
  use iso_fortran_env, only: real64
  implicit none
  integer :: nsubs, nsubsmax, npos
  real(real64) :: x, pbinomial, momentaa0
  momentaa0 = 0.0
  do nsubs = nsubsmax + 1, npos
    momentaa0 = momentaa0 + pbinomial(nsubs, npos, x)
  end do
end function momentaa0

function momentaa1(nsubsmax, npos, x)
  use iso_fortran_env, only: real64
  implicit none
  integer :: nsubs, nsubsmax, npos
  real(real64) :: x, pbinomial, momentaa1
  momentaa1 = 0.0
  do nsubs = nsubsmax + 1, npos
    momentaa1 = momentaa1 + pbinomial(nsubs, npos, x)*nsubs
  end do
end function momentaa1

function momentaa2(nsubsmax, npos, x)
  use iso_fortran_env, only: real64
  implicit none
  integer :: nsubs, nsubsmax, npos
  real(real64) :: x, pbinomial, momentaa2
  momentaa2 = 0.0
  do nsubs = nsubsmax + 1, npos
    momentaa2 = momentaa2 + pbinomial(nsubs, npos, x)*(nsubs**2)
  end do
end function momentaa2

function momentab0(nsubsmin, npos, x)
  use iso_fortran_env, only: real64
  implicit none
  integer :: nsubs, nsubsmin, npos
  real(real64) :: x, pbinomial, momentab0
  momentab0 = 0.0
  do nsubs = 0, nsubsmin - 1
    momentab0 = momentab0 + pbinomial(nsubs, npos, x)
  end do
end function momentab0

function momentab1(nsubsmin, npos, x)
  use iso_fortran_env, only: real64
  implicit none
  integer :: nsubs, nsubsmin, npos
  real(real64) :: x, pbinomial, momentab1
  momentab1 = 0.0
  do nsubs = 0, nsubsmin - 1
    momentab1 = momentab1 + pbinomial(nsubs, npos, x)*nsubs
  end do
end function momentab1

function momentab2(nsubsmin, npos, x)
  use iso_fortran_env, only: real64
  implicit none
  integer :: nsubs, nsubsmin, npos
  real(real64) :: x, pbinomial, momentab2
  momentab2 = 0.0
  do nsubs = 0, nsubsmin - 1
    momentab2 = momentab2 + pbinomial(nsubs, npos, x)*(nsubs**2)
  end do
end function momentab2

