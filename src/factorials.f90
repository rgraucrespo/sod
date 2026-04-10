function pbinomial(nsubs, npos, x)
  use iso_fortran_env, only: int64, real64
  implicit none
  integer :: nsubs, npos
  integer(int64) :: combinations
  real(real64) :: x, pbinomial
  pbinomial = combinations(nsubs, npos)*(x**nsubs)*((1 - x)**(npos - nsubs))
end function pbinomial

function combinations(nsubs, npos)
  use iso_fortran_env, only: int64, real64
  implicit none
  integer :: nsubs, npos, k, p
  integer(int64)::  factorial, combinations, ratio

  if (nsubs <= npos/2) then
    k = nsubs
  else
    k = npos - nsubs
  end if
  ratio = 1
  do p = npos - k + 1, npos
    ratio = ratio*p
  end do
  combinations = ratio/factorial(k)
end function combinations

recursive function factorial(n) result(aux)
  use iso_fortran_env, only: int64, real64
  implicit none
  integer, intent(in) :: n
  integer(int64):: aux
  if (n == 0) then
    aux = 1
  else
    aux = n*factorial(n - 1)
  end if
end function

