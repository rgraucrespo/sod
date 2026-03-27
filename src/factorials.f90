function pbinomial(nsubs, npos, x)
  implicit none
  integer :: nsubs, npos
  integer(kind=8) :: combinations
  real*8 :: x, pbinomial
  pbinomial = combinations(nsubs, npos)*(x**nsubs)*((1 - x)**(npos - nsubs))
end function pbinomial

function combinations(nsubs, npos)
  implicit none
  integer :: nsubs, npos, k, p
  integer(kind=8)::  factorial, combinations, ratio

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
  implicit none
  integer, intent(in) :: n
  integer(kind=8):: aux
  if (n == 0) then
    aux = 1
  else
    aux = n*factorial(n - 1)
  end if
end function

