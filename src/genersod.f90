!*******************************************************************************
!    Copyright (c) 2018 Ricardo Grau-Crespo, Said Hamad
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

program genersod
  implicit none

  integer, parameter :: nspmax = 10, natmax = 10000, nlineamax = 200

  integer :: i, l, m
  integer :: ifound
  integer :: sp, nsp, sptarget, ssp
  integer :: at0, nat0, at, nat, att, filer, genxtl, genarc, ndigits
  character :: outfilename*20, fmtstr*20
  integer :: na, nb, nc, nsubs, atini, atfin, cumnatsp
  integer :: npos, nic, indcount
  integer, dimension(:), allocatable :: newconf, degen
  integer, dimension(2) :: newshell
  integer, dimension(nspmax) :: natsp0, natsp, snatsp, ishell
  integer, dimension(natmax) :: spat
  real, dimension(natmax, 3) :: coords0, coords
  real, dimension(3, 3) :: cellvector
  integer, dimension(:, :), allocatable :: indconf
  real :: a1, b1, c1, alpha, beta, gamma, a, b, c
  character, dimension(nspmax) :: symbol*3, ssymbol*3
  character, dimension(2) :: newsymbol*3
  character :: linea*85, title*15, runtitle*40, trashtext*20, symboltrash*3

! Input files

  open (unit=9, file="INSOD")
  open (unit=30, file="OUTSOD")
  open (unit=31, file="SUPERCELL")

! Output files

  open (unit=43, file="filer")

!
! DEFINITION OF VARIABLES:
!
! sp                  Index for the species
! nsp                 Total number of species
! sptarget            Number of the species to be substituted
! at0,at1,at          Indexes for the atoms in the asymmetric unit, unit cell and supercell
! nat0,nat1,nat       Total numbers of atoms in the asymmetric unit, unit cell and supercell
! at1r,nat1r            Idem for the atoms in the redundant cell (with repeated positions)
! atini,atfin           Initial and final atom indexes of the species to be substituted
! pos                 Index for atomic positions of the target species
! npos                      Number of atoms of the target species
! conf                      List of all configurations (each configuration is a list of the substituted positions)
! count                     Index for the configurations (conf)
! ntc                 Total number of configurations in conf (count=1,ntc)
! indconf             List of independent configurations
! indcount            Index for the independent configurations
! nic                 Total number of independent configurations in indconf (indcount=1,nic)
! equivconf           Temporary list containing the equivalent configurations at every step of the algorithm
! equivcount          Index for the equivalent configurations (equivconf)
! tol0                      General tolerance
! tol1                      Tolerance used for correcting the x-FLOOR(x) function
!
!
!

  write (*, *) "Reading INSOD, SUPERCELL, and OUTSOD to generate calculation input files"
  write (*, *) " "

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the INSOD file
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  read (9, *)

  read (9, *) runtitle
  read (9, *)
  read (9, *)
  read (9, *) a1, b1, c1, alpha, beta, gamma
  read (9, *)
  read (9, *)
  read (9, *) nsp
  read (9, *)
  read (9, *)
  read (9, *) (symbol(sp), sp=1, nsp)
  read (9, *)
  read (9, *)
  read (9, *) (natsp0(sp), sp=1, nsp)
  read (9, *)
  read (9, *)

  nat0 = 0
  do sp = 1, nsp
    nat0 = nat0 + natsp0(sp)
  end do

  do at0 = 1, nat0
    read (9, *) (coords0(at0, i), i=1, 3)
  end do
  read (9, *)
  read (9, *)
  read (9, *) na, nb, nc
  read (9, *)
  read (9, *)
  read (9, *) sptarget
  read (9, *)
  read (9, *)
  read (9, *) nsubs
  if (nsubs == 0) then
    write (*, *) "Illegal number of substitutions"
    stop
  end if
  read (9, *)
  read (9, *)
  read (9, *)
  read (9, *) (newsymbol(i), i=1, 2)
  read (9, *)
  read (9, *)
  read (9, *) filer
  read (9, *)

  if ((filer > 0) .and. (filer < 10)) then
    read (9, *)
    read (9, *)
    read (9, *) (ishell(sp), sp=1, nsp)
    read (9, *)
    read (9, *) (newshell(i), i=1, 2)
    read (9, *)
    read (9, *) genxtl, genarc
  end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the OUTSOD file
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  read (30, *) nsubs, trashtext, trashtext, npos
  read (30, *) nic

  ndigits = max(1, int(log10(real(nic))) + 1)
  write (fmtstr, '(a,i0,a,i0,a)') '(a,i', ndigits, '.', ndigits, ')'

!!!!!!!!Allocating array sizes

  allocate (degen(1:nic))
  allocate (newconf(1:nsubs))
  allocate (indconf(1:nic, 1:nsubs))

  do indcount = 1, nic

    read (30, *) m, degen(indcount), indconf(indcount, 1:nsubs)
    if (m /= indcount) then
      write (*, *) "Error in configuration numbering in OUTSOD. Aborting..."
      stop
    end if
  end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the SUPERCELL file
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  read (31, *) a, b, c, alpha, beta, gamma
  read (31, *) natsp(1:nsp)

  nat = sum(natsp(1:nsp))
  do at = 1, nat
    read (31, *) symboltrash, coords(at, 1), coords(at, 2), coords(at, 3)
  end do

!cccccccccccccccccccccccccccccccccccc
! Generating spat array
!cccccccccccccccccccccccccccccccccccc

  do at = 1, natsp(1)
    spat(at) = 1
  end do
  cumnatsp = natsp(1)
  do sp = 2, nsp
    do at = cumnatsp + 1, cumnatsp + natsp(sp)
      spat(at) = sp
    end do
    cumnatsp = cumnatsp + natsp(sp)
  end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       Calculate the initial and final nat of the target species
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if (sptarget == 1) then
    atini = 1
  else
    atini = 1
    do sp = 1, sptarget - 1
      atini = atini + natsp(sp)
    end do
  end if

  atfin = atini + natsp(sptarget) - 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!    GENERATE INPUT FILES !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!! Write a temporary file with FILER, to be read later by the script
  write (43, *) filer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  select case (filer)

  case (0)
    write (*, *) "Calculation files not created.&
               & change filer value in insod file if you want to create calculation files."
    write (*, *) ""

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!  GULP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case (1)

    write (*, *) " "
    write (*, *) "Creating input files for GULP in ./CALCS... "
    write (*, *) " "

    open (unit=41, file="top.gulp")
    open (unit=42, file="bottom.gulp")

    do indcount = 1, nic
      newconf(1:nsubs) = indconf(indcount, 1:nsubs)

      do l = 1, nlineamax
        read (41, 332, end=99) linea
332     format(a85)
        write (indcount + 100000, 332) linea
      end do
99    continue
      rewind (41)

      write (indcount + 100000, 211) a, b, c, alpha, beta, gamma
      write (indcount + 100000, *) "frac"
211   format(6(f10.4, 2x))
      do at = 1, nat
        sp = spat(at)
        if (sp /= sptarget) then
          write (indcount + 100000, 331) symbol(sp), "core", coords(at, 1), coords(at, 2), coords(at, 3)
          if (ishell(sp) == 1) write (indcount + 100000, 331) symbol(sp), "shel", coords(at, 1), coords(at, 2), coords(at, 3)
        else
          att = at - atini + 1
          call member(nsubs, newconf, att, ifound)
          if (ifound == 1) then
            write (indcount + 100000, 331) newsymbol(1), "core", coords(at, 1), coords(at, 2), coords(at, 3)
            if (newshell(1) == 1) write (indcount + 100000, 331) newsymbol(1), "shel", coords(at, 1), coords(at, 2), coords(at, 3)
          else
            write (indcount + 100000, 331) newsymbol(2), "core", coords(at, 1), coords(at, 2), coords(at, 3)
            if (newshell(2) == 1) write (indcount + 100000, 331) newsymbol(2), "shel", coords(at, 1), coords(at, 2), coords(at, 3)
          end if
        end if
331     format(a3, 2x, a4, 2x, 3(f11.7, 2x))
      end do

      do l = 1, nlineamax
        read (42, 332, end=199) linea
        write (indcount + 100000, 333) linea
333     format(a85)
      end do
199   continue
      rewind (42)

      if (genxtl == 1 .or. genarc == 1) then
        write (outfilename, fmtstr) 'c', indcount
        if (genxtl == 1) write (indcount + 100000, *) 'output xtl '//trim(outfilename)
        if (genarc == 1) write (indcount + 100000, *) 'output arc '//trim(outfilename)
      end if

      close (unit=indcount + 100000)

    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! METADISE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case (2)

    write (*, *) " "
    write (*, *) "Creating input files for METADISE in ./CALCS... "
    write (*, *) " "

    open (unit=51, file="top.metadise")
    open (unit=52, file="bottom.metadise")

    do indcount = 1, nic
      newconf(1:nsubs) = indconf(indcount, 1:nsubs)

      do l = 1, nlineamax
        read (51, 332, end=299) linea
        write (indcount + 100000, 332) linea
      end do
299   continue
      rewind (51)

      write (indcount + 100000, 212) "CELL ", a, b, c, alpha, beta, gamma
212   format(a5, 6(f10.4, 2x))
      write (indcount + 100000, *) "SPACE FULL P1 1 1"
      write (indcount + 100000, *) "FRAC"
      do at = 1, nat
        sp = spat(at)
        if (sp /= sptarget) then
          write (indcount + 100000, 331) symbol(sp), "CORE", coords(at, 1), coords(at, 2), coords(at, 3)
          if (ishell(sp) == 1) write (indcount + 100000, 331) symbol(sp), "SHEL", coords(at, 1), coords(at, 2), coords(at, 3)
        else
          att = at - atini + 1
          call member(nsubs, newconf, att, ifound)
          if (ifound == 1) then
            write (indcount + 100000, 331) newsymbol(1), "CORE", coords(at, 1), coords(at, 2), coords(at, 3)
            if (newshell(1) == 1) write (indcount + 100000, 331) newsymbol(1), "SHEL", coords(at, 1), coords(at, 2), coords(at, 3)
          else
            write (indcount + 100000, 331) newsymbol(2), "CORE", coords(at, 1), coords(at, 2), coords(at, 3)
            if (newshell(2) == 1) write (indcount + 100000, 331) newsymbol(2), "SHEL", coords(at, 1), coords(at, 2), coords(at, 3)
          end if
        end if
      end do
      do l = 1, nlineamax
        read (52, 332, end=399) linea
        write (indcount + 100000, 333) linea
      end do
399   continue
      rewind (52)

      close (unit=indcount + 100000)

    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!  VASP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case (11)

    write (*, *) " "
    write (*, *) "Creating input files for VASP in ./CALCS... "
    write (*, *) " "

    title = 'vasp'

    do indcount = 1, nic
      newconf(1:nsubs) = indconf(indcount, 1:nsubs)

      write (indcount + 100000, *) title
      write (indcount + 100000, *) '1.00000000'
      call cell(cellvector, a, b, c, alpha, beta, gamma)
      write (indcount + 100000, 335) cellvector(1, 1), cellvector(2, 1), cellvector(3, 1)
      write (indcount + 100000, 335) cellvector(1, 2), cellvector(2, 2), cellvector(3, 2)
      write (indcount + 100000, 335) cellvector(1, 3), cellvector(2, 3), cellvector(3, 3)
335   format(3(f10.6, 2x))

      do ssp = 1, nsp
        if (ssp < sptarget) then
          snatsp(ssp) = natsp(ssp)
          ssymbol(ssp) = symbol(ssp)
        end if
        if (ssp == sptarget) then
          snatsp(ssp) = nsubs
          snatsp(ssp + 1) = npos - nsubs
          ssymbol(ssp) = newsymbol(1)
          ssymbol(ssp + 1) = newsymbol(2)
        end if
        if (ssp > sptarget) then
          snatsp(ssp + 1) = natsp(ssp)
          ssymbol(ssp + 1) = symbol(ssp)
        end if
      end do

      write (indcount + 100000, *) ssymbol(1:nsp + 1)
      write (indcount + 100000, 336) snatsp(1:nsp + 1)
336   format(10(i4, 1x))
      write (indcount + 100000, *) 'Direct'
      do at = 1, atini - 1
        sp = spat(at)
        write (indcount + 100000, 335) coords(at, 1), coords(at, 2), coords(at, 3)
      end do
      do at = atini, atfin
        sp = spat(at)
        att = at - atini + 1
        call member(nsubs, newconf, att, ifound)
        if (ifound == 1) then
          write (indcount + 100000, 335) coords(at, 1), coords(at, 2), coords(at, 3)
        end if
      end do
      do at = atini, atfin
        sp = spat(at)
        att = at - atini + 1
        call member(nsubs, newconf, att, ifound)
        if (ifound == 0) then
          write (indcount + 100000, 335) coords(at, 1), coords(at, 2), coords(at, 3)
        end if
      end do
      do at = atfin + 1, nat
        sp = spat(at)
        write (indcount + 100000, 335) coords(at, 1), coords(at, 2), coords(at, 3)
      end do

      close (indcount + 100000)

    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!  CASTEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case (12)

    write (*, *) " "
    write (*, *) "Creating input files for CASTEP in ./CALCS... "
    write (*, *) " "

    open (unit=62, file="bottom.castep")

    do indcount = 1, nic
      newconf(1:nsubs) = indconf(indcount, 1:nsubs)

      write (indcount + 100000, '(a19)') "%BLOCK lattice_cart"

      call cell(cellvector, a, b, c, alpha, beta, gamma)
      write (indcount + 100000, 335) cellvector(1, 1), cellvector(2, 1), cellvector(3, 1)
      write (indcount + 100000, 335) cellvector(1, 2), cellvector(2, 2), cellvector(3, 2)
      write (indcount + 100000, 335) cellvector(1, 3), cellvector(2, 3), cellvector(3, 3)

      write (indcount + 100000, '(a22)') "%ENDBLOCK lattice_cart"
      write (indcount + 100000, '(a21)') "%BLOCK positions_frac"

      do at = 1, nat
        sp = spat(at)
        if (sp /= sptarget) then
          write (indcount + 100000, 337) symbol(sp), coords(at, 1), coords(at, 2), coords(at, 3)
        else
          att = at - atini + 1
          call member(nsubs, newconf, att, ifound)
          if (ifound == 1) then
            write (indcount + 100000, 337) newsymbol(1), coords(at, 1), coords(at, 2), coords(at, 3)
          else
            write (indcount + 100000, 337) newsymbol(2), coords(at, 1), coords(at, 2), coords(at, 3)
          end if
        end if
337     format(a3, 3(f11.7, 2x))
      end do

      write (indcount + 100000, '(a24)') "%ENDBLOCK positions_frac"

      if (indcount > 99999) then
        write (indcount + 100000, *) "Error, too many configurations (>99999)! Calculation files not written!"
      end if

      do l = 1, nlineamax
        read (62, 332, end=499) linea
        write (indcount + 100000, 333) linea
      end do
499   continue
      rewind (62)

      if (indcount > 99999) then
        write (indcount + 100000, *) "Error, too many configurations (>99999)! Calculation files not written!"
      end if

      close (unit=indcount + 100000)

    end do

    close (unit=62)

  end select

!!!!!!!Deallocating arrays
  deallocate (newconf)
  deallocate (degen)
  deallocate (indconf)

  close (43)

!!!!!!!Reporting the end
  write (*, *) "Done!!!"
  write (*, *) ""
  write (*, *) ""

end program genersod

