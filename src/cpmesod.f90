!*******************************************************************************
!    Copyright (c) 2018 Ricardo Grau-Crespo and co-authors
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
!*******************************************************************************

program cpmesod
  use iso_fortran_env, only: error_unit
  use cpmemod
  implicit none

  integer :: target_level, i, nargs
  character(len=256) :: arg, model_filename

  write (*, '(A)') "SOD (Site-Occupancy Disorder) version 0.91 - cpmesod"
  write (*, *) " > Constrained Periodic Motif Expansion from low/high effective Hamiltonians..."
  write (*, *) ""

  nargs = command_argument_count()
  i = 1
  do while (i <= nargs)
    call get_command_argument(i, arg)
    if (trim(arg) == '-model') then
      if (i + 1 > nargs) then
        write(error_unit,'(A)') ' Error: -model requires a filename argument.'
        stop 1
      end if
      i = i + 1
      call get_command_argument(i, model_filename)
      call cpme_set_model_filename(trim(model_filename))
      write(*,'(A,A)') ' > Using model file: ', trim(model_filename)
    end if
    i = i + 1
  end do

  call cpme_get_target_level_from_insod(target_level)
  call cpme_initialize_model(target_level)
  call cpme_print_model_summary()
  call cpme_write_level_outputs(target_level)
  call cpme_finalize_model()

  write (*, *) ""
  write (*, *) " > CPME extrapolation completed."
  write (*, *) ""

end program cpmesod
