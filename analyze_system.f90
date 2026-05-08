! MLKMC - program for analyzing and debugging a .system file
!
! Copyright (C) 2026 Jyri Kimari
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

program analyze_system
  use mlkmc_atoms_module

  type(t_system) :: system
  integer :: error,fileunit
  logical :: exist
  character(len=80) :: filename

  call read_cmd()

  inquire(file=filename,exist=exist)
  if (exist) then
    write(6,"(a)") "# Reading "//trim(filename)//"..."
    call read_system(system,filename,error=error)
    if (error/=0) then
      write(0,"(a)") "*** Error reading system file"
      stop
    endif
  else
    write(0,"(a)") "*** Error: filename "//trim(filename)//" not found"
    stop
  endif

  call check_validity(system,error=error)
  if (error==0) then
    write(6,"(a)") "# System seems valid!"
  else
    write(6,"(a,i0,a)") "# ",error," error(s) detected in system"
  endif
  open(newunit=fileunit,file="validate.xyz",action="write")
  call print_xyz_frame(system,0.0d0,0,fileunit,hard_debug=.true.)
  close(fileunit)

contains

  subroutine read_cmd()
    integer :: iarg
    character(len=80) :: arg

    iarg=command_argument_count()
    ! File
    if (iarg<1) then
      call get_command_argument(0,arg)
      write(0,"(a)") "Usage: "//trim(arg)//" file"
      stop
    endif

    call get_command_argument(1,filename)

  end subroutine read_cmd

end program analyze_system
