! MLKMC - program for converting GAP predictor data to FANN format
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

program convert_to_FANN_format
  use libAtoms_module
  use ml_module
  use utils_module
  implicit none

  type(t_predictor) :: predictor

  integer :: i,ios
  character(len=80) :: filename,new_filename

  call read_cmd()

  call system_initialise()
  call verbosity_push(PRINT_SILENT)

  call read_predictor_xml(predictor,filename)
  write(6,"(a,i0,a)") "Read predictor with ",predictor%n_y," data points from "//trim(filename)
  open(unit=10,file=new_filename,action="write",iostat=ios)
  if (ios/=0) then
    write(0,"(a)") "*** Error: can't open file",trim(new_filename)
    stop
  endif
  write(10,"(i0,x,i0,x,i0)",iostat=ios) predictor%n_y,predictor%d,1
  if (ios/=0) then
    write(0,"(a)") "*** Error: can't write header to file",trim(new_filename)
    stop
  endif
  do i=1,predictor%n_y
    write(10,*,iostat=ios) predictor%x(:,i)
    if (ios/=0) then
      write(0,"(a)") "*** Error: can't write x data of point ",i," to file",trim(new_filename)
      stop
    endif
    write(10,*,iostat=ios) predictor%y(i)
    if (ios/=0) then
      write(0,"(a)") "*** Error: can't write y data of point ",i," to file",trim(new_filename)
      stop
    endif
  end do
  close(10)

contains

  subroutine read_cmd()
    integer :: iarg
    character(len=80) :: arg

    iarg=command_argument_count()
    ! File, new file
    if (iarg<2) then
      call get_command_argument(0,arg)
      write(0,"(a)") "Usage: "//trim(arg)//" file new_file"
      stop
    endif

    call get_command_argument(1,filename)
    call get_command_argument(2,new_filename)

  end subroutine read_cmd
  
end program convert_to_FANN_format
