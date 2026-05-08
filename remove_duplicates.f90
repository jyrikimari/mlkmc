! MLKMC - program for removing duplicates in predictor training data
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

program remove_duplicates
  use libAtoms_module
  use ml_module
  use utils_module
  implicit none

  type(t_predictor) :: predictor,new_predictor

  real(kind=rk) :: tol,dist2
  integer :: i,j,error
  logical :: duplicate
  character(len=80) :: filename,new_filename

  call read_cmd()

  call system_initialise()
  call verbosity_push(PRINT_SILENT)

  call read_predictor_xml(predictor,filename)
  write(6,"(a,i0,a)") "Read predictor with ",predictor%n_y," data points from "//trim(filename)
  call initialise_predictor(new_predictor,predictor%n_y_min,predictor%n_y_max,predictor%n_sparse_max,&
                            predictor%desc_str,delta=predictor%delta,zeta=predictor%zeta,&
                            covariance_type=predictor%covariance_type,regularisation=predictor%regularisation)
  do i=1,predictor%n_y
    duplicate=.false.
    do j=1,new_predictor%n_y
      dist2=normsq(predictor%x(:,i)-new_predictor%x(:,j))
      if (dist2<tol) then
        write(6,"(a,i0,a,f10.6,a,f10.6,a)") "Data point ",i," has a duplicate: ",predictor%y(i)," eV vs. ",new_predictor%y(j)," eV"
        duplicate=.true.
        exit
      endif
    end do
    if (.not. duplicate) then
      call add_data_point(new_predictor,predictor%x(:,i),predictor%y(i),&
                          predictor%y_error(i),predictor%covariance_cutoff(i),error=error)
      if (error/=0) then
        if (error/=2 .and. error/=3) then
          write(0,"(a,i0,a)") "*** Error adding data point ",j," to the new predictor"
          stop
        endif
      endif
    endif
  end do

  write(6,"(a,i0,a)") "Found ",predictor%n_y-new_predictor%n_y," duplicates"
  write(6,"(a,i0,a)") "Writing new predictor with ",new_predictor%n_y," data points to file "//trim(new_filename)
  call print_predictor_xml(new_predictor,new_filename,error=error)
  if (error/=0) then
    write(0,"(a)") "*** Error writing the new predictor to file "//trim(new_filename)
    stop
  endif

contains

  subroutine read_cmd()
    integer :: iarg,ios
    character(len=80) :: arg

    iarg=command_argument_count()
    ! File, new file, tolerance
    if (iarg<3) then
      call get_command_argument(0,arg)
      write(0,"(a)") "Usage: "//trim(arg)//" file new_file tol"
      stop
    endif

    call get_command_argument(1,filename)
    call get_command_argument(2,new_filename)
    call get_command_argument(3,arg)
    read(arg,*,iostat=ios) tol
    if (ios/=0) then
      write(0,"(a)") "*** Error: invalid integer "//trim(arg)//" for tol, the tolerance for detecting duplicates"
      stop
    endif

  end subroutine read_cmd
  
end program remove_duplicates
