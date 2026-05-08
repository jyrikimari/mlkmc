! MLKMC - program for refitting a predictor with new parameters
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

program refit_ml
  use libAtoms_module
  use ml_module
  implicit none

  type(t_predictor) :: predictor

  integer :: i
  integer :: n_y_min
  real(kind=rk) :: y_error,delta,zeta,regularisation,covariance_cutoff
  character(len=80) :: filename,new_filename
  integer :: error

  call read_cmd()

  call system_initialise()
  call verbosity_push(PRINT_SILENT)

  write(6,"(a)") "# Reading predictor at "//trim(filename)
  call read_predictor_xml(predictor,filename)
  if (n_y_min>=0) then
    predictor%n_y_min=n_y_min
  endif
  if (y_error>=0) then
    predictor%y_error=y_error
  endif
  if (delta>=0) then
    predictor%delta=delta
  endif
  if (zeta>=0) then
    predictor%zeta=zeta
  endif
  if (regularisation>=0) then
    predictor%regularisation=regularisation
  endif
  if (covariance_cutoff>=0) then
    do i=1,predictor%n_y
      predictor%covariance_cutoff(i)=covariance_cutoff
    end do
  endif
  predictor%fitted=.false.
  if (predictor%n_y>predictor%n_y_min) then
    write(6,"(a,i0,a,i0,a)") "# Predictor n_y=",predictor%n_y," > n_y_min=",predictor%n_y_min,"; fitting"
    call fit_predictor(predictor,error=error)
    if (error/=0) then
      write(0,"(a)") "*** Error fitting predictor"
      stop
    endif
  else
    predictor%fitted=.false.
  endif
  predictor%changed=.true.
  write(6,"(a)") "# Writing predictor at "//trim(new_filename)
  call print_predictor_xml(predictor,new_filename,error=error)
  if (error/=0) then
    write(0,"(a)") "*** Error printing predictor XML to file "//trim(new_filename)
  endif

contains

  subroutine read_cmd()
    integer :: iarg,ios
    character(len=80) :: arg

    iarg=command_argument_count()
    ! File, new file, y_error, delta, zeta, regularisation, covariance cutoff
    if (iarg<2) then
      call get_command_argument(0,arg)
      write(0,"(a)") "Usage: "//trim(arg)//" file new_file [n_y_min y_error delta zeta regularisation covariance_cutoff]"
      stop
    endif

    call get_command_argument(1,filename)
    call get_command_argument(2,new_filename)

    n_y_min=-1
    if (iarg>2) then
      call get_command_argument(3,arg)
      read(arg,*,iostat=ios) n_y_min
      if (ios/=0) then
        write(0,"(a)") "*** Error: invalid integer "//trim(arg)//" for n_y_min"
        stop
      endif
      if (n_y_min<0) then
        write(0,"(a,i0,a)") "*** Error: n_y_min ",n_y_min,"<0"
        stop
      endif
    endif

    y_error=-1
    if (iarg>3) then
      call get_command_argument(5,arg)
      read(arg,*,iostat=ios) y_error
      if (ios/=0) then
        write(0,"(a)") "*** Error: invalid real "//trim(arg)//" for y_error"
        stop
      endif
      if (y_error<0) then
        write(0,"(a,f10.3,a)") "*** Error: y_error ",y_error,"<0"
        stop
      endif
    endif
  
    delta=-1
    if (iarg>4) then
      call get_command_argument(6,arg)
      read(arg,*,iostat=ios) delta
      if (ios/=0) then
        write(0,"(a)") "*** Error: invalid real "//trim(arg)//" for delta"
        stop
      endif
      if (delta<0) then
        write(0,"(a,f10.3,a)") "*** Error: delta ",delta,"<0"
        stop
      endif
    endif
    
    zeta=-1
    if (iarg>5) then
      call get_command_argument(7,arg)
      read(arg,*,iostat=ios) zeta
      if (ios/=0) then
        write(0,"(a)") "*** Error: invalid real "//trim(arg)//" for zeta"
        stop
      endif
      if (zeta<0) then
        write(0,"(a,f10.3,a)") "*** Error: zeta ",zeta,"<0"
        stop
      endif
    endif
    
    regularisation=-1
    if (iarg>6) then
      call get_command_argument(7,arg)
      read(arg,*,iostat=ios) regularisation
      if (ios/=0) then
        write(0,"(a)") "*** Error: invalid real "//trim(arg)//" for regularisation"
        stop
      endif
      if (regularisation<0) then
        write(0,"(a,f10.3,a)") "*** Error: regularisation ",regularisation,"<0"
        stop
      endif
    endif
    
    covariance_cutoff=-1
    if (iarg>7) then
      call get_command_argument(8,arg)
      read(arg,*,iostat=ios) covariance_cutoff
      if (ios/=0) then
        write(0,"(a)") "*** Error: invalid real "//trim(arg)//" for covariance_cutoff"
        stop
      endif
      if (covariance_cutoff<0) then
        write(0,"(a,f10.3,a)") "*** Error: covariance_cutoff ",covariance_cutoff,"<0"
        stop
      endif
    endif
    
  end subroutine read_cmd
  
end program refit_ml
