! MLKMC - program for merging predictor directories and md stats
! The merge result can contain duplicates, run remove_duplicates after
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

program merge_predictors
  use libAtoms_module
  use kmc_module
  implicit none

  type(t_predictor) :: source1_predictor,source2_predictor,target_predictor

  integer :: i,shell,coord_target,coord_jumping,coord_middle
  character(len=80) :: source1_dir,source2_dir,target_dir,predictor_file
  integer :: error
  logical :: source1_exist,source1_pred_exist,source2_exist,source2_pred_exist,target_exist,target_pred_exist
  logical :: exist

  error=0
  call read_cmd()

  source1_pred_exist=.false.
  source2_pred_exist=.false.
  target_pred_exist=.false.
#ifdef INTEL_COMPILER
  inquire(directory=source1_dir,exist=source1_exist)
  if (source1_exist) then
    inquire(directory=trim(source1_dir)//"/predictors",exist=source1_pred_exist)
  endif
  inquire(directory=source2_dir,exist=source2_exist)
  if (source2_exist) then
    inquire(directory=trim(source2_dir)//"/predictors",exist=source2_pred_exist)
  endif
  inquire(directory=target_dir,exist=target_exist)
  if (target_exist) then
    inquire(directory=trim(target_dir)//"/predictors",exist=target_pred_exist)
  endif
#else
  inquire(file=source1_dir,exist=source1_exist)
  if (source1_exist) then
    inquire(file=trim(source1_dir)//"/predictors",exist=source1_pred_exist)
  endif
  inquire(file=source2_dir,exist=source2_exist)
  if (source2_exist) then
    inquire(file=trim(source2_dir)//"/predictors",exist=source2_pred_exist)
  endif
  inquire(file=target_dir,exist=target_exist)
  if (target_exist) then
    inquire(file=trim(target_dir)//"/predictors",exist=target_pred_exist)
  endif
#endif

  if (.not. source1_pred_exist) then
    write(0,"(a)") "*** Error: source directory 1 "//trim(source1_dir)//"/predictors doesn't exist"
    error=1
  endif
  if (.not. source2_pred_exist) then
    write(0,"(a)") "*** Error: source directory 2 "//trim(source2_dir)//"/predictors doesn't exist"
    error=1
  endif
  if (target_pred_exist) then
    write(0,"(a)") "*** Error: target directory "//trim(target_dir)//"/predictors DOES already exist, unsafe to merge"
    error=1
  endif
  if (error/=0) then
    stop
  endif

  if (.not. target_exist) then
    call execute_command_line("mkdir "//trim(target_dir),exitstat=error)
    if (error/=0) then
      write(0,"(a)") "*** Error making target directory"
      stop
    endif
  endif
  call execute_command_line("mkdir "//trim(target_dir)//"/predictors",exitstat=error)
  if (error/=0) then
    write(0,"(a)") "*** Error making target directory/predictors"
    stop
  endif

  call system_initialise()
  call verbosity_push(PRINT_SILENT)

  do shell=1,N_JUMP_SHELLS
    do coord_target=3,N_NN_MAX
      do coord_jumping=3,LOCK_THRESHOLD
        predictor_file="predictors/hop_predictor_"//coord_jumping//"_"//coord_target//"_"//shell//".xml"

        inquire(file=trim(source1_dir)//"/"//trim(predictor_file),exist=exist)
        if (exist) then
          call read_predictor_xml(source1_predictor,trim(source1_dir)//"/"//trim(predictor_file))
          ! Source 1 takes precedence in parameters
          call initialise_predictor(target_predictor,source1_predictor%n_y_min,source1_predictor%n_y_max,&
                                    source1_predictor%n_sparse_max,source1_predictor%desc_str,&
                                    source1_predictor%delta,source1_predictor%zeta,&
                                    source1_predictor%covariance_type,source1_predictor%regularisation,error=error)
          if (error/=0) then
            write(0,"(a)") "*** Error initializing target_predictor based on source 1 "//trim(predictor_file)
            stop
          endif
        endif
        inquire(file=trim(source2_dir)//"/"//trim(predictor_file),exist=exist)
        if (exist) then
          call read_predictor_xml(source2_predictor,trim(source2_dir)//"/"//trim(predictor_file))
          if (.not. target_predictor%initialised) then
            ! Source 2 parameters are used if source 1 didn't have any
            call initialise_predictor(target_predictor,source2_predictor%n_y_min,source2_predictor%n_y_max,&
                                      source2_predictor%n_sparse_max,source2_predictor%desc_str,&
                                      source2_predictor%delta,source2_predictor%zeta,&
                                      source2_predictor%covariance_type,source2_predictor%regularisation,error=error)
            if (error/=0) then
              write(0,"(a)") "*** Error initializing target_predictor based on source 2 "//trim(predictor_file)
              stop
            endif
          endif
        endif
        if (target_predictor%initialised) then
          if (source1_predictor%initialised) then
            do i=1,source1_predictor%n_y
              call add_data_point(target_predictor,source1_predictor%x(:,i),source1_predictor%y(i),source1_predictor%y_error(i),&
                                  source1_predictor%covariance_cutoff(i),error=error)
              if (error/=0 .and. error/=2) then
                write(0,"(a,i0,a)") "*** Error adding data point ",i," from source 1 "//trim(predictor_file)
                stop
              endif
            end do
          endif
          if (source2_predictor%initialised) then
            do i=1,source2_predictor%n_y
              call add_data_point(target_predictor,source2_predictor%x(:,i),source2_predictor%y(i),source2_predictor%y_error(i),&
                                  source2_predictor%covariance_cutoff(i),error=error)
              if (error/=0 .and. error/=2) then
                write(0,"(a,i0,a)") "*** Error adding data point ",i," from source 2 "//trim(predictor_file)
                stop
              endif
            end do
          endif
          call print_predictor_xml(target_predictor,trim(target_dir)//"/"//trim(predictor_file),error=error)
          if (error/=0) then
            write(0,"(a)") "*** Error printing target "//trim(predictor_file)
            stop
          endif
        endif
        call finalise_predictor(target_predictor)
        call finalise_predictor(source1_predictor)
        call finalise_predictor(source2_predictor)

        if (shell==1) then
          do coord_middle=3,N_NN_MAX-1
            predictor_file="predictors/exchange_predictor_"//coord_middle//"_"//coord_jumping//"_"//coord_target//".xml"

            inquire(file=trim(source1_dir)//"/"//trim(predictor_file),exist=exist)
            if (exist) then
              call read_predictor_xml(source1_predictor,trim(source1_dir)//"/"//trim(predictor_file))
              ! Source 1 takes precedence in parameters
              call initialise_predictor(target_predictor,source1_predictor%n_y_min,source1_predictor%n_y_max,&
                                        source1_predictor%n_sparse_max,source1_predictor%desc_str,source1_predictor%delta,&
                                        source1_predictor%zeta,source1_predictor%covariance_type,&
                                        source1_predictor%regularisation,error=error)
              if (error/=0) then
                write(0,"(a)") "*** Error initializing target_predictor based on source 1 "//trim(predictor_file)
                stop
              endif
            endif
            inquire(file=trim(source2_dir)//"/"//trim(predictor_file),exist=exist)
            if (exist) then
              call read_predictor_xml(source2_predictor,trim(source2_dir)//"/"//trim(predictor_file))
              if (.not. target_predictor%initialised) then
                ! Source 2 parameters are used if source 1 didn't have any
                call initialise_predictor(target_predictor,source2_predictor%n_y_min,source2_predictor%n_y_max,&
                                          source2_predictor%n_sparse_max,source2_predictor%desc_str,source2_predictor%delta,&
                                          source2_predictor%zeta,source2_predictor%covariance_type,&
                                          source2_predictor%regularisation,error=error)
                if (error/=0) then
                  write(0,"(a)") "*** Error initializing target_predictor based on source 2 "//trim(predictor_file)
                  stop
                endif
              endif
            endif
            if (target_predictor%initialised) then
              if (source1_predictor%initialised) then
                do i=1,source1_predictor%n_y
                  call add_data_point(target_predictor,source1_predictor%x(:,i),source1_predictor%y(i),&
                                      source1_predictor%y_error(i),source1_predictor%covariance_cutoff(i),error=error)
                  if (error/=0 .and. error/=2) then
                    write(0,"(a,i0,a)") "*** Error adding data point ",i," from source 1 "//trim(predictor_file)
                    stop
                  endif
                end do
              endif
              if (source2_predictor%initialised) then
                do i=1,source2_predictor%n_y
                  call add_data_point(target_predictor,source2_predictor%x(:,i),source2_predictor%y(i),&
                                      source2_predictor%y_error(i),source2_predictor%covariance_cutoff(i),error=error)
                  if (error/=0 .and. error/=2) then
                    write(0,"(a,i0,a)") "*** Error adding data point ",i," from source 2 "//trim(predictor_file)
                    stop
                  endif
                end do
              endif
              call print_predictor_xml(target_predictor,trim(target_dir)//"/"//trim(predictor_file),error=error)
              if (error/=0) then
                write(0,"(a)") "*** Error printing target "//trim(predictor_file)
                stop
              endif
            endif
            call finalise_predictor(target_predictor)
            call finalise_predictor(source1_predictor)
            call finalise_predictor(source2_predictor)

          end do
        endif
      end do
    end do
  end do

contains

  subroutine read_cmd()
    integer :: iarg
    character(len=80) :: arg

    iarg=command_argument_count()
    if (iarg<3) then
      call get_command_argument(0,arg)
      write(0,"(a)") "Usage: "//trim(arg)//" source-dir1 source-dir2 target-dir"
      stop
    endif

    call get_command_argument(1,source1_dir)
    call get_command_argument(2,source2_dir)
    call get_command_argument(3,target_dir)
    
  end subroutine read_cmd
  
  integer function get_max_jump_shell(filename,error)
    ! Linus Torvalds forgive me
    character(len=*),intent(in) :: filename
    integer,optional,intent(out) :: error

    integer :: n,ios,fileunit

    if (present(error)) then
      error=0
    endif

    open(newunit=fileunit,file=filename,action="read",iostat=ios)
    if (ios/=0) then
      write(0,"(a)") "*** Error opening file "//trim(filename)
      if (present(error)) then
        error=ios
      endif
      return
    endif

    n=0
    do
      read(fileunit,*,iostat=ios)
      if (ios/=0) then
        exit
      endif
      n=n+1
    end do

    get_max_jump_shell=(n-(N_NN_MAX-3+1)*(LOCK_THRESHOLD-3+1)*(N_NN_MAX-3))/(LOCK_THRESHOLD-3+1)/(N_NN_MAX-3+1)

    close(fileunit)
  end function get_max_jump_shell

end program merge_predictors
