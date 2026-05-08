! MLKMC - program for k-fold training and validation of a predictor
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

program train_ml
  use libAtoms_module
  use ml_module
  use utils_module
  implicit none

  type(t_predictor) :: predictor

  double precision :: frac,y_error,delta,zeta,regularisation,covariance_cutoff
  double precision :: y_pred,variance,sumsq
  double precision,allocatable :: y(:),x(:,:),y_k_train(:,:),y_k_test(:,:),x_k_train(:,:,:),x_k_test(:,:,:)
  integer,allocatable :: indices(:)
  integer :: i,j,k,seed,n_sparse,n_y,d,n_y_train,n_y_test,train_count,test_count,covariance_type,error,calib_count,n_sparse_max
  character(len=80) :: filename
  type(extendable_str) :: desc_str

  seed=123456
  call read_cmd()

  call system_initialise()
  call verbosity_push(PRINT_SILENT)

  call seed_mt(seed)

  call read_predictor_xml(predictor,filename)
  n_y=predictor%n_y 
  d=predictor%d
  allocate(y(n_y))
  allocate(indices(n_y))
  allocate(x(d,n_y))
  y=predictor%y
  x=predictor%x
  covariance_type=predictor%covariance_type
  desc_str=predictor%desc_str

  n_y_test=n_y/k
  n_y_train=n_y-n_y_test
  if (n_sparse_max>0) then
    n_sparse=min(int(n_y_train*frac),n_sparse_max)
  else
    n_sparse=int(n_y_train*frac)
  endif
  
  allocate(y_k_train(n_y_train,k))
  allocate(y_k_test(n_y_test,k))
  allocate(x_k_train(d,n_y_train,k))
  allocate(x_k_test(d,n_y_test,k))

  call shuffle(y,indices=indices)
  do i=1,k
    test_count=0
    train_count=0
    do j=1,n_y
      if (j>(i-1)*n_y_test .and. j<=i*n_y_test) then
        test_count=test_count+1
        y_k_test(test_count,i)=y(j)
        x_k_test(:,test_count,i)=x(:,indices(j))
      else
        train_count=train_count+1
        y_k_train(train_count,i)=y(j)
        x_k_train(:,train_count,i)=x(:,indices(j))
      end if
    end do
  end do

  write(6,"(a5,a10,a10)") "i","RMSE(eV)","calib(%)"
  do i=1,k
    open(unit=10,file="result_"//i//".dat",action="write")
    write(10,"(a10,a10,a10,a10)") "i","y","y_pred","stdev"
    call initialise_predictor(predictor,n_y_train,n_y_train,n_sparse,desc_str,&
                              delta=delta,zeta=zeta,&
                              covariance_type=covariance_type,regularisation=regularisation,error=error)
    if (error/=0) then
      write(0,"(a,i0)") "*** Error: couldn't initialise predictor ",i
      stop
    endif
    do j=1,n_y_train
      call add_data_point(predictor,x_k_train(:,j,i),y_k_train(j,i),y_error,covariance_cutoff,error=error)
      if (error/=0) then
        if (error/=2 .and. error/=3) then
          write(0,"(a,i0,a,i0)") "*** Error: couldn't add data point ",j," to predictor ",i
          stop
        endif
      endif
    end do
    call fit_predictor(predictor,error=error)
    if (error/=0) then
      write(0,"(a,i0)") "*** Error: couldn't fit predictor ",i
      stop
    endif
    sumsq=0.0d0
    calib_count=0
    do j=1,n_y_test
      call predict(predictor,x_k_test(:,j,i),y_pred,variance,error=error)
      if (error/=0) then
        write(0,"(a,i0,a,i0)") "*** Error: couldn't get prediction for data point ",j," from predictor ",i
        stop
      endif
      sumsq=sumsq+(y_k_test(j,i)-y_pred)**2 
      if (abs(y_k_test(j,i)-y_pred)<1.96*sqrt(variance)) then
        calib_count=calib_count+1
      endif
      write(10,"(i10,f10.6,f10.6,f10.6)") j,y_k_test(j,i),y_pred,sqrt(variance)
    end do
    write(6,"(i5,f10.6,f10.2)") i,sqrt(sumsq/n_y_test),100.0d0*calib_count/n_y_test
    close(10)
  end do

contains

  subroutine read_cmd()
    integer :: iarg,ios
    character(len=80) :: arg

    iarg=command_argument_count()
    ! File, k (of k-fold x-validation), frac_sparse, y_error, delta, zeta, regularisation, covariance cutoff
    if (iarg<8) then
      call get_command_argument(0,arg)
      write(0,"(a)") "Usage: "//trim(arg)//" file k frac y_error delta zeta regularisation covariance_cutoff [seed [n_sparse_max]]"
      stop
    endif

    call get_command_argument(1,filename)
    call get_command_argument(2,arg)
    read(arg,*,iostat=ios) k
    if (ios/=0) then
      write(0,"(a)") "*** Error: invalid integer "//trim(arg)//" for k of k-fold cross-validation"
      stop
    endif

    call get_command_argument(3,arg)
    read(arg,*,iostat=ios) frac
    if (ios/=0) then
      write(0,"(a)") "*** Error: invalid real "//trim(arg)//" for fraction of sparse points"
      stop
    endif

    call get_command_argument(4,arg)
    read(arg,*,iostat=ios) y_error
    if (ios/=0) then
      write(0,"(a)") "*** Error: invalid real "//trim(arg)//" for y_error"
      stop
    endif
  
    call get_command_argument(5,arg)
    read(arg,*,iostat=ios) delta
    if (ios/=0) then
      write(0,"(a)") "*** Error: invalid real "//trim(arg)//" for delta"
      stop
    endif
    
    call get_command_argument(6,arg)
    read(arg,*,iostat=ios) zeta
    if (ios/=0) then
      write(0,"(a)") "*** Error: invalid real "//trim(arg)//" for zeta"
      stop
    endif
    
    call get_command_argument(7,arg)
    read(arg,*,iostat=ios) regularisation
    if (ios/=0) then
      write(0,"(a)") "*** Error: invalid real "//trim(arg)//" for regularisation"
      stop
    endif
    
    call get_command_argument(8,arg)
    read(arg,*,iostat=ios) covariance_cutoff
    if (ios/=0) then
      write(0,"(a)") "*** Error: invalid real "//trim(arg)//" for covariance_cutoff"
      stop
    endif
    
    if (iarg>8) then
      call get_command_argument(9,arg)
      read(arg,*,iostat=ios) seed
      if (ios/=0) then
        write(0,"(a)") "*** Error: invalid integer "//trim(arg)//" for seed"
        stop
      endif
    endif

    n_sparse_max=-1
    if (iarg>9) then
      call get_command_argument(10,arg)
      read(arg,*,iostat=ios) n_sparse_max
      if (ios/=0) then
        write(0,"(a)") "*** Error: invalid integer "//trim(arg)//" for n_sparse_max"
        stop
      endif
      if (n_sparse_max<1) then
        write(0,"(a,i0,a)") "*** Error: n_sparse_max ",n_sparse_max," must be larger than 0"
        stop
      endif
    endif

  end subroutine read_cmd
  
end program train_ml
