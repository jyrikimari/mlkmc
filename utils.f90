! MLKMC - module for generic utilities
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

module utils_module
  use constants_module
  use mlkmc_atoms_types_module
  use mt19937_64
  implicit none

  interface mergesort
    module procedure mergesort_int,mergesort_real
  end interface

  interface shuffle
    module procedure shuffle_int,shuffle_real,shuffle_double
  end interface

  interface size_up
    module procedure size_up_int,size_up_int_2D,size_up_int_4D,size_up_real,&
                     size_up_real_2D,size_up_real_3D,size_up_logical,size_up_sites,&
                     size_up_double,size_up_double_2D
  end interface

  interface seed_mt
    module procedure seed_mt4,seed_mt8
  end interface

contains

  recursive subroutine mergesort_int(x,indices)
    integer,intent(inout) :: x(:)
    integer,optional,intent(out) :: indices(size(x))

    integer,allocatable :: sorted(:),ind(:)
    integer :: pivot,i,ii,ri,li,s

    s=size(x)
    if (s==1) then
      if (present(indices)) then
        indices(1)=1
      endif
      return
    endif
    allocate(sorted(s))

    pivot=s/2

    if (present(indices)) then
      call mergesort_int(x(:pivot),indices=indices(:pivot))
      call mergesort_int(x(pivot+1:),indices=indices(pivot+1:))
      indices(pivot+1:)=indices(pivot+1:)+pivot
      allocate(ind(s))
    else
      call mergesort_int(x(:pivot))
      call mergesort_int(x(pivot+1:))
    endif
    
    li=1
    ri=pivot+1
    do i=1,s
      if (li>pivot) then
        ii=ri
        ri=ri+1
      elseif (ri>s) then
        ii=li
        li=li+1
      elseif (x(li)<x(ri)) then
        ii=li
        li=li+1
      else
        ii=ri
        ri=ri+1
      endif
      sorted(i)=x(ii)
      if (present(indices)) then
        ind(i)=indices(ii)
      endif
    end do
    x=sorted
    deallocate(sorted)
    if (present(indices)) then
      indices=ind
      deallocate(ind)
    endif
  end subroutine mergesort_int

  recursive subroutine mergesort_real(x,indices)
    real(kind=rk),intent(inout) :: x(:)
    integer,optional,intent(out) :: indices(size(x))

    integer,allocatable :: ind(:)
    real(kind=rk),allocatable :: sorted(:)
    integer :: pivot,i,ii,ri,li,s

    s=size(x)
    if (s==1) then
      if (present(indices)) then
        indices(1)=1
      endif
      return
    endif
    allocate(sorted(s))

    pivot=s/2
    if (present(indices)) then
      call mergesort_real(x(:pivot),indices=indices(:pivot))
      call mergesort_real(x(pivot+1:),indices=indices(pivot+1:))
      indices(pivot+1:)=indices(pivot+1:)+pivot
      allocate(ind(s))
    else
      call mergesort_real(x(:pivot))
      call mergesort_real(x(pivot+1:))
    endif

    li=1
    ri=pivot+1
    do i=1,s
      if (li>pivot) then
        ii=ri
        ri=ri+1
      elseif (ri>s) then
        ii=li
        li=li+1
      elseif (x(li)<x(ri)) then
        ii=li
        li=li+1
      else
        ii=ri
        ri=ri+1
      endif
      sorted(i)=x(ii)
      if (present(indices)) then
        ind(i)=indices(ii)
      endif
    end do
    x=sorted
    deallocate(sorted)
    if (present(indices)) then
      indices=ind
      deallocate(ind)
    endif
  end subroutine mergesort_real

  subroutine shuffle_int(x,indices)
    integer,intent(inout) :: x(:)
    integer,optional,intent(out) :: indices(size(x))

    integer :: i,j,tmp,itmp
    real(kind=rk) :: u

    if (present(indices)) then
      indices=[(i,i=1,size(x))]
    endif

    do i=size(x),2,-1
      u=genrand64_real2()
      j=int(i*u)+1
      tmp=x(j)
      x(j)=x(i)
      x(i)=tmp
      if (present(indices)) then
        itmp=indices(j)
        indices(j)=indices(i)
        indices(i)=itmp
      endif
    end do
    
  end subroutine shuffle_int

  subroutine shuffle_real(x,indices)
    real(kind=sp),intent(inout) :: x(:)
    integer,optional,intent(out) :: indices(size(x))

    integer :: i,j,itmp
    real(kind=rk) :: tmp,u

    if (present(indices)) then
      indices=[(i,i=1,size(x))]
    endif

    do i=size(x),2,-1
      u=genrand64_real2()
      j=int(i*u)+1
      tmp=x(j)
      x(j)=x(i)
      x(i)=tmp
      if (present(indices)) then
        itmp=indices(j)
        indices(j)=indices(i)
        indices(i)=itmp
      endif
    end do
    
  end subroutine shuffle_real

  subroutine shuffle_double(x,indices)
    double precision,intent(inout) :: x(:)
    integer,optional,intent(out) :: indices(size(x))

    integer :: i,j,itmp
    double precision :: tmp
    real(kind=rk) :: u

    if (present(indices)) then
      indices=[(i,i=1,size(x))]
    endif

    do i=size(x),2,-1
      u=genrand64_real2()
      j=int(i*u)+1
      tmp=x(j)
      x(j)=x(i)
      x(i)=tmp
      if (present(indices)) then
        itmp=indices(j)
        indices(j)=indices(i)
        indices(i)=itmp
      endif
    end do
    
  end subroutine shuffle_double

  subroutine size_up_int(array,default_value,custom_buffer_size,error)
    integer,allocatable,intent(inout) :: array(:)
    integer,optional,intent(in) :: default_value
    integer,optional,intent(in) :: custom_buffer_size
    integer,optional,intent(out) :: error

    integer,allocatable :: array_tmp(:)
    integer :: l(1),u(1),my_buffer_size

    if (present(error)) then
      error=0
    endif
    if (present(custom_buffer_size)) then
      my_buffer_size=custom_buffer_size
    else
      my_buffer_size=BUFFER_SIZE
    endif

    if (.not. allocated(array)) then
      write(0,"(a)") "*** Error: array not allocated in size_up_int. Bug?"
      if (present(error)) then
        error=1
      endif
      return
    endif

    l=lbound(array)
    u=ubound(array)

    allocate(array_tmp(l(1):u(1)))
    array_tmp=array
    deallocate(array)
    allocate(array(l(1):u(1)+my_buffer_size))
    if (present(default_value)) then
      array=default_value
    endif
    array(l(1):u(1))=array_tmp
    deallocate(array_tmp)

  end subroutine size_up_int

  subroutine size_up_sites(array,error)
    type(t_site),allocatable,intent(inout) :: array(:)
    integer,optional,intent(out) :: error

    type(t_site),allocatable :: array_tmp(:)
    integer :: s

    if (present(error)) then
      error=0
    endif

    if (.not. allocated(array)) then
      write(0,"(a)") "*** Error: array not allocated in size_up_sites. Bug?"
      if (present(error)) then
        error=1
      endif
      return
    endif

    s=size(array)

    allocate(array_tmp(s))
    array_tmp=array
    deallocate(array)
    allocate(array(s+BUFFER_SIZE))
    array(1:s)=array_tmp
    deallocate(array_tmp)

  end subroutine size_up_sites

  subroutine size_up_int_2D(array,rank,default_value,custom_buffer_size,error)
    integer,allocatable,intent(inout) :: array(:,:)
    integer,optional,intent(in) :: rank
    integer,optional,intent(in) :: default_value
    integer,optional,intent(in) :: custom_buffer_size
    integer,optional,intent(out) :: error

    integer,allocatable :: array_tmp(:,:)
    integer :: l(2),u(2),my_buffer_size

    if (present(error)) then
      error=0
    endif
    if (present(custom_buffer_size)) then
      my_buffer_size=custom_buffer_size
    else
      my_buffer_size=BUFFER_SIZE
    endif

    if (.not. allocated(array)) then
      write(0,"(a)") "*** Error: array not allocated in size_up_int_2D. Bug?"
      if (present(error)) then
        error=1
      endif
      return
    endif

    l=lbound(array)
    u=ubound(array)

    allocate(array_tmp(l(1):u(1),l(2):u(2)))
    array_tmp=array
    deallocate(array)
    if (present(rank)) then
      select case(rank)
        case(RANK_FIRST)
          allocate(array(l(1):u(1)+my_buffer_size,l(2):u(2)))
        case(RANK_SECOND)
          allocate(array(l(1):u(1),l(2):u(2)+my_buffer_size))
        case(RANK_ALL)
          allocate(array(l(1):u(1)+my_buffer_size,l(2):u(2)+my_buffer_size))
        case default
          write(0,"(a,i0,a)") "*** Error: unknown rank ",rank," in size_up_int_2D"
          if (present(error)) then
            error=1
          endif
          return
      end select
    else
      allocate(array(l(1):u(1),l(2):u(2)+my_buffer_size))
    endif
    if (present(default_value)) then
      array=default_value
    endif
    array(l(1):u(1),l(2):u(2))=array_tmp
    deallocate(array_tmp)

  end subroutine size_up_int_2D

  subroutine size_up_real(array,default_value,custom_buffer_size,error)
    real(kind=sp),allocatable,intent(inout) :: array(:)
    real(kind=sp),optional,intent(in) :: default_value
    integer,optional,intent(in) :: custom_buffer_size
    integer,optional,intent(out) :: error

    real(kind=rk),allocatable :: array_tmp(:)
    integer :: l(1),u(1),my_buffer_size

    if (present(error)) then
      error=0
    endif
    if (present(custom_buffer_size)) then
      my_buffer_size=custom_buffer_size
    else
      my_buffer_size=BUFFER_SIZE
    endif

    if (.not. allocated(array)) then
      write(0,"(a)") "*** Error: array not allocated in size_up_real. Bug?"
      if (present(error)) then
        error=1
      endif
      return
    endif

    l=lbound(array)
    u=ubound(array)

    allocate(array_tmp(l(1):u(1)))
    array_tmp=array
    deallocate(array)
    allocate(array(l(1):u(1)+my_buffer_size))
    if (present(default_value)) then
      array=default_value
    endif
    array(l(1):u(1))=array_tmp
    deallocate(array_tmp)

  end subroutine size_up_real

  subroutine size_up_double(array,default_value,custom_buffer_size,error)
    double precision,allocatable,intent(inout) :: array(:)
    double precision,optional,intent(in) :: default_value
    integer,optional,intent(in) :: custom_buffer_size
    integer,optional,intent(out) :: error

    double precision,allocatable :: array_tmp(:)
    integer :: l(1),u(1),my_buffer_size

    if (present(error)) then
      error=0
    endif
    if (present(custom_buffer_size)) then
      my_buffer_size=custom_buffer_size
    else
      my_buffer_size=BUFFER_SIZE
    endif

    if (.not. allocated(array)) then
      write(0,"(a)") "*** Error: array not allocated in size_up_double. Bug?"
      if (present(error)) then
        error=1
      endif
      return
    endif

    l=lbound(array)
    u=ubound(array)

    allocate(array_tmp(l(1):u(1)))
    array_tmp=array
    deallocate(array)
    allocate(array(l(1):u(1)+my_buffer_size))
    if (present(default_value)) then
      array=default_value
    endif
    array(l(1):u(1))=array_tmp
    deallocate(array_tmp)

  end subroutine size_up_double

  subroutine size_up_real_2D(array,rank,default_value,custom_buffer_size,error)
    real(kind=sp),allocatable,intent(inout) :: array(:,:)
    integer,optional,intent(in) :: rank
    real(kind=sp),optional,intent(in) :: default_value
    integer,optional,intent(in) :: custom_buffer_size
    integer,optional,intent(out) :: error

    real(kind=rk),allocatable :: array_tmp(:,:)
    integer :: l(2),u(2),my_buffer_size

    if (present(error)) then
      error=0
    endif

    if (present(custom_buffer_size)) then
      my_buffer_size=custom_buffer_size
    else
      my_buffer_size=BUFFER_SIZE
    endif

    if (.not. allocated(array)) then
      write(0,"(a)") "*** Error: array not allocated in size_up_real_2D. Bug?"
      if (present(error)) then
        error=1
      endif
      return
    endif

    l=lbound(array)
    u=ubound(array)

    allocate(array_tmp(l(1):u(1),l(2):u(2)))
    array_tmp=array
    deallocate(array)
    if (present(rank)) then
      select case(rank)
        case(RANK_FIRST)
          allocate(array(l(1):u(1)+my_buffer_size,l(2):u(2)))
        case(RANK_SECOND)
          allocate(array(l(1):u(1),l(2):u(2)+my_buffer_size))
        case(RANK_ALL)
          allocate(array(l(1):u(1)+my_buffer_size,l(2):u(2)+my_buffer_size))
        case default
          write(0,"(a,i0,a)") "*** Error: unknown rank ",rank," in size_up_real_2D"
          if (present(error)) then
            error=1
          endif
          return
      end select
    else
      allocate(array(l(1):u(1),l(2):u(2)+my_buffer_size))
    endif
    if (present(default_value)) then
      array=default_value
    endif
    array(l(1):u(1),l(2):u(2))=array_tmp
    deallocate(array_tmp)

  end subroutine size_up_real_2D

  subroutine size_up_double_2D(array,rank,default_value,custom_buffer_size,error)
    double precision,allocatable,intent(inout) :: array(:,:)
    integer,optional,intent(in) :: rank
    double precision,optional,intent(in) :: default_value
    integer,optional,intent(in) :: custom_buffer_size
    integer,optional,intent(out) :: error

    double precision,allocatable :: array_tmp(:,:)
    integer :: l(2),u(2),my_buffer_size

    if (present(error)) then
      error=0
    endif

    if (present(custom_buffer_size)) then
      my_buffer_size=custom_buffer_size
    else
      my_buffer_size=BUFFER_SIZE
    endif

    if (.not. allocated(array)) then
      write(0,"(a)") "*** Error: array not allocated in size_up_double_2D. Bug?"
      if (present(error)) then
        error=1
      endif
      return
    endif

    l=lbound(array)
    u=ubound(array)

    allocate(array_tmp(l(1):u(1),l(2):u(2)))
    array_tmp=array
    deallocate(array)
    if (present(rank)) then
      select case(rank)
        case(RANK_FIRST)
          allocate(array(l(1):u(1)+my_buffer_size,l(2):u(2)))
        case(RANK_SECOND)
          allocate(array(l(1):u(1),l(2):u(2)+my_buffer_size))
        case(RANK_ALL)
          allocate(array(l(1):u(1)+my_buffer_size,l(2):u(2)+my_buffer_size))
        case default
          write(0,"(a,i0,a)") "*** Error: unknown rank ",rank," in size_up_double_2D"
          if (present(error)) then
            error=1
          endif
          return
      end select
    else
      allocate(array(l(1):u(1),l(2):u(2)+my_buffer_size))
    endif
    if (present(default_value)) then
      array=default_value
    endif
    array(l(1):u(1),l(2):u(2))=array_tmp
    deallocate(array_tmp)

  end subroutine size_up_double_2D

  subroutine size_up_real_3D(array,rank,default_value,custom_buffer_size,error)
    real(kind=rk),allocatable,intent(inout) :: array(:,:,:)
    integer,optional,intent(in) :: rank
    real(kind=rk),optional,intent(in) :: default_value
    integer,optional,intent(in) :: custom_buffer_size
    integer,optional,intent(out) :: error

    real(kind=rk),allocatable :: array_tmp(:,:,:)
    integer :: l(3),u(3),my_buffer_size

    if (present(error)) then
      error=0
    endif

    if (present(custom_buffer_size)) then
      my_buffer_size=custom_buffer_size
    else
      my_buffer_size=BUFFER_SIZE
    endif

    if (.not. allocated(array)) then
      write(0,"(a)") "*** Error: array not allocated in size_up_real_3D. Bug?"
      if (present(error)) then
        error=1
      endif
      return
    endif

    l=lbound(array)
    u=ubound(array)

    allocate(array_tmp(l(1):u(1),l(2):u(2),l(3):u(3)))
    array_tmp=array
    deallocate(array)
    if (present(rank)) then
      select case(rank)
        case(RANK_FIRST)
          allocate(array(l(1):u(1)+my_buffer_size,l(2):u(2),l(3):u(3)))
        case(RANK_SECOND)
          allocate(array(l(1):u(1),l(2):u(2)+my_buffer_size,l(3):u(3)))
        case(RANK_THIRD)
          allocate(array(l(1):u(1),l(2):u(2),l(3):u(3)+my_buffer_size))
        case(RANK_LAST_TWO)
          allocate(array(l(1):u(1),l(2):u(2)+my_buffer_size,l(3):u(3)+my_buffer_size))
        case(RANK_ALL)
          allocate(array(l(1):u(1)+my_buffer_size,l(2):u(2)+my_buffer_size,l(3):u(3)+my_buffer_size))
        case default
          write(0,"(a,i0,a)") "*** Error: unknown rank ",rank," in size_up_real_3D"
          if (present(error)) then
            error=1
          endif
          return
      end select
    else
      allocate(array(l(1):u(1),l(2):u(2),l(3):u(3)+my_buffer_size))
    endif
    if (present(default_value)) then
      array=default_value
    endif
    array(l(1):u(1),l(2):u(2),l(3):u(3))=array_tmp
    deallocate(array_tmp)

  end subroutine size_up_real_3D

  subroutine size_up_int_4D(array,default_value,error)
    integer,allocatable,intent(inout) :: array(:,:,:,:)
    integer,optional,intent(in) :: default_value
    integer,optional,intent(out) :: error

    integer,allocatable :: array_tmp(:,:,:,:)
    integer :: l(4),u(4)

    if (present(error)) then
      error=0
    endif

    if (.not. allocated(array)) then
      write(0,"(a)") "*** Error: array not allocated in size_up_int_4D. Bug?"
      if (present(error)) then
        error=1
      endif
      return
    endif

    l=lbound(array)
    u=ubound(array)

    allocate(array_tmp(l(1):u(1),l(2):u(2),l(3):u(3),l(4):u(4)))
    array_tmp=array
    deallocate(array)
    allocate(array(l(1):u(1)+BUFFER_SIZE,l(2):u(2),l(3):u(3),l(4):u(4)))
    if (present(default_value)) then
      array=default_value
    endif
    array(l(1):u(1),l(2):u(2),l(3):u(3),l(4):u(4))=array_tmp
    deallocate(array_tmp)

  end subroutine size_up_int_4D

  subroutine seed_mt4(seed)
    integer,intent(in) :: seed

    call init_genrand64(int(seed,kind=8))
  end subroutine seed_mt4

  subroutine seed_mt8(seed)
    integer(kind=8),intent(in) :: seed

    call init_genrand64(seed)
  end subroutine seed_mt8

  subroutine size_up_logical(array,default_value,error)
    logical,allocatable,intent(inout) :: array(:)
    logical,optional,intent(in) :: default_value
    integer,optional,intent(out) :: error

    logical,allocatable :: array_tmp(:)
    integer :: l(1),u(1)

    if (present(error)) then
      error=0
    endif

    if (.not. allocated(array)) then
      write(0,"(a)") "*** Error: array not allocated in size_up_logical. Bug?"
      if (present(error)) then
        error=1
      endif
      return
    endif

    l=lbound(array)
    u=ubound(array)

    allocate(array_tmp(l(1):u(1)))
    array_tmp=array
    deallocate(array)
    allocate(array(l(1):u(1)+BUFFER_SIZE))
    if (present(default_value)) then
      array=default_value
    endif
    array(l(1):u(1))=array_tmp
    deallocate(array_tmp)

  end subroutine size_up_logical

  subroutine print_progress(progress)
    integer,intent(in) :: progress

    write(6,"(12a)",advance="no") char(8),char(8),char(8),char(8),char(8),char(8),&
                                  char(8),char(8),char(8),char(8),char(8),char(8)
    write(6,"(a,i0,a)",advance="no") "# ",progress," % done"
  end subroutine print_progress

  character(len=80) function print_time(t)
    double precision,intent(in) :: t
    
    double precision :: t_remaining
    integer :: days,hours,minutes

    if (t<0.001d0) then
      print_time="  < 1 ms"
      return
    endif

    t_remaining=t
    days=int(t_remaining/24/3600)
    t_remaining=t_remaining-days*24*3600
    hours=int(t_remaining/3600)
    t_remaining=t_remaining-hours*3600
    minutes=int(t_remaining/60)
    t_remaining=t_remaining-minutes*60

    if (days>0) then
      write(print_time,"(a,i0,a)") "  ",days," d,"
    else
      print_time=""
    endif
    if (days>0 .or. hours>0) then
      write(print_time,"(a,i0,a)") " "//trim(print_time)//" ",hours," h,"
    endif
    if (days>0 .or. hours>0 .or. minutes>0) then
      write(print_time,"(a,i0,a)") " "//trim(print_time)//" ",minutes," min,"
    endif
    if (t_remaining>=10.0d0) then
      write(print_time,"(a,f6.3,a)") " "//trim(print_time)//" ",t_remaining," s"
    else
      write(print_time,"(a,f5.3,a)") " "//trim(print_time)//" ",t_remaining," s"
    endif
  end function print_time

end module utils_module
