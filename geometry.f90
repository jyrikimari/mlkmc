! MLKMC - module for geometry functions and subroutines
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

module geometry_module
  use constants_module
  use linearalgebra_module
  implicit none

contains

  function cross_product(u,v) result(c)
    real(kind=rk),intent(in),dimension(3) :: u,v
    real(kind=rk),dimension(3) :: c

    c=[ u(2)*v(3)-u(3)*v(2),&
        u(3)*v(1)-u(1)*v(3),&
        u(1)*v(2)-u(2)*v(1) ]
  end function cross_product

  function denormal_surface_plane(vectors)
    real(kind=rk),intent(in),dimension(3,3) :: vectors
    real(kind=rk),dimension(3) :: denormal_surface_plane

    real(kind=rk),dimension(3) :: v1,v2

    v1=vectors(:,2)-vectors(:,1)
    v2=vectors(:,3)-vectors(:,1)

    denormal_surface_plane=cross_product(v1,v2)

  end function denormal_surface_plane

  function geometric_centre(vectors) result(centre)
    real(kind=rk),intent(in),dimension(3,3) :: vectors
    real(kind=rk),dimension(3) :: centre

    centre=sum(vectors,dim=2)/3
  end function geometric_centre

  function circumcentre(vectors) result(centre)
    real(kind=rk),intent(in),dimension(3,3) :: vectors
    real(kind=rk),dimension(3) :: centre,a,b,c,crossbaca

    a=vectors(:,1)
    b=vectors(:,2)
    c=vectors(:,3)

    crossbaca=cross_product(b-a,c-a)

    centre=a+( normsq(real(c-a,kind=my_dp))*cross_product( crossbaca,b-a ) &
              +normsq(real(b-a,kind=my_dp))*cross_product( c-a,crossbaca ) )&
             / ( 2*normsq(real(crossbaca,kind=my_dp)) )

  end function circumcentre

  real(kind=rk) function dist2_to_line(vectors)
    real(kind=rk),intent(in),dimension(3,3) :: vectors
    real(kind=rk),dimension(3) :: a,b,c

    a=vectors(:,1)
    b=vectors(:,2)
    c=vectors(:,3)

    dist2_to_line=normsq(real(cross_product(a-c,b-c),kind=my_dp))/normsq(real(a-b,kind=my_dp))

  end function dist2_to_line

  logical function are_parallel(u,v)
    real(kind=rk),intent(in),dimension(3) :: u,v
    real(kind=rk),dimension(3) :: w
    real(kind=rk),parameter :: eps=1e-5

    are_parallel=.false.
    w=u/v
    if (abs(w(1)-w(2))<eps .and. abs(w(2)-w(3))<eps) then
      are_parallel=.true.
    endif

  end function are_parallel

end module geometry_module
