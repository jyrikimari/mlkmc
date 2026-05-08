! MLKMC - module describing the atoms derived type
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

module mlkmc_atoms_types_module
  use constants_module
  use geometry_module
  implicit none

  type t_atoms
    real(kind=rk) :: lattice(3,3)
    integer :: first_id=0,last_id=0
    real(kind=rk),allocatable :: pos(:,:)
    integer,allocatable :: Z(:)
    integer,allocatable :: atom_type(:)
    logical :: initialised=.false.
    logical :: periodic(3)=.false.
  end type t_atoms

  type t_neighbour_shell
    logical :: initialised=.false.
    integer :: n_neighbour_sites=0
    integer,allocatable :: neighbour_sites(:)
    integer,allocatable :: reverse_indices(:)
  end type t_neighbour_shell

  type t_site
    integer :: Z=-1,nn=0,n_near=0,atom_type=0,id=0,max_nn=0,nn_substrate=0,inactive_counter=0
    logical :: locked=.false.,lock_checked=.false.,jumps_checked=.false.,substrate=.false.,initial_site=.false.
    real(kind=rk),dimension(3) :: cartesian_coords

    integer :: my_box(3),my_index_in_box

    type(t_neighbour_shell) :: neighbour_shells(NN_SHELLS)

    integer,allocatable :: hop_states(:,:)
    real(kind=rk),allocatable :: hop_barriers(:,:,:)
    double precision,allocatable :: hop_rates(:,:)

    integer,allocatable :: exchange_states(:,:)
    real(kind=rk),allocatable :: exchange_barriers(:,:,:)
    double precision,allocatable :: exchange_rates(:,:)

    double precision :: hop_rate=0.0d0,exchange_rate=0.0d0

  end type t_site

  type t_system
    integer :: n_sites=0,n_substrate=0,n_deposited=0
    integer,allocatable :: deposited_atoms(:)
    real(kind=rk) :: lattice_constant,max_jump_dist,lae_cutoff_dist,height
    real(kind=rk),dimension(3,3) :: lattice
    type(t_site),allocatable :: sites(:)
    type(t_atoms) :: substrate,frozen_substrate

    integer :: max_nn=0,max_1nn=0,max_jump_shell,lae_cutoff_shell
    integer :: substrate_type=SUBSTRATE_TYPE_REGULAR

    real(kind=rk) :: first_atom_initial_coords(3)
    real(kind=rk) :: first_atom_unwrapped_coords(3)

    real(kind=rk) :: nn_dist2,relax_displacement_tol2
    real(kind=rk) :: layer_separation_111,layer_separation_111_squared

    real(kind=rk) :: one_box_dim(3)
    integer :: n_boxes(3)
    integer,allocatable :: sites_in_boxes(:,:,:,:),n_sites_in_boxes(:,:,:)

    logical :: allow_exchange=.true.
    logical :: periodic(3)=.false.
  end type t_system

end module mlkmc_atoms_types_module
