! MLKMC - module for atom-related functions and subroutines
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

module mlkmc_atoms_module
  use libAtoms_module
  use utils_module
  use mlkmc_atoms_types_module
  implicit none

  interface are_1nn_neighbours
    module procedure are_1nn_neighbours_dimer,are_1nn_neighbours_trimer
  end interface 

contains

  subroutine create_system(system,xmax,ymax,zmax,lattice_constant,substrate_layers,&
                           species,substrate_type,max_jump_nth,lae_cutoff_nth,allow_exchange,periodic,&
                           initial_atoms,initial_atoms_layers,initial_atoms_termination,initial_atoms_coords,error)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: xmax,ymax,zmax,substrate_layers,species,substrate_type,max_jump_nth,lae_cutoff_nth
    integer,intent(in) :: initial_atoms(:),initial_atoms_layers(:),initial_atoms_termination(:)
    double precision,intent(in) :: initial_atoms_coords(:,:)
    real(kind=rk),intent(in) :: lattice_constant
    logical,optional,intent(in) :: allow_exchange,periodic(3)
    integer,optional,intent(out) :: error

    integer :: i,j,x,y,z,p,internal_error,n_substrate,n_frozen_substrate,n_sites_old,deptype
    real(kind=rk) :: furthest_nn_dist

    if (present(error)) then
      error=0
    endif

    system%lattice_constant=lattice_constant
    system%nn_dist2=(NEIGH_DIST_FACT(1)*lattice_constant)**2
    system%relax_displacement_tol2=0.25_rk*system%nn_dist2
    system%layer_separation_111=lattice_constant/sqrt(3.0_rk)
    system%layer_separation_111_squared=system%layer_separation_111**2
    system%lattice=0
    system%lattice(1,1)=xmax*UNIT_CELL(1)*lattice_constant
    system%lattice(2,2)=ymax*UNIT_CELL(2)*lattice_constant
    system%lattice(3,3)=zmax*UNIT_CELL(3)*lattice_constant
    write(6,"(a,3(f8.4,a))") "# Lattice: ",system%lattice(1,1)," x ",system%lattice(2,2)," x ",system%lattice(3,3)," Å^3"
    write(6,"(a,i0,a)") "# Monolayer size: ",2*xmax*ymax," atoms"
    select case(max_jump_nth)
      case(1)
        system%max_jump_shell=TRUE_1NN_JUMP_SHELL
      case(2)
        system%max_jump_shell=TRUE_2NN_JUMP_SHELL
      case(3)
        system%max_jump_shell=TRUE_3NN_JUMP_SHELL
      case default
        write(0,"(a,i0)") "*** Error: invalid max_jump_nth ",max_jump_nth
        if (present(error)) then
          error=1
        endif
        return
    end select

    select case(lae_cutoff_nth)
      case(1)
        system%lae_cutoff_shell=TRUE_1NN_SHELL
      case(2)
        system%lae_cutoff_shell=TRUE_2NN_SHELL
      case(3)
        system%lae_cutoff_shell=TRUE_3NN_SHELL
      case(4)
        system%lae_cutoff_shell=TRUE_4NN_SHELL
      case default
        write(0,"(a,i0)") "*** Error: invalid lae_cutoff_nth ",lae_cutoff_nth
        if (present(error)) then
          error=1
        endif
        return
    end select

    system%lae_cutoff_dist=get_Xnn_dist(lattice_constant,lae_cutoff_nth,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error: couldn't set lae_cutoff_dist based on lae_cutoff_nth"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

    system%max_jump_dist=get_Xnn_dist(lattice_constant,max_jump_nth,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error: couldn't set max_jump_dist based on max_jump_nth"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

    if (present(allow_exchange)) then
      system%allow_exchange=allow_exchange
    else
      system%allow_exchange=.true.
    endif
    if (present(periodic)) then
      system%periodic=periodic
    else
      system%periodic=.true.
    endif
    allocate(system%sites(BUFFER_SIZE))

    system%frozen_substrate%first_id=1
    system%frozen_substrate%last_id=2*xmax*ymax
    allocate(system%frozen_substrate%pos(3,system%frozen_substrate%first_id:system%frozen_substrate%last_id))
    system%frozen_substrate%pos=0.0_rk
    allocate(system%frozen_substrate%atom_type(system%frozen_substrate%first_id:system%frozen_substrate%last_id))
    system%frozen_substrate%atom_type=1
    allocate(system%frozen_substrate%Z(system%frozen_substrate%first_id:system%frozen_substrate%last_id))
    system%frozen_substrate%Z=species

    system%substrate%first_id=system%frozen_substrate%last_id+1
    system%substrate%last_id=6*xmax*ymax*substrate_layers
    allocate(system%substrate%pos(3,system%substrate%first_id:system%substrate%last_id))
    system%substrate%pos=0.0d0
    allocate(system%substrate%atom_type(system%substrate%first_id:system%substrate%last_id))
    system%substrate%atom_type=1
    allocate(system%substrate%Z(system%substrate%first_id:system%substrate%last_id))
    system%substrate%Z=species

    furthest_nn_dist=get_Xnn_dist(lattice_constant,lae_cutoff_nth,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error: couldn't get ",lae_cutoff_nth,"-nearest neighbour distance"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    system%n_boxes(1)=floor(system%lattice(1,1)/furthest_nn_dist)
    system%n_boxes(2)=floor(system%lattice(2,2)/furthest_nn_dist)
    system%n_boxes(3)=floor(system%lattice(3,3)/furthest_nn_dist)
    system%one_box_dim(1)=system%lattice(1,1)/system%n_boxes(1)
    system%one_box_dim(2)=system%lattice(2,2)/system%n_boxes(2)
    system%one_box_dim(3)=system%lattice(3,3)/system%n_boxes(3)
    allocate(system%sites_in_boxes(BUFFER_SIZE,system%n_boxes(1),system%n_boxes(2),system%n_boxes(3)))
    allocate(system%n_sites_in_boxes(system%n_boxes(1),system%n_boxes(2),system%n_boxes(3)))
    system%sites_in_boxes=0
    system%n_sites_in_boxes=0

    system%n_sites=0
    n_frozen_substrate=0
    n_substrate=0
    system%height=0.0_rk
    write(6,"(a)",advance="no") "# Generating initial fcc sites..."
    do z=0,substrate_layers-1
      do y=0,ymax-1
        do x=0,xmax-1
          do p=1,6
            call add_site(system,get_cartesian(x,y,z,p)*lattice_constant,error=internal_error)
            if (internal_error/=0) then
              write(0,"(a)") "*** Error adding site while generating substrate in create_system()"
              if (present(error)) then
                error=internal_error
              endif
              return
            endif
            system%n_substrate=system%n_substrate+1
            system%sites(system%n_sites)%Z=species
            system%sites(system%n_sites)%atom_type=1
            system%sites(system%n_sites)%id=system%n_substrate
            system%sites(system%n_sites)%substrate=.true.
            if (system%sites(system%n_sites)%cartesian_coords(3)>system%height) then
              system%height=system%sites(system%n_sites)%cartesian_coords(3)
            endif
            if (z==0 .and. p<=2) then
              n_frozen_substrate=n_frozen_substrate+1
              system%frozen_substrate%pos(:,n_frozen_substrate)=system%sites(system%n_sites)%cartesian_coords(:)
            else
              n_substrate=n_substrate+1
              system%substrate%pos(:,system%frozen_substrate%last_id+n_substrate)=system%sites(system%n_sites)%cartesian_coords(:)
            endif
            if (z==substrate_layers-1 .and. p>=5) then
              call create_neighbour_list(system,system%n_sites)
              system%sites(system%n_sites)%nn=N_NN_MAX-3
            else
              system%sites(system%n_sites)%nn=N_NN_MAX
            endif
          end do
        end do
      end do
    end do
    write(6,"(a)") " done"
    write(6,"(a,i0,a)") "# Generated ",system%n_sites," initial fcc sites"

    write(6,"(a)",advance="no") "# Adding initial adsorption sites..."
    system%substrate_type=substrate_type
    n_sites_old=system%n_sites
    call add_adsorption_sites(system,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error adding adsorption sites in create_system()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    write(6,"(a)") " done"
    write(6,"(a,i0,a)") "# Added ",system%n_sites-n_sites_old," adsorption sites"

    select case(substrate_type)
      case(SUBSTRATE_TYPE_REGULAR)
        deptype=1
      case(SUBSTRATE_TYPE_WEAK)
        deptype=2
    end select

    allocate(system%deposited_atoms(BUFFER_SIZE))
    system%deposited_atoms=0
    do i=1,size(initial_atoms)
      select case(initial_atoms(i))
        case(INITIAL_ATOMS_RANDOM)
          do j=1,initial_atoms_layers(i)
            call add_atom_random(system,species,deptype,error)
            if (internal_error/=0) then
              write(0,"(a,i0,a)") "*** Error adding initial random atom ",j," in create_system()"
              if (present(error)) then
                error=internal_error
              endif
              return
            endif
          end do
        case(INITIAL_ATOMS_HEXAGON)
          call add_hexagon(system,initial_atoms_layers(i),species,deptype,error=internal_error)
          if (internal_error/=0) then
            write(0,"(a)") "*** Error adding initial hexagon in create_system()"
            if (present(error)) then
              error=internal_error
            endif
            return
          endif
        case(INITIAL_ATOMS_HEXAGON_PLUS_ADATOM)
          call add_hexagon(system,initial_atoms_layers(i),species,deptype,error=internal_error)
          if (internal_error/=0) then
            write(0,"(a)") "*** Error adding initial hexagon in create_system()"
            if (present(error)) then
              error=internal_error
            endif
            return
          endif
          call add_hexagon(system,1,species,deptype,error=internal_error)
          if (internal_error/=0) then
            write(0,"(a)") "*** Error adding initial adatom on top of hexagon in create_system()"
            if (present(error)) then
              error=internal_error
            endif
            return
          endif
        case(INITIAL_ATOMS_TRIANGLE)
          call add_triangle(system,initial_atoms_layers(i),species,deptype,initial_atoms_termination(i),error=internal_error)
          if (internal_error/=0) then
            write(0,"(a)") "*** Error adding initial triangle in create_system()"
            if (present(error)) then
              error=internal_error
            endif
            return
          endif
        case(INITIAL_ATOMS_TRIANGLE_PLUS_ADATOM)
          call add_triangle(system,initial_atoms_layers(i),species,deptype,initial_atoms_termination(i),error=internal_error)
          if (internal_error/=0) then
            write(0,"(a)") "*** Error adding initial triangle in create_system()"
            if (present(error)) then
              error=internal_error
            endif
            return
          endif
          call add_hexagon(system,1,species,deptype,error=internal_error)
          if (internal_error/=0) then
            write(0,"(a)") "*** Error adding initial adatom on top of triangle in create_system()"
            if (present(error)) then
              error=internal_error
            endif
            return
          endif
        case(INITIAL_ATOMS_SINGLE)
          call add_atom_xyz(system,initial_atoms_coords(1,i),initial_atoms_coords(2,i),initial_atoms_coords(3,i),&
                            species,deptype,error=internal_error)
          if (internal_error/=0) then
            write(0,"(a)") "*** Error adding initial single adatom in create_system()"
            if (present(error)) then
              error=internal_error
            endif
            return
          endif
        case(INITIAL_ATOMS_NONE)
          ! Do nothing
        case default
          write(0,"(a,i0,a)") "*** Error: invalid initial_atoms argument ",initial_atoms(i)," in create_system()"
          if (present(error)) then
            error=1
          endif
          return
      end select
    end do
    write(6,"(a,i0,a)") "# Generated ",system%n_deposited," initial atoms on substrate"

  end subroutine create_system

  subroutine add_site(system,coords,error)
    type(t_system),intent(inout) :: system
    real(kind=rk),intent(in) :: coords(3)
    integer,optional,intent(out) :: error

    integer :: internal_error,my_box(3),c

    if (present(error)) then
      error=0
    endif

    if (system%n_sites>=size(system%sites)) then
      call size_up(system%sites,error=internal_error)
      if (internal_error/=0) then
        write(0,"(a)") "*** Error sizing up sites in add_site()"
        if (present(error)) then
          error=internal_error
        endif
        return
      endif
    endif
    system%n_sites=system%n_sites+1
    system%sites(system%n_sites)%cartesian_coords=coords

    my_box(:)=int(coords(:)/system%one_box_dim(:))+1
    do c=1,3
      if (my_box(c)>system%n_boxes(c)) then
        my_box(c)=system%n_boxes(c)
      endif
    end do
    if (system%n_sites_in_boxes(my_box(1),my_box(2),my_box(3))>=size(system%sites_in_boxes,1)) then
      call size_up(system%sites_in_boxes,default_value=0,error=internal_error)
      if (internal_error/=0) then
        write(0,"(a)") "*** Error sizing up sites_in_boxes in add_site()"
        if (present(error)) then
          error=internal_error
        endif
        return
      endif
    endif
    system%n_sites_in_boxes(my_box(1),my_box(2),my_box(3))=system%n_sites_in_boxes(my_box(1),my_box(2),my_box(3))+1
    system%sites_in_boxes(system%n_sites_in_boxes(my_box(1),my_box(2),my_box(3)),my_box(1),my_box(2),my_box(3))=system%n_sites
    system%sites(system%n_sites)%my_box=my_box
    system%sites(system%n_sites)%my_index_in_box=system%n_sites_in_boxes(my_box(1),my_box(2),my_box(3))

  end subroutine add_site

  subroutine create_neighbour_list(system,i,error)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: i
    integer,optional,intent(out) :: error

    integer :: j,id,dist_index,internal_error,my_box(3),ii,jj,kk,wrapped_ii,wrapped_jj,wrapped_kk

    if (present(error)) then
      error=0
    endif

    my_box=system%sites(i)%my_box

    do kk=my_box(3)-1,my_box(3)+1
      wrapped_kk=kk
      if (wrapped_kk<1) then
        if (system%periodic(3)) then
          wrapped_kk=wrapped_kk+system%n_boxes(3)
        else
          cycle
        endif
      endif
      if (wrapped_kk>system%n_boxes(3)) then
        if (system%periodic(3)) then
          wrapped_kk=wrapped_kk-system%n_boxes(3)
        else
          cycle
        endif
      endif
      do jj=my_box(2)-1,my_box(2)+1
        wrapped_jj=jj
        if (wrapped_jj<1) then
          if (system%periodic(2)) then
            wrapped_jj=wrapped_jj+system%n_boxes(2)
          else
            cycle
          endif
        endif
        if (wrapped_jj>system%n_boxes(2)) then
          if (system%periodic(2)) then
            wrapped_jj=wrapped_jj-system%n_boxes(2)
          else
            cycle
          endif
        endif
        do ii=my_box(1)-1,my_box(1)+1
          wrapped_ii=ii
          if (wrapped_ii<1) then
            if (system%periodic(1)) then
              wrapped_ii=wrapped_ii+system%n_boxes(1)
            else
              cycle
            endif
          endif
          if (wrapped_ii>system%n_boxes(1)) then
            if (system%periodic(1)) then
              wrapped_ii=wrapped_ii-system%n_boxes(1)
            else
              cycle
            endif
          endif
          do j=1,system%n_sites_in_boxes(wrapped_ii,wrapped_jj,wrapped_kk)
            id=system%sites_in_boxes(j,wrapped_ii,wrapped_jj,wrapped_kk)
            if (id==i) then
              cycle
            endif
            dist_index=get_dist_index(system,id,i,max_shell=system%lae_cutoff_shell,error=internal_error)
            if (internal_error/=0) then
              write(0,"(a,i0,a,i0,a)") "*** Error: couldn't get distance index for pair ",i,", ",id," in create_neighbour_list()"
              if (present(error)) then
                error=internal_error
              endif
              return
            endif
            if (dist_index>0) then
              call add_neighbours(system,i,id,dist_index,error=internal_error)
              if (internal_error/=0) then
                write(0,"(a,i0,a,i0,a,i0,a)") "*** Error connecting sites ",i," and ",id," at distance index ",dist_index,&
                                              " in create_neighbour_list()"
                if (present(error)) then
                  error=internal_error
                endif
                return
              endif
            endif
          end do
        end do
      end do
    end do

  end subroutine create_neighbour_list

  subroutine update_nn(system,i)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: i

    integer :: j,shell,nn_site

    system%sites(i)%n_near=0
    do shell=1,TRUE_1NN_SHELL-1
      do j=1,system%sites(i)%neighbour_shells(shell)%n_neighbour_sites
        nn_site=system%sites(i)%neighbour_shells(shell)%neighbour_sites(j)
        if (system%sites(nn_site)%Z>0) then
          system%sites(i)%n_near=system%sites(i)%n_near+1
        endif
      end do
    end do

    system%sites(i)%nn=0
    do j=1,system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
      nn_site=system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(j)
      if (system%sites(nn_site)%Z>0) then
        system%sites(i)%nn=system%sites(i)%nn+1
        if (system%sites(nn_site)%substrate) then
          system%sites(i)%nn_substrate=system%sites(i)%nn_substrate+1
        endif
      endif
    end do

  end subroutine update_nn

  subroutine add_neighbours(system,i,j,dist_index,error)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: i,j,dist_index
    integer,optional,intent(out) :: error

    integer :: internal_error,jump_shell

    if (present(error)) then
      error=0
    endif

    call initialise_neighbour_shell(system,i,dist_index)
    call initialise_neighbour_shell(system,j,dist_index)

    if (system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites>=&
        size(system%sites(i)%neighbour_shells(dist_index)%neighbour_sites)) then
      call size_up_neighbour_shell(system,i,dist_index,error=internal_error)
      if (internal_error/=0) then
        write(0,"(a,i0,a,i0,a)") "*** Error sizing up neighbour data of shell ",dist_index," of site ",i," in add_neighbours()"
        if (present(error)) then
          error=internal_error
        endif
        return
      endif
    endif

    if (system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites>=&
        size(system%sites(j)%neighbour_shells(dist_index)%neighbour_sites)) then
      call size_up_neighbour_shell(system,j,dist_index,error=internal_error)
      if (internal_error/=0) then
        write(0,"(a,i0,a,i0,a)") "*** Error sizing up neighbour data of shell ",dist_index," of site ",j," in add_neighbours()"
        if (present(error)) then
          error=internal_error
        endif
        return
      endif
    endif

    system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites=&
      system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites+1
    system%sites(i)%neighbour_shells(dist_index)%neighbour_sites(&
      system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites)=j
    if (.not. system%sites(i)%substrate) then
      jump_shell=NN_SHELL_TO_JUMP_SHELL(dist_index)
      if (jump_shell>0) then
        if (system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites>system%sites(i)%max_nn) then
          system%sites(i)%max_nn=system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites
          if (system%sites(i)%max_nn>system%max_nn) then
            system%max_nn=system%sites(i)%max_nn
          endif
          if (allocated(system%sites(i)%hop_states)) then
            if (system%sites(i)%max_nn>size(system%sites(i)%hop_states,1)) then
              call size_up_hop_events(system,i,error=internal_error)
              if (internal_error/=0) then
                write(0,"(a,i0,a)") "*** Error sizing up hop events of site ",i," in add_neighbours()"
                if (present(error)) then
                  error=internal_error
                endif
                return
              endif
            endif
          endif
        endif
        if (jump_shell==TRUE_1NN_JUMP_SHELL) then
          if (system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites>system%max_1nn) then
            system%max_1nn=system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites
          endif
          if (system%allow_exchange) then
            if (allocated(system%sites(i)%exchange_states)) then
              if (system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites>size(system%sites(i)%exchange_states,1)) then
                call size_up_exchange_events(system,i,error=internal_error)
                if (internal_error/=0) then
                  write(0,"(a,i0,a)") "*** Error sizing up exchange events of site ",i," in add_neighbours()"
                  if (present(error)) then
                    error=internal_error
                  endif
                  return
                endif
              endif
            endif
          endif
        endif
        system%sites(i)%jumps_checked=.false.
        if (allocated(system%sites(i)%hop_states)) then
          system%sites(i)%hop_states(system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites,jump_shell)=&
            EVENT_LAE_CHANGED
          system%sites(i)%hop_barriers(:,system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites,jump_shell)=&
            100.0_rk
          system%sites(i)%hop_rates(system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites,jump_shell)=&
            0.0_rk
        endif
        if (system%allow_exchange) then
          if (jump_shell==TRUE_1NN_JUMP_SHELL) then
            if (allocated(system%sites(i)%exchange_states)) then
              system%sites(i)%exchange_states(1:system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites,&
                system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites)=EVENT_LAE_CHANGED
              system%sites(i)%exchange_barriers(:,:,system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites)=&
                100.0_rk
              system%sites(i)%exchange_rates(:,system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites)=&
                0.0_rk
              system%sites(i)%exchange_states(system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites,&
                1:system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites)=EVENT_LAE_CHANGED
              system%sites(i)%exchange_barriers(:,system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites,:)=&
                100.0_rk
              system%sites(i)%exchange_rates(system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites,:)=&
                0.0_rk
            endif
          endif
        endif
      endif
    endif

    system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites=&
      system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites+1
    system%sites(j)%neighbour_shells(dist_index)%neighbour_sites(&
      system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites)=i
    if (.not. system%sites(j)%substrate) then
      jump_shell=NN_SHELL_TO_JUMP_SHELL(dist_index)
      if (jump_shell>0) then
        if (system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites>system%sites(j)%max_nn) then
          system%sites(j)%max_nn=system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites
          if (system%sites(j)%max_nn>system%max_nn) then
            system%max_nn=system%sites(j)%max_nn
          endif
          if (allocated(system%sites(j)%hop_states)) then
            if (system%sites(j)%max_nn>size(system%sites(j)%hop_states,1)) then
              call size_up_hop_events(system,j,error=internal_error)
              if (internal_error/=0) then
                write(0,"(a,i0,a)") "*** Error sizing up hop events of site ",j," in add_neighbours()"
                if (present(error)) then
                  error=internal_error
                endif
                return
              endif
            endif
          endif
        endif
        if (jump_shell==TRUE_1NN_JUMP_SHELL) then
          if (system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites>system%max_1nn) then
            system%max_1nn=system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites
          endif
          if (system%allow_exchange) then
            if (allocated(system%sites(j)%exchange_states)) then
              if (system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites>size(system%sites(j)%exchange_states,1)) then
                call size_up_exchange_events(system,j,error=internal_error)
                if (internal_error/=0) then
                  write(0,"(a,i0,a)") "*** Error sizing up exchange events of site ",j," in add_neighbours()"
                  if (present(error)) then
                    error=internal_error
                  endif
                  return
                endif
              endif
            endif
          endif
        endif
        system%sites(j)%jumps_checked=.false.
        if (allocated(system%sites(j)%hop_states)) then
          system%sites(j)%hop_states(system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites,jump_shell)=&
            EVENT_LAE_CHANGED
          system%sites(j)%hop_barriers(:,system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites,jump_shell)=&
            100.0_rk
          system%sites(j)%hop_rates(system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites,jump_shell)=&
            0.0_rk
        endif
        if (system%allow_exchange) then
          if (jump_shell==TRUE_1NN_JUMP_SHELL) then
            if (allocated(system%sites(j)%exchange_states)) then
              system%sites(j)%exchange_states(1:system%sites(j)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites,&
                system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites)=EVENT_LAE_CHANGED
              system%sites(j)%exchange_barriers(:,:,system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites)=&
                100.0_rk
              system%sites(j)%exchange_rates(:,system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites)=&
                0.0_rk
              system%sites(j)%exchange_states(system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites,&
                1:system%sites(j)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites)=EVENT_LAE_CHANGED
              system%sites(j)%exchange_barriers(:,system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites,:)=&
                100.0_rk
              system%sites(j)%exchange_rates(system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites,:)=&
                0.0_rk
            endif
          endif
        endif
      endif
    endif

    system%sites(i)%neighbour_shells(dist_index)%reverse_indices(&
      system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites)=&
      system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites
    system%sites(j)%neighbour_shells(dist_index)%reverse_indices(&
      system%sites(j)%neighbour_shells(dist_index)%n_neighbour_sites)=&
      system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites

  end subroutine add_neighbours

  subroutine size_up_neighbour_shell(system,i,dist_index,error)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: i,dist_index
    integer,optional,intent(out) :: error

    integer :: internal_error

    if (present(error)) then
      error=0
    endif

    call size_up(system%sites(i)%neighbour_shells(dist_index)%neighbour_sites,default_value=0,&
                 custom_buffer_size=NN_BUFFER_SIZE,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a,i0,a)") "*** Error sizing up the neighbour list shell ",dist_index,&
                               " of site ",i," in size_up_neighbour_shell()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    call size_up(system%sites(i)%neighbour_shells(dist_index)%reverse_indices,default_value=0,&
                 custom_buffer_size=NN_BUFFER_SIZE,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a,i0,a)") "*** Error sizing up the reverse index shell ",dist_index,&
                               " of site ",i," in size_up_neighbour_shell()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

  end subroutine size_up_neighbour_shell

  subroutine size_up_hop_events(system,i,error)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: i
    integer,optional,intent(out) :: error

    integer :: internal_error

    if (present(error)) then
      error=0
    endif
    
    call size_up(system%sites(i)%hop_states,rank=RANK_FIRST,default_value=EVENT_EMPTY,&
                 custom_buffer_size=NN_BUFFER_SIZE,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error sizing up the hop states of site ",i," in size_up_hop_events()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    call size_up(system%sites(i)%hop_barriers,rank=RANK_SECOND,default_value=100.0_rk,&
                 custom_buffer_size=NN_BUFFER_SIZE,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error sizing up the hop barriers of site ",i," in size_up_hop_events()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    call size_up(system%sites(i)%hop_rates,rank=RANK_FIRST,default_value=0.0d0,&
                 custom_buffer_size=NN_BUFFER_SIZE,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error sizing up the hop rates of site ",i," in size_up_hop_events()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

  end subroutine size_up_hop_events

  subroutine size_up_exchange_events(system,i,error)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: i
    integer,optional,intent(out) :: error

    integer :: internal_error

    if (present(error)) then
      error=0
    endif

    call size_up(system%sites(i)%exchange_states,rank=RANK_ALL,default_value=EVENT_EMPTY,&
                 custom_buffer_size=NN_BUFFER_SIZE,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error sizing up the hop states of site ",i," in size_up_exchange_events()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    call size_up(system%sites(i)%exchange_barriers,rank=RANK_LAST_TWO,default_value=100.0_rk,&
                 custom_buffer_size=NN_BUFFER_SIZE,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error sizing up the hop barriers of site ",i," in size_up_exchange_events()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    call size_up(system%sites(i)%exchange_rates,rank=RANK_ALL,default_value=0.0d0,&
                 custom_buffer_size=NN_BUFFER_SIZE,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error sizing up the hop rates of site ",i," in size_up_exchange_events()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

  end subroutine size_up_exchange_events

  subroutine initialise_neighbour_shell(system,i,dist_index)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: i,dist_index

    if (system%sites(i)%neighbour_shells(dist_index)%initialised) then
      return
    endif

    allocate(system%sites(i)%neighbour_shells(dist_index)%neighbour_sites(NN_BUFFER_SIZE))
    system%sites(i)%neighbour_shells(dist_index)%neighbour_sites=0
    allocate(system%sites(i)%neighbour_shells(dist_index)%reverse_indices(NN_BUFFER_SIZE))
    system%sites(i)%neighbour_shells(dist_index)%reverse_indices=0
    system%sites(i)%neighbour_shells(dist_index)%initialised=.true.
  end subroutine initialise_neighbour_shell

  subroutine add_adsorption_sites(system,error)
    type(t_system),intent(inout) :: system
    integer,optional,intent(out) :: error

    integer :: i,nn1,nn2,nn_site1,nn_site2,internal_error
    real(kind=rk) :: triplet(3,3),diff(3),centre(3)
    real(kind=rk) :: new_site(3)
    logical,allocatable :: checked(:)

    if (present(error)) then
      error=0
    endif

    allocate(checked(system%n_sites))
    checked=.false.
    do i=1,system%n_sites
      if (system%sites(i)%nn>=N_NN_MAX) then
        cycle
      endif
      triplet(:,1)=system%sites(i)%cartesian_coords(:)
      do nn1=1,system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites-1
        nn_site1=system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(nn1)
        if (system%sites(nn_site1)%nn>=N_NN_MAX) then
          cycle
        endif
        if (system%sites(nn_site1)%Z<=0) then
          cycle
        endif
        if (checked(nn_site1)) then
          cycle
        endif
        diff=system_diff_min_image(system,i,nn_site1)
        triplet(:,2)=system%sites(i)%cartesian_coords(:)+diff
        do nn2=nn1+1,system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
          nn_site2=system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(nn2)
          if (system%sites(nn_site2)%nn>=N_NN_MAX) then
            cycle
          endif
          if (system%sites(nn_site2)%Z<=0) then
            cycle
          endif
          if (checked(nn_site2)) then
            cycle
          endif
          if (.not. are_1nn_neighbours(system,nn_site1,nn_site2)) then
            cycle
          endif

          diff=system_diff_min_image(system,i,nn_site2)
          triplet(:,3)=system%sites(i)%cartesian_coords(:)+diff
          centre=geometric_centre(triplet)
          new_site=centre
          new_site(3)=new_site(3)+system%layer_separation_111
          new_site=wrap_to_cell(system%lattice,new_site,system%periodic)
          call add_site(system,new_site)
          system%sites(system%n_sites)%Z=0
          system%sites(system%n_sites)%initial_site=.true.
          call create_neighbour_list(system,system%n_sites,error=internal_error)
          if (internal_error/=0) then
            write(0,"(a,i0,a)") "*** Error creating neighbour list for site ",system%n_sites," in add_adsorption_sites()"
            if (present(error)) then
              error=internal_error
            endif
            return
          endif
          call update_nn(system,system%n_sites)
        end do
      end do
      checked(i)=.true.
    end do
    deallocate(checked)

  end subroutine add_adsorption_sites

  subroutine find_new_sites(system,i,error)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: i
    integer,optional,intent(out) :: error

    integer :: j,k,nn1,nn2,dir,nn_site,found_nn,internal_error
    integer :: n_triplets,nn_dist,found_sites_about_triplet
    integer :: triplets(2,(N_NN_MAX-1)*(N_NN_MAX-1))
    real(kind=rk) :: triplet(3,3),diff(3),centre(3),splane(3)
    real(kind=rk) :: new_site(3),mirror_coords(3)

    if (present(error)) then
      error=0
    endif

    if (system%sites(i)%nn>=N_NN_MAX) then
      return
    endif
    call find_triplets(system,i,triplets,n_triplets)
    do j=1,n_triplets
      nn1=triplets(1,j)
      if (system%sites(nn1)%nn>=N_NN_MAX) then
        cycle
      endif
      nn2=triplets(2,j)
      if (system%sites(nn2)%nn>=N_NN_MAX) then
        cycle
      endif
      if (are_1nn_neighbours(system,nn1,nn2)) then
        nn_dist=1
      elseif (are_2nn_neighbours(system,nn1,nn2)) then
        nn_dist=2
      else
        cycle
      endif
      found_sites_about_triplet=0
      do k=1,system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
        nn_site=system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(k)
        if (are_1nn_neighbours(system,nn_site,nn1)) then
          if (are_1nn_neighbours(system,nn_site,nn2)) then
            found_sites_about_triplet=found_sites_about_triplet+1
            found_nn=nn_site
            if (found_sites_about_triplet>=2) then
              exit
            endif
          endif
        endif
      end do
      if (found_sites_about_triplet>=2) then
        cycle
      endif

      triplet(:,1)=system%sites(i)%cartesian_coords(:)
      diff=system_diff_min_image(system,i,nn1)
      triplet(:,2)=system%sites(i)%cartesian_coords(:)+diff
      diff=system_diff_min_image(system,i,nn2)
      triplet(:,3)=system%sites(i)%cartesian_coords(:)+diff
      if (nn_dist==1) then
        centre=geometric_centre(triplet)
      else ! nn_dist==2
        centre=0.5_rk*(triplet(:,2)+triplet(:,3))
      endif
      if (found_sites_about_triplet==1) then
        mirror_coords=system%sites(i)%cartesian_coords(:)+system_diff_min_image(system,i,found_nn)
        diff=centre-mirror_coords
        new_site=wrap_to_cell(system%lattice,centre+diff,system%periodic)
        call add_site(system,new_site)
        call create_neighbour_list(system,system%n_sites,error=internal_error)
        if (internal_error/=0) then
          write(0,"(a,i0,a)") "*** Error creating neighbour list for site ",system%n_sites," in find_new_sites()"
          if (present(error)) then
            error=internal_error
          endif
          return
        endif
        call update_nn(system,system%n_sites)
        if (nn_dist==1) then
          system%sites(system%n_sites)%Z=0
        elseif (has_support(system,system%n_sites)) then
          system%sites(system%n_sites)%Z=0
        endif
      elseif (found_sites_about_triplet==0) then
        if (nn_dist==1) then
          splane=denormal_surface_plane(triplet)/system%lattice_constant/0.75_rk
        else ! nn_dist==2
          splane=denormal_surface_plane(triplet)/system%lattice_constant
        endif
        do dir=-1,1,2
          new_site=wrap_to_cell(system%lattice,centre+dir*splane,system%periodic)
          call add_site(system,new_site)
          call create_neighbour_list(system,system%n_sites,error=internal_error)
          if (internal_error/=0) then
            write(0,"(a,i0,a)") "*** Error creating neighbour list for site ",system%n_sites," in find_new_sites()"
            if (present(error)) then
              error=internal_error
            endif
            return
          endif
          call update_nn(system,system%n_sites)
          if (nn_dist==1) then
            system%sites(system%n_sites)%Z=0
          elseif (has_support(system,system%n_sites)) then
            system%sites(system%n_sites)%Z=0
          endif
        end do
      endif
    end do

  end subroutine find_new_sites

  subroutine find_triplets(system,i,triplets,n_triplets)
    type(t_system),intent(in) :: system
    integer,intent(in) :: i
    integer,intent(out) :: triplets(2,(N_NN_MAX-1)*(N_NN_MAX-1)),n_triplets

    integer :: nn1,nn2,nn_site1,nn_site2
    
    n_triplets=0
    do nn1=1,system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites-1
      nn_site1=system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(nn1)
      if (system%sites(nn_site1)%substrate) then
        cycle
      endif
      if (system%sites(nn_site1)%Z<=0) then
        cycle
      endif
      do nn2=nn1+1,system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
        nn_site2=system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(nn2)
        if (system%sites(nn_site2)%substrate) then
          cycle
        endif
        if (system%sites(nn_site2)%Z<=0) then
          cycle
        endif
        n_triplets=n_triplets+1
        triplets(1,n_triplets)=nn_site1
        triplets(2,n_triplets)=nn_site2
      end do
    end do
  end subroutine find_triplets

  subroutine print_xyz_frame(system,time,step,fileunit,debug,hard_debug,error)
    type(t_system),intent(in) :: system
    double precision,intent(in) :: time
    integer,intent(in) :: step,fileunit
    logical,optional,intent(in) :: debug,hard_debug
    integer,optional,intent(out) :: error

    integer :: i,site_i,ad_site_id,my_error
    logical :: my_debug,my_hard_debug
    
    if (present(error)) then
      error=0
    endif
    my_error=0

    my_debug=.false.
    if (present(debug)) then
      my_debug=debug
    endif

    my_hard_debug=.false.
    if (present(hard_debug)) then
      my_hard_debug=hard_debug
      if (my_hard_debug) then
        my_debug=.true.
      endif
    endif

    if (my_debug) then
      if (my_hard_debug) then
        write(fileunit,"(i0)",iostat=my_error) system%n_sites
      else
        write(fileunit,"(i0)",iostat=my_error) system%n_substrate+system%n_deposited+count(system%sites(1:system%n_sites)%Z==0)
      endif
    else
      write(fileunit,"(i0)",iostat=my_error) system%n_substrate+system%n_deposited
    endif
    if (my_error/=0) then
      if (present(error)) then
        error=my_error
      endif
      write(0,"(a)") "*** Error writing the number of particles in the xyz file in print_xyz_frame()"
      return
    endif
    if (my_debug) then
      if (my_hard_debug) then
        write(fileunit,"(a,f10.6,a,f10.6,a,f10.6,a,g20.10,a,i0,a)",iostat=my_error) "Lattice=""",&
                         system%lattice(1,1)," 0.0 0.0 0.0 ",&
                         system%lattice(2,2)," 0.0 0.0 0.0 ",system%lattice(3,3),&
                         """ Properties=species:S:1:pos:R:3:id:I:1:n_neighb:I:1:n_near:I:1:"//&
                         "site_id:I:1:Z:I:1:substrate:L:1:initial_site:L:1:inactive_counter:I:1:locked:L:1 Time=""",&
                         time,""" Step=""",step,""""
      else
        write(fileunit,"(a,f10.6,a,f10.6,a,f10.6,a,g20.10,a,i0,a)",iostat=my_error) "Lattice=""",&
                                   system%lattice(1,1)," 0.0 0.0 0.0 ",&
                                   system%lattice(2,2)," 0.0 0.0 0.0 ",system%lattice(3,3),&
                                   """ Properties=species:S:1:pos:R:3:id:I:1:n_neighb:I:1:n_near:I:1:site_id:I:1 Time=""",&
                                   time,""" Step=""",step,""""
      endif
    else
      write(fileunit,"(a,f10.6,a,f10.6,a,f10.6,a,g20.10,a,i0,a)",iostat=my_error) "Lattice=""",system%lattice(1,1)," 0.0 0.0 0.0 ",&
                                   system%lattice(2,2)," 0.0 0.0 0.0 ",system%lattice(3,3),&
                                   """ Properties=species:S:1:pos:R:3:id:I:1 Time=""",time,""" Step=""",step,""""
    endif
    if (my_error/=0) then
      if (present(error)) then
        error=my_error
      endif
      write(0,"(a)") "*** Error writing the comment line in the xyz file in print_xyz_frame()"
      return
    endif
    if (my_debug) then
      ad_site_id=system%n_substrate+system%n_deposited
      do i=1,system%n_sites
        if (my_hard_debug) then
          ad_site_id=ad_site_id+1
          write(fileunit,"(a,3(f15.6),5(1x,i0),2(1x,L1),1x,i0,1x,L1)",iostat=my_error) ElementName(max(0,system%sites(i)%Z)),&
                                                    system%sites(i)%cartesian_coords(:),&
                                                    ad_site_id,&
                                                    system%sites(i)%nn,&
                                                    system%sites(i)%n_near,&
                                                    i,system%sites(i)%Z,&
                                                    system%sites(i)%substrate,system%sites(i)%initial_site,&
                                                    system%sites(i)%inactive_counter,&
                                                    system%sites(i)%locked
        else
          if (system%sites(i)%Z>0) then
            write(fileunit,"(a,3(f15.6),4(1x,i0))",iostat=my_error) ElementName(system%sites(i)%Z),&
                                                      system%sites(i)%cartesian_coords(:),&
                                                      system%sites(i)%id,&
                                                      system%sites(i)%nn,&
                                                      system%sites(i)%n_near,&
                                                      i
          elseif (system%sites(i)%Z==0) then
            ad_site_id=ad_site_id+1
            write(fileunit,"(a,3(f15.6),4(1x,i0))",iostat=my_error) ElementName(max(0,system%sites(i)%Z)),&
                                                      system%sites(i)%cartesian_coords(:),&
                                                      ad_site_id,&
                                                      system%sites(i)%nn,&
                                                      system%sites(i)%n_near,&
                                                      i
          endif
        endif
      end do
    else
      do i=system%frozen_substrate%first_id,system%frozen_substrate%last_id
        write(fileunit,"(a,3(f15.6),1x,i0)",iostat=my_error) ElementName(system%frozen_substrate%Z(i)),&
                                                             system%frozen_substrate%pos(:,i),i
        if (my_error/=0) then
          if (present(error)) then
            error=my_error
          endif
          write(0,"(a,i0,a)") "*** Error writing frozen substrate atom ",i," in the xyz file in print_xyz_frame()"
          return
        endif
      end do
      do i=system%substrate%first_id,system%substrate%last_id
        write(fileunit,"(a,3(f15.6),1x,i0)",iostat=my_error) ElementName(system%substrate%Z(i)),&
                                                             system%substrate%pos(:,i),i
        if (my_error/=0) then
          if (present(error)) then
            error=my_error
          endif
          write(0,"(a,i0,a)") "*** Error writing substrate atom ",i," in the xyz file in print_xyz_frame()"
          return
        endif
      end do

      do i=1,system%n_deposited
        site_i=system%deposited_atoms(i)
        write(fileunit,"(a,3(f15.6),1x,i0)",iostat=my_error) ElementName(system%sites(site_i)%Z),&
                                                          system%sites(site_i)%cartesian_coords(:),&
                                                          system%sites(site_i)%id
        if (my_error/=0) then
          if (present(error)) then
            error=my_error
          endif
          write(0,"(a,i0,a)") "*** Error writing deposited atom ",i," in the xyz file in print_xyz_frame()"
          return
        endif
      end do
    endif

  end subroutine print_xyz_frame

  function get_cartesian(x,y,z,p)
    integer,intent(in) :: x,y,z,p
    real(kind=rk),dimension(3) :: get_cartesian

    get_cartesian(1)=x*UNIT_CELL(1)+P_BASE(1,p)
    get_cartesian(2)=y*UNIT_CELL(2)+P_BASE(2,p)
    get_cartesian(3)=z*UNIT_CELL(3)+P_BASE(3,p)

  end function get_cartesian

  type(t_atoms) function get_atoms(system)
    type(t_system),intent(in) ::system

    integer :: i,site_id,first_id,last_id

    get_atoms%lattice=system%lattice
    get_atoms%periodic=system%periodic
    first_id=system%n_substrate+1
    last_id=system%n_substrate+system%n_deposited
    get_atoms%first_id=first_id
    get_atoms%last_id=last_id
    if (last_id<first_id) then
      return
    endif
    allocate(get_atoms%pos(3,first_id:last_id))
    allocate(get_atoms%atom_type(first_id:last_id))
    allocate(get_atoms%Z(first_id:last_id))

    do i=1,system%n_deposited
      site_id=system%deposited_atoms(i)
      get_atoms%pos(:,system%sites(site_id)%id)=system%sites(site_id)%cartesian_coords(:)
      get_atoms%atom_type(system%sites(site_id)%id)=system%sites(site_id)%atom_type
      get_atoms%Z(system%sites(site_id)%id)=system%sites(site_id)%Z
    end do
    get_atoms%initialised=.true.
  end function get_atoms

  subroutine dealloc_atoms(at)
    type(t_atoms),intent(inout) :: at

    if (at%initialised) then
      at%first_id=0
      at%last_id=0
      deallocate(at%pos,at%atom_type,at%Z)
      at%initialised=.false.
    endif
  end subroutine dealloc_atoms

  real(kind=rk) function system_dist2_min_image(system,i,j)
    type(t_system),intent(in) :: system
    integer,intent(in) :: i,j

    system_dist2_min_image=normsq(real(system_diff_min_image(system,i,j),kind=dp))
  end function system_dist2_min_image

  function system_diff_min_image(system,i,j)
    type(t_system),intent(in) :: system
    integer,intent(in) :: i,j
    real(kind=rk) :: system_diff_min_image(3)

    system_diff_min_image=atoms_diff_min_image(system%lattice,system%sites(i)%cartesian_coords,&
                                               system%sites(j)%cartesian_coords,system%periodic)
  end function system_diff_min_image

  real(kind=rk) function atoms_dist2_min_image(lattice,x1,x2,periodic)
    real(kind=rk),intent(in) :: lattice(3,3),x1(3),x2(3)
    logical,intent(in) :: periodic(3)

    atoms_dist2_min_image=normsq(real(atoms_diff_min_image(lattice,x1,x2,periodic),kind=dp))
  end function atoms_dist2_min_image

  function atoms_diff_min_image(lattice,x1,x2,periodic)
    real(kind=rk),intent(in) :: lattice(3,3),x1(3),x2(3)
    logical,intent(in) :: periodic(3)
    real(kind=rk) :: atoms_diff_min_image(3)

    integer :: c

    atoms_diff_min_image=x2-x1
    do c=1,3
      if (periodic(c)) then
        if (atoms_diff_min_image(c)>0.5_rk*lattice(c,c)) then
          atoms_diff_min_image(c)=atoms_diff_min_image(c)-lattice(c,c)
        endif
        if (atoms_diff_min_image(c)<-0.5_rk*lattice(c,c)) then
          atoms_diff_min_image(c)=atoms_diff_min_image(c)+lattice(c,c)
        endif
      endif
    end do
  end function atoms_diff_min_image

  function wrap_to_cell(lattice,x,periodic)
    real(kind=rk),intent(in) :: lattice(3,3),x(3)
    logical,intent(in) :: periodic(3)
    real(kind=rk) :: wrap_to_cell(3)

    integer :: c

    wrap_to_cell=x
    do c=1,3
      if (periodic(c)) then
        if (wrap_to_cell(c)>lattice(c,c)) then
          wrap_to_cell(c)=wrap_to_cell(c)-lattice(c,c)
        endif
        if (wrap_to_cell(c)<0.0d0) then
          wrap_to_cell(c)=wrap_to_cell(c)+lattice(c,c)
        endif
      endif
    end do
  end function wrap_to_cell

  subroutine check_lockedness(system,i,error)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: i
    integer,optional,intent(out) :: error

    integer :: nn_atom,nn_nn_atom,nn_id,nn_nn_id,n_substrate_attached,internal_error

    if (present(error)) then
      error=0
    endif
    internal_error=0

    if (system%sites(i)%Z<=0 .or. system%sites(i)%substrate) then
      return
    endif

    system%sites(i)%locked=.false.
    if (system%sites(i)%nn>=lock_threshold) then
      system%sites(i)%locked=.true.
      return
    endif

    do nn_atom=1,system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
      nn_id=system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(nn_atom)
      if (system%sites(nn_id)%substrate) then
        cycle
      endif
      if (system%sites(nn_id)%Z>0) then
        if (system%sites(nn_id)%nn<=3) then
          system%sites(i)%locked=.true.
          return
        endif
        ! If site i is substrate-attached and it supports non-substrate-attached sites that have up to three
        ! substrate-attached supports, mark it locked
        if (system%sites(i)%nn_substrate>0) then
          if (system%sites(nn_id)%nn_substrate<=0) then
            n_substrate_attached=0
            do nn_nn_atom=1,system%sites(nn_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
              nn_nn_id=system%sites(nn_id)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(nn_nn_atom)
              if (system%sites(nn_nn_id)%Z<=0) then
                cycle
              endif
              if (system%sites(nn_nn_id)%nn_substrate>0) then
                n_substrate_attached=n_substrate_attached+1
              endif
            end do
            if (n_substrate_attached<=3) then
              system%sites(i)%locked=.true.
              return
            endif
          endif
        endif
      endif
    end do

  end subroutine check_lockedness

  subroutine do_lockedness_checks(system,i,error)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: i
    integer,optional,intent(out) :: error

    integer :: nn_atom,nn_id
    integer :: nn_nn_atom,nn_nn_id
    integer :: internal_error

    if (present(error)) then
      error=0
    endif

    do nn_atom=1,system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
      nn_id=system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(nn_atom)
      if (system%sites(nn_id)%Z>0) then
        if (.not. system%sites(nn_id)%lock_checked) then
          call check_lockedness(system,nn_id,error=internal_error)
          if (internal_error/=0) then
            write(0,"(a,i0,a,i0,a)") "*** Error checking lockedness of site ",&
                                      nn_id," in do_lockedness_checks()"
            if (present(error)) then
              error=internal_error
            endif
            return
          endif
          system%sites(nn_id)%lock_checked=.true.
        endif
        do nn_nn_atom=1,system%sites(nn_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
          nn_nn_id=system%sites(nn_id)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(nn_nn_atom)
          if (system%sites(nn_nn_id)%Z>0) then
            if (.not. system%sites(nn_nn_id)%lock_checked) then
              call check_lockedness(system,nn_nn_id,error=internal_error)
              if (internal_error/=0) then
                write(0,"(a,i0,a,i0,a)") "*** Error checking lockedness of site ",&
                                         nn_nn_id," in do_lockedness_checks()"
                if (present(error)) then
                  error=internal_error
                endif
                return
              endif
              system%sites(nn_nn_id)%lock_checked=.true.
            endif
          endif
        end do
      endif
    end do

  end subroutine do_lockedness_checks

  integer function get_dist_index(system,i,j,max_shell,error)
    type(t_system),intent(in) :: system
    integer,intent(in) :: i,j
    integer,optional,intent(in) :: max_shell
    integer,intent(out),optional :: error

    integer :: my_max_shell,internal_error
    real(kind=rk) :: d2,d2overl2

    if (present(error)) then
      error=0
    endif
    internal_error=0
    get_dist_index=0

    d2=system_dist2_min_image(system,i,j)
    if (d2<COLLISION2_TOL) then
      write(0,"(a,i0,a,i0,a)") "*** Error: sites ",i," and ",j," collide in get_dist_index()"
      if (present(error)) then
        error=1
      endif
      return
    endif
    d2overl2=d2/system%lattice_constant/system%lattice_constant

    if (present(max_shell)) then
      my_max_shell=max_shell
    else
      my_max_shell=NN_SHELLS
    endif

    ! Check if the neighbour is further than max distance
    if (d2overl2>JUMP_DIST2_INDEX(my_max_shell)+DIST2_TOL) then
      get_dist_index=0
      return
    endif

    internal_error=1
    do get_dist_index=2,my_max_shell,2
      if (d2overl2<=JUMP_DIST2_INDEX(get_dist_index)-DIST2_TOL) then
        ! Hits the previous mid-shell, handled at the end
        internal_error=0
        exit
      endif
      if (d2overl2>JUMP_DIST2_INDEX(get_dist_index)-DIST2_TOL .and. &
          d2overl2<=JUMP_DIST2_INDEX(get_dist_index)+DIST2_TOL) then
          return
      endif
    end do

    if (internal_error==0) then
      get_dist_index=get_dist_index-1
    else
      write(0,"(a,i0,a,i0,a)") "*** Error getting dist_index for pair ",i,", ",j," in get_dist_index()"
      if (present(error)) then
        error=internal_error
      endif
    endif

  end function get_dist_index

  logical function are_1nn_neighbours_trimer(system,i,j,k)
    type(t_system),intent(in) :: system
    integer,intent(in) :: i,j,k

    are_1nn_neighbours_trimer=.false.
    if (.not. are_1nn_neighbours_dimer(system,i,j)) then
      return
    endif
    if (.not. are_1nn_neighbours_dimer(system,i,k)) then
      return
    endif
    if (.not. are_1nn_neighbours_dimer(system,j,k)) then
      return
    endif
    are_1nn_neighbours_trimer=.true.

  end function are_1nn_neighbours_trimer

  logical function are_1nn_neighbours_dimer(system,i,j)
    type(t_system),intent(in) :: system
    integer,intent(in) :: i,j

    integer :: nn_atom

    are_1nn_neighbours_dimer=.false.
    do nn_atom=1,system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
      if (j==system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(nn_atom)) then
        are_1nn_neighbours_dimer=.true.
        return
      endif
    end do

  end function are_1nn_neighbours_dimer

  logical function are_2nn_neighbours(system,i,j)
    type(t_system),intent(in) :: system
    integer,intent(in) :: i,j

    integer :: nn_atom

    are_2nn_neighbours=.false.
    do nn_atom=1,system%sites(i)%neighbour_shells(TRUE_2NN_SHELL)%n_neighbour_sites
      if (j==system%sites(i)%neighbour_shells(TRUE_2NN_SHELL)%neighbour_sites(nn_atom)) then
        are_2nn_neighbours=.true.
        return
      endif
    end do

  end function are_2nn_neighbours

  logical function are_closer_than(system,i,j,max_shell)
    type(t_system),intent(in) :: system
    integer,intent(in) :: i,j,max_shell

    integer :: nn_atom,shell

    are_closer_than=.false.
    do shell=1,max_shell-1
      do nn_atom=1,system%sites(i)%neighbour_shells(shell)%n_neighbour_sites
        if (system%sites(i)%neighbour_shells(shell)%neighbour_sites(nn_atom)==j) then
          are_closer_than=.true.
          return
        endif
      end do
    end do

  end function are_closer_than

  subroutine print_system(system,filename,error)
    type(t_system),intent(in) :: system
    character(len=*),intent(in) :: filename
    integer,optional,intent(out) :: error

    integer :: ii,jj,kk,internal_error,fileunit

    if (present(error)) then
      error=0
    endif

    open(newunit=fileunit,file=filename,action="write",iostat=internal_error,RECL=1000000)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error opening file "//trim(filename)// " for writing in print_system()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    write(fileunit,*) system%n_sites
    write(fileunit,*) system%n_substrate
    write(fileunit,*) system%n_deposited
    if (system%n_deposited>0) then
      write(fileunit,*) system%deposited_atoms(1:system%n_deposited)
    endif
    write(fileunit,*) system%lattice_constant
    write(fileunit,*) system%nn_dist2
    write(fileunit,*) system%relax_displacement_tol2
    write(fileunit,*) system%layer_separation_111
    write(fileunit,*) system%layer_separation_111_squared
    write(fileunit,*) system%max_jump_dist
    write(fileunit,*) system%lae_cutoff_dist
    write(fileunit,*) system%height
    write(fileunit,*) system%lattice
    write(fileunit,*) system%max_nn
    write(fileunit,*) system%max_1nn
    write(fileunit,*) system%max_jump_shell
    write(fileunit,*) system%lae_cutoff_shell
    write(fileunit,*) system%substrate_type
    write(fileunit,*) system%one_box_dim
    write(fileunit,*) system%n_boxes
    do kk=1,system%n_boxes(3)
      do jj=1,system%n_boxes(2)
        write(fileunit,*) system%n_sites_in_boxes(:,jj,kk)
      end do
    end do
    do kk=1,system%n_boxes(3)
      do jj=1,system%n_boxes(2)
        do ii=1,system%n_boxes(1)
          if (system%n_sites_in_boxes(ii,jj,kk)>0) then
            write(fileunit,*) system%sites_in_boxes(1:system%n_sites_in_boxes(ii,jj,kk),ii,jj,kk)
          endif
        end do
      end do
    end do
    write(fileunit,*) system%allow_exchange
    write(fileunit,*) system%periodic
    call print_sites(system,fileunit,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error writing sites to file "//trim(filename)// " in print_system()"
      if (present(error)) then
        error=internal_error
      endif
      close(fileunit)
      return
    endif
    call print_atoms(system%frozen_substrate,fileunit,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error writing frozen substrate to file "//trim(filename)// " in print_system()"
      if (present(error)) then
        error=internal_error
      endif
      close(fileunit)
      return
    endif
    call print_atoms(system%substrate,fileunit,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error writing substrate to file "//trim(filename)// " in print_system()"
      if (present(error)) then
        error=internal_error
      endif
      close(fileunit)
      return
    endif
    close(fileunit)
    
  end subroutine print_system

  subroutine print_sites(system,fileunit,error)
    type(t_system),intent(in) :: system
    integer,intent(in) :: fileunit
    integer,optional,intent(out) :: error

    integer :: i,j,internal_error

    if (present(error)) then
      error=0
    endif

    do i=1,system%n_sites
      write(fileunit,*,iostat=internal_error) system%sites(i)%id
      if (internal_error/=0) then
        write(0,"(a,i0,a)") "*** Error writing id of site ",i," in print_sites()"
        if (present(error)) then
          error=internal_error
        endif
        return
      endif
      write(fileunit,*) system%sites(i)%Z
      write(fileunit,*) system%sites(i)%atom_type
      write(fileunit,*) system%sites(i)%nn
      write(fileunit,*) system%sites(i)%nn_substrate
      write(fileunit,*) system%sites(i)%n_near
      write(fileunit,*) system%sites(i)%max_nn
      write(fileunit,*) system%sites(i)%inactive_counter
      write(fileunit,*) system%sites(i)%locked
      write(fileunit,*) system%sites(i)%lock_checked
      write(fileunit,*) system%sites(i)%jumps_checked
      write(fileunit,*) system%sites(i)%substrate
      write(fileunit,*) system%sites(i)%initial_site
      write(fileunit,*) system%sites(i)%cartesian_coords
      write(fileunit,*) system%sites(i)%my_box
      write(fileunit,*) system%sites(i)%my_index_in_box
      write(fileunit,*) system%sites(i)%neighbour_shells(:)%initialised
      do j=1,NN_SHELLS
        if (.not. system%sites(i)%neighbour_shells(j)%initialised) then
          cycle
        endif
        write(fileunit,*) system%sites(i)%neighbour_shells(j)%n_neighbour_sites
        write(fileunit,*) system%sites(i)%neighbour_shells(j)%neighbour_sites&
                          (1:system%sites(i)%neighbour_shells(j)%n_neighbour_sites)
        write(fileunit,*) system%sites(i)%neighbour_shells(j)%reverse_indices&
                          (1:system%sites(i)%neighbour_shells(j)%n_neighbour_sites)
      end do
      write(fileunit,*) system%sites(i)%hop_rate
      write(fileunit,*) system%sites(i)%exchange_rate
      if (.not. system%sites(i)%substrate .and. system%sites(i)%Z>0) then
        do j=1,N_JUMP_SHELLS
          write(fileunit,*) system%sites(i)%hop_states(:,j)
          write(fileunit,*) system%sites(i)%hop_barriers(:,:,j)
          write(fileunit,*) system%sites(i)%hop_rates(:,j)
        end do
        if (system%allow_exchange) then
          do j=1,system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
            write(fileunit,*) system%sites(i)%exchange_states(:,j)
            write(fileunit,*) system%sites(i)%exchange_barriers(:,:,j)
            write(fileunit,*) system%sites(i)%exchange_rates(:,j)
          end do
        endif
      endif
    end do
  end subroutine print_sites

  subroutine print_atoms(atoms,fileunit,error)
    type(t_atoms),intent(in) :: atoms
    integer,intent(in) :: fileunit
    integer,optional,intent(out) :: error

    integer :: i,internal_error

    if (present(error)) then
      error=0
    endif
    internal_error=0

    write(fileunit,*,iostat=internal_error) atoms%initialised
    if (internal_error/=0) then
      write(0,"(a)") "*** Error writing the initialised-status in print_atoms()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    write(fileunit,*) atoms%lattice
    write(fileunit,*) atoms%periodic
    write(fileunit,*) atoms%first_id,atoms%last_id
    do i=atoms%first_id,atoms%last_id
      write(fileunit,*) atoms%Z(i)
      write(fileunit,*) atoms%atom_type(i)
      write(fileunit,*) atoms%pos(:,i)
    end do
  end subroutine print_atoms

  subroutine print_atoms_xyz(atoms,fileunit,error)
    type(t_atoms),intent(in) :: atoms
    integer,intent(in) :: fileunit
    integer,optional,intent(out) :: error

    integer :: i,internal_error

    if (present(error)) then
      error=0
    endif
    internal_error=0

    write(fileunit,*,iostat=internal_error) atoms%last_id-atoms%first_id+1
    if (internal_error/=0) then
      write(0,"(a)") "*** Error writing number of atoms in unit ",fileunit," in print_atoms_xyz()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    write(fileunit,"(a,f10.6,a,f10.6,a,f10.6,a)",iostat=internal_error) "Lattice=""",&
                                                                        atoms%lattice(1,1)," 0.0 0.0 0.0 ",&
                                                                        atoms%lattice(2,2)," 0.0 0.0 0.0 ",atoms%lattice(3,3),&
                                                                        """ Properties=species:S:1:pos:R:3:id:I:1"
    if (internal_error/=0) then
      write(0,"(a)") "*** Error writing comment in unit ",fileunit," in print_atoms_xyz()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    do i=atoms%first_id,atoms%last_id
      write(fileunit,*) ElementName(max(0,atoms%Z(i))),atoms%pos(:,i),i
    end do
  end subroutine print_atoms_xyz

  subroutine read_system(system,filename,error)
    type(t_system),intent(out) :: system
    character(len=*),intent(in) :: filename
    integer,optional,intent(out) :: error

    integer :: ii,jj,kk,internal_error,fileunit
    character(len=1000000) :: line
    
    if (present(error)) then
      error=0
    endif

    open(newunit=fileunit,file=filename,action="read",iostat=internal_error,RECL=1000000)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error opening file "//trim(filename)// " for reading in read_system()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    read(fileunit,*) system%n_sites
    allocate(system%sites(system%n_sites))
    read(fileunit,*) system%n_substrate
    read(fileunit,*) system%n_deposited
    if (system%n_deposited>0) then
      allocate(system%deposited_atoms(system%n_deposited))
      read(fileunit,*,iostat=internal_error) system%deposited_atoms(1:system%n_deposited)
      if (internal_error/=0) then
        write(0,"(a)") "*** Error reading deposited atoms indices in read_system()"
        if (present(error)) then
          error=internal_error
        endif
        close(fileunit)
        return
      endif
    else
      allocate(system%deposited_atoms(BUFFER_SIZE))
    endif
    read(fileunit,*) system%lattice_constant
    read(fileunit,*) system%nn_dist2
    read(fileunit,*) system%relax_displacement_tol2
    read(fileunit,*) system%layer_separation_111
    read(fileunit,*) system%layer_separation_111_squared
    read(fileunit,*) system%max_jump_dist
    read(fileunit,*) system%lae_cutoff_dist
    read(fileunit,*) system%height
    read(fileunit,*) system%lattice
    read(fileunit,*) system%max_nn
    read(fileunit,*) system%max_1nn
    read(fileunit,*) system%max_jump_shell
    read(fileunit,*) system%lae_cutoff_shell
    read(fileunit,*) system%substrate_type
    read(fileunit,*) system%one_box_dim
    read(fileunit,*) system%n_boxes
    allocate(system%n_sites_in_boxes(system%n_boxes(1),system%n_boxes(2),system%n_boxes(3)))
    system%n_sites_in_boxes=0
    do kk=1,system%n_boxes(3)
      do jj=1,system%n_boxes(2)
        read(fileunit,*,iostat=internal_error) system%n_sites_in_boxes(1:system%n_boxes(1),jj,kk)
        if (internal_error/=0) then
          write(0,"(a,i0,a,i0,a)") "*** Error reading sites in boxes ",jj,", ",kk," in read_system()"
          if (present(error)) then
            error=internal_error
          endif
          close(fileunit)
          return
        endif
      end do
    end do
    allocate(system%sites_in_boxes(maxval(system%n_sites_in_boxes),system%n_boxes(1),system%n_boxes(2),system%n_boxes(3)))
    system%sites_in_boxes=0
    do kk=1,system%n_boxes(3)
      do jj=1,system%n_boxes(2)
        do ii=1,system%n_boxes(1)
          if (system%n_sites_in_boxes(ii,jj,kk)>0) then
            read(fileunit,*,iostat=internal_error) system%sites_in_boxes(1:system%n_sites_in_boxes(ii,jj,kk),ii,jj,kk)
            if (internal_error/=0) then
              write(0,"(a,i0,a,i0,a,i0,a)") "*** Error reading sites in box ",ii,", ",jj,", ",kk," in read_system()"
              if (present(error)) then
                error=internal_error
              endif
              close(fileunit)
              return
            endif
          endif
        end do
      end do
    end do
    ! Backwards-compatibility, since the earlier version wrote a lot of zeros here...
    do
      read(fileunit,"(a)") line
      if (len_trim(adjustl(line))==1) then
        read(line,*) system%allow_exchange
        exit
      endif
    end do
    read(fileunit,*) system%periodic
    write(6,"(a)") "# Reading sites..."
    call read_sites(system,fileunit,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error reading sites from file "//trim(filename)// " in read_system()"
      if (present(error)) then
        error=internal_error
      endif
      close(fileunit)
      return
    endif
    write(6,"(a)") "# Reading frozen substrate atomic system..."
    call read_my_atoms(system%frozen_substrate,fileunit,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error reading frozen substrate from file "//trim(filename)// " in read_system()"
      if (present(error)) then
        error=internal_error
      endif
      close(fileunit)
      return
    endif
    write(6,"(a,i0,a)") "# Read ",system%frozen_substrate%last_id-system%frozen_substrate%first_id+1,&
                        " atoms in the frozen substrate."
    write(6,"(a)") "# Reading mobile substrate atomic system..."
    call read_my_atoms(system%substrate,fileunit,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error reading substrate from file "//trim(filename)// " in read_system()"
      if (present(error)) then
        error=internal_error
      endif
      close(fileunit)
      return
    endif
    write(6,"(a,i0,a)") "# Read ",system%substrate%last_id-system%substrate%first_id+1," atoms in the mobile substrate."
    close(fileunit)
  end subroutine read_system

  subroutine read_sites(system,fileunit,error)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: fileunit
    integer,optional,intent(out) :: error

    integer :: i,j,internal_error,print_interval,print_remainder

    if (present(error)) then
      error=0
    endif

    print_interval=system%n_sites/100
    print_remainder=system%n_sites-100*print_interval
    call print_progress(0)
    do i=1,system%n_sites
      if (i>print_remainder) then
        if (mod(i-print_remainder,print_interval)==0) then
          call print_progress((i-print_remainder)/print_interval)
        endif
      endif
      read(fileunit,*,iostat=internal_error) system%sites(i)%id
      if (internal_error/=0) then
        write(0,"(a,i0,a)") "*** Error: couldn't read id of site ",i," in read_sites()"
        if (present(error)) then
          error=internal_error
        endif
        return
      endif
      read(fileunit,*) system%sites(i)%Z
      read(fileunit,*) system%sites(i)%atom_type
      read(fileunit,*) system%sites(i)%nn
      read(fileunit,*) system%sites(i)%nn_substrate
      read(fileunit,*) system%sites(i)%n_near
      read(fileunit,*) system%sites(i)%max_nn
      read(fileunit,*) system%sites(i)%inactive_counter
      read(fileunit,*) system%sites(i)%locked
      read(fileunit,*) system%sites(i)%lock_checked
      read(fileunit,*) system%sites(i)%jumps_checked
      read(fileunit,*) system%sites(i)%substrate
      read(fileunit,*) system%sites(i)%initial_site
      read(fileunit,*) system%sites(i)%cartesian_coords
      read(fileunit,*) system%sites(i)%my_box
      read(fileunit,*) system%sites(i)%my_index_in_box
      read(fileunit,*) system%sites(i)%neighbour_shells(:)%initialised
      do j=1,NN_SHELLS
        if (.not. system%sites(i)%neighbour_shells(j)%initialised) then
          cycle
        endif
        read(fileunit,*) system%sites(i)%neighbour_shells(j)%n_neighbour_sites
        allocate(system%sites(i)%neighbour_shells(j)%neighbour_sites(system%sites(i)%neighbour_shells(j)%n_neighbour_sites))
        allocate(system%sites(i)%neighbour_shells(j)%reverse_indices(system%sites(i)%neighbour_shells(j)%n_neighbour_sites))
        read(fileunit,*) system%sites(i)%neighbour_shells(j)%neighbour_sites
        read(fileunit,*) system%sites(i)%neighbour_shells(j)%reverse_indices
      end do
      read(fileunit,*) system%sites(i)%hop_rate
      read(fileunit,*) system%sites(i)%exchange_rate
      if (.not. system%sites(i)%substrate .and. system%sites(i)%Z>0) then
        allocate(system%sites(i)%hop_states(system%sites(i)%max_nn,system%max_jump_shell))
        allocate(system%sites(i)%hop_barriers(2,system%sites(i)%max_nn,system%max_jump_shell))
        allocate(system%sites(i)%hop_rates(system%sites(i)%max_nn,system%max_jump_shell))
        do j=1,N_JUMP_SHELLS
          read(fileunit,*) system%sites(i)%hop_states(:,j)
          read(fileunit,*) system%sites(i)%hop_barriers(:,:,j)
          read(fileunit,*) system%sites(i)%hop_rates(:,j)
        end do
        if (system%allow_exchange) then
          allocate(system%sites(i)%exchange_states(system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites,&
                                                   system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites))
          allocate(system%sites(i)%exchange_barriers(2,system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites,&
                                                       system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites))
          allocate(system%sites(i)%exchange_rates(system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites,&
                                                system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites))
          do j=1,system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
            read(fileunit,*) system%sites(i)%exchange_states(:,j)
            read(fileunit,*) system%sites(i)%exchange_barriers(:,:,j)
            read(fileunit,*) system%sites(i)%exchange_rates(:,j)
          end do
        endif
      endif
    end do
    write(6,*)
  end subroutine read_sites

  subroutine read_my_atoms(atoms,fileunit,error)
    type(t_atoms),intent(out) :: atoms
    integer,intent(in) :: fileunit
    integer,optional,intent(out) :: error

    integer :: i,internal_error,print_interval,print_remainder

    if (present(error)) then
      error=0
    endif

    read(fileunit,*,iostat=internal_error) atoms%initialised
    if (internal_error/=0) then
      write(0,"(a)") "*** Error: couldn't read the initialised-status of atoms in read_my_atoms()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    read(fileunit,*) atoms%lattice
    read(fileunit,*) atoms%periodic
    read(fileunit,*) atoms%first_id,atoms%last_id
    allocate(atoms%Z(atoms%first_id:atoms%last_id))
    allocate(atoms%atom_type(atoms%first_id:atoms%last_id))
    allocate(atoms%pos(3,atoms%first_id:atoms%last_id))
    print_interval=(atoms%last_id-atoms%first_id+1)/100
    print_remainder=atoms%last_id-atoms%first_id+1-100*print_interval
    call print_progress(0)
    do i=atoms%first_id,atoms%last_id
      if (i-atoms%first_id+1>print_remainder) then
        if (mod(i-atoms%first_id+1-print_remainder,atoms%last_id-atoms%first_id+1)==0) then
          call print_progress((i-atoms%first_id+1-print_remainder)/(atoms%last_id-atoms%first_id+1))
        endif
      endif
      read(fileunit,*) atoms%Z(i)
      read(fileunit,*) atoms%atom_type(i)
      read(fileunit,*) atoms%pos(:,i)
    end do
    write(6,*)
  end subroutine read_my_atoms

  logical function has_support(system,i)
    type(t_system),intent(in) :: system
    integer,intent(in) :: i

    integer :: nn_atom,nn_id,nnn,nns(3)

    has_support=.false.
    if (system%sites(i)%nn>3) then
      has_support=.true.
      return
    endif

    nnn=0
    do nn_atom=1,system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
      nn_id=system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(nn_atom)
      if (system%sites(nn_id)%Z>0) then
        nnn=nnn+1
        nns(nnn)=nn_id
        if (nnn>=3) then
          exit
        endif
      endif
    end do

    if (are_1nn_neighbours(system,nns(1),nns(2),nns(3))) then
      has_support=.true.
    endif
  end function has_support

  subroutine alloc_jumps(system,id,error)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: id
    integer,optional,intent(out) :: error

    integer :: max_nn,nn

    if (present(error)) then
      error=0
    endif

    if (allocated(system%sites(id)%hop_states)) then
      return
    endif

    max_nn=system%sites(id)%max_nn
    nn=system%sites(id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites

    allocate(system%sites(id)%hop_states(max_nn,system%max_jump_shell))
    system%sites(id)%hop_states=EVENT_EMPTY
    allocate(system%sites(id)%hop_barriers(2,max_nn,system%max_jump_shell))
    system%sites(id)%hop_barriers=100.0_rk
    allocate(system%sites(id)%hop_rates(max_nn,system%max_jump_shell))
    system%sites(id)%hop_rates=0.0_rk

    if (system%allow_exchange) then
      allocate(system%sites(id)%exchange_states(nn,&
                                                nn))
      system%sites(id)%exchange_states=EVENT_EMPTY
      allocate(system%sites(id)%exchange_barriers(2,nn,&
                                                    nn))
      system%sites(id)%exchange_barriers=100.0_rk
      allocate(system%sites(id)%exchange_rates(nn,&
                                               nn))
      system%sites(id)%exchange_rates=0.0_rk
    endif

  end subroutine alloc_jumps

  real(kind=rk) function get_Xnn_dist(lattice_constant,x,error)
    real(kind=rk),intent(in) :: lattice_constant
    integer,intent(in) :: x
    integer,optional,intent(out) :: error

    if (present(error)) then
      error=0
    endif
    get_Xnn_dist=0

    if (x<1 .or. x>6) then
      write(0,"(a)") "*** Error: requesting NN distance for X<1 or X>6 nearest neighbour in get_Xnn_dist()"
      if (present(error)) then
        error=1
      endif
      return
    endif

    if (x<6) then
      get_Xnn_dist=0.5_rk*lattice_constant*(NEIGH_DIST_FACT(x)+NEIGH_DIST_FACT(x+1))
    else
      get_Xnn_dist=1.1_rk*lattice_constant*(NEIGH_DIST_FACT(x))
    endif
  end function get_Xnn_dist

  subroutine check_validity(system,verbosity,error)
    type(t_system),intent(in) :: system
    integer,optional,intent(in) :: verbosity
    integer,optional,intent(out) :: error

    integer :: error_count,old_error_count,internal_error,my_verbosity
    integer :: id,n_sites_in_boxes,dist_index,shell,nn_id,nn_atom,c,reverse_index
    integer :: i,j,ii,jj,kk,wrapped_ii,wrapped_jj,wrapped_kk,my_box(3),nn,n_near
    integer :: jump_shell,site_id,initial_site,final_site,initial_id,final_id
    double precision :: hop_rate,exchange_rate
    real(kind=rk) :: coords(3)
    logical :: found
    logical,allocatable :: checked(:)

    if (present(error)) then
      error=0
    endif
    internal_error=0
    error_count=0

    if (present(verbosity)) then
      my_verbosity=verbosity
    else
      my_verbosity=VERBOSITY_LOUD
    endif

    allocate(checked(system%n_sites))
    checked=.false.

    if (my_verbosity>=VERBOSITY_LOUD) then
      write(6,"(a)") "# Checking the number of sites in boxes..."
    endif
    old_error_count=error_count
    n_sites_in_boxes=sum(system%n_sites_in_boxes)
    if (n_sites_in_boxes/=system%n_sites) then
      error_count=error_count+1
      write(0,"(a,i0,a,i0,a)") "*** Error: total count of sites in boxes ",n_sites_in_boxes,&
                               " not equal to system%n_sites ",system%n_sites
    endif
    if (my_verbosity>=VERBOSITY_LOUD) then
      write(6,"(a,i0,a)") "# ",error_count-old_error_count," errors found"
    endif

    if (my_verbosity>=VERBOSITY_LOUD) then
      write(6,"(a)") "# Checking the cartesian coordinates correspond to boxes..."
    endif
    old_error_count=error_count
    do i=1,system%n_sites
      coords(:)=system%sites(i)%cartesian_coords
      my_box(:)=int(coords(:)/system%one_box_dim(:))+1
      do c=1,3
        if (my_box(c)>system%n_boxes(c)) then
          my_box(c)=system%n_boxes(c)
        endif
      end do
      if (any(my_box/=system%sites(i)%my_box)) then
        do c=1,3
          if (coords(c)<(system%sites(i)%my_box(c)-1)*system%one_box_dim(c)-1e-6 .or.&
              coords(c)>system%sites(i)%my_box(c)*system%one_box_dim(c)+1e-6) then
            error_count=error_count+1
            write(0,"(a,i0,a,i0,a,i0,a,i0)") "*** Error: site ",i," my_box(",c,") likely incorrect: should be ",&
                                             my_box(c),"; is saved as ",&
                                             system%sites(i)%my_box(c)
            write(0,"(a,f11.6)") "*** Coords(",c,"): ",coords(c)
            write(0,"(a,f11.6,a,f11.6)") "*** Saved box bounds: ",&
                                        (system%sites(i)%my_box(c)-1)*system%one_box_dim(c)," to ",&
                                        system%sites(i)%my_box(c)*system%one_box_dim(c)
          endif
        end do
      endif
    end do
    if (my_verbosity>=VERBOSITY_LOUD) then
      write(6,"(a,i0,a)") "# ",error_count-old_error_count," errors found"
    endif

    if (my_verbosity>=VERBOSITY_LOUD) then
      write(6,"(a)") "# Checking box-to-site mapping..."
    endif
    old_error_count=error_count
    do kk=1,system%n_boxes(3)
      do jj=1,system%n_boxes(2)
        do ii=1,system%n_boxes(1)
          do i=1,system%n_sites_in_boxes(ii,jj,kk)
            id=system%sites_in_boxes(i,ii,jj,kk)
            if (id==0) then
              error_count=error_count+1
              write(0,"(a,i0,a,i0,a,i0,a,i0,a)") "*** Error: site number ",i," in box ",ii,", ",jj,", ",kk," has index 0"
            elseif (id>system%n_sites) then
              error_count=error_count+1
              write(0,"(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)") "*** Error: site number ",i," in box ",ii,", ",jj,&
                ", ",kk," has index ",id,">n_sites=",system%n_sites
            else
              if (system%sites(id)%my_index_in_box/=i) then
                error_count=error_count+1
                write(0,"(a,i0,a,i0,a,i0)") "*** Error: site ",id," has wrong information on own index in box: my_index=",&
                  system%sites(id)%my_index_in_box,"!=",i
              endif
              if (any(system%sites(id)%my_box/=[ii,jj,kk])) then
                error_count=error_count+1
                write(0,"(a,i0,a)") "*** Error: site ",id," not in correct box"
                write(0,"(a,i0,a,i0,a,i0)") "*** Saved in system%sites_in_boxes: ",ii,", ",jj,", ",kk
                write(0,"(a,i0,a,i0,a,i0)") "*** Saved in system%sites(:)%my_box: ",&
                                            system%sites(id)%my_box(1),", ",system%sites(id)%my_box(2),&
                                            ", ",system%sites(id)%my_box(3)
              endif
              if (checked(id)) then
                error_count=error_count+1
                write(0,"(a,i0,a)") "*** Error: site ",id," in multiple boxes"
              else
                checked(id)=.true.
              endif
            endif
          end do
        end do
      end do
    end do
    if (.not. all(checked)) then
      error_count=error_count+1
      write(0,"(a)") "*** Error: not all sites properly in boxes. Missing sites:"
      do i=1,system%n_sites
        if (.not. checked(i)) then
          write(0,"(a,i0)") "*** ",i
        endif
      end do
    endif
    if (my_verbosity>=VERBOSITY_LOUD) then
      write(6,"(a,i0,a)") "# ",error_count-old_error_count," errors found"
    endif

    if (my_verbosity>=VERBOSITY_LOUD) then
      write(6,"(a)") "# Checking non-substrate site neighbour lists and collisions..."
    endif
    old_error_count=error_count
    do i=system%n_substrate+1,system%n_sites
      my_box=system%sites(i)%my_box
      do kk=my_box(3)-1,my_box(3)+1
        wrapped_kk=kk
        if (wrapped_kk<1) then
          if (system%periodic(3)) then
            wrapped_kk=wrapped_kk+system%n_boxes(3)
          else
            cycle
          endif
        endif
        if (wrapped_kk>system%n_boxes(3)) then
          if (system%periodic(3)) then
            wrapped_kk=wrapped_kk-system%n_boxes(3)
          else
            cycle
          endif
        endif
        do jj=my_box(2)-1,my_box(2)+1
          wrapped_jj=jj
          if (wrapped_jj<1) then
            if (system%periodic(2)) then
              wrapped_jj=wrapped_jj+system%n_boxes(2)
            else
              cycle
            endif
          endif
          if (wrapped_jj>system%n_boxes(2)) then
            if (system%periodic(2)) then
              wrapped_jj=wrapped_jj-system%n_boxes(2)
            else
              cycle
            endif
          endif
          do ii=my_box(1)-1,my_box(1)+1
            wrapped_ii=ii
            if (wrapped_ii<1) then
              if (system%periodic(1)) then
                wrapped_ii=wrapped_ii+system%n_boxes(1)
              else
                cycle
              endif
            endif
            if (wrapped_ii>system%n_boxes(1)) then
              if (system%periodic(1)) then
                wrapped_ii=wrapped_ii-system%n_boxes(1)
              else
                cycle
              endif
            endif
            do j=1,system%n_sites_in_boxes(wrapped_ii,wrapped_jj,wrapped_kk)
              id=system%sites_in_boxes(j,wrapped_ii,wrapped_jj,wrapped_kk)
              if (id>system%n_sites) then
                error_count=error_count+1
                write(0,"(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)") "*** Error: site ",j," in box ",wrapped_ii,", ",wrapped_jj,&
                                                             ", ",wrapped_kk,&
                                                             " has id ",id,">n_sites=",system%n_sites," in check_validity()"
                cycle
              endif
              if (id==i) then
                cycle
              endif
              dist_index=get_dist_index(system,id,i,max_shell=system%lae_cutoff_shell,error=internal_error)
              if (internal_error/=0) then
                error_count=error_count+1
                write(0,"(a,i0,a,i0,a)") "*** Error: couldn't get distance index for pair ",i,", ",id," in check_validity()"
                if (dist_index>0) then
                  write(0,"(a,f11.6,a,f11.6,a)") "*** Distance: ",sqrt(system_dist2_min_image(system,i,id))," Å vs expected ",&
                                                   sqrt(JUMP_DIST2_INDEX(dist_index))*system%lattice_constant," Å"
                else
                  write(0,"(a,f11.6,a)") "*** Distance: ",sqrt(system_dist2_min_image(system,i,id))," Å; dist_index=0"
                endif
                cycle
              endif
              if (dist_index>0) then
                found=.false.
                do nn_atom=1,system%sites(i)%neighbour_shells(dist_index)%n_neighbour_sites
                  nn_id=system%sites(i)%neighbour_shells(dist_index)%neighbour_sites(nn_atom)
                  if (nn_id<=0) then
                    error_count=error_count+1
                    write(0,"(a,i0,a,i0,a,i0,a,i0)") "*** Error: site ",i," neighbour number ",nn_atom,&
                                                     " has nn_id ",nn_id,"<=0 on shell ",dist_index
                  elseif (nn_id>system%n_sites) then
                    error_count=error_count+1
                    write(0,"(a,i0,a,i0,a,i0,a,i0,a,i0)") "*** Error: site ",i," neighbour number ",nn_atom,&
                                                          " has nn_id ",nn_id,">n_sites=",system%n_sites," on shell ",dist_index
                  else
                    if (nn_id==id) then
                      if (found) then
                        error_count=error_count+1
                        write(0,"(a,i0,a,i0,a,i0)") "*** Error: site ",i," neighbour site ",id,&
                                                    " duplicated on shell ",dist_index
                        write(0,"(a,f11.6,a,f11.6,a)") "*** Distance: ",sqrt(system_dist2_min_image(system,i,id)),&
                                                       " Å vs expected ",&
                                                       sqrt(JUMP_DIST2_INDEX(dist_index))*system%lattice_constant," Å"
                      else
                        found=.true.
                      endif
                    endif
                  endif
                end do
                if (.not. found) then
                  error_count=error_count+1
                  write(0,"(a,i0,a,i0,a,i0)") "*** Error: site ",i," neighbour site ",id," not found on shell ",dist_index
                  write(0,"(a,f11.6,a,f11.6,a)") "*** Distance: ",sqrt(system_dist2_min_image(system,i,id))," Å vs expected ",&
                                                 sqrt(JUMP_DIST2_INDEX(dist_index))*system%lattice_constant," Å"
                endif
              endif
              do shell=1,system%lae_cutoff_shell
                if (shell==dist_index) then
                  cycle
                endif
                do nn_atom=1,system%sites(i)%neighbour_shells(shell)%n_neighbour_sites
                  nn_id=system%sites(i)%neighbour_shells(shell)%neighbour_sites(nn_atom)
                  if (nn_id==id) then
                    error_count=error_count+1
                    if (dist_index>0) then
                      write(0,"(a,i0,a,i0,a,i0,a,i0)") "*** Error: site ",i," neighbour site ",&
                                                       id," on wrong shell ",shell,"; expected only on ",dist_index
                      write(0,"(a,f11.6,a,f11.6,a)") "*** Distance: ",sqrt(system_dist2_min_image(system,i,id))," Å vs expected ",&
                                                       sqrt(JUMP_DIST2_INDEX(dist_index))*system%lattice_constant," Å"
                    else
                      write(0,"(a,i0,a,i0,a,i0,a)") "*** Error: site ",i," neighbour site ",&
                                                       id," on shell ",shell,"; not expected"
                      write(0,"(a,f11.6,a)") "*** Distance: ",sqrt(system_dist2_min_image(system,i,id))," Å; dist_index=0"
                    endif
                  endif
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    if (my_verbosity>=VERBOSITY_LOUD) then
      write(6,"(a,i0,a)") "# ",error_count-old_error_count," errors found"
    endif

    if (my_verbosity>=VERBOSITY_LOUD) then
      write(6,"(a)") "# Checking coordination numbers..."
    endif
    old_error_count=error_count
    do i=system%n_substrate+1,system%n_sites
      nn=0
      n_near=0
      do shell=1,TRUE_1NN_SHELL-1
        do nn_atom=1,system%sites(i)%neighbour_shells(shell)%n_neighbour_sites
          nn_id=system%sites(i)%neighbour_shells(shell)%neighbour_sites(nn_atom)
          if (nn_id<=0) then
            error_count=error_count+1
            write(0,"(a,i0,a,i0,a,i0,a,i0)") "*** Error: site ",i," neighbour number ",nn_atom,&
                                             " has nn_id ",nn_id,"<=0 on shell ",shell
          elseif (nn_id>system%n_sites) then
            error_count=error_count+1
            write(0,"(a,i0,a,i0,a,i0,a,i0,a,i0)") "*** Error: site ",i," neighbour number ",nn_atom,&
                                             " has nn_id ",nn_id,">n_sites=",system%n_sites," on shell ",shell
          else
            if (system%sites(nn_id)%Z>0) then
              n_near=n_near+1
            endif
          endif
        end do
      end do
      if (n_near/=system%sites(i)%n_near) then
        error_count=error_count+1
        write(0,"(a,i0,a,i0,a,i0)") "*** Error: site ",i," n_near count incorrect: found ",&
                                    n_near,", expected ",system%sites(i)%n_near
      endif
      do nn_atom=1,system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
        nn_id=system%sites(i)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(nn_atom)
        if (nn_id<=0) then
          error_count=error_count+1
          write(0,"(a,i0,a,i0,a,i0,a,i0)") "*** Error: site ",i," neighbour number ",nn_atom,&
                                           " has nn_id ",nn_id,"<=0 on shell ",TRUE_1NN_SHELL
        elseif (nn_id>system%n_sites) then
          error_count=error_count+1
          write(0,"(a,i0,a,i0,a,i0,a,i0,a,i0)") "*** Error: site ",i," neighbour number ",nn_atom,&
                                           " has nn_id ",nn_id,">n_sites=",system%n_sites," on shell ",TRUE_1NN_SHELL
        else
          if (system%sites(nn_id)%Z>0) then
            nn=nn+1
          endif
        endif
      end do
      if (nn/=system%sites(i)%nn) then
        error_count=error_count+1
        write(0,"(a,i0,a,i0,a,i0)") "*** Error: site ",i," nn count incorrect: found ",nn,", expected ",system%sites(i)%nn
      endif
    end do
    if (my_verbosity>=VERBOSITY_LOUD) then
      write(6,"(a,i0,a)") "# ",error_count-old_error_count," errors found"
    endif

    if (my_verbosity>=VERBOSITY_LOUD) then
      write(6,"(a)") "# Checking neighbour lists and reverse index mapping..."
    endif
    old_error_count=error_count
    do i=system%n_substrate+1,system%n_sites
      do shell=1,NN_SHELLS
        if (any(system%sites(i)%neighbour_shells(shell)%neighbour_sites(&
          1:system%sites(i)%neighbour_shells(shell)%n_neighbour_sites)==0)) then
          error_count=error_count+1
          write(0,"(a,i0,a,i0,a)") "*** Error: site ",i," shell ",shell," has zero indexed sites in the neighbour list:"
          write(0,"(a,i0)") "*** n_neighbour_sites=",system%sites(i)%neighbour_shells(shell)%n_neighbour_sites
          write(0,*) "*** neighbour_sites:",system%sites(i)%neighbour_shells(shell)%neighbour_sites(&
            1:system%sites(i)%neighbour_shells(shell)%n_neighbour_sites)
        endif
        do nn_atom=1,system%sites(i)%neighbour_shells(shell)%n_neighbour_sites
          nn_id=system%sites(i)%neighbour_shells(shell)%neighbour_sites(nn_atom)
          if (nn_id<=0) then
            error_count=error_count+1
            write(0,"(a,i0,a,i0,a,i0,a,i0)") "*** Error: site ",i," neighbour number ",nn_atom,&
                                             " has nn_id ",nn_id,"<=0 on shell ",shell
          elseif (nn_id>system%n_sites) then
            error_count=error_count+1
            write(0,"(a,i0,a,i0,a,i0,a,i0,a,i0)") "*** Error: site ",i," neighbour number ",nn_atom,&
                                                  " has nn_id ",nn_id,">n_sites=",system%n_sites," on shell ",shell
          else
            reverse_index=system%sites(i)%neighbour_shells(shell)%reverse_indices(nn_atom)
            if (reverse_index>system%sites(nn_id)%neighbour_shells(shell)%n_neighbour_sites) then
              error_count=error_count+1
              write(0,"(a,i0,a,i0,a,i0,a,i0,a,i0)") "*** Error: site ",i," nn_atom ",nn_atom," on shell ",shell,&
                                         " reverse index ",reverse_index," exceeds n_neighbour_sites of site ",nn_id
            else
              if (system%sites(nn_id)%neighbour_shells(shell)%neighbour_sites(reverse_index)/=i) then
                error_count=error_count+1
                write(0,"(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)") "*** Error: site ",i," nn atom ",nn_atom," on shell ",shell,&
                                         " reverse index ",reverse_index," points to site ",&
                                         system%sites(nn_id)%neighbour_shells(shell)%neighbour_sites(reverse_index),&
                                         " on site ",nn_id,"'s neighbour list"
              endif
            endif
          endif
        end do
      end do
    end do
    if (my_verbosity>=VERBOSITY_LOUD) then
      write(6,"(a,i0,a)") "# ",error_count-old_error_count," errors found"
    endif

    if (my_verbosity>=VERBOSITY_LOUD) then
      write(6,"(a)") "# Checking deposited atoms indexing..."
    endif
    checked=.false.
    old_error_count=error_count
    do i=1,system%n_deposited
      id=system%deposited_atoms(i)
      if (id<=0) then
        error_count=error_count+1
        write(0,"(a,i0,a,i0,a)") "*** Error: site number ",i," in deposited_atoms array has index ",id,"<=0"
      elseif (id>system%n_sites) then
        error_count=error_count+1
        write(0,"(a,i0,a,i0,a,i0)") "*** Error: site number ",i," in deposited_atoms array index ",&
                                     id,">n_sites=",system%n_sites
      else
        if (checked(id)) then
          error_count=error_count+1
          write(0,"(a,i0,a,i0,a)") "*** Error: site number ",i," with id ",id," is duplicated in deposited_atoms array"
        elseif (system%sites(id)%Z<=0) then
          error_count=error_count+1
          write(0,"(a,i0,a,i0,a)") "*** Error: site number ",i," in deposited_atoms array with id ",id," points to an empty site"
        elseif (system%sites(id)%substrate) then
          error_count=error_count+1
          write(0,"(a,i0,a,i0,a)") "*** Error: site number ",i," in deposited_atoms array with id ",id," points to substrate"
        endif
        checked(id)=.true.
      endif
    end do
    do i=system%n_substrate+1,system%n_sites
      if (system%sites(i)%Z<=0) then
        cycle
      endif
      if (.not. checked(i)) then
        error_count=error_count+1
        write(0,"(a,i0,i0,a)") "*** site ",i,", with deposited atom id ",system%sites(i)%id," missing from deposited_atoms"
      endif
    end do
    if (my_verbosity>=VERBOSITY_LOUD) then
      write(6,"(a,i0,a)") "# ",error_count-old_error_count," errors found"
    endif
    deallocate(checked)

    if (my_verbosity>=VERBOSITY_LOUD) then
      write(6,"(a)") "# Checking jump events..."
    endif
    old_error_count=error_count
    do i=1,system%n_deposited
      site_id=system%deposited_atoms(i)
      if (.not. system%sites(site_id)%jumps_checked) then
        cycle
      endif
      hop_rate=0.0d0
      if (system%sites(site_id)%locked .and. system%sites(site_id)%hop_rate>0.0d0) then
        write(0,"(a,i0,a,g13.6,a)") "*** Error: site ",site_id,&
                                    " is marked locked but it has total hop rate ",system%sites(site_id)%hop_rate," Hz"
        error_count=error_count+1
      endif
      do jump_shell=1,system%max_jump_shell
        shell=JUMP_SHELLS(jump_shell)
        do nn_atom=1,system%sites(site_id)%neighbour_shells(shell)%n_neighbour_sites
          if (system%sites(site_id)%hop_states(nn_atom,jump_shell)==EVENT_VALID) then
            hop_rate=hop_rate+system%sites(site_id)%hop_rates(nn_atom,jump_shell)
            if (system%sites(site_id)%hop_rates(nn_atom,jump_shell)<1.0d-100) then
              nn_id=system%sites(site_id)%neighbour_shells(shell)%neighbour_sites(nn_atom)
              write(0,"(a,i0,a,i0,a,g13.6,a)") "*** Warning: Hop from site ",site_id," to site ",nn_id,&
                                               " marked valid but has rate ",&
                                               system%sites(site_id)%hop_rates(nn_atom,jump_shell)," Hz"
            endif
          endif
        end do
      end do
      if (abs(1.0d0-hop_rate/system%sites(site_id)%hop_rate)>1d-6 .and. system%sites(site_id)%hop_rate>1d-16) then
        write(0,"(a,i0,a,g13.6,a,g13.6,a)") "*** Calculated hop rate for site ",site_id," ",hop_rate,&
                                            " Hz doesn't match the saved total ",system%sites(site_id)%hop_rate," Hz"
        error_count=error_count+1
      endif
      if (system%allow_exchange) then
        exchange_rate=0.0d0
        do initial_site=1,system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
          do final_site=1,system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
            if (system%sites(site_id)%exchange_states(final_site,initial_site)==EVENT_VALID) then
              initial_id=system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(initial_site)
              if (system%sites(initial_id)%locked) then
                final_id=system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(final_site)
                write(0,"(a,i0,a,g13.6,a,i0)") "*** Error: site ",initial_id,&
                                               " is marked locked but it has a valid event with exchange rate ",&
                                               system%sites(site_id)%hop_rate,&
                                               " Hz via site ",site_id," to ",final_id
                error_count=error_count+1
              endif
              exchange_rate=exchange_rate+system%sites(site_id)%exchange_rates(final_site,initial_site)
              if (system%sites(site_id)%exchange_rates(final_site,initial_site)<1.0d-100) then
                final_id=system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(final_site)
                write(0,"(a,i0,a,i0,a,i0,a,g13.6,a)") "*** Warning: exchange from site ",initial_id," via ",site_id,&
                                                      " to ",final_id," marked valid but has rate ",&
                                                      system%sites(site_id)%exchange_rates(final_site,initial_site)," Hz"
              endif
            endif
          end do
        end do
        if (abs(1.0d0-exchange_rate/system%sites(site_id)%exchange_rate)>1d-6 .and. system%sites(site_id)%exchange_rate>1d-16) then
          write(0,"(a,i0,a,g13.6,a,g13.6,a)") "*** Calculated exchange rate for site ",site_id," ",exchange_rate,&
                                              " Hz doesn't match the saved total ",system%sites(site_id)%exchange_rate," Hz"
          error_count=error_count+1
        endif
      endif
    end do
    if (my_verbosity>=VERBOSITY_LOUD) then
      write(6,"(a,i0,a)") "# ",error_count-old_error_count," errors found"
    endif

    if (present(error)) then
      error=error_count
    endif

  end subroutine check_validity

  subroutine add_atom_random(system,species,atom_type,verbosity,error)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: species,atom_type
    integer,optional,intent(in) :: verbosity
    integer,optional,intent(out) :: error
    
    integer :: dep_id,nn_id,internal_error
    integer :: nn_atom,shell
    real(kind=rk) :: u,summ
    logical :: found_site
    
    if (present(error)) then
      error=0
    endif
    internal_error=0

    u=genrand64_real2()*count(system%sites(1:system%n_sites)%Z==0 .and. system%sites(1:system%n_sites)%n_near==0)
    summ=0.0_rk
    found_site=.false.
    do dep_id=system%n_substrate+1,system%n_sites
      if (system%sites(dep_id)%Z==0 .and. system%sites(dep_id)%n_near==0) then
        summ=summ+1
        if (summ>=u) then
          found_site=.true.
          exit
        endif
      endif
    end do
    if (.not. found_site) then
      write(0,"(a)") "*** Error: couldn't find adsorption site to deposit atom to in add_atom_random()"
      if (present(error)) then
        error=1
      endif
      return
    endif

    ! If deposition site is 3-coordinated, check if there is a more stable site at shortest jump range
    found_site=.false.
    if (system%sites(dep_id)%nn==3) then
      do shell=1,TRUE_1NN_SHELL
        do nn_atom=1,system%sites(dep_id)%neighbour_shells(shell)%n_neighbour_sites
          nn_id=system%sites(dep_id)%neighbour_shells(shell)%neighbour_sites(nn_atom)
          if (system%sites(nn_id)%Z==0) then
            if (system%sites(nn_id)%n_near==0) then
              if (shell<TRUE_1NN_SHELL) then
                if (system%sites(nn_id)%nn>system%sites(dep_id)%nn) then
                  dep_id=nn_id
                  found_site=.true.
                  exit
                endif
              else
                if (system%sites(nn_id)%nn>system%sites(dep_id)%nn+1) then
                  dep_id=nn_id
                  found_site=.true.
                  exit
                endif
              endif
            endif
          endif
        end do
        if (found_site) then
          exit
        endif
      end do
    endif
    if (present(verbosity)) then
      if (verbosity>=VERBOSITY_LOUD) then
        write(6,"(a,i0,a,i0)") "# Depositing atom with Z=",species," at site ",dep_id
        write(6,"(a,3(f9.3),a)") "# Coordinates: ",system%sites(dep_id)%cartesian_coords," Å"
      endif
    endif

    call add_atom_site(system,dep_id,species,atom_type,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error adding atom on site ",dep_id," in add_atom_random()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

  end subroutine add_atom_random

  subroutine add_atom_site(system,dep_id,species,atom_type,error)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: dep_id,species,atom_type
    integer,optional,intent(out) :: error

    integer :: internal_error
    integer :: i,id,nn_id,nn_nn_id,reverse_nn_ind
    integer :: nn_atom,nn_nn_atom,shell,nn_shell,jump_shell

    if (present(error)) then
      error=0
    endif

    system%n_deposited=system%n_deposited+1
    if (system%n_deposited>size(system%deposited_atoms)) then
      call size_up(system%deposited_atoms,default_value=0,error=internal_error)
      if (internal_error/=0) then
        write(0,"(a)") "*** Error sizing up system%deposited_atoms in add_atom_site()"
        if (present(error)) then
          error=internal_error
        endif
        return
      endif
    endif
    system%deposited_atoms(system%n_deposited)=dep_id
    system%sites(dep_id)%id=system%n_substrate+system%n_deposited
    system%sites(dep_id)%Z=species
    system%sites(dep_id)%atom_type=atom_type
    system%sites(dep_id)%inactive_counter=0
    call alloc_jumps(system,dep_id,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error allocating jumps for site ",dep_id," in add_atom_site()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

    if (system%sites(dep_id)%cartesian_coords(3)>system%height) then
      system%height=system%sites(dep_id)%cartesian_coords(3)
    endif

    do jump_shell=1,system%max_jump_shell
      shell=JUMP_SHELLS(jump_shell)
      system%sites(dep_id)%hop_states(1:system%sites(dep_id)%neighbour_shells(shell)%n_neighbour_sites,&
        jump_shell)=EVENT_LAE_CHANGED
    end do
    if (system%allow_exchange) then
      system%sites(dep_id)%exchange_states(1:system%sites(dep_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites,&
                                             1:system%sites(dep_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites)&
        =EVENT_LAE_CHANGED
    endif
    system%sites(dep_id)%jumps_checked=.false.

    ! If this is the first atom, record initial position for mean square displacement tracking
    if (system%n_deposited==1) then
      system%first_atom_initial_coords(:)=system%sites(dep_id)%cartesian_coords(:)
      system%first_atom_unwrapped_coords(:)=system%sites(dep_id)%cartesian_coords(:)
    endif

    ! If this is the second, the third or the fourth atom, release the first deposited atoms to move again
    if (system%n_deposited<=4) then
      do i=1,system%n_deposited-1
        id=system%deposited_atoms(i)
        call release_atom(system,id)
      end do
    endif

    do shell=1,system%lae_cutoff_shell
      do nn_atom=1,system%sites(dep_id)%neighbour_shells(shell)%n_neighbour_sites
        nn_id=system%sites(dep_id)%neighbour_shells(shell)%neighbour_sites(nn_atom)
        if (shell<TRUE_1NN_SHELL) then
          system%sites(nn_id)%n_near=system%sites(nn_id)%n_near+1
        elseif (shell==TRUE_1NN_SHELL) then
          system%sites(nn_id)%nn=system%sites(nn_id)%nn+1
          if (.not. system%sites(nn_id)%substrate) then
            if (system%sites(nn_id)%Z<0) then
              if (system%sites(nn_id)%nn>3) then
                system%sites(nn_id)%Z=0
                system%sites(nn_id)%inactive_counter=0
              elseif (system%sites(nn_id)%nn==3) then
                if (has_support(system,nn_id)) then
                  system%sites(nn_id)%Z=0
                  system%sites(nn_id)%inactive_counter=0
                endif
              endif
            endif
          endif
        endif
        if (system%sites(nn_id)%Z>0 .and. .not. system%sites(nn_id)%substrate) then
          do jump_shell=1,system%max_jump_shell
            nn_shell=JUMP_SHELLS(jump_shell)
            system%sites(nn_id)%hop_states(1:system%sites(nn_id)%neighbour_shells(nn_shell)%n_neighbour_sites,&
              jump_shell)=EVENT_LAE_CHANGED
          end do
          if (system%allow_exchange) then
            system%sites(nn_id)%exchange_states(&
              1:system%sites(nn_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites,&
              1:system%sites(nn_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites)=EVENT_LAE_CHANGED
          endif
          system%sites(nn_id)%jumps_checked=.false.
          jump_shell=NN_SHELL_TO_JUMP_SHELL(shell)
          if (jump_shell>0) then
            reverse_nn_ind=system%sites(dep_id)%neighbour_shells(shell)%reverse_indices(nn_atom)
            system%sites(nn_id)%hop_states(reverse_nn_ind,jump_shell)=EVENT_FORBIDDEN
            if (system%allow_exchange) then
              if (shell==TRUE_1NN_SHELL) then
                system%sites(nn_id)%exchange_states(reverse_nn_ind,:)=EVENT_FORBIDDEN
              endif
            endif
          elseif (system%sites(nn_id)%Z==0) then
            do nn_shell=1,system%max_jump_shell
              do nn_nn_atom=1,system%sites(nn_id)%neighbour_shells(nn_shell)%n_neighbour_sites
                nn_nn_id=system%sites(nn_id)%neighbour_shells(nn_shell)%neighbour_sites(nn_nn_atom)
                system%sites(nn_nn_id)%jumps_checked=.false.
                if (system%allow_exchange) then
                  if (nn_shell==TRUE_1NN_SHELL) then
                    reverse_nn_ind=system%sites(nn_id)%neighbour_shells(nn_shell)%reverse_indices(nn_nn_atom)
                    system%sites(nn_nn_id)%exchange_states(reverse_nn_ind,&
                      1:system%sites(nn_nn_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites&
                      )=EVENT_LAE_CHANGED
                  endif
                endif
              end do
            end do
          endif
        endif
      end do
    end do

    ! Find new adsorption sites around the deposited
    call find_new_sites(system,dep_id,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error finding new sites around site ",&
                          dep_id," in add_atom_site()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

    system%sites(1:system%n_sites)%lock_checked=.false.
    call check_lockedness(system,dep_id,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error checking lockedness of site ",&
                          dep_id," in add_atom_site()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    system%sites(dep_id)%lock_checked=.true.

    call do_lockedness_checks(system,dep_id,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error doing lockedness checks to site ",dep_id,",in add_atom_site()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

  end subroutine add_atom_site

  subroutine release_atom(system,id)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: id

    integer :: jump_shell,shell

    do jump_shell=1,system%max_jump_shell
      shell=JUMP_SHELLS(jump_shell)
      system%sites(id)%hop_states(1:system%sites(id)%neighbour_shells(shell)%n_neighbour_sites,&
        jump_shell)=EVENT_LAE_CHANGED
    end do
    if (system%allow_exchange) then
      system%sites(id)%exchange_states(1:system%sites(id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites,&
                                       1:system%sites(id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites)&
        =EVENT_LAE_CHANGED
    endif
    system%sites(id)%jumps_checked=.false.
  end subroutine release_atom

  subroutine add_atom_xyz(system,x,y,z,species,atom_type,zmax,force_fcc,error)
    type(t_system),intent(inout) :: system
    real(kind=rk),intent(in) :: x,y,z
    integer,intent(in) :: species,atom_type
    real(kind=rk),optional,intent(in) :: zmax
    logical,optional,intent(in) :: force_fcc
    integer,optional,intent(out) :: error

    integer :: internal_error
    integer :: my_box(3),c,ii,jj,kk,wrapped_ii,wrapped_jj,wrapped_kk
    integer :: i,site_id,dep_id,shell,nn_atom,nn_id
    real(kind=rk) :: d2,d2_min,my_zmax
    logical :: my_force_fcc,not_fcc

    if (present(error)) then
      error=0
    endif

    my_zmax=system%lattice(3,3)
    if (present(zmax)) then
      my_zmax=zmax
    endif

    my_force_fcc=.false.
    if (present(force_fcc)) then
      my_force_fcc=force_fcc
    endif

    my_box(:)=int([x,y,z]/system%one_box_dim(:))+1
    do c=1,3
      if (my_box(c)>system%n_boxes(c)) then
        my_box(c)=system%n_boxes(c)
      endif
    end do

    d2_min=huge(d2_min)
    dep_id=0
    do kk=my_box(3)-1,my_box(3)+1
      wrapped_kk=kk
      if (wrapped_kk<1) then
        if (system%periodic(3)) then
          wrapped_kk=wrapped_kk+system%n_boxes(3)
        else
          cycle
        endif
      endif
      if (wrapped_kk>system%n_boxes(3)) then
        if (system%periodic(3)) then
          wrapped_kk=wrapped_kk-system%n_boxes(3)
        else
          cycle
        endif
      endif
      do jj=my_box(2)-1,my_box(2)+1
        wrapped_jj=jj
        if (wrapped_jj<1) then
          if (system%periodic(2)) then
            wrapped_jj=wrapped_jj+system%n_boxes(2)
          else
            cycle
          endif
        endif
        if (wrapped_jj>system%n_boxes(2)) then
          if (system%periodic(2)) then
            wrapped_jj=wrapped_jj-system%n_boxes(2)
          else
            cycle
          endif
        endif
        do ii=my_box(1)-1,my_box(1)+1
          wrapped_ii=ii
          if (wrapped_ii<1) then
            if (system%periodic(1)) then
              wrapped_ii=wrapped_ii+system%n_boxes(1)
            else
              cycle
            endif
          endif
          if (wrapped_ii>system%n_boxes(1)) then
            if (system%periodic(1)) then
              wrapped_ii=wrapped_ii-system%n_boxes(1)
            else
              cycle
            endif
          endif
          do i=1,system%n_sites_in_boxes(wrapped_ii,wrapped_jj,wrapped_kk)
            site_id=system%sites_in_boxes(i,wrapped_ii,wrapped_jj,wrapped_kk)
            if (system%sites(site_id)%Z/=0) then
              cycle
            endif
            if (system%sites(site_id)%n_near>0) then
              cycle
            endif
            if (system%sites(site_id)%cartesian_coords(3)>my_zmax) then
              cycle
            endif
            if (my_force_fcc) then
              ! Look for atoms between true fcc 2nn and 3nn shells, disregard site if any are found
              not_fcc=.false.
              do shell=TRUE_2NN_SHELL+1,TRUE_3NN_SHELL-1
                do nn_atom=1,system%sites(site_id)%neighbour_shells(shell)%n_neighbour_sites
                  nn_id=system%sites(site_id)%neighbour_shells(shell)%neighbour_sites(nn_atom)
                  if (system%sites(nn_id)%Z>0) then
                    not_fcc=.true.
                    exit
                  endif
                end do
                if (not_fcc) then
                  exit
                endif
              end do
              if (not_fcc) then
                cycle
              endif
            endif
            d2=atoms_dist2_min_image(system%lattice,system%sites(site_id)%cartesian_coords,[x,y,z],system%periodic)
            if (d2<d2_min) then
              d2_min=d2
              dep_id=site_id
            endif
          end do
        end do
      end do
    end do

    if (dep_id<=0) then
      write(0,"(a)") "*** Error finding the nearest site in add_atom_xyz()"
      if (present(error)) then
        error=1
      endif
      return
    endif

    call add_atom_site(system,dep_id,species,atom_type,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error adding atom on site ",dep_id," in add_atom_xyz()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

  end subroutine add_atom_xyz

  subroutine add_hexagon(system,max_layers,species,atom_type,error)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: max_layers,species,atom_type
    integer,optional,intent(out) :: error

    integer :: internal_error
    integer :: layers,n_deposited,i,id
    real(kind=rk) :: x,y,z

    if (present(error)) then
      error=0
    endif

    if (max_layers<1) then
      return
    endif

    x=system%lattice(1,1)/2
    y=system%lattice(2,2)/2
    z=system%height

    call add_atom_xyz(system,x,y,z,species,atom_type,force_fcc=.true.,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error adding the first atom in add_hexagon()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    layers=1

    do
      if (layers>=max_layers) then
        exit
      endif
      n_deposited=system%n_deposited
      do i=1,n_deposited
        id=system%deposited_atoms(i)
        x=system%sites(id)%cartesian_coords(1)
        y=system%sites(id)%cartesian_coords(2)
        z=system%sites(id)%cartesian_coords(3)
        do
          if (system%sites(id)%nn>=9) then
            exit
          endif
          call add_atom_xyz(system,x,y,z,species,atom_type,zmax=z+DIST2_TOL,error=internal_error)
          if (internal_error/=0) then
            write(0,"(a,i0,a)") "*** Error adding atom in layer ",layers+1," in add_hexagon()"
            if (present(error)) then
              error=internal_error
            endif
            return
          endif
        end do
      end do
      layers=layers+1
    end do

  end subroutine add_hexagon

  subroutine add_triangle(system,max_layers,species,atom_type,termination,error)
    type(t_system),intent(inout) :: system
    integer,intent(in) :: max_layers,species,atom_type,termination
    integer,optional,intent(out) :: error

    integer :: internal_error
    integer :: layers,n_deposited,i,id
    real(kind=rk) :: x,y,z
    real(kind=rk) :: xcenter,ycenter,xmin,xmax,ymin,ymax

    if (present(error)) then
      error=0
    endif

    if (max_layers<1) then
      return
    endif

    x=system%lattice(1,1)/2
    y=system%lattice(2,2)/2
    z=system%height

    call add_atom_xyz(system,x,y,z,species,atom_type,force_fcc=.true.,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error adding the first atom in add_triangle()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    xcenter=system%sites(system%deposited_atoms(1))%cartesian_coords(1)
    ycenter=system%sites(system%deposited_atoms(1))%cartesian_coords(2)
    layers=1

    do
      if (layers>=max_layers) then
        exit
      endif
      n_deposited=system%n_deposited
      do i=1,n_deposited
        id=system%deposited_atoms(i)
        x=system%sites(id)%cartesian_coords(1)
        y=system%sites(id)%cartesian_coords(2)
        z=system%sites(id)%cartesian_coords(3)
        do
          if (system%sites(id)%nn>=9) then
            exit
          endif
          call add_atom_xyz(system,x,y,z,species,atom_type,zmax=z+DIST2_TOL,error=internal_error)
          if (internal_error/=0) then
            write(0,"(a,i0,a)") "*** Error adding atom in layer ",layers+1," in add_triangle()"
            if (present(error)) then
              error=internal_error
            endif
            return
          endif
        end do
      end do
      ! Find correct corners for the vertices, based on desired termination
      ! Loop over newly deposited atoms
      xmin=xcenter
      xmax=xcenter
      ymax=ycenter
      ymin=ycenter
      z=system%height
      do i=n_deposited,system%n_deposited
        id=system%deposited_atoms(i)
        x=system%sites(id)%cartesian_coords(1)
        y=system%sites(id)%cartesian_coords(2)
        if (x>xmax) then
          xmax=x
        endif
        if (x<xmin) then
          xmin=x
        endif
        if (y>ymax) then
          ymax=y
        endif
        if (y<ymin) then
          ymin=y
        endif
      end do
      if (termination==TERMINATION_A) then
        call add_atom_xyz(system,xmax,ymax,z,species,atom_type,zmax=z+DIST2_TOL,error=internal_error)
        call add_atom_xyz(system,xmax,ymin,z,species,atom_type,zmax=z+DIST2_TOL,error=internal_error)
        call add_atom_xyz(system,xmin,ycenter,z,species,atom_type,zmax=z+DIST2_TOL,error=internal_error)
      else ! TERMINATION_B
        call add_atom_xyz(system,xmin,ymax,z,species,atom_type,zmax=z+DIST2_TOL,error=internal_error)
        call add_atom_xyz(system,xmin,ymin,z,species,atom_type,zmax=z+DIST2_TOL,error=internal_error)
        call add_atom_xyz(system,xmax,ycenter,z,species,atom_type,zmax=z+DIST2_TOL,error=internal_error)
      endif

      layers=layers+1
    end do

  end subroutine add_triangle

  subroutine delete_site(system,id,error)
    ! Remove site id by swapping it with the
    ! last site, and decrementing the site count by 1.
    ! Remaps all of the last site's neighbour lists
    ! to the site id.
    type(t_system),intent(inout) :: system
    integer,intent(in) :: id
    integer,optional,intent(out) :: error

    integer :: i,shell,jump_shell,nn_atom,nn_id,nn_nn_id,reverse_index,nn_reverse_index,n_neighbour_sites
    integer :: n_sites_in_box,my_box(3),my_index_in_box

    if (present(error)) then
      error=0
    endif

    if (id>system%n_sites) then
      write(0,"(a,i0,a,i0,a)") "*** Error deleting site id=",id,">system%n_sites=",system%n_sites," in delete_site()"
      if (present(error)) then
        error=1
      endif
      return
    endif

    if (system%sites(id)%Z>0) then
      write(0,"(a,i0,a,i0)") "*** Error: attempting to delete site ",id," occupied with atom Z=",system%sites(id)%Z
      if (present(error)) then
        error=2
      endif
      return
    endif

    if (system%sites(id)%Z==0 .and. system%sites(id)%n_near<2) then
      write(0,"(a,i0)") "*** Error: attempting to delete active site ",id
      if (present(error)) then
        error=3
      endif
      return
    endif

    do shell=1,NN_SHELLS
      jump_shell=NN_SHELL_TO_JUMP_SHELL(shell)
      do nn_atom=1,system%sites(id)%neighbour_shells(shell)%n_neighbour_sites
        nn_id=system%sites(id)%neighbour_shells(shell)%neighbour_sites(nn_atom)
        n_neighbour_sites=system%sites(nn_id)%neighbour_shells(shell)%n_neighbour_sites
        reverse_index=system%sites(id)%neighbour_shells(shell)%reverse_indices(nn_atom)

        system%sites(nn_id)%neighbour_shells(shell)%neighbour_sites(reverse_index)=&
          system%sites(nn_id)%neighbour_shells(shell)%neighbour_sites(n_neighbour_sites)
        system%sites(nn_id)%neighbour_shells(shell)%reverse_indices(reverse_index)=&
          system%sites(nn_id)%neighbour_shells(shell)%reverse_indices(n_neighbour_sites)

        if (jump_shell>0) then
          if (allocated(system%sites(nn_id)%hop_states)) then
            system%sites(nn_id)%hop_states(reverse_index,jump_shell)=&
              system%sites(nn_id)%hop_states(n_neighbour_sites,jump_shell)
            system%sites(nn_id)%hop_barriers(:,reverse_index,jump_shell)=&
              system%sites(nn_id)%hop_barriers(:,n_neighbour_sites,jump_shell)
            system%sites(nn_id)%hop_rates(reverse_index,jump_shell)=&
              system%sites(nn_id)%hop_rates(n_neighbour_sites,jump_shell)

            system%sites(nn_id)%hop_states(n_neighbour_sites,jump_shell)=EVENT_EMPTY
            if (jump_shell==TRUE_1NN_JUMP_SHELL) then
              if (system%allow_exchange) then
                system%sites(nn_id)%exchange_states(reverse_index,:)=&
                  system%sites(nn_id)%exchange_states(n_neighbour_sites,:)
                system%sites(nn_id)%exchange_states(:,reverse_index)=&
                  system%sites(nn_id)%exchange_states(:,n_neighbour_sites)
                system%sites(nn_id)%exchange_barriers(:,reverse_index,:)=&
                  system%sites(nn_id)%exchange_barriers(:,n_neighbour_sites,:)
                system%sites(nn_id)%exchange_barriers(:,:,reverse_index)=&
                  system%sites(nn_id)%exchange_barriers(:,:,n_neighbour_sites)
                system%sites(nn_id)%exchange_rates(reverse_index,:)=&
                  system%sites(nn_id)%exchange_rates(n_neighbour_sites,:)
                system%sites(nn_id)%exchange_rates(:,reverse_index)=&
                  system%sites(nn_id)%exchange_rates(:,n_neighbour_sites)

                system%sites(nn_id)%exchange_states(:,n_neighbour_sites)=EVENT_EMPTY
                system%sites(nn_id)%exchange_states(n_neighbour_sites,:)=EVENT_EMPTY
              endif
            endif
          endif
        endif

        nn_nn_id=system%sites(nn_id)%neighbour_shells(shell)%neighbour_sites(reverse_index)
        nn_reverse_index=system%sites(nn_id)%neighbour_shells(shell)%reverse_indices(reverse_index)
        system%sites(nn_nn_id)%neighbour_shells(shell)%reverse_indices(nn_reverse_index)=reverse_index

        system%sites(nn_id)%neighbour_shells(shell)%neighbour_sites(n_neighbour_sites)=0
        system%sites(nn_id)%neighbour_shells(shell)%reverse_indices(n_neighbour_sites)=0
        system%sites(nn_id)%neighbour_shells(shell)%n_neighbour_sites=n_neighbour_sites-1

      end do
    end do

    if (id/=system%n_sites) then
      do shell=1,NN_SHELLS
        do nn_atom=1,system%sites(system%n_sites)%neighbour_shells(shell)%n_neighbour_sites
          nn_id=system%sites(system%n_sites)%neighbour_shells(shell)%neighbour_sites(nn_atom)
          reverse_index=system%sites(system%n_sites)%neighbour_shells(shell)%reverse_indices(nn_atom)
          system%sites(nn_id)%neighbour_shells(shell)%neighbour_sites(reverse_index)=id
        end do
      end do
    endif

    ! Remove the site from its box
    my_box(:)=system%sites(id)%my_box(:)
    my_index_in_box=system%sites(id)%my_index_in_box
    n_sites_in_box=system%n_sites_in_boxes(my_box(1),my_box(2),my_box(3))
    system%sites_in_boxes(my_index_in_box,my_box(1),my_box(2),my_box(3))=system%sites_in_boxes(n_sites_in_box,&
      my_box(1),my_box(2),my_box(3))
    i=system%sites_in_boxes(my_index_in_box,my_box(1),my_box(2),my_box(3))
    system%sites(i)%my_index_in_box=my_index_in_box
    system%sites_in_boxes(n_sites_in_box,my_box(1),my_box(2),my_box(3))=0
    system%n_sites_in_boxes(my_box(1),my_box(2),my_box(3))=n_sites_in_box-1

    if (id/=system%n_sites) then
      ! Remap last site's box to have to correct index
      my_box(:)=system%sites(system%n_sites)%my_box(:)
      my_index_in_box=system%sites(system%n_sites)%my_index_in_box
      system%sites_in_boxes(my_index_in_box,my_box(1),my_box(2),my_box(3))=id
    endif

    if (.not. system%sites(system%n_sites)%substrate .and. system%sites(system%n_sites)%Z>0) then
      ! Update deposited_atoms index
      i=system%sites(system%n_sites)%id-system%n_substrate
      system%deposited_atoms(i)=id
    endif

    system%sites(id)=system%sites(system%n_sites)

    ! Deallocate memory and de-initialise the site
    system%sites(system%n_sites)%Z=-1
    system%sites(system%n_sites)%nn=0
    system%sites(system%n_sites)%n_near=0
    system%sites(system%n_sites)%atom_type=0
    system%sites(system%n_sites)%id=0
    system%sites(system%n_sites)%max_nn=0
    system%sites(system%n_sites)%nn_substrate=0
    system%sites(system%n_sites)%inactive_counter=0
    system%sites(system%n_sites)%locked=.false.
    system%sites(system%n_sites)%lock_checked=.false.
    system%sites(system%n_sites)%jumps_checked=.false.
    system%sites(system%n_sites)%substrate=.false.
    system%sites(system%n_sites)%initial_site=.false.
    if (allocated(system%sites(system%n_sites)%hop_states)) then
      deallocate(system%sites(system%n_sites)%hop_states)
      deallocate(system%sites(system%n_sites)%hop_barriers)
      deallocate(system%sites(system%n_sites)%hop_rates)
      deallocate(system%sites(system%n_sites)%exchange_states)
      deallocate(system%sites(system%n_sites)%exchange_barriers)
      deallocate(system%sites(system%n_sites)%exchange_rates)
    endif
    system%sites(system%n_sites)%hop_rate=0.0d0
    system%sites(system%n_sites)%exchange_rate=0.0d0
    do shell=1,NN_SHELLS
      if (system%sites(system%n_sites)%neighbour_shells(shell)%initialised) then
        system%sites(system%n_sites)%neighbour_shells(shell)%initialised=.false.
        system%sites(system%n_sites)%neighbour_shells(shell)%n_neighbour_sites=0
        deallocate(system%sites(system%n_sites)%neighbour_shells(shell)%neighbour_sites)
        deallocate(system%sites(system%n_sites)%neighbour_shells(shell)%reverse_indices)
      endif
    end do

    system%n_sites=system%n_sites-1

  end subroutine delete_site

  integer function delete_inactive_sites(system,error)
    type(t_system),intent(inout) :: system
    integer,optional,intent(out) :: error

    integer :: id,internal_error

    if (present(error)) then
      error=0
    endif
    delete_inactive_sites=0

    do id=system%n_sites,system%n_substrate+1,-1
      if (system%sites(id)%inactive_counter>INACTIVE_DELETION_THRESHOLD) then
        call delete_site(system,id,error=internal_error)
        if (internal_error/=0) then
          write(0,"(a,i0,a)") "*** Error deleting site ",id," in delete_inactive_sites()"
          if (present(error)) then
            error=internal_error
          endif
          return
        endif
        delete_inactive_sites=delete_inactive_sites+1
      endif
    end do

  end function delete_inactive_sites

end module mlkmc_atoms_module
