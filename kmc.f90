! MLKMC - module for all KMC simulation functions and subroutines
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

module kmc_module
  use constants_module
  use mlkmc_atoms_module
  use liblammps
  use libAtoms_module
  use descriptors_module
  use ml_module
  implicit none

  type t_kmc
    integer :: relax_iter,neb_images,neb_steps
    double precision :: relax_tol,ml_tol,neb_spring,neb_error
    type(t_predictor),allocatable :: hop_predictor(:,:,:),exchange_predictor(:,:,:)
    integer,allocatable :: minim_cache(:,:,:)
    integer :: minim_cache_results(MINIM_CACHE_SIZE)=MINIM_RESULT_UNINIT
    integer :: minim_cache_index=0
    integer :: neb_calculations=0,minimizations=0

    integer :: verbosity

    integer :: steps_since_deposition_or_rate_lowering=0
    integer :: rate_lowered_times=0
    integer :: rate_lowering_threshold=0

    integer :: nZ
    integer,allocatable :: Z(:)
    double precision :: temp,attempt_frequency,kT
    double precision,allocatable :: dep_rate(:)

    type(t_system) :: system

    real(kind=rk) :: substrate_fcc_hcp_barrier=-1.0_rk
    real(kind=rk) :: substrate_hcp_fcc_barrier=-1.0_rk
    real(kind=rk) :: last_barrier=0.0_rk

    logical :: debug,track_MSD=.false.

    integer :: n_data_min,n_data_max,n_sparse_max,covariance_type
    double precision :: delta,zeta,regularisation
    type(extendable_str) :: desc_str
    type(descriptor) :: desc

    integer :: step
    double precision :: time,total_rate

    double precision :: kmc_run_time
    ! Timers for largest functions
    double precision :: get_rates_time,pick_event_time,execute_jump_event_time,deposit_atom_time
    ! Timers for bigger components within functions
    double precision :: get_rates_overhead_time,event_finding_time,lae_time
    double precision :: lae_sub_time,lae_init_time,lae_add_time,lae_connect_time
    double precision :: desc_time,get_barrier_time,assigning_time
    double precision :: rate_fingerprint_time,rate_cache_time,event_fingerprint_time,event_cache_time
    double precision :: minimize_time,neb_barrier_time
    double precision :: predict_time,fit_time,write_time

    logical :: initialised=.false.
  end type t_kmc

contains

  subroutine run_kmc(kmc,lmp,temp,dep_rate,max_steps,max_time,max_wall_time,min_height,max_height,&
                     start_time,fileunit,print_interval,write_interval,myrank,error)
    type(t_kmc),intent(inout) :: kmc
    type(lammps),intent(inout) :: lmp
    integer,intent(in) :: max_steps,print_interval
    real(kind=rk),intent(in) :: max_time,write_interval,min_height,max_height
    double precision,intent(in) :: dep_rate(kmc%nZ),max_wall_time,temp
    double precision,intent(in) :: start_time
    integer,intent(in) :: fileunit
    integer,intent(in) :: myrank
    integer,optional,intent(out) :: error

    integer :: first_step,step,i,a,b,event_type,n_deleted
    double precision :: t1,t2,dump_t1
    integer :: final_site,initial_site,nn_atom,shell,jump_shell,atom_type
    double precision :: u,dt,elapsed_time,wall_time_interval
    integer :: mpi_error,internal_error,write_error
    logical :: end_loop

    if (present(error)) then
      error=0
    endif
    internal_error=0
    
    call mpi_bcast(max_steps,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)

    if (myrank==0) then
      ! Frozen substrate is one monolayer
      kmc%dep_rate=dep_rate*(kmc%system%frozen_substrate%last_id-kmc%system%frozen_substrate%first_id+1)

      ! Set diffusion rate lowering threshold to be equal to ten times the monolayer size
      kmc%rate_lowering_threshold=10*(kmc%system%frozen_substrate%last_id-kmc%system%frozen_substrate%first_id+1)
  
      kmc%temp=temp
      kmc%kT=kB*temp

      first_step=kmc%step
      if (kmc%verbosity>=VERBOSITY_QUIET) then
        write(6,"(a,i0,a)") "# Starting KMC run up to ",max_steps," steps"
      endif

    endif

    call mpi_bcast(first_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)

    if (myrank==0) then
      wall_time_interval=max_wall_time/10.0d0
      ! Always printed:
      call write_output_header(kmc)
      call write_output_line(kmc,error=write_error)
      if (write_error/=0) then
        write(0,"(a)") "*** Warning: error writing output line"
      else
        t1=mpi_wtime()
        call print_xyz_frame(kmc%system,kmc%time,kmc%step,fileunit,debug=kmc%debug,error=internal_error)
        t2=mpi_wtime()
        kmc%write_time=kmc%write_time+t2-t1
        if (internal_error/=0) then
          write(0,"(a)") "*** Error writing atoms xyz in run_kmc()"
        endif
        dump_t1=mpi_wtime()
      endif
    endif
    call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
    if (internal_error/=0) then
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

    elapsed_time=0.0d0
    end_loop=.false.
    do step=first_step+1,max_steps
      if (myrank==0) then
        kmc%step=step
        t1=mpi_wtime()
      endif
      call get_rates(kmc,lmp,myrank,error=internal_error)
      if (myrank==0) then
        t2=mpi_wtime()
        kmc%get_rates_time=kmc%get_rates_time+t2-t1
        if (internal_error/=0) then
          write(0,"(a)") "*** Error when calculating rates in run_kmc()"
        endif
      endif
      call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
      if (internal_error/=0) then
        exit
      endif
      if (myrank==0) then
        if (kmc%total_rate<=0.0_rk) then
          write(0,"(a,i0,a)") "*** Warning: no valid events found at step ",step,"!"
          internal_error=1
        endif
      endif
      call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
      if (internal_error/=0) then
        internal_error=0
        exit
      endif
      if (myrank==0) then
        t1=mpi_wtime()
        call pick_event(kmc,i,a,b,event_type,error=internal_error)
        t2=mpi_wtime()
        kmc%pick_event_time=kmc%pick_event_time+t2-t1
        if (internal_error/=0) then
          write(0,"(a,i0,a)") "*** Error: failed to pick event at step ",step," in run_kmc()"
          call write_output_line(kmc,error=write_error)
          if (write_error/=0) then
            write(0,"(a)") "*** Warning: error writing output line"
          endif
        endif
      endif
      call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
      if (internal_error/=0) then
        exit
      endif
      if (myrank==0) then
        u=genrand64_real2()
        dt=-log(u)/kmc%total_rate
        kmc%time=kmc%time+dt
        elapsed_time=elapsed_time+dt

        if (event_type==EVENT_TYPE_HOP) then
          nn_atom=a
          shell=b
          jump_shell=NN_SHELL_TO_JUMP_SHELL(shell)
          if (jump_shell<=0) then
            write(0,"(a,i0,a)") "*** Error: invalid shell ",shell," doesn't map to an allowed jump_shell in run_kmc()"
            internal_error=1
          endif
        endif
        if (internal_error==0) then
          t1=mpi_wtime()
          if (event_type==EVENT_TYPE_HOP) then
            kmc%last_barrier=kmc%system%sites(i)%hop_barriers(1,nn_atom,jump_shell)
            call execute_jump_event(kmc,i,nn_atom,shell,event_type,error=internal_error)
            t2=mpi_wtime()
            kmc%execute_jump_event_time=kmc%execute_jump_event_time+t2-t1
          elseif (event_type==EVENT_TYPE_EXCHANGE) then
            final_site=a
            initial_site=b
            kmc%last_barrier=kmc%system%sites(i)%exchange_barriers(1,final_site,initial_site)
            call execute_jump_event(kmc,i,final_site,initial_site,event_type,error=internal_error)
            t2=mpi_wtime()
            kmc%execute_jump_event_time=kmc%execute_jump_event_time+t2-t1
          elseif (event_type==EVENT_TYPE_DEPOSITION) then
            atom_type=a
            kmc%steps_since_deposition_or_rate_lowering=0
            kmc%rate_lowered_times=0
            call add_atom_random(kmc%system,kmc%Z(atom_type),atom_type,verbosity=kmc%verbosity,error=internal_error)
            t2=mpi_wtime()
            kmc%deposit_atom_time=kmc%deposit_atom_time+t2-t1
          else ! Event type none; error
            write(0,"(a)") "*** Error: no event type to execute in run_kmc()"
            internal_error=1
          endif
          if (internal_error==0) then
            where(.not. kmc%system%sites(1:kmc%system%n_sites)%initial_site .and. &
                (kmc%system%sites(1:kmc%system%n_sites)%Z<0 .or. kmc%system%sites(1:kmc%system%n_sites)%n_near>1))
              kmc%system%sites(1:kmc%system%n_sites)%inactive_counter=kmc%system%sites(1:kmc%system%n_sites)%inactive_counter+1
            end where
          endif
        else
          if (event_type==EVENT_TYPE_HOP) then
            write(0,"(a)") "*** Error: could not execute hop event from site ",i," in run_kmc()"
          elseif (event_type==EVENT_TYPE_EXCHANGE) then
            write(0,"(a)") "*** Error: could not execute exchange event about site ",i," in run_kmc()"
          elseif (event_type==EVENT_TYPE_DEPOSITION) then
            write(0,"(a)") "*** Error: could not execute deposition event in run_kmc()"
          else
            write(0,"(a)") "*** Error: no event type to execute in run_kmc()"
            internal_error=1
          endif
          call write_output_line(kmc,error=write_error)
          if (write_error/=0) then
            write(0,"(a)") "*** Warning: error writing output line"
          endif
        endif
      endif
      call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
      if (internal_error/=0) then
        exit
      endif
      if (myrank==0) then
        if (mod(step,INACTIVE_DELETION_INTERVAL)==0) then
          n_deleted=delete_inactive_sites(kmc%system,error=internal_error)
          if (internal_error==0) then
            if (kmc%verbosity>VERBOSITY_NORMAL) then
              write(6,"(a,i0,a,i0)") "# Deleted ",n_deleted," inactive sites at step ",step
            endif
            if (kmc%debug) then
              call check_validity(kmc%system,verbosity=kmc%verbosity,error=internal_error)
              if (internal_error/=0) then
                write(0,"(a,i0,a,i0,a)") "*** ",internal_error," errors detected in check_validity() at step ",step," in run_kmc()"
              endif
            endif
          else
            write(0,"(a,i0,a)") "*** Error deleting inactive sites at step ",step," in run_kmc()"
          endif
        endif
        if (internal_error==0) then
          if (print_interval>0) then
            if (mod(step,print_interval)==0) then
              call write_output_line(kmc,error=write_error)
              if (write_error/=0) then
                write(0,"(a)") "*** Warning: error writing output line"
              endif
            endif
          endif
        endif
      endif
      call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
      if (internal_error/=0) then
        exit
      endif
      if (myrank==0) then
        if (write_interval>=0) then
          if (elapsed_time>=write_interval) then
            if (kmc%verbosity>=VERBOSITY_QUIET) then
              write(6,"(a)",advance="no") "# Writing xyz frame..."
            endif
            elapsed_time=0.0_rk
            t1=mpi_wtime()
            call print_xyz_frame(kmc%system,kmc%time,kmc%step,fileunit,debug=kmc%debug,error=internal_error)
            if (internal_error/=0) then
              write(0,"(a)") "*** Error writing atoms xyz in run_kmc()"
            endif
            t2=mpi_wtime()
            kmc%write_time=kmc%write_time+t2-t1
            if (kmc%verbosity>=VERBOSITY_QUIET) then
              write(6,"(a)") " done"
            endif
          endif
        endif
      endif
      call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
      if (internal_error/=0) then
        exit
      endif

      if (myrank==0) then
        if (max_time>0) then
          if (kmc%time>=max_time) then
            if (kmc%verbosity>=VERBOSITY_QUIET) then
              write(6,"(a)") "# Reached the maximum allowed KMC time"
            endif
            call write_output_line(kmc,error=write_error)
            if (write_error/=0) then
              write(0,"(a)") "*** Warning: error writing output line"
            endif
            end_loop=.true.
          endif
        endif
        if (internal_error==0 .and. .not. end_loop) then
          if (kmc%system%height>=max_height .or. kmc%system%height<=min_height) then
            if (kmc%verbosity>=VERBOSITY_QUIET) then
              write(6,"(a)") "# Reached the minimum/maximum system height"
            endif
            call write_output_line(kmc,error=write_error)
            if (write_error/=0) then
              write(0,"(a)") "*** Warning: error writing output line"
            endif
            end_loop=.true.
          endif
        endif
        if (internal_error==0 .and. .not. end_loop .and. max_wall_time>0) then
          t2=mpi_wtime()
          if (t2-dump_t1>wall_time_interval) then
            write(6,"(a)",advance="no") "# Reached (another) tenth of maximum wall time, dumping intermediate data..."
            call dump_data_intermediate(kmc,error=internal_error)
            write(6,"(a)") " done"
            dump_t1=mpi_wtime()
          endif
          if (internal_error==0) then
            if (t2-start_time>max_wall_time) then
              if (kmc%verbosity>=VERBOSITY_QUIET) then
                write(6,"(a)") "# Reached the maximum allowed wall time"
              endif
              call write_output_line(kmc,error=write_error)
              if (write_error/=0) then
                write(0,"(a)") "*** Warning: error writing output line"
              endif
              end_loop=.true.
            endif
          endif
        endif
      endif
      call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
      if (internal_error/=0) then
        exit
      endif
      call mpi_bcast(end_loop,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
      if (end_loop) then
        exit
      endif
    end do

    if (internal_error/=0) then
      if (present(error)) then
        error=internal_error
      endif
      if (myrank==0) then
        ! Dump data
        t1=mpi_wtime()
        call dump_data_error(kmc)
        call print_xyz_frame(kmc%system,kmc%time,kmc%step,fileunit,debug=.true.)
        t2=mpi_wtime()
        kmc%write_time=kmc%write_time+t2-t1
      endif
      call mpi_barrier(mpi_comm_world,mpi_error)
      return
    endif

    if (myrank==0) then
      if (kmc%verbosity>=VERBOSITY_QUIET) then
        write(6,"(a)") "# KMC finished!"
        write(6,"(a)") "# Writing files..."
      endif
      t1=mpi_wtime()
      call print_xyz_frame(kmc%system,kmc%time,kmc%step,fileunit,debug=kmc%debug,error=internal_error)
      if (internal_error/=0) then
        write(0,"(a)") "*** Error writing atoms xyz in run_kmc()"
      endif
      t2=mpi_wtime()
      kmc%write_time=kmc%write_time+t2-t1
    endif
    call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
    if (internal_error/=0) then
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

  end subroutine run_kmc

  subroutine initialise_kmc(kmc,relax_iter,neb_images,neb_steps,neb_spring,&
                            neb_error,relax_tol,ml_tol,&
                            nZ,Z,attempt_frequency,system,verbosity,&
                            debug,track_MSD,n_data_min,n_data_max,n_sparse_max,desc_str,&
                            delta,zeta,covariance_type,regularisation,&
                            error)
    type(t_kmc),intent(inout) :: kmc
    integer,intent(in) :: relax_iter,neb_images,neb_steps
    double precision,intent(in) :: relax_tol,ml_tol,neb_spring,neb_error
    integer,intent(in) :: nZ,Z(nZ)
    double precision,intent(in) :: attempt_frequency
    type(t_system),intent(in) :: system
    integer,intent(in) :: verbosity
    logical,intent(in) :: debug,track_MSD
    integer,intent(in) :: n_data_min,n_data_max,n_sparse_max,covariance_type
    type(extendable_str),intent(in) :: desc_str
    double precision,intent(in) :: delta,zeta,regularisation
    integer,optional,intent(out) :: error
  
    integer :: internal_error
  
    if (present(error)) then
      error=0
    endif
    internal_error=0
  
    kmc%relax_iter=relax_iter
    kmc%neb_images=neb_images
    kmc%neb_steps=neb_steps
    kmc%neb_spring=neb_spring
    kmc%neb_error=neb_error
  
    kmc%relax_tol=relax_tol
    kmc%ml_tol=ml_tol
  
    kmc%nZ=nZ
    allocate(kmc%Z(nZ))
    kmc%Z=Z
  
    kmc%system=system
  
    allocate(kmc%dep_rate(nZ))
    kmc%attempt_frequency=attempt_frequency
  
    kmc%step=0
    kmc%time=0.0d0
    kmc%total_rate=0.0d0
  
    kmc%n_data_min=n_data_min
    kmc%n_data_max=n_data_max
    kmc%n_sparse_max=n_sparse_max
    kmc%desc_str=desc_str
    call initialise(kmc%desc,string(desc_str),error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error initialising descriptor in initialise_kmc()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    kmc%delta=delta
    kmc%zeta=zeta
    kmc%covariance_type=covariance_type
    kmc%regularisation=regularisation
  
    kmc%verbosity=verbosity
    kmc%debug=debug
    kmc%track_MSD=track_MSD
  
    kmc%kmc_run_time=0.0d0
    kmc%get_rates_time=0.0d0
    kmc%pick_event_time=0.0d0
    kmc%execute_jump_event_time=0.0d0
    kmc%deposit_atom_time=0.0d0

    kmc%get_rates_overhead_time=0.0d0
    kmc%get_barrier_time=0.0d0
    kmc%minimize_time=0.0d0
    kmc%neb_barrier_time=0.0d0
    kmc%lae_time=0.0d0
    kmc%lae_sub_time=0.0d0
    kmc%lae_init_time=0.0d0
    kmc%lae_add_time=0.0d0
    kmc%lae_connect_time=0.0d0
    kmc%event_finding_time=0.0d0
    kmc%rate_fingerprint_time=0.0d0
    kmc%rate_cache_time=0.0d0
    kmc%event_fingerprint_time=0.0d0
    kmc%event_cache_time=0.0d0
    kmc%assigning_time=0.0d0
    kmc%desc_time=0.0d0
    kmc%predict_time=0.0d0
    kmc%fit_time=0.0d0
    kmc%write_time=0.0d0

    kmc%initialised=.true.
  end subroutine initialise_kmc
  
  subroutine get_rates(kmc,lmp,myrank,error)
    type(t_kmc),intent(inout) :: kmc
    type(lammps),intent(inout) :: lmp
    integer,intent(in) :: myrank
    integer,optional,intent(out) :: error
  
    integer :: internal_error
    double precision :: t1,t2
    integer,save :: lae_index=1
    real(kind=rk) :: barrier,reverse_barrier
    logical :: event_ok,nan_detected,obstructed,locked
    integer :: i,site_id,n_deposited,final_id,initial_id
    integer :: shell,jump_shell,nn_atom,nn_nn_atom,initial_site,final_site
    integer :: n_near_final,nn_final,nn_initial,nn_mid
    integer :: only_substrate_jump
    integer :: max_jump_shell,max_nn,max_1nn,nn_id,nn_nn_id
    integer :: n_neighbour_sites(NN_SHELLS),fingerprint(kmc%nZ,0:NN_SHELLS),minim_result
    integer,allocatable :: nn_tmp(:)
    integer,allocatable :: n_near_tmp(:)
    logical,allocatable :: site_blocked(:)
    logical :: jumps_checked,allow_exchange
    integer,allocatable :: hop_event_states(:,:),exchange_event_states(:,:)
    type(t_atoms) :: initial_config,final_config
    logical,pointer :: jumping_ptr(:)
  
    type(Atoms) :: lae
    type(descriptor_data) :: desc_forward,desc_reverse
    integer :: mpi_error
  
    if (present(error)) then
      error=0
    endif
    internal_error=0
  
    if (myrank==0) then
      if (.not. kmc%initialised) then
        write(0,"(a)") "*** Error: kmc not initialised before get_rates()"
      endif
      t1=mpi_wtime()
    endif
    call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
    if (internal_error/=0) then
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
  
    if (myrank==0) then
      allocate(nn_tmp(kmc%system%n_sites))
      allocate(n_near_tmp(kmc%system%n_sites))
      allocate(site_blocked(kmc%system%n_sites))

      initial_config=get_atoms(kmc%system)

      max_nn=kmc%system%max_nn
      max_1nn=kmc%system%max_1nn
      max_jump_shell=kmc%system%max_jump_shell
      allow_exchange=kmc%system%allow_exchange
    endif
    call mpi_bcast(max_nn,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
    call mpi_bcast(max_1nn,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
    call mpi_bcast(max_jump_shell,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
    call mpi_bcast(allow_exchange,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
    allocate(hop_event_states(max_nn,max_jump_shell))
    if (allow_exchange) then
      allocate(exchange_event_states(max_1nn,max_1nn))
    endif
    if (myrank==0) then
      n_deposited=kmc%system%n_deposited
    endif
    call mpi_bcast(n_deposited,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
    if (myrank==0) then
      t2=mpi_wtime()
      kmc%get_rates_overhead_time=kmc%get_rates_overhead_time+t2-t1
    endif

    do i=1,n_deposited
      if (myrank==0) then
        t1=mpi_wtime()
        site_id=kmc%system%deposited_atoms(i)
        jumps_checked=kmc%system%sites(site_id)%jumps_checked
        locked=kmc%system%sites(site_id)%locked
        if (locked) then
          kmc%system%sites(site_id)%hop_states(1:kmc%system%sites(site_id)%max_nn,1:max_jump_shell)=EVENT_FORBIDDEN
          kmc%system%sites(site_id)%hop_rate=0.0_rk
        endif
      endif
      call mpi_bcast(jumps_checked,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
      if (jumps_checked) then
        if (myrank==0) then
          t2=mpi_wtime()
          kmc%get_rates_overhead_time=kmc%get_rates_overhead_time+t2-t1
        endif
        cycle
      endif
      if (myrank==0) then
        hop_event_states=EVENT_FORBIDDEN
        if (.not. locked) then
          hop_event_states(1:kmc%system%sites(site_id)%max_nn,1:max_jump_shell)=&
            kmc%system%sites(site_id)%hop_states(1:kmc%system%sites(site_id)%max_nn,1:max_jump_shell)
        endif
        if (allow_exchange) then
          exchange_event_states=EVENT_FORBIDDEN
          exchange_event_states(1:kmc%system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites,&
                                1:kmc%system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites)=&
            kmc%system%sites(site_id)%exchange_states(&
              1:kmc%system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites,&
              1:kmc%system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites)
        endif
        n_neighbour_sites=kmc%system%sites(site_id)%neighbour_shells(:)%n_neighbour_sites
      endif
      call mpi_bcast(hop_event_states,max_nn*max_jump_shell,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
      if (allow_exchange) then
        call mpi_bcast(exchange_event_states,max_1nn*max_1nn,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
      endif
      call mpi_bcast(n_neighbour_sites,NN_SHELLS,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
      call mpi_bcast(locked,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
  
      if (myrank==0) then
        nn_tmp=kmc%system%sites(1:kmc%system%n_sites)%nn
        n_near_tmp=kmc%system%sites(1:kmc%system%n_sites)%n_near
        do shell=1,TRUE_1NN_SHELL-1
          do nn_atom=1,kmc%system%sites(site_id)%neighbour_shells(shell)%n_neighbour_sites
            nn_id=kmc%system%sites(site_id)%neighbour_shells(shell)%neighbour_sites(nn_atom)
            n_near_tmp(nn_id)=n_near_tmp(nn_id)-1
          end do
        end do
        do nn_atom=1,kmc%system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
          nn_id=kmc%system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(nn_atom)
          nn_tmp(nn_id)=nn_tmp(nn_id)-1
        end do
      endif

      if (myrank==0) then
        site_blocked=.false.
        t2=mpi_wtime()
        kmc%get_rates_overhead_time=kmc%get_rates_overhead_time+t2-t1
      endif

      ! Hop events of site site_id
      if (.not. locked) then
        do jump_shell=1,max_jump_shell
          shell=JUMP_SHELLS(jump_shell)
          do nn_atom=1,n_neighbour_sites(shell)
            obstructed=.false.
            if (myrank==0) then
              t1=mpi_wtime()
              final_id=kmc%system%sites(site_id)%neighbour_shells(shell)%neighbour_sites(nn_atom)
              if (hop_event_states(nn_atom,jump_shell)==EVENT_VALID) then
                site_blocked(final_id)=.true.
              endif
              if (site_blocked(final_id)) then
                obstructed=.true.
                do nn_nn_atom=1,kmc%system%sites(final_id)%neighbour_shells(JUMP_SHELLS(1))%n_neighbour_sites
                  nn_nn_id=kmc%system%sites(final_id)%neighbour_shells(JUMP_SHELLS(1))%neighbour_sites(nn_nn_atom)
                  if (kmc%system%sites(nn_nn_id)%Z==0) then
                    site_blocked(nn_nn_id)=.true.
                  endif
                end do
              endif
            endif
            call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
            if (internal_error/=0) then
              if (present(error)) then
                error=internal_error
              endif
              return
            endif
            if (hop_event_states(nn_atom,jump_shell)==EVENT_VALID) then
              if (myrank==0) then
                t2=mpi_wtime()
                kmc%event_finding_time=kmc%event_finding_time+t2-t1
              endif
              cycle
            endif
            if (hop_event_states(nn_atom,jump_shell)==EVENT_FORBIDDEN) then
              if (myrank==0) then
                t2=mpi_wtime()
                kmc%event_finding_time=kmc%event_finding_time+t2-t1
              endif
              cycle
            endif
            if (myrank==0) then
              if (kmc%verbosity>=VERBOSITY_LOUD) then
                write(6,"(2(a,i0))") "# Getting rate for hop ",site_id,&
                                     " -> ",final_id
                t2=mpi_wtime()
                kmc%event_finding_time=kmc%event_finding_time+t2-t1
              endif
            endif
            call mpi_bcast(obstructed,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
            if (.not. obstructed) then
              if (myrank==0) then
                if (kmc%system%sites(final_id)%Z/=0) then
                  obstructed=.true.
                endif
                if (.not. obstructed) then
                  t1=mpi_wtime()
                  n_near_final=n_near_tmp(final_id)
                  nn_final=nn_tmp(final_id)
                  nn_initial=kmc%system%sites(site_id)%nn
                  if (nn_final<3 .or. &
                      n_near_final>0) then
                    obstructed=.true.
                  endif
                  ! Check minimization cache for a hit before doing a new minimization
                  if (.not. obstructed) then
                    t1=mpi_wtime()
                    call create_fingerprint(kmc,site_id,final_id,fingerprint,error=internal_error)
                    t2=mpi_wtime()
                    kmc%rate_fingerprint_time=kmc%rate_fingerprint_time+t2-t1
                    if (internal_error/=0) then
                      write(0,"(a)") "*** Error getting fingerprint in get_rates()"
                    endif
                    if (internal_error==0) then
                      t1=mpi_wtime()
                      minim_result=check_minim_cache(kmc,fingerprint,error=internal_error)
                      t2=mpi_wtime()
                      kmc%rate_cache_time=kmc%rate_cache_time+t2-t1
                      if (internal_error/=0) then
                        write(0,"(a)") "*** Error checking minimization cache in get_rates()"
                      endif
                      if (minim_result==MINIM_RESULT_FAIL) then
                        obstructed=.true.
                        if (kmc%verbosity>=VERBOSITY_LOUD) then
                          write(6,"(a)") "# Found a failed minimization in minim_cache"
                        endif
                      endif
                    endif
                  endif
                endif
              endif
            endif
            call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
            if (internal_error/=0) then
              if (present(error)) then
                error=internal_error
              endif
              return
            endif
            call mpi_bcast(obstructed,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
            if (.not. obstructed) then
              if (hop_event_states(nn_atom,jump_shell)==EVENT_LAE_CHANGED) then
                if (myrank==0) then
                  t1=mpi_wtime()
                  call only_substrate_lae(kmc,site_id,final_id,only_substrate_jump,error=internal_error)
                  t2=mpi_wtime()
                  kmc%lae_sub_time=kmc%lae_sub_time+t2-t1
                  if (internal_error/=0) then
                    write(0,"(a)") "*** Error determining only-substrate lae case in get_rates()"
                  else
                    if (only_substrate_jump==SUBSTRATE_NOT .or. kmc%substrate_fcc_hcp_barrier<0.0_rk &
                                                           .or. kmc%substrate_hcp_fcc_barrier<0.0_rk) then
                      t1=mpi_wtime()
                      call get_jump_lae(kmc,site_id,final_id,jump_type=EVENT_TYPE_HOP,&
                                        lae=lae,error=internal_error)
                      t2=mpi_wtime()
                      kmc%lae_time=kmc%lae_time+t2-t1
                      if (internal_error/=0) then
                        write(0,"(a)") "*** Error getting lae in get_rates()"
                      endif
                    endif
                    event_ok=.false.
                    if (only_substrate_jump==SUBSTRATE_FCC_HCP) then
                      if (kmc%substrate_fcc_hcp_barrier>=0.0_rk .and. kmc%substrate_hcp_fcc_barrier>=0.0_rk) then
                        barrier=kmc%substrate_fcc_hcp_barrier
                        reverse_barrier=kmc%substrate_hcp_fcc_barrier
                        event_ok=.true.
                      endif
                    else if (only_substrate_jump==SUBSTRATE_HCP_FCC) then
                      if (kmc%substrate_fcc_hcp_barrier>=0.0_rk .and. kmc%substrate_hcp_fcc_barrier>=0.0_rk) then
                        barrier=kmc%substrate_hcp_fcc_barrier
                        reverse_barrier=kmc%substrate_fcc_hcp_barrier
                        event_ok=.true.
                      endif
                    endif
                    nan_detected=.false.
                    if (.not. event_ok) then
                      t1=mpi_wtime()
                      call calc(kmc%desc,lae,desc_forward,do_descriptor=.true.,&
                                args_str="atom_mask_name=jumping",error=internal_error)
                      t2=mpi_wtime()
                      kmc%desc_time=kmc%desc_time+t2-t1
                      if (internal_error/=0) then
                        write(0,"(a)") "*** Error: couldn't calculate forward descriptor in get_rates()"
                        call write(lae,"lae_"//lae_index//"_1.xyz")
                      else
                        if (any(desc_forward%x(1)%data(:)/=desc_forward%x(1)%data(:))) then
                          write(0,"(a,i0)") "*** Warning: NaN values in forward descriptor of LAE ",lae_index
                          nan_detected=.true.
                          call write(lae,"lae_"//lae_index//"_1.xyz")
                        endif
                        t1=mpi_wtime()
                        lae%pos(:,1)=kmc%system%sites(final_id)%cartesian_coords(:)
                        call calc_connect(lae)
                        t2=mpi_wtime()
                        kmc%lae_time=kmc%lae_time+t2-t1

                        t1=mpi_wtime()
                        call calc(kmc%desc,lae,desc_reverse,do_descriptor=.true.,&
                                  args_str="atom_mask_name=jumping",error=internal_error)
                        t2=mpi_wtime()
                        kmc%desc_time=kmc%desc_time+t2-t1
                        if (internal_error/=0) then
                          write(0,"(a)") "*** Error: couldn't calculate reverse descriptor in get_rates()"
                          call write(lae,"lae_"//lae_index//"_2.xyz")
                        else
                          if (any(desc_reverse%x(1)%data(:)/=desc_reverse%x(1)%data(:))) then
                            write(0,"(a,i0)") "*** Warning: NaN values in reverse descriptor of LAE index ",lae_index
                            nan_detected=.true.
                            call write(lae,"lae_"//lae_index//"_2.xyz")
                          endif
                        endif
                      endif
                      call finalise(lae)
                      lae_index=lae_index+1
                    endif
                  endif
                endif
                call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
                if (internal_error/=0) then
                  if (present(error)) then
                    error=internal_error
                  endif
                  return
                endif
                call mpi_bcast(nan_detected,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
                call mpi_bcast(event_ok,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
                if (.not. nan_detected .and. .not. event_ok) then
                  if (myrank==0) then
                    t1=mpi_wtime()
                    final_config=initial_config
                    final_config%pos(:,kmc%system%sites(site_id)%id)=&
                      kmc%system%sites(final_id)%cartesian_coords(:)
                  endif
                  call get_barrier(kmc,lmp,site_id,final_id,nn_tmp,fingerprint,minim_result,&
                                   desc_forward,desc_reverse,initial_config,&
                                   final_config,barrier,reverse_barrier,event_ok,myrank,shell=shell,error=internal_error)
                  if (myrank==0) then
                    call finalise(desc_forward)
                    call finalise(desc_reverse)
                    call dealloc_atoms(final_config)
                    t2=mpi_wtime()
                    kmc%get_barrier_time=kmc%get_barrier_time+t2-t1
                  endif
                  if (internal_error/=0) then
                    if (myrank==0) then
                      write(0,"(a)") "*** Error: couldn't get barrier in get_rates()"
                    endif
                    if (present(error)) then
                      error=internal_error
                    endif
                    return
                  endif
                endif
              else
                nan_detected=.false.
                event_ok=.true.
              endif
            endif
            if (myrank==0) then
              t1=mpi_wtime()
              if (obstructed .or. nan_detected .or. .not. event_ok) then
                kmc%system%sites(site_id)%hop_states(nn_atom,jump_shell)=EVENT_FORBIDDEN
                if (kmc%verbosity>=VERBOSITY_LOUD) then
                  write(6,"(a)") "# Event marked forbidden"
                  if (kmc%debug .and. kmc%verbosity>=VERBOSITY_ABSURD) then
                    if (obstructed) then
                      write(6,"(a)") "# Reason: obstructed site"
                    elseif (nan_detected) then
                      write(6,"(a)") "# Reason: NaN detected"
                    elseif (.not. event_ok) then
                      write(6,"(a)") "# Reason: get_barrier() returned event_ok=FALSE"
                    endif
                  endif
                endif
              else
                do nn_nn_atom=1,kmc%system%sites(final_id)%neighbour_shells(JUMP_SHELLS(1))%n_neighbour_sites
                  nn_id=kmc%system%sites(final_id)%neighbour_shells(JUMP_SHELLS(1))%neighbour_sites(nn_nn_atom)
                  if (kmc%system%sites(nn_id)%Z==0) then
                    site_blocked(nn_id)=.true.
                  endif
                end do
                if (kmc%system%sites(site_id)%hop_states(nn_atom,jump_shell)==EVENT_LAE_CHANGED) then
                  barrier=max(0.0_rk,barrier)
                  if (barrier>0) then
                    reverse_barrier=max(0.0_rk,reverse_barrier)
                  else
                    reverse_barrier=100.0_rk
                  endif
                  kmc%system%sites(site_id)%hop_barriers(1,nn_atom,jump_shell)=barrier
                  kmc%system%sites(site_id)%hop_barriers(2,nn_atom,jump_shell)=reverse_barrier
                  kmc%system%sites(site_id)%hop_rates(nn_atom,jump_shell)=&
                    kmc%attempt_frequency*exp(-real(barrier,kind=dp)/kmc%kT)
                endif
                kmc%system%sites(site_id)%hop_states(nn_atom,jump_shell)=EVENT_VALID
  
                if (kmc%verbosity>=VERBOSITY_LOUD) then
                  write(6,"(a,f10.6,a)") "# E_{m} = ",&
                    kmc%system%sites(site_id)%hop_barriers(1,nn_atom,jump_shell)," eV"
                  write(6,"(a,f10.6,a)") "# E_{m,reverse} = ",&
                    kmc%system%sites(site_id)%hop_barriers(2,nn_atom,jump_shell)," eV"
                  write(6,"(a,g10.3,a)") "# rate = ",kmc%system%sites(site_id)%hop_rates(nn_atom,jump_shell)," Hz"
                endif
                if (only_substrate_jump==SUBSTRATE_FCC_HCP) then
                  if (kmc%substrate_fcc_hcp_barrier<0.0_rk .and. kmc%substrate_hcp_fcc_barrier<0.0_rk) then
                    kmc%substrate_fcc_hcp_barrier=barrier
                    kmc%substrate_hcp_fcc_barrier=reverse_barrier
                  endif
                elseif (only_substrate_jump==SUBSTRATE_HCP_FCC) then
                  if (kmc%substrate_fcc_hcp_barrier<0.0_rk .and. kmc%substrate_hcp_fcc_barrier<0.0_rk) then
                    kmc%substrate_hcp_fcc_barrier=barrier
                    kmc%substrate_fcc_hcp_barrier=reverse_barrier
                  endif
                endif
              endif
              t2=mpi_wtime()
              kmc%assigning_time=kmc%assigning_time+t2-t1
            endif
            call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
            if (internal_error/=0) then
              if (present(error)) then
                error=internal_error
              endif
              return
            endif
          end do
        end do
      endif
      if (myrank==0) then
        if (any(kmc%system%sites(site_id)%hop_states(:,:)==EVENT_LAE_CHANGED)) then
          write(0,"(a,i0,a)") "*** Error: site ",site_id," has unchecked hop events in get_rates()"
          internal_error=1
        endif
        if (internal_error==0) then
          kmc%system%sites(site_id)%hop_rate=sum(&
            kmc%system%sites(site_id)%hop_rates(:,:),&
            mask=kmc%system%sites(site_id)%hop_states(:,:)==EVENT_VALID)
        endif
      endif
      call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
      if (internal_error/=0) then
        if (present(error)) then
          error=internal_error
        endif
        return
      endif

      call mpi_bcast(allow_exchange,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
      if (.not. allow_exchange) then
        cycle
      endif

      ! Exchange events
      do initial_site=1,n_neighbour_sites(TRUE_1NN_SHELL)
        obstructed=.false.
        if (myrank==0) then
          t1=mpi_wtime()
          initial_id=kmc%system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(initial_site)
          if (kmc%system%sites(initial_id)%substrate) then
            obstructed=.true.
          elseif (kmc%system%sites(initial_id)%Z<=0) then
            obstructed=.true.
          elseif (kmc%system%sites(initial_id)%locked) then
            obstructed=.true.
          endif
          if (.not. obstructed) then
            nn_tmp=kmc%system%sites(1:kmc%system%n_sites)%nn
            n_near_tmp=kmc%system%sites(1:kmc%system%n_sites)%n_near
            do shell=1,TRUE_1NN_SHELL-1
              do nn_atom=1,kmc%system%sites(initial_id)%neighbour_shells(shell)%n_neighbour_sites
              n_near_tmp(kmc%system%sites(initial_id)%neighbour_shells(shell)%neighbour_sites(nn_atom)) &
                         =n_near_tmp(kmc%system%sites(initial_id)%neighbour_shells(shell)%neighbour_sites(nn_atom))-1
              end do
            end do
            do nn_atom=1,kmc%system%sites(initial_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
              nn_tmp(kmc%system%sites(initial_id)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(nn_atom)) &
                     =nn_tmp(kmc%system%sites(initial_id)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(nn_atom))-1
            end do
          endif
        endif
        call mpi_bcast(obstructed,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
        if (obstructed) then
          if (myrank==0) then
            if (kmc%debug .and. kmc%verbosity>=VERBOSITY_ABSURD) then
              write(6,"(a,i0,a)") "# All exchange events from initial_site ",&
                                   initial_id," marked forbidden. Reason:"
              if (kmc%system%sites(initial_id)%substrate) then
                write(6,"(a)") "# initial_site is substrate"
              elseif (kmc%system%sites(initial_id)%Z<=0) then
                write(6,"(a)") "# initial_site is not an atom"
              elseif (kmc%system%sites(initial_id)%locked) then
                write(6,"(a)") "# initial_site is locked"
              else
                write(6,"(a)") "# other. Bug?"
              endif
            endif
            kmc%system%sites(site_id)%exchange_states(:,initial_site)=EVENT_FORBIDDEN
          endif
          exchange_event_states(:,initial_site)=EVENT_FORBIDDEN
          if (myrank==0) then
            t2=mpi_wtime()
            kmc%event_finding_time=kmc%event_finding_time+t2-t1
          endif
          cycle
        endif
        do final_site=1,n_neighbour_sites(TRUE_1NN_SHELL)
          if (myrank==0) then
            t1=mpi_wtime()
          endif
          if (exchange_event_states(final_site,initial_site)==EVENT_VALID) then
            if (myrank==0) then
              t2=mpi_wtime()
              kmc%event_finding_time=kmc%event_finding_time+t2-t1
            endif
            cycle
          endif
          if (exchange_event_states(final_site,initial_site)==EVENT_FORBIDDEN) then
            if (myrank==0) then
              t2=mpi_wtime()
              kmc%event_finding_time=kmc%event_finding_time+t2-t1
            endif
            cycle
          endif
          t2=mpi_wtime()
          kmc%event_finding_time=kmc%event_finding_time+t2-t1
          obstructed=.false.
          if (myrank==0) then
            t1=mpi_wtime()
            final_id=kmc%system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(final_site)
            if (kmc%verbosity>=VERBOSITY_LOUD) then
              write(6,"(3(a,i0))") "# Getting rate for exchange ",&
                                   initial_id," -> ",site_id," -> ",final_id
            endif
            nn_final=nn_tmp(final_id)
            nn_initial=kmc%system%sites(initial_id)%nn
            nn_mid=kmc%system%sites(site_id)%nn
            if (kmc%system%sites(final_id)%Z/=0) then
              if (kmc%debug .and. kmc%verbosity>=VERBOSITY_ABSURD) then
                write(6,"(a)") "# Not a valid adsorption site"
              endif
              obstructed=.true.
            elseif (nn_final<3) then
              if (kmc%debug .and. kmc%verbosity>=VERBOSITY_ABSURD) then
                write(6,"(a)") "# Not a stable adsorption site - nn<3"
              endif
              obstructed=.true.
            elseif (n_near_tmp(final_id)>0) then
              if (kmc%debug .and. kmc%verbosity>=VERBOSITY_ABSURD) then
                write(6,"(a)") "# Not a stable adsorption site - n_near>0"
              endif
              obstructed=.true.
            endif
            ! Check minimization cache for a hit before doing a new minimization
            if (.not. obstructed) then
              t1=mpi_wtime()
              call create_fingerprint(kmc,site_id,final_id,fingerprint,initial_id=initial_id,error=internal_error)
              t2=mpi_wtime()
              kmc%rate_fingerprint_time=kmc%rate_fingerprint_time+t2-t1
              if (internal_error/=0) then
                write(0,"(a)") "*** Error getting fingerprint in get_rates()"
              endif
              if (internal_error==0) then
                t1=mpi_wtime()
                minim_result=check_minim_cache(kmc,fingerprint,error=internal_error)
                t2=mpi_wtime()
                kmc%rate_cache_time=kmc%rate_cache_time+t2-t1
                if (internal_error/=0) then
                  write(0,"(a)") "*** Error checking minimization cache in get_rates()"
                endif
                if (minim_result==MINIM_RESULT_FAIL) then
                  obstructed=.true.
                  if (kmc%verbosity>=VERBOSITY_LOUD) then
                    write(6,"(a)") "# Found a failed minimization in minim_cache"
                  endif
                endif
              endif
            endif
          endif
          call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
          if (internal_error/=0) then
            if (present(error)) then
              error=internal_error
            endif
            return
          endif
          call mpi_bcast(obstructed,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
          if (.not. obstructed) then
            if (exchange_event_states(final_site,initial_site)==EVENT_LAE_CHANGED) then
              if (myrank==0) then
                t1=mpi_wtime()
                call get_jump_lae(kmc,site_id,final_id,jump_type=EVENT_TYPE_EXCHANGE,&
                                  lae=lae,jumping_ptr=jumping_ptr,&
                                  initial_id=initial_id,&
                                  error=internal_error)
                t2=mpi_wtime()
                kmc%lae_time=kmc%lae_time+t2-t1
                if (internal_error/=0) then
                  write(0,"(a)") "*** Error getting lae in get_rates()"
                  if (present(error)) then
                    error=internal_error
                  endif
                endif
                nan_detected=.false.
                t1=mpi_wtime()
                call calc(kmc%desc,lae,desc_forward,do_descriptor=.true.,&
                          args_str="atom_mask_name=jumping",error=internal_error)
                t2=mpi_wtime()
                kmc%desc_time=kmc%desc_time+t2-t1
                if (internal_error/=0) then
                  write(0,"(a)") "*** Error: couldn't calculate forward descriptor in get_rates()"
                  if (present(error)) then
                    error=internal_error
                  endif
                else
                  if (any(desc_forward%x(1)%data(:)/=desc_forward%x(1)%data(:))) then
                    write(0,"(a,i0)") "*** Warning: NaN values in forward descriptor of LAE index ",lae_index
                    call write(lae,"lae_"//lae_index//"_1.xyz")
                    nan_detected=.true.
                  endif
                  t1=mpi_wtime()
                  lae%pos(:,1)=kmc%system%sites(site_id)%cartesian_coords(:)
                  lae%pos(:,2)=kmc%system%sites(final_id)%cartesian_coords(:)
                  jumping_ptr(1)=.false.
                  jumping_ptr(2)=.true.
                  call calc_connect(lae)
                  t2=mpi_wtime()
                  kmc%lae_time=kmc%lae_time+t2-t1

                  t1=mpi_wtime()
                  call calc(kmc%desc,lae,desc_reverse,do_descriptor=.true.,args_str="atom_mask_name=jumping",error=internal_error)
                  t2=mpi_wtime()
                  kmc%desc_time=kmc%desc_time+t2-t1
                  if (internal_error/=0) then
                    write(0,"(a)") "*** Error: couldn't calculate reverse descriptor in get_rates()"
                    if (present(error)) then
                      error=internal_error
                    endif
                  else
                    if (any(desc_reverse%x(1)%data(:)/=desc_reverse%x(1)%data(:))) then
                      write(0,"(a,i0)") "*** Warning: NaN values in reverse descriptor of LAE index ",lae_index
                      call write(lae,"lae_"//lae_index//"_2.xyz")
                      nan_detected=.true.
                    endif
                  endif
                endif
                call finalise(lae)
                nullify(jumping_ptr)
                lae_index=lae_index+1
              endif
              call mpi_bcast(nan_detected,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
              if (.not. nan_detected) then
                if (myrank==0) then
                  t1=mpi_wtime()
                  final_config=initial_config
                  final_config%pos(:,kmc%system%sites(site_id)%id)=&
                    kmc%system%sites(final_id)%cartesian_coords(:)
                  final_config%pos(:,kmc%system%sites(initial_id)%id)=&
                    kmc%system%sites(site_id)%cartesian_coords(:)
                endif
                call get_barrier(kmc,lmp,site_id,final_id,nn_tmp,fingerprint,minim_result,&
                                 desc_forward,desc_reverse,initial_config,&
                                 final_config,barrier,reverse_barrier,event_ok,myrank,&
                                 initial_id=initial_id,&
                                 error=internal_error)
                if (myrank==0) then
                  call finalise(desc_forward)
                  call finalise(desc_reverse)
                  call dealloc_atoms(final_config)
                  t2=mpi_wtime()
                  kmc%get_barrier_time=kmc%get_barrier_time+t2-t1
                endif
                if (internal_error/=0) then
                  if (myrank==0) then
                    write(0,"(a)") "*** Error: couldn't get barrier in get_rates()"
                  endif
                  if (present(error)) then
                    error=internal_error
                  endif
                  return
                endif
              endif
            else
              event_ok=.true.
              nan_detected=.false.
            endif
          endif
          if (myrank==0) then
            t1=mpi_wtime()
            if (obstructed .or. nan_detected .or. .not. event_ok) then
              kmc%system%sites(site_id)%exchange_states(final_site,initial_site)=EVENT_FORBIDDEN
              if (kmc%verbosity>=VERBOSITY_LOUD) then
                write(6,"(a)") "# Event marked forbidden"
                if (kmc%debug .and. kmc%verbosity>=VERBOSITY_ABSURD) then
                  if (obstructed) then
                    write(6,"(a)") "# Reason: obstructed site"
                  elseif (nan_detected) then
                    write(6,"(a)") "# Reason: NaN detected"
                  elseif (.not. event_ok) then
                    write(6,"(a)") "# Reason: get_barrier() returned event_ok=FALSE"
                  endif
                endif
              endif
            else
              if (kmc%system%sites(site_id)%exchange_states(final_site,initial_site)==EVENT_LAE_CHANGED) then
                barrier=max(0.0_rk,barrier)
                if (barrier>0) then
                  reverse_barrier=max(0.0_rk,reverse_barrier)
                else
                  reverse_barrier=100.0_rk
                endif
                kmc%system%sites(site_id)%exchange_barriers(1,final_site,initial_site)=barrier
                kmc%system%sites(site_id)%exchange_barriers(2,final_site,initial_site)=reverse_barrier
                kmc%system%sites(site_id)%exchange_rates(final_site,initial_site)=&
                  kmc%attempt_frequency*exp(-real(barrier,kind=dp)/kmc%kT)
              endif
              kmc%system%sites(site_id)%exchange_states(final_site,initial_site)=EVENT_VALID

              if (kmc%verbosity>=VERBOSITY_LOUD) then
                write(6,"(a,f10.6,a)") "# E_{m} = ",&
                  kmc%system%sites(site_id)%exchange_barriers(1,final_site,initial_site)," eV"
                write(6,"(a,f10.6,a)") "# E_{m,reverse} = ",&
                  kmc%system%sites(site_id)%exchange_barriers(2,final_site,initial_site)," eV"
                write(6,"(a,g10.3,a)") "# rate = ",kmc%system%sites(site_id)%exchange_rates(final_site,initial_site)," Hz"
              endif
            endif
            t2=mpi_wtime()
            kmc%assigning_time=kmc%assigning_time+t2-t1
          endif
        end do
      end do

      if (myrank==0) then
        if (any(kmc%system%sites(site_id)%exchange_states(:,:)==EVENT_LAE_CHANGED)) then
          write(0,"(a,i0,a)") "*** Error: site ",site_id," has unchecked exchange events in get_rates()"
          internal_error=1
        endif
        if (internal_error==0) then
          kmc%system%sites(site_id)%exchange_rate=sum(&
            kmc%system%sites(site_id)%exchange_rates(:,:),&
            mask=kmc%system%sites(site_id)%exchange_states(:,:)==EVENT_VALID)
          kmc%system%sites(site_id)%jumps_checked=.true.
        endif
      endif
      call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
      if (internal_error/=0) then
        if (present(error)) then
          error=internal_error
        endif
        return
      endif

    end do

    if (myrank==0) then
      t1=mpi_wtime()
    endif

    deallocate(hop_event_states)
    if (allow_exchange) then
      deallocate(exchange_event_states)
    endif
    if (myrank==0) then
      deallocate(nn_tmp)
      deallocate(n_near_tmp)
      deallocate(site_blocked)
      kmc%total_rate=RATE_LOWERING_FACTOR**kmc%rate_lowered_times&
                       *sum(kmc%system%sites(1:kmc%system%n_sites)%hop_rate)+sum(kmc%dep_rate)
      if (kmc%system%allow_exchange) then
        kmc%total_rate=kmc%total_rate+RATE_LOWERING_FACTOR**kmc%rate_lowered_times&
                         *sum(kmc%system%sites(1:kmc%system%n_sites)%exchange_rate)
      endif
      call dealloc_atoms(initial_config)
    endif
    if (present(error)) then
      call mpi_bcast(error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
    endif

    if (myrank==0) then
      t2=mpi_wtime()
      kmc%get_rates_overhead_time=kmc%get_rates_overhead_time+t2-t1
    endif

  end subroutine get_rates

  subroutine get_barrier(kmc,lmp,site_id,final_id,&
                         nn_tmp,fingerprint,minim_result,desc_forward,desc_reverse,init,&
                         final,forward_barrier,reverse_barrier,event_ok,myrank,&
                         initial_id,shell,error)
    type(t_kmc),intent(inout) :: kmc
    type(lammps),intent(inout) :: lmp
    integer,intent(in) :: site_id,final_id
    integer,intent(in) :: nn_tmp(kmc%system%n_sites)
    integer,intent(in) :: fingerprint(kmc%nZ,0:NN_SHELLS),minim_result
    type(descriptor_data) :: desc_forward,desc_reverse
    type(t_atoms),intent(inout) :: init,final
    real(kind=rk),intent(out) :: forward_barrier,reverse_barrier
    logical,intent(out) :: event_ok
    integer,intent(in) :: myrank
    integer,optional,intent(in) :: initial_id,shell
    integer,optional,intent(out) :: error

    type(t_atoms) :: final_config_relaxed
    double precision :: ml_tol2,zero_tol
    double precision :: neb_forward_barrier,neb_reverse_barrier
    double precision :: ml_forward_barrier,ml_reverse_barrier
    double precision :: forward_variance,reverse_variance
    real(kind=rk) :: displacement2
    double precision,allocatable :: xdata(:)
    integer :: i,internal_error
    double precision :: t1,t2
    integer :: coord_jumping_forward,coord_target_forward,coord_middle
    integer :: coord_jumping_reverse,coord_target_reverse
    integer :: coord_jumping_forward_substrate,coord_target_forward_substrate,coord_middle_substrate
    integer :: coord_jumping_reverse_substrate,coord_target_reverse_substrate
    logical :: final_relaxed,ml_forward_called,ml_reverse_called,forward_ml_full,reverse_ml_full
    integer :: jump_shell,nn_shell
    integer :: check_list(25),n_check_list,nn_atom,nn_id,id
    integer :: mpi_error
    logical :: forward_data_exists,reverse_data_exists,in_check_list

    if (present(error)) then
      error=0
    endif
    internal_error=0
    event_ok=.false.

    if (myrank==0) then
      ml_tol2=kmc%ml_tol*kmc%ml_tol
      zero_tol=kmc%neb_error
      ml_forward_called=.false.
      ml_reverse_called=.false.
      forward_data_exists=.false.
      reverse_data_exists=.false.
      forward_ml_full=.false.
      reverse_ml_full=.false.
      final_relaxed=.false.
      if (present(shell)) then
        jump_shell=NN_SHELL_TO_JUMP_SHELL(shell)
        if (jump_shell<=0) then
          write(0,"(a,i0,a)") "*** Error: invalid shell ",shell," doesn't map to an allowed jump shell in get_barrier()"
          internal_error=1
        endif
      endif
      coord_target_forward=nn_tmp(final_id)
      coord_target_forward_substrate=kmc%system%sites(final_id)%nn_substrate
      if (present(initial_id)) then ! Exchange event
        coord_jumping_forward=kmc%system%sites(initial_id)%nn
        coord_middle=kmc%system%sites(site_id)%nn
        coord_target_reverse=nn_tmp(initial_id)

        coord_jumping_forward_substrate=kmc%system%sites(initial_id)%nn_substrate
        coord_middle_substrate=kmc%system%sites(site_id)%nn_substrate
        coord_target_reverse_substrate=kmc%system%sites(initial_id)%nn_substrate
      else ! Hop event
        coord_jumping_forward=kmc%system%sites(site_id)%nn
        coord_target_reverse=nn_tmp(site_id)

        coord_jumping_forward_substrate=kmc%system%sites(site_id)%nn_substrate
        coord_target_reverse_substrate=kmc%system%sites(site_id)%nn_substrate
      endif
      coord_jumping_reverse=nn_tmp(final_id)
      coord_jumping_reverse_substrate=kmc%system%sites(final_id)%nn_substrate

      if (kmc%verbosity>=VERBOSITY_ABSURD) then
        write(6,"(a,i0)") "# Coordination of the jumping atom: ",coord_jumping_forward
        write(6,"(a,i0)") "#                 for reverse jump: ",coord_jumping_reverse
        write(6,"(a,i0)") "# Coordination of the target site: ",coord_target_forward
        write(6,"(a,i0)") "#                for reverse jump: ",coord_target_reverse
        if (present(initial_id)) then
          write(6,"(a,i0)") "# Coordination of the middle atom: ",coord_middle
        else
          write(6,"(a,i0)") "# Jump distance shell: ",shell
        endif
        if (kmc%system%substrate_type==SUBSTRATE_TYPE_WEAK) then
          write(6,"(a,i0)") "# Substrate coordination of the jumping atom: ",coord_jumping_forward_substrate
          write(6,"(a,i0)") "#                           for reverse jump: ",coord_jumping_reverse_substrate
          write(6,"(a,i0)") "# Substrate coordination of the target site: ",coord_target_forward_substrate
          write(6,"(a,i0)") "#                          for reverse jump: ",coord_target_reverse_substrate
          if (present(initial_id)) then
            write(6,"(a,i0)") "# Substrate coordination of the middle atom: ",coord_middle_substrate
          endif
        endif
      endif
    endif
    call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
    if (internal_error/=0) then
      if (present(error)) then
        error=1
      endif
      return
    endif

    if (myrank==0) then
      if (.not. present(initial_id)) then
        if (.not. kmc%hop_predictor(coord_jumping_forward,coord_target_forward,jump_shell)%initialised) then
          call initialise_predictor(kmc%hop_predictor(coord_jumping_forward,coord_target_forward,jump_shell),&
                          kmc%n_data_min,kmc%n_data_max,kmc%n_sparse_max,kmc%desc_str,delta=kmc%delta,zeta=kmc%zeta,&
                          covariance_type=kmc%covariance_type,regularisation=kmc%regularisation)
        endif
        if (.not. kmc%hop_predictor(coord_jumping_reverse,coord_target_reverse,jump_shell)%initialised) then
          call initialise_predictor(kmc%hop_predictor(coord_jumping_reverse,coord_target_reverse,jump_shell),&
                          kmc%n_data_min,kmc%n_data_max,kmc%n_sparse_max,kmc%desc_str,delta=kmc%delta,zeta=kmc%zeta,&
                          covariance_type=kmc%covariance_type,regularisation=kmc%regularisation)
        endif
      else
        if (.not. kmc%exchange_predictor(coord_middle,coord_jumping_forward,&
                                         coord_target_forward)%initialised) then
          call initialise_predictor(kmc%exchange_predictor(coord_middle,coord_jumping_forward,&
                                                           coord_target_forward),&
                          kmc%n_data_min,kmc%n_data_max,kmc%n_sparse_max,kmc%desc_str,delta=kmc%delta,zeta=kmc%zeta,&
                          covariance_type=kmc%covariance_type,regularisation=kmc%regularisation)
        endif
        if (.not. kmc%exchange_predictor(coord_middle,coord_jumping_reverse,&
                                         coord_target_reverse)%initialised) then
          call initialise_predictor(kmc%exchange_predictor(coord_middle,coord_jumping_reverse,&
                                                           coord_target_reverse),&
                          kmc%n_data_min,kmc%n_data_max,kmc%n_sparse_max,kmc%desc_str,delta=kmc%delta,zeta=kmc%zeta,&
                          covariance_type=kmc%covariance_type,regularisation=kmc%regularisation)
        endif
      endif
    endif

    call mpi_bcast(ml_tol2,1,MPI_DOUBLE,0,MPI_COMM_WORLD,mpi_error)
    call mpi_bcast(zero_tol,1,MPI_DOUBLE,0,MPI_COMM_WORLD,mpi_error)

    if (myrank==0) then
      if (.not. present(initial_id)) then ! hop event
        if (kmc%hop_predictor(coord_jumping_forward,coord_target_forward,jump_shell)%fitted) then
          t1=mpi_wtime()
          call predict(kmc%hop_predictor(coord_jumping_forward,coord_target_forward,jump_shell),&
                       desc_forward%x(1)%data(:),ml_forward_barrier,forward_variance,error=internal_error)
          t2=mpi_wtime()
          kmc%predict_time=kmc%predict_time+t2-t1
          if (internal_error/=0) then
            write(0,"(a)") "*** Error: couldn't get forward barrier prediction in get_barrier()"
          endif
          ml_forward_called=.true.
        endif
        if (internal_error==0) then
          if (kmc%hop_predictor(coord_jumping_reverse,coord_target_reverse,jump_shell)%fitted) then
            t1=mpi_wtime()
            call predict(kmc%hop_predictor(coord_jumping_reverse,coord_target_reverse,jump_shell),&
                         desc_reverse%x(1)%data(:),ml_reverse_barrier,reverse_variance,error=internal_error)
            t2=mpi_wtime()
            kmc%predict_time=kmc%predict_time+t2-t1
            if (internal_error/=0) then
              write(0,"(a)") "*** Error: couldn't get reverse barrier prediction in get_barrier()"
            endif
            ml_reverse_called=.true.
          endif
        endif
      else ! exchange event
        if (kmc%exchange_predictor(coord_middle,coord_jumping_forward,coord_target_forward)%fitted) then
          t1=mpi_wtime()
          call predict(kmc%exchange_predictor(coord_middle,coord_jumping_forward,&
                                              coord_target_forward),desc_forward%x(1)%data(:),&
                       ml_forward_barrier,forward_variance,error=internal_error)
          t2=mpi_wtime()
          kmc%predict_time=kmc%predict_time+t2-t1
          if (internal_error/=0) then
            write(0,"(a)") "*** Error: couldn't get forward barrier prediction in get_barrier()"
          endif
          ml_forward_called=.true.
        endif
        if (internal_error==0) then
          if (kmc%exchange_predictor(coord_middle,coord_jumping_reverse,coord_target_reverse)%fitted) then
            t1=mpi_wtime()
            call predict(kmc%exchange_predictor(coord_middle,coord_jumping_reverse,&
                                                coord_target_reverse),desc_reverse%x(1)%data(:),&
                         ml_reverse_barrier,reverse_variance,error=internal_error)
            t2=mpi_wtime()
            kmc%predict_time=kmc%predict_time+t2-t1
            if (internal_error/=0) then
              write(0,"(a)") "*** Error: couldn't get reverse barrier prediction in get_barrier()"
            endif
            ml_reverse_called=.true.
          endif
        endif
      endif
      if (internal_error==0) then
        if (ml_forward_called .and. kmc%verbosity>=VERBOSITY_ABSURD) then
          write(6,"(a,f10.6,a,f10.6,a)") "# E_{ML} = ",ml_forward_barrier," eV, stdev = ",sqrt(forward_variance)," eV"
        endif
        if (ml_reverse_called .and. kmc%verbosity>=VERBOSITY_ABSURD) then
          write(6,"(a,f10.6,a,f10.6,a)") "# E_{ML,reverse} = ",ml_reverse_barrier," eV, stdev = ",sqrt(reverse_variance)," eV"
        endif
      endif
    endif
    call mpi_bcast(internal_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
    if (internal_error/=0) then
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

    call mpi_bcast(ml_forward_called,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
    call mpi_bcast(ml_reverse_called,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
    call mpi_bcast(ml_forward_barrier,1,MPI_DOUBLE,0,MPI_COMM_WORLD,mpi_error)
    call mpi_bcast(forward_variance,1,MPI_DOUBLE,0,MPI_COMM_WORLD,mpi_error)
    call mpi_bcast(ml_reverse_barrier,1,MPI_DOUBLE,0,MPI_COMM_WORLD,mpi_error)
    call mpi_bcast(reverse_variance,1,MPI_DOUBLE,0,MPI_COMM_WORLD,mpi_error)

    if (ml_forward_called .and. ml_reverse_called &
        .and. forward_variance<ml_tol2 .and. ml_forward_barrier==ml_forward_barrier .and. ml_forward_barrier>=-zero_tol &
        .and. reverse_variance<ml_tol2 .and. ml_reverse_barrier==ml_reverse_barrier .and. ml_reverse_barrier>=-zero_tol) then
      if (myrank==0) then
        forward_barrier=ml_forward_barrier
        reverse_barrier=ml_reverse_barrier
        if (reverse_barrier>0.0_rk) then
          event_ok=.true.
        endif
      endif
    else
      if (myrank==0) then
        if (.not. present(initial_id)) then ! hop event
          forward_data_exists=data_point_exists(kmc%hop_predictor(coord_jumping_forward,coord_target_forward,jump_shell),&
                                                desc_forward%x(1)%data(:),y=neb_forward_barrier)
          reverse_data_exists=data_point_exists(kmc%hop_predictor(coord_jumping_reverse,coord_target_reverse,jump_shell),&
                                                desc_reverse%x(1)%data(:),y=neb_reverse_barrier)
          forward_ml_full=kmc%hop_predictor(coord_jumping_forward,coord_target_forward,jump_shell)%n_y&
                            >=kmc%hop_predictor(coord_jumping_forward,coord_target_forward,jump_shell)%n_y_max
          reverse_ml_full=kmc%hop_predictor(coord_jumping_reverse,coord_target_reverse,jump_shell)%n_y&
                            >=kmc%hop_predictor(coord_jumping_reverse,coord_target_reverse,jump_shell)%n_y_max
        else
          forward_data_exists=data_point_exists(kmc%exchange_predictor(coord_middle,coord_jumping_forward,&
                                                                       coord_target_forward),&
                                                desc_forward%x(1)%data(:),y=neb_forward_barrier)
          reverse_data_exists=data_point_exists(kmc%exchange_predictor(coord_middle,coord_jumping_reverse,&
                                                                       coord_target_reverse),&
                                                desc_reverse%x(1)%data(:),y=neb_reverse_barrier)
          forward_ml_full=kmc%exchange_predictor(coord_middle,coord_jumping_forward,coord_target_forward)%n_y&
                            >=kmc%exchange_predictor(coord_middle,coord_jumping_forward,coord_target_forward)%n_y_max
          reverse_ml_full=kmc%exchange_predictor(coord_middle,coord_jumping_reverse,coord_target_reverse)%n_y&
                            >=kmc%exchange_predictor(coord_middle,coord_jumping_reverse,coord_target_reverse)%n_y
        endif
        if (forward_data_exists .and. reverse_data_exists) then
          if (kmc%verbosity>=VERBOSITY_ABSURD) then
            write(6,"(a)") "# ML is uncertain, but both forward and the reverse descriptor are already in data set"
            write(6,"(a)") "# Taking barriers directly from the data set"
          endif
        elseif (forward_ml_full .and. reverse_ml_full) then
          if (kmc%verbosity>=VERBOSITY_ABSURD) then
            write(6,"(a)") "# ML is uncertain, no data in data set but both forward and the reverse models are full"
          endif
          if (ml_forward_called .and. ml_reverse_called &
                .and. ml_forward_barrier==ml_forward_barrier &
                .and. ml_reverse_barrier==ml_reverse_barrier) then
            forward_barrier=ml_forward_barrier
            reverse_barrier=ml_reverse_barrier
            if (reverse_barrier>0.0_rk) then
              event_ok=.true.
            endif
            if (kmc%verbosity>=VERBOSITY_ABSURD) then
              write(6,"(a)") "# Accepting uncertain data point"
            endif
          else
            forward_barrier=100.0_rk
            reverse_barrier=100.0_rk
            event_ok=.false.
            if (kmc%verbosity>=VERBOSITY_ABSURD) then
              write(6,"(a)") "# ML wasn't called or were called with NaN result(s), forbidding jump"
            endif
          endif
        endif
        if ((.not. forward_data_exists .and. .not. forward_ml_full) .or. &
            (.not. reverse_data_exists .and. .not. reverse_ml_full)) then
          if (kmc%verbosity>=VERBOSITY_NORMAL) then
            write(6,"(a)") "# Getting barrier from NEB..."
          endif
          if (minim_result==MINIM_RESULT_SUCCESS) then
            final_relaxed=.true.
          elseif (minim_result==MINIM_RESULT_FAIL) then
            forward_barrier=100.0_rk
            reverse_barrier=100.0_rk
            write (0,"(a)") "*** Warning: reached minimization with a failed minim_cache result...?"
          else
            if (kmc%verbosity>=VERBOSITY_NORMAL) then
              write(6,"(a)") "# Minimizing final configuration..."
            endif
            kmc%minim_cache_index=kmc%minim_cache_index+1
            if (kmc%minim_cache_index>MINIM_CACHE_SIZE) then
              kmc%minim_cache_index=1
            endif
            kmc%minim_cache(:,:,kmc%minim_cache_index)=fingerprint
            call lmp%create_atoms(id=[(i,i=kmc%system%substrate%first_id,kmc%system%substrate%last_id)],&
                                  type=kmc%system%substrate%atom_type,&
                                  x=reshape(real(kmc%system%substrate%pos,kind=dp),&
                                  [3*(kmc%system%substrate%last_id-kmc%system%substrate%first_id+1)]))
            call lmp%create_atoms(id=[(i,i=final%first_id,final%last_id)],&
                                  type=final%atom_type,&
                                  x=reshape(real(final%pos,kind=dp),[3*(final%last_id-final%first_id+1)]))
            call lmp%command("group mobile subtract all frozen")
            t1=mpi_wtime()
            call lmp%command("minimize "//kmc%relax_tol//" "//0.1d0*kmc%relax_tol//" "//kmc%neb_steps//" "//10*kmc%neb_steps)
            ! Increment minimizations counter by one
            kmc%minimizations=kmc%minimizations+1
            t2=mpi_wtime()
            kmc%minimize_time=kmc%minimize_time+t2-t1
            call lmp%gather_atoms_subset("x",3,[(i,i=final%first_id,final%last_id)],xdata)
            call lmp%command("delete_atoms group mobile")
            final_config_relaxed=final
            final_config_relaxed%pos=reshape(xdata,shape(final_config_relaxed%pos))
            deallocate(xdata)
            n_check_list=1
            check_list(1)=kmc%system%sites(site_id)%id
            if (present(initial_id)) then
              n_check_list=n_check_list+1
              check_list(n_check_list)=kmc%system%sites(initial_id)%id
            endif
            do nn_shell=TRUE_1NN_SHELL,CLOSE_NEIGH_SHELL
              do nn_atom=1,kmc%system%sites(site_id)%neighbour_shells(nn_shell)%n_neighbour_sites
                nn_id=kmc%system%sites(site_id)%neighbour_shells(nn_shell)%neighbour_sites(nn_atom)
                if (kmc%system%sites(nn_id)%Z>0 .and. .not. kmc%system%sites(nn_id)%substrate) then
                  id=kmc%system%sites(nn_id)%id
                  in_check_list=.false.
                  do i=1,n_check_list
                    if (check_list(i)==id) then
                      in_check_list=.true.
                      exit
                    endif
                  end do
                  if (.not. in_check_list) then
                    n_check_list=n_check_list+1
                    check_list(n_check_list)=id
                  endif
                endif
              end do
              do nn_atom=1,kmc%system%sites(final_id)%neighbour_shells(nn_shell)%n_neighbour_sites
                nn_id=kmc%system%sites(final_id)%neighbour_shells(nn_shell)%neighbour_sites(nn_atom)
                if (kmc%system%sites(nn_id)%Z>0 .and. .not. kmc%system%sites(nn_id)%substrate) then
                  id=kmc%system%sites(nn_id)%id
                  in_check_list=.false.
                  do i=1,n_check_list
                    if (check_list(i)==id) then
                      in_check_list=.true.
                      exit
                    endif
                  end do
                  if (.not. in_check_list) then
                    n_check_list=n_check_list+1
                    check_list(n_check_list)=id
                  endif
                endif
              end do
              if (present(initial_id)) then
                do nn_atom=1,kmc%system%sites(initial_id)%neighbour_shells(nn_shell)%n_neighbour_sites
                  nn_id=kmc%system%sites(initial_id)%neighbour_shells(nn_shell)%neighbour_sites(nn_atom)
                  if (kmc%system%sites(nn_id)%Z>0 .and. .not. kmc%system%sites(nn_id)%substrate) then
                    id=kmc%system%sites(nn_id)%id
                    in_check_list=.false.
                    do i=1,n_check_list
                      if (check_list(i)==id) then
                        in_check_list=.true.
                        exit
                      endif
                    end do
                    if (.not. in_check_list) then
                      n_check_list=n_check_list+1
                      check_list(n_check_list)=id
                    endif
                  endif
                end do
              endif
            end do
            displacement2=max_displacement2(final,final_config_relaxed,check_list,n_check_list)
            call dealloc_atoms(final_config_relaxed)
            if (displacement2<kmc%system%relax_displacement_tol2) then
              final_relaxed=.true.
              kmc%minim_cache_results(kmc%minim_cache_index)=MINIM_RESULT_SUCCESS
              if (kmc%verbosity>=VERBOSITY_NORMAL) then
                write(6,"(a)") "# Final configuration minimized"
              endif
            else
              forward_barrier=100.0_rk
              reverse_barrier=100.0_rk
              kmc%minim_cache_results(kmc%minim_cache_index)=MINIM_RESULT_FAIL
              if (kmc%verbosity>=VERBOSITY_NORMAL) then
                write(6,"(a)") "# Final configuration minimization failed:"
                write(6,"(a)") "# too large displacement when relaxing the final configuration."
              endif
            endif
          endif
        endif
      endif
      call mpi_bcast(final_relaxed,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
      call mpi_bcast(forward_data_exists,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
      call mpi_bcast(reverse_data_exists,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
      call mpi_bcast(forward_ml_full,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
      call mpi_bcast(reverse_ml_full,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error)
      if ((.not. forward_data_exists .and. .not. forward_ml_full) .or. &
          (.not. reverse_data_exists .and. .not. reverse_ml_full)) then
        if (.not. final_relaxed) then
          return
        endif
        if (myrank==0) then
          if (kmc%verbosity>=VERBOSITY_NORMAL) then
            write(6,"(a)") "# Running NEB..."
          endif
          t1=mpi_wtime()
        endif
        if (.not. present(initial_id)) then
          call neb_barrier(kmc,lmp,site_id,final_id,init,&
                           neb_forward_barrier,neb_reverse_barrier,myrank)
        else
          call neb_barrier(kmc,lmp,site_id,final_id,init,&
                           neb_forward_barrier,neb_reverse_barrier,myrank,&
                           initial_id=initial_id)
        endif
        if (myrank==0) then
          t2=mpi_wtime()
          kmc%neb_barrier_time=kmc%neb_barrier_time+t2-t1
        endif
      endif
      if (myrank==0) then
        if (.not. forward_ml_full .or. .not. reverse_ml_full) then
          if (kmc%verbosity>=VERBOSITY_ABSURD) then
            write(6,"(a,f10.6,a)") "# E_{NEB} = ",neb_forward_barrier," eV"
            write(6,"(a,f10.6,a)") "# E_{NEB,reverse} = ",neb_reverse_barrier," eV"
          endif
          forward_barrier=neb_forward_barrier
          reverse_barrier=neb_reverse_barrier
          if (neb_reverse_barrier>0.0_rk) then
            event_ok=.true.
          endif
          if (.not. present(initial_id)) then ! hop event
            if (kmc%hop_predictor(coord_jumping_forward,coord_target_forward,jump_shell)%n_y<&
                kmc%hop_predictor(coord_jumping_forward,coord_target_forward,jump_shell)%n_y_max &
                .and. .not. forward_data_exists) then
              call add_data_point(kmc%hop_predictor(coord_jumping_forward,coord_target_forward,jump_shell),&
                                  desc_forward%x(1)%data(:),&
                                  neb_forward_barrier,kmc%neb_error,desc_forward%x(1)%covariance_cutoff,error=internal_error)
              if (internal_error/=0) then
                if (internal_error/=2) then
                  write(0,"(a)") "*** Error: couldn't add new data point to hop_predictor in get_barrier()"
                  if (present(error)) then
                    error=internal_error
                  endif
                endif
              else
                if (kmc%hop_predictor(coord_jumping_forward,coord_target_forward,jump_shell)%n_y>=&
                    kmc%hop_predictor(coord_jumping_forward,coord_target_forward,jump_shell)%n_y_min) then
                  t1=mpi_wtime()
                  call fit_predictor(kmc%hop_predictor(coord_jumping_forward,coord_target_forward,jump_shell),error=internal_error)
                  t2=mpi_wtime()
                  kmc%fit_time=kmc%fit_time+t2-t1
                  if (internal_error/=0) then
                    write(0,"(a)") "*** Error: couldn't fit hop_predictor in get_barrier()"
                    if (present(error)) then
                      error=internal_error
                    endif
                  endif
                endif
              endif
            endif
            reverse_data_exists=data_point_exists(kmc%hop_predictor(coord_jumping_reverse,coord_target_reverse,jump_shell),&
                                                  desc_reverse%x(1)%data(:))
            if (kmc%hop_predictor(coord_jumping_reverse,coord_target_reverse,jump_shell)%n_y<&
                kmc%hop_predictor(coord_jumping_reverse,coord_target_reverse,jump_shell)%n_y_max &
                .and. .not. reverse_data_exists) then
              call add_data_point(kmc%hop_predictor(coord_jumping_reverse,coord_target_reverse,jump_shell),&
                                  desc_reverse%x(1)%data(:),&
                                  neb_reverse_barrier,kmc%neb_error,desc_reverse%x(1)%covariance_cutoff,error=internal_error)
              if (internal_error/=0) then
                if (internal_error/=2) then
                  write(0,"(a)") "*** Error: couldn't add new data point to hop_predictor in get_barrier()"
                  if (present(error)) then
                    error=internal_error
                  endif
                endif
              else
                if (kmc%hop_predictor(coord_jumping_reverse,coord_target_reverse,jump_shell)%n_y>=&
                    kmc%hop_predictor(coord_jumping_reverse,coord_target_reverse,jump_shell)%n_y_min) then
                  t1=mpi_wtime()
                  call fit_predictor(kmc%hop_predictor(coord_jumping_reverse,coord_target_reverse,jump_shell),error=internal_error)
                  t2=mpi_wtime()
                  kmc%fit_time=kmc%fit_time+t2-t1
                  if (internal_error/=0) then
                    write(0,"(a)") "*** Error: couldn't fit hop_predictor in get_barrier()"
                    if (present(error)) then
                      error=internal_error
                    endif
                  endif
                endif
              endif
            endif
          else ! exchange event
            if (kmc%exchange_predictor(coord_middle,coord_jumping_forward,coord_target_forward)%n_y<&
                kmc%exchange_predictor(coord_middle,coord_jumping_forward,coord_target_forward)%n_y_max &
                .and. .not. forward_data_exists) then
              call add_data_point(kmc%exchange_predictor(coord_middle,coord_jumping_forward,&
                                                         coord_target_forward),desc_forward%x(1)%data(:),&
                                  neb_forward_barrier,kmc%neb_error,desc_forward%x(1)%covariance_cutoff,error=internal_error)
              if (internal_error/=0) then
                if (internal_error/=2) then
                  write(0,"(a)") "*** Error: couldn't add new data point to hop_predictor in get_barrier()"
                  if (present(error)) then
                    error=internal_error
                  endif
                endif
              else
                if (kmc%exchange_predictor(coord_middle,coord_jumping_forward,coord_target_forward)%n_y>=&
                    kmc%exchange_predictor(coord_middle,coord_jumping_forward,coord_target_forward)%n_y_min) then
                  t1=mpi_wtime()
                  call fit_predictor(kmc%exchange_predictor(coord_middle,coord_jumping_forward,&
                                                            coord_target_forward),error=internal_error)
                  t2=mpi_wtime()
                  kmc%fit_time=kmc%fit_time+t2-t1
                  if (internal_error/=0) then
                    write(0,"(a)") "*** Error: couldn't fit exchange_predictor in get_barrier()"
                    if (present(error)) then
                      error=internal_error
                    endif
                  endif
                endif
              endif
            endif
            reverse_data_exists=data_point_exists(kmc%exchange_predictor(coord_middle,coord_jumping_reverse,&
                                                                         coord_target_reverse),&
                                                  desc_reverse%x(1)%data(:))
            if (kmc%exchange_predictor(coord_middle,coord_jumping_reverse,coord_target_reverse)%n_y<&
                kmc%exchange_predictor(coord_middle,coord_jumping_reverse,coord_target_reverse)%n_y_max &
                .and. .not. reverse_data_exists) then
              call add_data_point(kmc%exchange_predictor(coord_middle,coord_jumping_reverse,&
                                                         coord_target_reverse),desc_reverse%x(1)%data(:),&
                                  neb_reverse_barrier,kmc%neb_error,desc_reverse%x(1)%covariance_cutoff,error=internal_error)
              if (internal_error/=0) then
                if (internal_error/=2) then
                  write(0,"(a)") "*** Error: couldn't add new data point to hop_predictor in get_barrier()"
                  if (present(error)) then
                    error=internal_error
                  endif
                endif
              else
                if (kmc%exchange_predictor(coord_middle,coord_jumping_reverse,coord_target_reverse)%n_y>=&
                    kmc%exchange_predictor(coord_middle,coord_jumping_reverse,coord_target_reverse)%n_y_min) then
                  t1=mpi_wtime()
                  call fit_predictor(kmc%exchange_predictor(coord_middle,coord_jumping_reverse,&
                                                            coord_target_reverse),error=internal_error)
                  t2=mpi_wtime()
                  kmc%fit_time=kmc%fit_time+t2-t1
                  if (internal_error/=0) then
                    write(0,"(a)") "*** Error: couldn't fit exchange_predictor in get_barrier()"
                    if (present(error)) then
                      error=internal_error
                    endif
                  endif
                endif
              endif
            endif
          endif
        endif
      endif
    endif

    if (present(error)) then
      call mpi_bcast(error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
    endif

  end subroutine get_barrier

  subroutine neb_barrier(kmc,lmp,site_id,final_id,init,&
                         barrier,reverse_barrier,myrank,&
                         initial_id)
    type(t_kmc),intent(inout) :: kmc
    type(lammps),intent(inout) :: lmp
    integer,intent(in) :: site_id,final_id
    type(t_atoms),intent(in) :: init
    double precision,intent(out) :: barrier,reverse_barrier
    integer,optional,intent(in) :: initial_id
    integer,intent(in) :: myrank

    integer :: i,first_id,last_id,nproc,jumping_atom_i,middle_atom_i
    double precision :: pe
    integer,allocatable :: types(:)
    real(kind=rk),allocatable :: x(:,:)
    real(kind=rk),dimension(3) :: diff_init_target,diff_init_middle,diff_middle_target
    double precision,allocatable :: pe_all(:)
    character(len=512) :: cmds
    integer :: mpi_error
    
    if (myrank==0) then
      middle_atom_i=0
      if (present(initial_id)) then
        diff_init_middle=system_diff_min_image(kmc%system,initial_id,site_id)
        diff_middle_target=system_diff_min_image(kmc%system,site_id,final_id)
        jumping_atom_i=kmc%system%sites(initial_id)%id
        middle_atom_i=kmc%system%sites(site_id)%id
      else
        diff_init_target=system_diff_min_image(kmc%system,site_id,final_id)
        jumping_atom_i=kmc%system%sites(site_id)%id
      endif

      nproc=kmc%neb_images
      first_id=kmc%system%substrate%first_id
      last_id=kmc%system%n_substrate+kmc%system%n_deposited
      allocate(types(first_id:last_id))
      allocate(x(3,first_id:last_id))
      types(kmc%system%substrate%first_id:kmc%system%substrate%last_id)=kmc%system%substrate%atom_type(:)
      x(:,kmc%system%substrate%first_id:kmc%system%substrate%last_id)=kmc%system%substrate%pos(:,:)
      types(kmc%system%substrate%last_id+1:last_id)=init%atom_type(:)
      x(:,kmc%system%substrate%last_id+1:last_id)=init%pos(:,:)
    endif
    call mpi_bcast(jumping_atom_i,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
    call mpi_bcast(middle_atom_i,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)

    call mpi_bcast(nproc,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
    call mpi_bcast(first_id,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
    call mpi_bcast(last_id,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
    if (myrank/=0) then
      allocate(types(first_id:last_id))
      allocate(x(3,first_id:last_id))
    endif
    call mpi_bcast(types,last_id-first_id+1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_error)
    call mpi_bcast(x,3*(last_id-first_id+1),MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,mpi_error)
    call mpi_bcast(diff_init_target,3,MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,mpi_error)
    call mpi_bcast(diff_init_middle,3,MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,mpi_error)
    call mpi_bcast(diff_middle_target,3,MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,mpi_error)

    if (middle_atom_i==0) then
      x(:,jumping_atom_i)=x(:,jumping_atom_i)+myrank*diff_init_target/(nproc-1)
    else
      x(:,jumping_atom_i)=x(:,jumping_atom_i)+myrank*diff_init_middle/(nproc-1)
      x(:,middle_atom_i)=x(:,middle_atom_i)+myrank*diff_middle_target/(nproc-1)
    endif
    call lmp%create_atoms(id=[(i,i=first_id,last_id)],type=types,x=reshape(real(x,kind=dp),[3*(last_id-first_id+1)]))
    call lmp%command("group mobile subtract all frozen")

    if (myrank==0) then
      if (kmc%verbosity>VERBOSITY_QUIET) then
        cmds="fix neb_fix all neb "//kmc%neb_spring//new_line("a")//&
             "min_style quickmin"//new_line("a")//&
             "min_modify dmax 0.01"//new_line("a")//&
             "neb "//kmc%relax_tol//" "//0.1d0*kmc%relax_tol//" "//kmc%neb_steps//" "//kmc%neb_steps//" "//&
              kmc%neb_steps/10//" none verbosity terse"
      else
        cmds="fix neb_fix all neb "//kmc%neb_spring//new_line("a")//&
             "min_style quickmin"//new_line("a")//&
             "min_modify dmax 0.01"//new_line("a")//&
             "neb "//kmc%relax_tol//" "//0.1d0*kmc%relax_tol//" "//kmc%neb_steps//" "//kmc%neb_steps//" "//&
              kmc%neb_steps/10//" none verbosity silent"
      endif
    endif
    call mpi_bcast(cmds,512,MPI_CHAR,0,MPI_COMM_WORLD,mpi_error)
    call lmp%commands_string(cmds)
    call lmp%command("unfix neb_fix")
    call lmp%command("min_style cg")
    call lmp%command("min_modify dmax 0.1")

    if (myrank==0) then
      allocate(pe_all(kmc%neb_images))
    endif
    pe=lmp%get_thermo("pe")
    call mpi_gather(pe,1,MPI_DOUBLE,pe_all,1,MPI_DOUBLE,0,MPI_COMM_WORLD,mpi_error)
    if (myrank==0) then
      barrier=maxval(pe_all)-pe_all(1)
      reverse_barrier=maxval(pe_all)-pe_all(kmc%neb_images)
    endif
    deallocate(x,types)
    if (myrank==0) then
      deallocate(pe_all)
      ! Increment neb_calculations counter
      kmc%neb_calculations=kmc%neb_calculations+1
    endif
    call lmp%command("delete_atoms group mobile")

  end subroutine neb_barrier

  subroutine pick_event(kmc,site_id,a,b,event_type,error)
    ! Pick an event
    type(t_kmc),intent(in) :: kmc
    integer,intent(out) :: site_id,a,b,event_type
    integer,optional,intent(out) :: error

    double precision :: u,summ,subsum,rate_lowering_f
    integer :: i,nn_atom,shell,jump_shell,initial_site,final_site

    if (present(error)) then
      error=0
    endif

    rate_lowering_f=RATE_LOWERING_FACTOR**kmc%rate_lowered_times
    a=0
    b=0
    event_type=EVENT_TYPE_NONE

    u=kmc%total_rate*genrand64_real2()
    summ=0.0d0
    
    do i=1,kmc%nZ ! Deposition event
      if (kmc%dep_rate(i)<=0.0_rk) then
        cycle
      endif
      summ=summ+kmc%dep_rate(i)
      if (summ>u) then
        a=i
        event_type=EVENT_TYPE_DEPOSITION
        return
      endif
    end do
    
    do i=1,kmc%system%n_deposited ! Jump event
      site_id=kmc%system%deposited_atoms(i)
      ! Check hop events
      if (kmc%system%sites(site_id)%hop_rate>0.0_rk) then
        subsum=summ
        summ=summ+rate_lowering_f*kmc%system%sites(site_id)%hop_rate
        if (summ>u) then
          do jump_shell=1,kmc%system%max_jump_shell
            shell=JUMP_SHELLS(jump_shell)
            do nn_atom=1,kmc%system%sites(site_id)%neighbour_shells(shell)%n_neighbour_sites
              if (kmc%system%sites(site_id)%hop_states(nn_atom,jump_shell)==EVENT_VALID) then
                subsum=subsum+rate_lowering_f*kmc%system%sites(site_id)%hop_rates(nn_atom,jump_shell)
                if (subsum>u) then
                  a=nn_atom
                  b=shell
                  event_type=EVENT_TYPE_HOP
                  return
                endif
              endif
            end do
          end do
          ! If we reach here, it's error time
          write (0,"(a)") "*** Error: failed to pick hop event for site ",site_id," in pick_event()."
          write (0,"(a,g20.10)")   "***          u = ",u
          write (0,"(a,g20.10,a)") "***  Summation = ",summ," Hz."
          write (0,"(a,g20.10,a)") "***     Subsum = ",subsum," Hz."
          if (present(error)) then
            error=1
          endif
          return
        endif
      endif
      ! Check exchange events
      if (kmc%system%allow_exchange) then
        if (kmc%system%sites(site_id)%exchange_rate>0.0_rk) then
          subsum=summ
          summ=summ+rate_lowering_f*kmc%system%sites(site_id)%exchange_rate
          if (summ>u) then
            do initial_site=1,kmc%system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
              do final_site=1,kmc%system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
                if (kmc%system%sites(site_id)%exchange_states(final_site,initial_site)==EVENT_VALID) then
                  subsum=subsum+rate_lowering_f*kmc%system%sites(site_id)%exchange_rates(final_site,initial_site)
                  if (subsum>u) then
                    a=final_site
                    b=initial_site
                    event_type=EVENT_TYPE_EXCHANGE
                    return
                  endif
                endif
              end do
            end do
            ! If we reach here, it's error time
            write (0,"(a)") "*** Error: failed to pick exchange event for site ",site_id," in pick_event()."
            write (0,"(a,g20.10)")   "***          u = ",u
            write (0,"(a,g20.10,a)") "***  Summation = ",summ," Hz."
            write (0,"(a,g20.10,a)") "***     Subsum = ",subsum," Hz."
            if (present(error)) then
              error=1
            endif
            return
          endif
        endif
      endif
    end do
    
    ! If we reach this point, it's error time
    write (0,"(a)") "*** Error: failed to pick event in pick_event()."
    write (0,"(a,g20.10)")   "***          u = ",u
    write (0,"(a,g20.10,a)") "*** Total rate = ",kmc%total_rate," Hz."
    write (0,"(a,g20.10,a)") "***  Summation = ",summ," Hz."
    if (present(error)) then
      error=1
    endif
  end subroutine pick_event

  subroutine execute_jump_event(kmc,site_id,a,b,event_type,error)
    type(t_kmc),intent(inout) :: kmc
    integer,intent(in) :: site_id,a,b,event_type
    integer,optional,intent(out) :: error

    integer :: i,nn_atom,shell,final_site,initial_site,id,reverse_nn_atom,saved_reverse_nn_atom
    integer :: final_id,initial_id,nn_id,nn_nn_id
    integer :: nn_shell,nn_nn_atom,jump_shell
    integer :: nn_initial,nn_mid,nn_final,fingerprint(kmc%nZ,0:NN_SHELLS),minim_result
    real(kind=rk) :: barrier,reverse_barrier
    double precision :: t1,t2
    integer :: internal_error

    if (present(error)) then
      error=0
    endif
    internal_error=0

    if (any(kmc%dep_rate>0.0d0)) then
      kmc%steps_since_deposition_or_rate_lowering=kmc%steps_since_deposition_or_rate_lowering+1
      if (kmc%steps_since_deposition_or_rate_lowering>=kmc%rate_lowering_threshold) then
        kmc%rate_lowered_times=kmc%rate_lowered_times+1
        kmc%steps_since_deposition_or_rate_lowering=0
      endif
    endif

    if (event_type==EVENT_TYPE_HOP) then
      nn_atom=a
      shell=b
      jump_shell=NN_SHELL_TO_JUMP_SHELL(shell)
      if (jump_shell<=0) then
        write(0,"(a,i0,a)") "*** Error: shell ",shell," doesn't map to an allowed jump shell in execute_jump_event()"
        if (present(error)) then
          error=1
        endif
        return
      endif

      final_id=kmc%system%sites(site_id)%neighbour_shells(shell)%neighbour_sites(nn_atom)

      if (kmc%verbosity>=VERBOSITY_LOUD) then
        write(6,"(2(a,i0))") "# Executing hop from ",site_id," to ",final_id
      endif

      if (kmc%system%sites(site_id)%Z<=0) then
        write(0,"(a,i0,a)") "*** Error: initial site ",site_id," does not hold an atom in execute_jump_event()"
        if (present(error)) then
          error=1
        endif
        return
      endif

      if (kmc%system%sites(final_id)%Z/=0) then
        write(0,"(a,i0,a)") "*** Error: final site ",final_id," is not a valid target site in execute_jump_event()"
        if (present(error)) then
          error=1
        endif
        return
      endif

      id=kmc%system%sites(site_id)%id
      kmc%system%deposited_atoms(id-kmc%system%n_substrate)=final_id

      kmc%system%sites(final_id)%atom_type=kmc%system%sites(site_id)%atom_type
      kmc%system%sites(site_id)%atom_type=0
      kmc%system%sites(final_id)%Z=kmc%system%sites(site_id)%Z
      kmc%system%sites(site_id)%Z=0
      kmc%system%sites(final_id)%id=id
      kmc%system%sites(site_id)%id=0
      call alloc_jumps(kmc%system,final_id,error=internal_error)
      if (internal_error/=0) then
        write(0,"(a,i0,a)") "*** Error allocating jumps for site ",final_id," in execute_jump_event()"
        if (present(error)) then
          error=internal_error
        endif
        return
      endif

      if (kmc%track_MSD) then
        if (id==kmc%system%n_substrate+1) then
          kmc%system%first_atom_unwrapped_coords(:)=kmc%system%first_atom_unwrapped_coords+&
                                                    system_diff_min_image(kmc%system,final_id,site_id)
        endif
      endif

      saved_reverse_nn_atom=kmc%system%sites(site_id)%neighbour_shells(shell)%reverse_indices(nn_atom)

      barrier=kmc%system%sites(site_id)%hop_barriers(1,nn_atom,jump_shell)
      reverse_barrier=kmc%system%sites(site_id)%hop_barriers(2,nn_atom,jump_shell)

      if (kmc%system%sites(final_id)%cartesian_coords(3)>kmc%system%height) then
        kmc%system%height=kmc%system%sites(final_id)%cartesian_coords(3)
      elseif (kmc%system%sites(final_id)%cartesian_coords(3)<kmc%system%sites(site_id)%cartesian_coords(3)) then
        if ((kmc%system%sites(final_id)%cartesian_coords(3)-kmc%system%sites(site_id)%cartesian_coords(3))**2>&
             kmc%system%layer_separation_111_squared-DIST2_TOL) then
          if (kmc%debug .and. kmc%verbosity>=VERBOSITY_LOUD) then
            write(6,"(a,f10.4,a)") "# Descension by hopping, barrier ",barrier," eV"
          endif
          kmc%system%height=0.0_rk
          do i=1,kmc%system%n_deposited
            id=kmc%system%deposited_atoms(i)
            if (kmc%system%sites(id)%cartesian_coords(3)>kmc%system%height) then
              kmc%system%height=kmc%system%sites(id)%cartesian_coords(3)
            endif
          end do
        endif
      endif

    elseif (event_type==EVENT_TYPE_EXCHANGE) then
      final_site=a
      initial_site=b

      initial_id=kmc%system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(initial_site)
      final_id=kmc%system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(final_site)

      if (kmc%verbosity>=VERBOSITY_LOUD) then
        write(6,"(3(a,i0))") "# Executing exchange from ",initial_id,&
                             " via ",site_id," to ",final_id
      endif
      
      if (kmc%system%sites(site_id)%Z<=0) then
        write(0,"(a)") "*** Error: middle site does not hold an atom in execute_jump_event()"
        if (present(error)) then
          error=1
        endif
        return
      endif

      if (kmc%system%sites(initial_id)%Z<=0) then
        write(0,"(a,i0,a)") "*** Error: initial site ",initial_id," does not hold an atom in execute_jump_event()"
        if (present(error)) then
          error=1
        endif
        return
      endif

      if (kmc%system%sites(final_id)%Z/=0) then
        write(0,"(a,i0,a)") "*** Error: final site ",final_id," is not a valid target site in execute_jump_event()"
        if (present(error)) then
          error=1
        endif
        return
      endif

      id=kmc%system%sites(site_id)%id

      kmc%system%deposited_atoms(id-kmc%system%n_substrate)=final_id

      kmc%system%sites(final_id)%atom_type=kmc%system%sites(site_id)%atom_type
      kmc%system%sites(final_id)%Z=kmc%system%sites(site_id)%Z
      kmc%system%sites(final_id)%id=id
      call alloc_jumps(kmc%system,final_id,error=internal_error)
      if (internal_error/=0) then
        write(0,"(a,i0,a)") "*** Error allocating jumps for site ",final_id," in execute_jump_event()"
        if (present(error)) then
          error=internal_error
        endif
        return
      endif

      if (kmc%track_MSD) then
        if (id==kmc%system%n_substrate+1) then
          kmc%system%first_atom_unwrapped_coords(:)=kmc%system%first_atom_unwrapped_coords+&
                                                    system_diff_min_image(kmc%system,final_id,site_id)
        endif
      endif

      id=kmc%system%sites(initial_id)%id

      kmc%system%deposited_atoms(id-kmc%system%n_substrate)=site_id

      kmc%system%sites(site_id)%atom_type=kmc%system%sites(initial_id)%atom_type
      kmc%system%sites(initial_id)%atom_type=0
      kmc%system%sites(site_id)%Z=kmc%system%sites(initial_id)%Z
      kmc%system%sites(initial_id)%Z=0
      kmc%system%sites(site_id)%id=id
      kmc%system%sites(initial_id)%id=0

      if (kmc%track_MSD) then
        if (id==kmc%system%n_substrate+1) then
          kmc%system%first_atom_unwrapped_coords(:)=kmc%system%first_atom_unwrapped_coords+&
                                                    system_diff_min_image(kmc%system,site_id,initial_id)
        endif
      endif

      barrier=kmc%system%sites(site_id)%exchange_barriers(1,final_site,initial_site)
      reverse_barrier=kmc%system%sites(site_id)%exchange_barriers(2,final_site,initial_site)

      if (kmc%system%sites(final_id)%cartesian_coords(3)>kmc%system%height) then
        kmc%system%height=kmc%system%sites(final_id)%cartesian_coords(3)
      elseif (kmc%system%sites(final_id)%cartesian_coords(3)<kmc%system%sites(initial_id)%cartesian_coords(3)) then
        if ((kmc%system%sites(final_id)%cartesian_coords(3)-kmc%system%sites(initial_id)%cartesian_coords(3))**2>&
            kmc%system%layer_separation_111_squared-DIST2_TOL) then
          if (kmc%debug .and. kmc%verbosity>=VERBOSITY_LOUD) then
            write(6,"(a,f10.4,a)") "# Descension by exchange, barrier ",barrier," eV"
          endif
          kmc%system%height=0.0_rk
          do i=1,kmc%system%n_deposited
            id=kmc%system%deposited_atoms(i)
            if (kmc%system%sites(id)%cartesian_coords(3)>kmc%system%height) then
              kmc%system%height=kmc%system%sites(id)%cartesian_coords(3)
            endif
          end do
        endif
      endif

    else
      write(0,"(a)") "*** Error: invalid event_type in execute_jump_event(). Bug?"
      if (present(error)) then
        error=1
      endif
      return
    endif

    ! Remove jumps to final position
    ! Update nn and n_near counts
    ! Mark jumps from and to LAE distance EVENT_LAE_CHANGED and associated sites to have jumps checked
    do shell=1,kmc%system%lae_cutoff_shell
      do nn_atom=1,kmc%system%sites(final_id)%neighbour_shells(shell)%n_neighbour_sites
        nn_id=kmc%system%sites(final_id)%neighbour_shells(shell)%neighbour_sites(nn_atom)
        if (shell<TRUE_1NN_SHELL) then
          kmc%system%sites(nn_id)%n_near=kmc%system%sites(nn_id)%n_near+1
        elseif (shell==TRUE_1NN_SHELL) then
          kmc%system%sites(nn_id)%nn=kmc%system%sites(nn_id)%nn+1
          if (.not. kmc%system%sites(nn_id)%substrate) then
            if (kmc%system%sites(nn_id)%Z<0) then
              if (kmc%system%sites(nn_id)%nn>3) then
                kmc%system%sites(nn_id)%Z=0
                kmc%system%sites(nn_id)%inactive_counter=0
              elseif (kmc%system%sites(nn_id)%nn==3) then
                if (has_support(kmc%system,nn_id)) then
                  kmc%system%sites(nn_id)%Z=0
                  kmc%system%sites(nn_id)%inactive_counter=0
                endif
              endif
            endif
          endif
        endif
        if (kmc%system%sites(nn_id)%Z>0 .and. .not. kmc%system%sites(nn_id)%substrate) then
          do jump_shell=1,kmc%system%max_jump_shell
            nn_shell=JUMP_SHELLS(jump_shell)
            kmc%system%sites(nn_id)%hop_states(1:kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%n_neighbour_sites,&
              jump_shell)=EVENT_LAE_CHANGED
          end do
          if (kmc%system%allow_exchange) then
            kmc%system%sites(nn_id)%exchange_states(&
              1:kmc%system%sites(nn_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites,&
              1:kmc%system%sites(nn_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites)=EVENT_LAE_CHANGED
          endif
          kmc%system%sites(nn_id)%jumps_checked=.false.

          jump_shell=NN_SHELL_TO_JUMP_SHELL(shell)
          if (jump_shell>0) then
            reverse_nn_atom=kmc%system%sites(final_id)%neighbour_shells(shell)%reverse_indices(nn_atom)
            kmc%system%sites(nn_id)%hop_states(reverse_nn_atom,jump_shell)=EVENT_FORBIDDEN
            if (kmc%system%allow_exchange) then
              if (shell==TRUE_1NN_SHELL) then
                kmc%system%sites(nn_id)%exchange_states(reverse_nn_atom,:)=EVENT_FORBIDDEN
              endif
            endif
          endif
        elseif (kmc%system%sites(nn_id)%Z==0) then
          do jump_shell=1,kmc%system%max_jump_shell
            nn_shell=JUMP_SHELLS(jump_shell)
            do nn_nn_atom=1,kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%n_neighbour_sites
              nn_nn_id=kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%neighbour_sites(nn_nn_atom)
              if (kmc%system%sites(nn_nn_id)%Z>0 .and. .not. kmc%system%sites(nn_nn_id)%substrate) then
                kmc%system%sites(nn_nn_id)%jumps_checked=.false.
                reverse_nn_atom=kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%reverse_indices(nn_nn_atom)
                kmc%system%sites(nn_nn_id)%hop_states(reverse_nn_atom,jump_shell)=EVENT_LAE_CHANGED
                if (kmc%system%allow_exchange) then
                  if (nn_shell==TRUE_1NN_SHELL) then
                    kmc%system%sites(nn_nn_id)%exchange_states(&
                      reverse_nn_atom,&
                      1:kmc%system%sites(nn_nn_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites&
                      )=EVENT_LAE_CHANGED
                  endif
                endif
              endif
            end do
          end do
        endif
      end do

      do nn_atom=1,kmc%system%sites(site_id)%neighbour_shells(shell)%n_neighbour_sites
        nn_id=kmc%system%sites(site_id)%neighbour_shells(shell)%neighbour_sites(nn_atom)
        if (event_type==EVENT_TYPE_HOP) then
          if (shell<TRUE_1NN_SHELL) then
            kmc%system%sites(nn_id)%n_near=kmc%system%sites(nn_id)%n_near-1
            if (kmc%system%sites(nn_id)%n_near<2) then
              kmc%system%sites(nn_id)%inactive_counter=0
            endif
          elseif (shell==TRUE_1NN_SHELL) then
            kmc%system%sites(nn_id)%nn=kmc%system%sites(nn_id)%nn-1
            if (kmc%system%sites(nn_id)%Z==0) then
              if (kmc%system%sites(nn_id)%nn<3) then
                kmc%system%sites(nn_id)%Z=-1
              elseif (kmc%system%sites(nn_id)%nn==3) then
                if (.not. has_support(kmc%system,nn_id)) then
                  kmc%system%sites(nn_id)%Z=-1
                endif
              endif
            endif
          endif
        endif
        if (kmc%system%sites(nn_id)%Z>0 .and. .not. kmc%system%sites(nn_id)%substrate) then
          do jump_shell=1,kmc%system%max_jump_shell
            nn_shell=JUMP_SHELLS(jump_shell)
            kmc%system%sites(nn_id)%hop_states(1:kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%n_neighbour_sites,&
              jump_shell)=EVENT_LAE_CHANGED
          end do
          if (kmc%system%allow_exchange) then
            kmc%system%sites(nn_id)%exchange_states(&
              1:kmc%system%sites(nn_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites,&
              1:kmc%system%sites(nn_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites)=EVENT_LAE_CHANGED
          endif
          kmc%system%sites(nn_id)%jumps_checked=.false.

          if (event_type==EVENT_TYPE_HOP) then
            if (shell==TRUE_1NN_SHELL) then
              reverse_nn_atom=kmc%system%sites(site_id)%neighbour_shells(shell)%reverse_indices(nn_atom)
              if (kmc%system%allow_exchange) then
                kmc%system%sites(nn_id)%exchange_states(:,reverse_nn_atom)=EVENT_FORBIDDEN
              endif
            endif
          endif
        elseif (kmc%system%sites(nn_id)%Z==0) then
          do jump_shell=1,kmc%system%max_jump_shell
            nn_shell=JUMP_SHELLS(jump_shell)
            do nn_nn_atom=1,kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%n_neighbour_sites
              nn_nn_id=kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%neighbour_sites(nn_nn_atom)
              if (kmc%system%sites(nn_nn_id)%Z>0 .and. .not. kmc%system%sites(nn_nn_id)%substrate) then
                kmc%system%sites(nn_nn_id)%jumps_checked=.false.
                reverse_nn_atom=kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%reverse_indices(nn_nn_atom)
                kmc%system%sites(nn_nn_id)%hop_states(reverse_nn_atom,jump_shell)=EVENT_LAE_CHANGED
                if (kmc%system%allow_exchange) then
                  if (nn_shell==TRUE_1NN_SHELL) then
                    kmc%system%sites(nn_nn_id)%exchange_states(&
                      reverse_nn_atom,&
                      1:kmc%system%sites(nn_nn_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites&
                      )=EVENT_LAE_CHANGED
                  endif
                endif
              endif
            end do
          end do
        else ! Z<0
          do jump_shell=1,kmc%system%max_jump_shell
            nn_shell=JUMP_SHELLS(jump_shell)
            do nn_nn_atom=1,kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%n_neighbour_sites
              nn_nn_id=kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%neighbour_sites(nn_nn_atom)
              if (kmc%system%sites(nn_nn_id)%Z>0 .and. .not. kmc%system%sites(nn_nn_id)%substrate) then
                kmc%system%sites(nn_nn_id)%jumps_checked=.false.
                reverse_nn_atom=kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%reverse_indices(nn_nn_atom)
                kmc%system%sites(nn_nn_id)%hop_states(reverse_nn_atom,jump_shell)=EVENT_FORBIDDEN
                if (kmc%system%allow_exchange) then
                  if (nn_shell==TRUE_1NN_SHELL) then
                    kmc%system%sites(nn_nn_id)%exchange_states(&
                      reverse_nn_atom,&
                      1:kmc%system%sites(nn_nn_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites&
                      )=EVENT_FORBIDDEN
                  endif
                endif
              endif
            end do
          end do
        endif
      end do

      ! If this is an exchange event, decrement the nn counts and mark events to be updated at the initial site
      if (event_type==EVENT_TYPE_EXCHANGE) then
        do nn_atom=1,kmc%system%sites(initial_id)%neighbour_shells(shell)%n_neighbour_sites
          nn_id=kmc%system%sites(initial_id)%neighbour_shells(shell)%neighbour_sites(nn_atom)
          if (shell<TRUE_1NN_SHELL) then
            kmc%system%sites(nn_id)%n_near=kmc%system%sites(nn_id)%n_near-1
            if (kmc%system%sites(nn_id)%n_near<2) then
              kmc%system%sites(nn_id)%inactive_counter=0
            endif
          elseif (shell==TRUE_1NN_SHELL) then
            kmc%system%sites(nn_id)%nn=kmc%system%sites(nn_id)%nn-1
            if (kmc%system%sites(nn_id)%Z==0) then
              if (kmc%system%sites(nn_id)%nn<3) then
                kmc%system%sites(nn_id)%Z=-1
              elseif (kmc%system%sites(nn_id)%nn==3) then
                if (.not. has_support(kmc%system,nn_id)) then
                  kmc%system%sites(nn_id)%Z=-1
                endif
              endif
            endif
          endif
          if (kmc%system%sites(nn_id)%Z>0 .and. .not. kmc%system%sites(nn_id)%substrate) then
            do jump_shell=1,kmc%system%max_jump_shell
              nn_shell=JUMP_SHELLS(jump_shell)
              kmc%system%sites(nn_id)%hop_states(1:kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%n_neighbour_sites,&
                jump_shell)=EVENT_LAE_CHANGED
            end do
            kmc%system%sites(nn_id)%exchange_states(&
              1:kmc%system%sites(nn_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites,&
              1:kmc%system%sites(nn_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites)=EVENT_LAE_CHANGED
            kmc%system%sites(nn_id)%jumps_checked=.false.

            if (shell==TRUE_1NN_SHELL) then
              reverse_nn_atom=kmc%system%sites(initial_id)%neighbour_shells(shell)%reverse_indices(nn_atom)
              kmc%system%sites(nn_id)%exchange_states(:,reverse_nn_atom)=EVENT_FORBIDDEN
            endif
          elseif (kmc%system%sites(nn_id)%Z==0) then
            do jump_shell=1,kmc%system%max_jump_shell
              nn_shell=JUMP_SHELLS(jump_shell)
              do nn_nn_atom=1,kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%n_neighbour_sites
                nn_nn_id=kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%neighbour_sites(nn_nn_atom)
                if (kmc%system%sites(nn_nn_id)%Z>0 .and. .not. kmc%system%sites(nn_nn_id)%substrate) then
                  kmc%system%sites(nn_nn_id)%jumps_checked=.false.
                  reverse_nn_atom=kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%reverse_indices(nn_nn_atom)
                  kmc%system%sites(nn_nn_id)%hop_states(reverse_nn_atom,jump_shell)=EVENT_LAE_CHANGED
                  if (nn_shell==TRUE_1NN_SHELL) then
                    kmc%system%sites(nn_nn_id)%exchange_states(&
                      reverse_nn_atom,&
                      1:kmc%system%sites(nn_nn_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites&
                      )=EVENT_LAE_CHANGED
                  endif
                endif
              end do
            end do
          else ! Z<0
            do jump_shell=1,kmc%system%max_jump_shell
              nn_shell=JUMP_SHELLS(jump_shell)
              do nn_nn_atom=1,kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%n_neighbour_sites
                nn_nn_id=kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%neighbour_sites(nn_nn_atom)
                if (kmc%system%sites(nn_nn_id)%Z>0 .and. .not. kmc%system%sites(nn_nn_id)%substrate) then
                  kmc%system%sites(nn_nn_id)%jumps_checked=.false.
                  reverse_nn_atom=kmc%system%sites(nn_id)%neighbour_shells(nn_shell)%reverse_indices(nn_nn_atom)
                  kmc%system%sites(nn_nn_id)%hop_states(reverse_nn_atom,jump_shell)=EVENT_FORBIDDEN
                  if (nn_shell==TRUE_1NN_SHELL) then
                    kmc%system%sites(nn_nn_id)%exchange_states(&
                      reverse_nn_atom,&
                      1:kmc%system%sites(nn_nn_id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites&
                      )=EVENT_FORBIDDEN
                  endif
                endif
              end do
            end do
          endif
        end do
      endif

    end do

    if (event_type==EVENT_TYPE_HOP) then
      nn_atom=a
      shell=b
      jump_shell=NN_SHELL_TO_JUMP_SHELL(shell)

      ! Reverse jump validity
      ! Check minimization cache for a hit of initial state minimization
      t1=mpi_wtime()
      call create_fingerprint(kmc,final_id,site_id,fingerprint,error=internal_error)
      t2=mpi_wtime()
      kmc%event_fingerprint_time=kmc%event_fingerprint_time+t2-t1
      if (internal_error/=0) then
        write(0,"(a)") "*** Error getting fingerprint in execute_jump_event()"
        return
      endif
      t1=mpi_wtime()
      minim_result=check_minim_cache(kmc,fingerprint,error=internal_error)
      t2=mpi_wtime()
      kmc%event_cache_time=kmc%event_cache_time+t2-t1
      if (internal_error/=0) then
        write(0,"(a)") "*** Error getting fingerprint in execute_jump_event()"
        return
      endif
      if (minim_result/=MINIM_RESULT_FAIL) then
        kmc%system%sites(final_id)%hop_states(saved_reverse_nn_atom,jump_shell)=EVENT_VALID
        kmc%system%sites(final_id)%hop_barriers(1,saved_reverse_nn_atom,jump_shell)=reverse_barrier
        kmc%system%sites(final_id)%hop_barriers(2,saved_reverse_nn_atom,jump_shell)=barrier
        kmc%system%sites(final_id)%hop_rates(saved_reverse_nn_atom,jump_shell)=&
          kmc%attempt_frequency*exp(-real(reverse_barrier,kind=dp)/kmc%kT)
      else
        if (kmc%debug .and. kmc%verbosity>=VERBOSITY_ABSURD) then
          write(6,"(a)") "# Reverse jump not to a stable adsorption site according to minim_cache"
        endif
        kmc%system%sites(final_id)%hop_states(saved_reverse_nn_atom,jump_shell)=EVENT_FORBIDDEN
        kmc%system%sites(final_id)%hop_barriers(:,saved_reverse_nn_atom,jump_shell)=100.0_rk
        kmc%system%sites(final_id)%hop_rates(saved_reverse_nn_atom,jump_shell)=0.0_rk
      endif
      kmc%system%sites(site_id)%hop_rate=0.0_rk
      kmc%system%sites(site_id)%exchange_rate=0.0_rk

    elseif (event_type==EVENT_TYPE_EXCHANGE) then
      final_site=a
      initial_site=b

      initial_id=kmc%system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(initial_site)
      final_id=kmc%system%sites(site_id)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(final_site)

      ! Reverse jump validity
      nn_initial=kmc%system%sites(final_id)%nn
      nn_final=kmc%system%sites(initial_id)%nn
      nn_mid=kmc%system%sites(site_id)%nn
      ! Check minimization cache for a hit of initial state minimization
      t1=mpi_wtime()
      call create_fingerprint(kmc,site_id,initial_id,fingerprint,initial_id=final_id,error=internal_error)
      t2=mpi_wtime()
      kmc%event_fingerprint_time=kmc%event_fingerprint_time+t2-t1
      if (internal_error/=0) then
        write(0,"(a)") "*** Error getting fingerprint in execute_jump_event()"
        return
      endif
      t1=mpi_wtime()
      minim_result=check_minim_cache(kmc,fingerprint,error=internal_error)
      t2=mpi_wtime()
      kmc%event_cache_time=kmc%event_cache_time+t2-t1
      if (internal_error/=0) then
        write(0,"(a)") "*** Error checking minimization cache in execute_jump_event()"
        return
      endif
      if (minim_result/=MINIM_RESULT_FAIL) then
        kmc%system%sites(site_id)%exchange_states(initial_site,final_site)=EVENT_VALID
        kmc%system%sites(site_id)%exchange_barriers(1,initial_site,final_site)=reverse_barrier
        kmc%system%sites(site_id)%exchange_barriers(2,initial_site,final_site)=barrier
        kmc%system%sites(site_id)%exchange_rates(initial_site,final_site)=&
          kmc%attempt_frequency*exp(-real(reverse_barrier,kind=dp)/kmc%kT)
        kmc%system%sites(site_id)%exchange_rate=&
          kmc%system%sites(site_id)%exchange_rates(initial_site,final_site)
      else
        if (kmc%debug .and. kmc%verbosity>=VERBOSITY_ABSURD) then
          write(6,"(a)") "# Reverse jump not to a stable adsorption site according to minim_cache"
        endif
        kmc%system%sites(site_id)%exchange_states(initial_site,final_site)=EVENT_FORBIDDEN
        kmc%system%sites(site_id)%exchange_barriers(:,initial_site,final_site)=100.0_rk
        kmc%system%sites(site_id)%exchange_rates(initial_site,final_site)=0.0_rk
      endif
      kmc%system%sites(initial_id)%hop_rate=0.0_rk
      kmc%system%sites(initial_id)%exchange_rate=0.0_rk

    endif

    ! Find new adsorption sites around the final site
    call find_new_sites(kmc%system,final_id,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error finding new sites around site ",&
                          final_id," in execute_jump_event()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

    ! Update lockedness
    kmc%system%sites(1:kmc%system%n_sites)%lock_checked=.false.

    call check_lockedness(kmc%system,final_id,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error checking lockedness of site ",&
                          final_id," in execute_jump_event()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    kmc%system%sites(final_id)%lock_checked=.true.

    call do_lockedness_checks(kmc%system,final_id,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error doing lockedness checks to site ",&
                          final_id," in execute_jump_event()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

    call do_lockedness_checks(kmc%system,site_id,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error doing lockedness checks to site ",&
                          site_id," in execute_jump_event()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

    if (event_type==EVENT_TYPE_EXCHANGE) then
      call do_lockedness_checks(kmc%system,initial_id,error=internal_error)
      if (internal_error/=0) then
        write(0,"(a,i0,a)") "*** Error doing lockedness checks to site ",&
                            initial_id," in execute_jump_event()"
        if (present(error)) then
          error=internal_error
        endif
        return
      endif

    endif

  end subroutine execute_jump_event

  function print_event(event_type) result(event_str)
    integer,intent(in) :: event_type
    character(len=80) :: event_str

    select case(event_type)
      case(EVENT_VALID)
        event_str="valid"
      case(EVENT_EMPTY)
        event_str="empty"
      case(EVENT_LAE_CHANGED)
        event_str="LAE changed"
      case(EVENT_FORBIDDEN)
        event_str="forbidden"
      case default
        event_str="Unknown! Bug?"
    end select
  end function print_event

  subroutine get_jump_lae(kmc,site_id,final_id,jump_type,&
                         lae,initial_id,jumping_ptr,error)
    type(t_kmc),intent(inout) :: kmc
    integer,intent(in) :: site_id,final_id
    integer,intent(in) :: jump_type
    type(Atoms),intent(out) :: lae
    integer,optional,intent(in) :: initial_id
    logical,pointer,optional,intent(out) :: jumping_ptr(:)
    integer,optional,intent(out) :: error

    double precision :: t1,t2
    integer :: internal_error
    integer :: nn_atom,shell
    integer :: nn_id,Z
    logical,allocatable :: jumping(:)
    logical,allocatable :: lae_mask(:)

    if (present(error)) then
      error=0
    endif
    internal_error=0

    allocate(lae_mask(kmc%system%n_sites))
    lae_mask=.false.

    t1=mpi_wtime()
    call initialise(lae,0,real(kmc%system%lattice,kind=dp))
    t2=mpi_wtime()
    kmc%lae_init_time=kmc%lae_init_time+t2-t1

    t1=mpi_wtime()
    if (jump_type==EVENT_TYPE_EXCHANGE) then
      call add_atoms(lae,real(kmc%system%sites(initial_id)%cartesian_coords(:),kind=dp),&
                     kmc%system%sites(initial_id)%Z)
      lae_mask(initial_id)=.true.
    endif

    call add_atoms(lae,real(kmc%system%sites(site_id)%cartesian_coords(:),kind=dp),kmc%system%sites(site_id)%Z)
    lae_mask(site_id)=.true.

    ! Neighbourhood of site_id
    t1=mpi_wtime()
    do shell=TRUE_1NN_SHELL,kmc%system%lae_cutoff_shell
      do nn_atom=1,kmc%system%sites(site_id)%neighbour_shells(shell)%n_neighbour_sites
        nn_id=kmc%system%sites(site_id)%neighbour_shells(shell)%neighbour_sites(nn_atom)
        if (kmc%system%sites(nn_id)%Z>0) then
          if (.not. lae_mask(nn_id)) then
            ! If modelling weakly-interacting substrate, add substrate atoms with dummy atomic number
            if (kmc%system%sites(nn_id)%substrate .and. kmc%system%substrate_type==SUBSTRATE_TYPE_WEAK) then
              Z=SUBSTRATE_Z
            else
              Z=kmc%system%sites(nn_id)%Z
            endif
            call add_atoms(lae,real(kmc%system%sites(nn_id)%cartesian_coords(:),kind=dp),Z)
            lae_mask(nn_id)=.true.
          endif
        endif
      end do
    end do

    ! Neighbourhood of final_id
    do shell=TRUE_1NN_SHELL,kmc%system%lae_cutoff_shell
      do nn_atom=1,kmc%system%sites(final_id)%neighbour_shells(shell)%n_neighbour_sites
        nn_id=kmc%system%sites(final_id)%neighbour_shells(shell)%neighbour_sites(nn_atom)
        if (kmc%system%sites(nn_id)%Z>0) then
          if (.not. lae_mask(nn_id)) then
            ! If modelling weakly-interacting substrate, add substrate atoms with dummy atomic number
            if (kmc%system%sites(nn_id)%substrate .and. kmc%system%substrate_type==SUBSTRATE_TYPE_WEAK) then
              Z=SUBSTRATE_Z
            else
              Z=kmc%system%sites(nn_id)%Z
            endif
            call add_atoms(lae,real(kmc%system%sites(nn_id)%cartesian_coords(:),kind=dp),Z)
            lae_mask(nn_id)=.true.
          endif
        endif
      end do
    end do

    if (jump_type==EVENT_TYPE_EXCHANGE) then
      ! Neighbourhood of initial_id
      do shell=TRUE_1NN_SHELL,kmc%system%lae_cutoff_shell
        do nn_atom=1,kmc%system%sites(initial_id)%neighbour_shells(shell)%n_neighbour_sites
          nn_id=kmc%system%sites(initial_id)%neighbour_shells(shell)%neighbour_sites(nn_atom)
          if (kmc%system%sites(nn_id)%Z>0) then
            if (.not. lae_mask(nn_id)) then
              ! If modelling weakly-interacting substrate, add substrate atoms with dummy atomic number
              if (kmc%system%sites(nn_id)%substrate .and. kmc%system%substrate_type==SUBSTRATE_TYPE_WEAK) then
                Z=SUBSTRATE_Z
              else
                Z=kmc%system%sites(nn_id)%Z
              endif
              call add_atoms(lae,real(kmc%system%sites(nn_id)%cartesian_coords(:),kind=dp),Z)
              lae_mask(nn_id)=.true.
            endif
          endif
        end do
      end do
    endif
    t2=mpi_wtime()
    kmc%lae_add_time=kmc%lae_add_time+t2-t1

    allocate(jumping(lae%N))
    jumping=.false.
    jumping(1)=.true.
    if (present(jumping_ptr)) then
      call add_property(lae,"jumping",jumping,ptr=jumping_ptr)
    else
      call add_property(lae,"jumping",jumping)
    endif
    deallocate(jumping)
    deallocate(lae_mask)
    call set_cutoff(lae,cutoff(kmc%desc))
    t1=mpi_wtime()
    call calc_connect(lae)
    t2=mpi_wtime()
    kmc%lae_connect_time=kmc%lae_connect_time+t2-t1

  end subroutine get_jump_lae

  subroutine only_substrate_lae(kmc,site_id,final_id,jump_type,error)
    type(t_kmc),intent(inout) :: kmc
    integer,intent(in) :: site_id,final_id
    integer,intent(out) :: jump_type
    integer,optional,intent(out) :: error

    integer :: internal_error
    integer :: nn_atom,nn_id,shell

    ! Default "only-substrate jump". Changed to hcp-fcc if off-lattice is detected
    jump_type=SUBSTRATE_FCC_HCP

    do shell=TRUE_1NN_SHELL,kmc%system%lae_cutoff_shell
      do nn_atom=1,kmc%system%sites(site_id)%neighbour_shells(shell)%n_neighbour_sites
        nn_id=kmc%system%sites(site_id)%neighbour_shells(shell)%neighbour_sites(nn_atom)
        if (kmc%system%sites(nn_id)%Z>0) then
          ! Not-substrate detected
          if (.not. kmc%system%sites(nn_id)%substrate) then
            jump_type=SUBSTRATE_NOT
            return
          ! Off-lattice site:
          elseif (shell>TRUE_2NN_SHELL .and. shell<TRUE_3NN_SHELL) then
            jump_type=SUBSTRATE_HCP_FCC
          endif
        endif
      end do
      do nn_atom=1,kmc%system%sites(final_id)%neighbour_shells(shell)%n_neighbour_sites
        nn_id=kmc%system%sites(final_id)%neighbour_shells(shell)%neighbour_sites(nn_atom)
        if (kmc%system%sites(nn_id)%Z>0) then
          ! Not-substrate detected
          if (.not. kmc%system%sites(nn_id)%substrate) then
            jump_type=SUBSTRATE_NOT
            return
          endif
        endif
      end do
    end do

  end subroutine only_substrate_lae

  real(kind=rk) function max_displacement2(at1,at2,check_list,n_check_list)
    type(t_atoms),intent(in) :: at1,at2
    integer,intent(in) :: check_list(:),n_check_list

    real(kind=rk) :: d2
    integer :: i,id

    max_displacement2=0.0_rk

    do i=1,n_check_list
      id=check_list(i)
      d2=atoms_dist2_min_image(at1%lattice,at2%pos(:,id),at1%pos(:,id),at1%periodic)
      if (d2>max_displacement2) then
        max_displacement2=d2
      endif
    end do
  end function max_displacement2

  subroutine print_predictors(kmc,error)
    type(t_kmc),intent(inout) :: kmc
    integer,optional,intent(out) :: error

    integer :: shell,coord_target,coord_jumping,coord_middle
    integer :: internal_error

    if (present(error)) then
      error=0
    endif
    internal_error=0

    do shell=1,kmc%system%max_jump_shell
      do coord_target=3,N_NN_MAX
        do coord_jumping=3,LOCK_THRESHOLD
          if (coord_jumping<N_NN_MAX) then
            if (kmc%hop_predictor(coord_jumping,coord_target,shell)%initialised) then
              if (kmc%hop_predictor(coord_jumping,coord_target,shell)%changed) then
                call print_predictor_xml(kmc%hop_predictor(coord_jumping,coord_target,shell),&
                      "predictors/hop_predictor_"//coord_jumping//"_"//coord_target//"_"//shell//".xml",&
                      error=internal_error)
                if (internal_error/=0) then
                  write(0,"(a)") "*** Error writing predictors/hop_predictor_"//coord_jumping&
                                  //"_"//coord_target//"_"//shell//".xml"
                  if (present(error)) then
                    error=internal_error
                  endif
                  return
                endif
              endif
            endif
          endif
          if (kmc%system%allow_exchange) then
            if (shell==1) then
              do coord_middle=3,N_NN_MAX-1
                if (kmc%exchange_predictor(coord_middle,coord_jumping,coord_target)%initialised) then
                  if (kmc%exchange_predictor(coord_middle,coord_jumping,coord_target)%changed) then
                    call print_predictor_xml(kmc%exchange_predictor(coord_middle,coord_jumping,&
                                                                    coord_target),&
                                             "predictors/exchange_predictor_"//coord_middle//"_"//coord_jumping&
                                             //"_"//coord_target//".xml",error=internal_error)
                    if (internal_error/=0) then
                      write(0,"(a)") "*** Error writing predictors/exchange_predictor_"//coord_middle//"_"//coord_jumping&
                                     //"_"//coord_target//".xml"
                      if (present(error)) then
                        error=internal_error
                      endif
                      return
                    endif
                  endif
                endif
              end do
            endif
          endif
        end do
      end do
    end do
  end subroutine print_predictors

  integer function get_island_count(system,error)
    type(t_system),intent(in) :: system
    integer,optional,intent(out) :: error

    integer :: i,id,internal_error
    integer,allocatable :: island_ids(:)

    if (present(error)) then
      error=0
    endif
    internal_error=0
    
    allocate(island_ids(system%n_sites))
    get_island_count=0
    island_ids=0

    do i=1,system%n_deposited
      id=system%deposited_atoms(i)
      if (id<=0 .or. id>system%n_sites) then
        write(0,"(a,i0,a,i0,a)") "*** Error: site ",i," in deposited_atoms array has invalid id ",&
                                      id," in get_island_count()"
        if (present(error)) then
          error=1
        endif
        get_island_count=0
        return
      endif

      if (island_ids(id)==0) then
        ! If the atom has non-substrate neighbours, i.e., it is not a monomer:
        if (system%sites(id)%nn>system%sites(id)%nn_substrate) then
          get_island_count=get_island_count+1
          island_ids(id)=get_island_count
          call color_islands(system,id,get_island_count,island_ids,error=internal_error)
          if (internal_error/=0) then
            write(0,"(a)") "*** Error from color_islands() in get_island_count()"
            if (present(error)) then
              error=internal_error
            endif
            get_island_count=0
            return
          endif
        endif
      endif
    end do
    deallocate(island_ids)

  end function get_island_count

  recursive subroutine color_islands(system,id,island_color,island_ids,error)
    type(t_system),intent(in) :: system
    integer,intent(in) :: id,island_color
    integer,intent(inout) :: island_ids(system%n_sites)
    integer,optional,intent(out) :: error

    integer :: i,nn_id,internal_error

    if (present(error)) then
      error=0
    endif
    internal_error=0

    do i=1,system%sites(id)%neighbour_shells(TRUE_1NN_SHELL)%n_neighbour_sites
      nn_id=system%sites(id)%neighbour_shells(TRUE_1NN_SHELL)%neighbour_sites(i)
      if (nn_id<=0 .or. nn_id>system%n_sites) then
        write(0,"(a,i0,a,i0,a,i0,a)") "*** Error: nn_atom ",i," of site ",id," on shell ",TRUE_1NN_SHELL,&
                                      " has invalid id ",nn_id," in color_islands()"
        if (present(error)) then
          error=1
        endif
        return
      endif
      if (system%sites(id)%Z<=0) then
        cycle
      endif
      if (system%sites(id)%substrate) then
        cycle
      endif
      if (island_ids(nn_id)==0) then
        island_ids(nn_id)=island_color
        call color_islands(system,nn_id,island_color,island_ids,error=internal_error)
        if (internal_error/=0) then
          if (present(error)) then
            error=internal_error
          endif
          return
        endif
      endif
    end do

  end subroutine color_islands

  real(kind=rk) function get_square_displacement(system)
    type(t_system),intent(in) :: system

    if (system%n_deposited<=0) then
      get_square_displacement=0.0_rk
      return
    endif

    ! Return square displacement in square centimeters
    get_square_displacement=sum(1d-16*(system%first_atom_unwrapped_coords(:)-system%first_atom_initial_coords(:))**2)
  end function get_square_displacement

  subroutine dump_data_final(kmc,error)
    type(t_kmc),intent(inout) :: kmc
    integer,optional,intent(out) :: error

    integer :: internal_error

    if (present(error)) then
      error=0
    endif
    internal_error=0

    call dump_data(kmc,"final",error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error dumping data in dump_data_final()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

    call print_predictors(kmc,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error printing predictors in dump_data_final()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

  end subroutine dump_data_final

  subroutine dump_data_intermediate(kmc,error)
    type(t_kmc),intent(inout) :: kmc
    integer,optional,intent(out) :: error

    integer :: internal_error

    if (present(error)) then
      error=0
    endif
    internal_error=0

    call dump_data(kmc,"intermediate",error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error dumping data in dump_data_intermediate()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

    call print_predictors(kmc,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error printing predictors in dump_data_intermediate()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

  end subroutine dump_data_intermediate

  subroutine dump_data_error(kmc,error)
    type(t_kmc),intent(inout) :: kmc
    integer,optional,intent(out) :: error

    integer :: internal_error

    if (present(error)) then
      error=0
    endif
    internal_error=0

    call dump_data(kmc,"error",error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error dumping data in dump_data_error()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

    ! Don't print predictors

  end subroutine dump_data_error

  subroutine dump_data(kmc,filename_stub,error)
    type(t_kmc),intent(inout) :: kmc
    character(len=*),intent(in) :: filename_stub
    integer,optional,intent(out) :: error

    character(len=80) :: filename
    integer :: internal_error

    if (present(error)) then
      error=0
    endif
    internal_error=0

    filename=trim(filename_stub)//".restart"
    call write_my_restart(kmc,filename,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error writing KMC restart info to "//trim(filename)//" in dump_data()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

    filename=trim(filename_stub)//".mt.state"
    call write_mt_state(filename,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error writing Mersenne Twister state to "//trim(filename)//" in dump_data()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

    filename=trim(filename_stub)//".system"
    call print_system(kmc%system,filename,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error writing restart system to "//trim(filename)//" in dump_data()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

  end subroutine dump_data

  subroutine write_my_restart(kmc,filename,error)
    type(t_kmc),intent(in) :: kmc
    character(len=*),intent(in) :: filename
    integer,optional,intent(out) :: error

    integer :: j,k,fileunit,ios

    if (present(error)) then
      error=0
    endif

    open(newunit=fileunit,file=filename,action="write",iostat=ios,RECL=20000)
    if (ios/=0) then
      write(0,"(a)") "*** Error opening file "//trim(filename)//" for writing in write_my_restart()"
      if (present(error)) then
        error=ios
      endif
      return
    endif

    ! Timers not written on purpose!
    ! Verbosity and debug flag and deposition rate also now written to
    ! be able to modify more easily between runs. They are read from in.mlkmc
    ! Same goes for deposition rates and the temperature

    write(fileunit,*) kmc%relax_iter
    write(fileunit,*) kmc%neb_images
    write(fileunit,*) kmc%neb_steps
    write(fileunit,*) kmc%relax_tol
    write(fileunit,*) kmc%ml_tol
    write(fileunit,*) kmc%neb_spring
    write(fileunit,*) kmc%neb_error

    write(fileunit,*) kmc%nZ
    write(fileunit,*) kmc%Z

    write(fileunit,*) kmc%attempt_frequency

    write(fileunit,*) kmc%substrate_fcc_hcp_barrier
    write(fileunit,*) kmc%substrate_hcp_fcc_barrier
    write(fileunit,*) kmc%last_barrier

    write(fileunit,*) kmc%steps_since_deposition_or_rate_lowering
    write(fileunit,*) kmc%rate_lowered_times
    write(fileunit,*) kmc%rate_lowering_threshold

    write(fileunit,*) kmc%n_data_min
    write(fileunit,*) kmc%n_data_max
    write(fileunit,*) kmc%n_sparse_max
    write(fileunit,*) kmc%covariance_type
    write(fileunit,*) kmc%delta
    write(fileunit,*) kmc%zeta
    write(fileunit,*) kmc%regularisation
    write(fileunit,*) string(kmc%desc_str)

    write(fileunit,*) kmc%step
    write(fileunit,*) kmc%time
    write(fileunit,*) kmc%total_rate

    write(fileunit,*) kmc%initialised

    write(fileunit,*) kmc%minim_cache_index
    do k=1,MINIM_CACHE_SIZE
      write(fileunit,*) kmc%minim_cache_results(k)
      do j=0,NN_SHELLS
        write(fileunit,*) kmc%minim_cache(:,j,k)
      end do
    end do

    write(fileunit,*) kmc%neb_calculations
    write(fileunit,*) kmc%minimizations

    close(fileunit)

  end subroutine write_my_restart

  subroutine read_my_restart(kmc,filename,error)
    type(t_kmc),intent(inout) :: kmc
    character(len=*),intent(in) :: filename
    integer,optional,intent(out) :: error

    integer :: j,k,fileunit,ios,internal_error
    character(len=512) :: desc_str

    if (present(error)) then
      error=0
    endif
    internal_error=0

    open(newunit=fileunit,file=filename,action="read",iostat=ios,RECL=20000)
    if (ios/=0) then
      write(0,"(a)") "*** Error opening file "//trim(filename)//" for writing in read_my_restart()"
      if (present(error)) then
        error=ios
      endif
      return
    endif

    read(fileunit,*) kmc%relax_iter
    read(fileunit,*) kmc%neb_images
    read(fileunit,*) kmc%neb_steps
    read(fileunit,*) kmc%relax_tol
    read(fileunit,*) kmc%ml_tol
    read(fileunit,*) kmc%neb_spring
    read(fileunit,*) kmc%neb_error

    read(fileunit,*) kmc%nZ
    allocate(kmc%Z(kmc%nZ))
    allocate(kmc%dep_rate(kmc%nZ))
    read(fileunit,*) kmc%Z

    read(fileunit,*) kmc%attempt_frequency

    read(fileunit,*) kmc%substrate_fcc_hcp_barrier
    read(fileunit,*) kmc%substrate_hcp_fcc_barrier
    read(fileunit,*) kmc%last_barrier

    read(fileunit,*) kmc%steps_since_deposition_or_rate_lowering
    read(fileunit,*) kmc%rate_lowered_times
    read(fileunit,*) kmc%rate_lowering_threshold

    read(fileunit,*) kmc%n_data_min
    read(fileunit,*) kmc%n_data_max
    read(fileunit,*) kmc%n_sparse_max
    read(fileunit,*) kmc%covariance_type
    read(fileunit,*) kmc%delta
    read(fileunit,*) kmc%zeta
    read(fileunit,*) kmc%regularisation
    read(fileunit,"(a)") desc_str
    kmc%desc_str=desc_str
    call initialise(kmc%desc,desc_str,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error initialising descriptor in read_my_restart()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

    read(fileunit,*) kmc%step
    read(fileunit,*) kmc%time
    read(fileunit,*) kmc%total_rate

    read(fileunit,*) kmc%initialised

    allocate(kmc%minim_cache(kmc%nZ,0:NN_SHELLS,MINIM_CACHE_SIZE))
    kmc%minim_cache=0
    read(fileunit,*) kmc%minim_cache_index
    do k=1,MINIM_CACHE_SIZE
      read(fileunit,*) kmc%minim_cache_results(k)
      do j=0,NN_SHELLS
        read(fileunit,*) kmc%minim_cache(:,j,k)
      end do
    end do

    read(fileunit,*) kmc%neb_calculations
    read(fileunit,*) kmc%minimizations
    
    close(fileunit)

  end subroutine read_my_restart

  subroutine create_fingerprint(kmc,site_id,final_id,fingerprint,initial_id,error)
    type(t_kmc),intent(in) :: kmc
    integer,intent(in) :: site_id,final_id
    integer,intent(out) :: fingerprint(kmc%nZ,0:NN_SHELLS)
    integer,optional,intent(in) :: initial_id
    integer,optional,intent(out) :: error
    
    integer :: nn_atom,nn_id,shell,internal_error

    if (present(error)) then
      error=0
    endif
    internal_error=0

    fingerprint=0

    ! The zeroth identifier is the species of the jumping atom, now at final_id
    fingerprint(kmc%system%sites(site_id)%atom_type,0)=1

    do shell=1,NN_SHELLS
      do nn_atom=1,kmc%system%sites(final_id)%neighbour_shells(shell)%n_neighbour_sites
        nn_id=kmc%system%sites(final_id)%neighbour_shells(shell)%neighbour_sites(nn_atom)
        if (kmc%system%sites(nn_id)%Z<=0) then
          cycle
        endif
        if (nn_id==site_id) then
          ! If this is an exchange event, the atom at initial_id is actually at site_id
          if (present(initial_id)) then
            fingerprint(kmc%system%sites(initial_id)%atom_type,shell)=fingerprint(kmc%system%sites(initial_id)%atom_type,shell)+1
          endif
          ! The atom at site_id has jumped to final_id
          cycle
        endif
        if (present(initial_id)) then
          ! The atom at initial_id has jumped to site_id
          if (nn_id==initial_id) then
            cycle
          endif
        endif
        fingerprint(kmc%system%sites(nn_id)%atom_type,shell)=fingerprint(kmc%system%sites(nn_id)%atom_type,shell)+1
      end do
    end do

  end subroutine create_fingerprint

  integer function check_minim_cache(kmc,fingerprint,error)
    type(t_kmc),intent(in) :: kmc
    integer,intent(in) :: fingerprint(kmc%nZ,0:NN_SHELLS)
    integer,optional,intent(out) :: error

    integer :: i,j,k
    logical :: match

    check_minim_cache=MINIM_RESULT_UNINIT

    do k=1,MINIM_CACHE_SIZE
      match=.true.
      do j=0,NN_SHELLS
        do i=1,kmc%nZ
          if (fingerprint(i,j)/=kmc%minim_cache(i,j,k)) then
            match=.false.
            exit
          endif
        end do
        if (.not. match) then
          exit
        endif
      end do
      if (match) then
        check_minim_cache=kmc%minim_cache_results(k)
        return
      endif
    end do

  end function check_minim_cache

  subroutine write_output_line(kmc,error)
    type(t_kmc),intent(in) :: kmc
    integer,optional,intent(out) :: error

    integer :: island_count,internal_error

    if (present(error)) then
      error=0
    endif

    island_count=get_island_count(kmc%system,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Warning: error getting island_count in write_output_line()"
      if (present(error)) then
        error=internal_error
      endif
    endif

    write(6,"(i12,g20.6,g20.6,f20.6,i12,i12,f12.3,g14.6)",advance="no") kmc%step,kmc%time,kmc%total_rate,kmc%last_barrier,&
                                                                        kmc%system%n_deposited,&
                                                                        island_count,kmc%system%height
    if (kmc%track_MSD) then
      write(6,"(g14.6)",advance="no") get_square_displacement(kmc%system)
    endif
    if (kmc%debug) then
      write(6,"(i12,i14)",advance="no") kmc%neb_calculations,kmc%minimizations
    endif
    write(6,*)

  end subroutine write_output_line

  subroutine write_output_header(kmc)
    type(t_kmc),intent(in) :: kmc

    write(6,"(a12,a20,a20,a20,a12,a12,a13)",advance="no") "step","time(s)","total_rate(Hz)","last_barrier(eV)",&
                                                          "n_deposited","n_islands","height(Å)"
    if (kmc%track_MSD) then
      write(6,"(a12)",advance="no") "SD(cm^2)"
    endif
    if (kmc%debug) then
      write(6,"(a12,a14)",advance="no") "neb_runs","minimizations"
    endif
    write(6,*)
  
  end subroutine write_output_header

end module kmc_module
