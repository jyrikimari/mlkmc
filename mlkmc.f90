! MLKMC - main program
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

program mlkmc
  use mpi
  use liblammps
  use libAtoms_module
  use system_module, only : OUTPUT
  use potential_module
  use constants_module
  use ml_module
  use kmc_module
  implicit none

  type(t_kmc) :: kmc
  type(lammps) :: lmp
  integer :: relax_iter,neb_images,neb_steps,max_jump_nth,lae_cutoff_nth
  double precision :: relax_tol,ml_tol,neb_spring
  real(kind=rk) :: lattice_constant
  integer :: nZ,n_soap_Z,n_soap_species
  integer,allocatable :: Z(:),Z_tmp(:),soap_Z(:),soap_species_Z(:)
  double precision :: temp,attempt_frequency
  type(t_system) :: system

  integer :: xmax,ymax,zmax,substrate_layers,substrate_species

  double precision,allocatable :: dep_rate(:),dep_rate_tmp(:)
  integer,allocatable :: initial_atoms(:),initial_atoms_layers(:),initial_atoms_termination(:)
  double precision,allocatable :: initial_atoms_coords(:,:)
  integer :: substrate_type

  type(extendable_str) :: desc_str
  integer :: n_data_min,n_data_max,n_sparse_max,covariance_type,l_max,n_max
  double precision :: delta,zeta,regularisation
  double precision :: atomic_gaussian_width
  double precision :: neb_error

  integer(8) :: seed
  integer :: max_steps,verbosity
  double precision :: t1,t2
  real(kind=rk) :: max_time,write_interval,min_height,max_height
  double precision :: max_wall_time
  double precision :: soap_cutoff
  integer :: print_interval
  character(len=80) :: xyz_out,pot_name,pot_style,predictor_file,restart_file
  character(len=80) :: zstring,depstring
  logical :: exist,debug,restart_read,allow_exchange,periodic(3),track_MSD
  integer :: i,error,internal_error,coord_jumping,coord_middle,coord_target,shell,ml_read,outfileunit

  integer :: first_id,last_id
  integer,allocatable :: atom_type(:)
  real(kind=rk),allocatable :: pos(:,:)

  character(len=512) :: cmds
  character(len=12),allocatable :: args(:)
  integer :: myrank,mpierror,ios
  double precision :: start_time

  error=0
  call read_input()

  allocate(args(9))
  args(1)="liblammps"
  args(2)="-log"
  args(3)="none"
  args(4)="-screen"
  args(5)="none"
  args(6)="-partition"
  args(7)=neb_images//"x1"
  args(8)="-in"
  args(9)="none"

  lmp=lammps(args=args)
  start_time=mpi_wtime()

  call mpi_comm_rank(MPI_COMM_WORLD,myrank,mpierror)

  call system_initialise()
  ! This silences most libAtoms prints:
  call verbosity_push(PRINT_SILENT)
  
  if (myrank==0) then
    inquire(file=pot_name,exist=exist)
    if (.not. exist) then
      write(0,"(a)") "*** Error: pot_name file "//trim(pot_name)//" doesn't exist"
      error=1
    endif
    if (error==0) then
      n_soap_Z=nZ
      n_soap_species=nZ
      allocate(soap_Z(n_soap_Z))
      allocate(soap_species_Z(n_soap_species))
      soap_Z(:)=Z(:)
      soap_species_Z(:)=Z(:)
      if (substrate_type==SUBSTRATE_TYPE_WEAK) then
        allocate(Z_tmp(nZ))
        Z_tmp=Z
        allocate(dep_rate_tmp(nZ))
        dep_rate_tmp=dep_rate
        nZ=nZ+1
        deallocate(Z,dep_rate)
        allocate(Z(nZ))
        allocate(dep_rate(nZ))
        Z(1)=Z_tmp(1)
        Z(2:)=Z_tmp(1:)
        deallocate(Z_tmp)
        dep_rate(1)=0.0_rk
        dep_rate(2:)=dep_rate_tmp(1:)
        deallocate(dep_rate_tmp)
        deallocate(soap_species_Z)
        n_soap_species=n_soap_species+1
        allocate(soap_species_Z(n_soap_species))
        soap_species_Z(2:)=Z(2:)
        soap_species_Z(1)=SUBSTRATE_Z
      endif
      restart_read=.false.
      inquire(file=restart_file,exist=exist)
      if (exist) then
        call read_my_restart(kmc,restart_file,error=internal_error)
        if (internal_error/=0) then
          write(0,"(a)") "*** Error reading "//trim(restart_file)
          error=internal_error
        endif
        if (error==0) then
          restart_read=.true.
          inquire(file="mt.state",exist=exist)
          if (exist) then
            call read_mt_state("mt.state",error=internal_error)
            if (internal_error/=0) then
              write(0,"(a)") "*** Error reading Mersenne Twister state from mt.state"
              error=internal_error
            else
              write(6,"(a)") "# Read Mersenne Twister state from mt.state, ignoring the seed parameter"
            endif
          else
            write(6,"(a)") "# Read kmc.restart but can't find mt.state; reseeding Mersenne Twister"
            call seed_mt(seed)
          endif

          nZ=kmc%nZ
          if (allocated(Z)) then
            deallocate(Z)
          endif
          allocate(Z(nZ))
          Z=kmc%Z
        endif
      else
        call seed_mt(seed)
      endif
    endif

    if (error==0) then
      write(zstring,*) nZ,"(x,i0)"
      write(depstring,*) nZ,"(g10.3)"
      if (verbosity>=VERBOSITY_QUIET) then
        if (restart_read) then
          write(6,"(a)") "# Parameters read from restart file, most of in.mlkmc ignored"
          kmc%debug=debug
          kmc%verbosity=verbosity
          kmc%dep_rate=dep_rate
          kmc%temp=temp
          kmc%kT=kB*temp
        else
          call print_params()
        endif
      endif
  
      if (restart_read) then
        inquire(file="final.system",exist=exist)
        if (exist) then
          write(6,"(a)") "# Reading final.system..."
          call read_system(system,"final.system",error=internal_error)
        else
          write(0,"(a)") "*** Error: read restart file but couldn't find final.system!"
          error=1
        endif
      else
        inquire(file="initial.system",exist=exist)
        if (exist) then
          write(6,"(a)") "# Reading initial.system..."
          call read_system(system,"initial.system",error=internal_error)
          if (internal_error/=0) then
            write(0,"(a)") "*** Error reading initial.system"
            error=internal_error
          endif
        else
          call create_system(system,xmax,ymax,zmax,lattice_constant,substrate_layers,&
                             substrate_species,substrate_type,max_jump_nth,lae_cutoff_nth,allow_exchange=allow_exchange,&
                             periodic=periodic,initial_atoms=initial_atoms,initial_atoms_layers=initial_atoms_layers,&
                             initial_atoms_termination=initial_atoms_termination,initial_atoms_coords=initial_atoms_coords,&
                             error=internal_error)
          if (internal_error/=0) then
            write(0,"(a)") "*** Error creating system"
            call print_system(system,"error.system")
            error=internal_error
          else
            write(6,"(a)") "# Created system"
            write(6,"(a)",advance="no") "# Writing initial.system..."
            call print_system(system,"initial.system",error=internal_error)
            if (internal_error/=0) then
              write(6,*)
              write(0,"(a)") "*** Error writing initial system to initial.system"
              error=internal_error
            else
              write(6,"(a)") " done!"
            endif
          endif
        endif
      endif
    endif
  endif
  call mpi_bcast(error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierror)
  if (error/=0) then
    call lmp%close()
    call mpi_finalize(mpierror)
    stop
  endif
  if (myrank==0) then
    cmds="units metal"//new_line("a")//&
         "dimension 3"//new_line("a")//&
         "atom_style atomic"//new_line("a")//&
         "atom_modify map array sort 0 0.0"//new_line("a")//&
         "region box block 0 "//real(system%lattice(1,1),kind=dp)//&
         " 0 "//real(system%lattice(2,2),kind=dp)//&
         " 0 "//real(system%lattice(3,3),kind=dp)//new_line("a")//&
         "create_box "//nZ//" box"//new_line("a")//&
         "pair_style "//pot_style//new_line("a")//&
         "pair_coeff * * "//pot_name
    do i=1,nZ
      if (substrate_TYPE==SUBSTRATE_TYPE_WEAK) then
        if (i==1) then
          cmds=trim(cmds)//" "//trim(ElementName(Z(i)))//"1"
        elseif (i==2) then
          cmds=trim(cmds)//" "//trim(ElementName(Z(i)))//"2"
        endif
      else
        cmds=trim(cmds)//" "//ElementName(Z(i))
      endif
    end do
  endif
  call mpi_bcast(error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierror)
  if (error/=0) then
    call lmp%close()
    call mpi_finalize(mpierror)
    stop
  endif

  call mpi_bcast(cmds,512,MPI_CHAR,0,MPI_COMM_WORLD,mpierror)
  call lmp%commands_string(cmds)

  ! Add frozen atoms to the LAMMPS object, since they will never move
  if (myrank==0) then
    first_id=system%frozen_substrate%first_id
    last_id=system%frozen_substrate%last_id
    allocate(atom_type(first_id:last_id))
    allocate(pos(3,first_id:last_id))
    atom_type=system%frozen_substrate%atom_type
    pos=system%frozen_substrate%pos
  endif
  call mpi_bcast(first_id,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierror)
  call mpi_bcast(last_id,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierror)
  if (myrank/=0) then
    allocate(atom_type(first_id:last_id))
    allocate(pos(3,first_id:last_id))
  endif
  call mpi_bcast(atom_type,last_id-first_id+1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierror)
  call mpi_bcast(pos,3*(last_id-first_id+1),MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,mpierror)
  call lmp%create_atoms(id=[(i,i=first_id,last_id)],&
                        type=atom_type,&
                        x=reshape(real(pos,kind=dp),[3*(last_id-first_id+1)]))
  deallocate(atom_type,pos)
  call lmp%command("group frozen union all")
  call lmp%command("fix freeze frozen setforce 0.0 0.0 0.0")
  
  if (myrank==0) then   
    if (restart_read) then
      kmc%system=system
    else
      soap_cutoff=get_Xnn_dist(lattice_constant,lae_cutoff_nth,error=internal_error)
      if (internal_error/=0) then
        write(0,"(a)") "*** Error setting LAE cutoff dist based on lae_cutoff_nth"
        error=internal_error
      endif
      soap_cutoff=soap_cutoff+get_Xnn_dist(lattice_constant,max_jump_nth,error=internal_error)
      if (internal_error/=0) then
        write(0,"(a)") "*** Error setting max jump dist based on max_jump_nth"
        error=internal_error
      endif
      desc_str="soap cutoff="//soap_cutoff//&
               " l_max="//l_max//" n_max="//n_max//&
               " cutoff_transition_width="//0.5d0//&
               " atom_gaussian_width="//atomic_gaussian_width//" n_Z="//n_soap_Z//" n_species="//n_soap_species//&
               " species_Z={"//soap_species_Z//"} Z={"//soap_Z//"}"
      call initialise_kmc(kmc,relax_iter,neb_images,neb_steps,neb_spring,&
                          neb_error,relax_tol,ml_tol,&
                          nZ,Z,attempt_frequency,system,verbosity,&
                          debug,track_MSD,n_data_min,n_data_max,n_sparse_max,desc_str,&
                          delta,zeta,covariance_type,regularisation,&
                          error=internal_error)
      if (internal_error/=0) then
        write(0,"(a)") "*** Error initialising KMC"
        error=internal_error
      endif
    endif

    allocate(kmc%hop_predictor(3:LOCK_THRESHOLD,3:N_NN_MAX,kmc%system%max_jump_shell))
    if (allow_exchange) then
      allocate(kmc%exchange_predictor(3:N_NN_MAX-1,3:LOCK_THRESHOLD,3:N_NN_MAX))
    endif

    write(6,"(a)",advance="no") "# Reading predictors..."

#ifdef INTEL_COMPILER
    inquire(directory="predictors",exist=exist)
#else
    inquire(file="predictors/",exist=exist)
#endif
    if (.not. exist) then
      inquire(file="predictors",exist=exist)
      if (.not. exist) then
        call execute_command_line("mkdir predictors",exitstat=error)
      else
        write(0,"(a)") "*** Error: can't mkdir predictors, because a file with that name exists"
        error=1
      endif
    endif

    if (error==0) then
      ml_read=0
      do shell=1,kmc%system%max_jump_shell
        do coord_target=3,N_NN_MAX
          do coord_jumping=3,LOCK_THRESHOLD
            if (coord_jumping<N_NN_MAX) then
              predictor_file="predictors/hop_predictor_"//coord_jumping//"_"//coord_target//"_"//shell//".xml"
              inquire(file=predictor_file,exist=exist)
              if (exist) then
                call read_predictor_xml(kmc%hop_predictor(coord_jumping,coord_target,shell),predictor_file)
                ml_read=ml_read+1
              endif
            endif
            if (allow_exchange) then
              if (shell==1) then
                do coord_middle=3,N_NN_MAX-1
                  predictor_file="predictors/exchange_predictor_"//coord_middle//"_"//coord_jumping//"_"//coord_target//".xml"
                  inquire(file=predictor_file,exist=exist)
                  if (exist) then
                    call read_predictor_xml(kmc%exchange_predictor(coord_middle,coord_jumping,&
                                                                          coord_target),predictor_file)
                    ml_read=ml_read+1
                  endif
                end do
              endif
            endif
          end do
        end do
      end do
      write(6,"(a,i0,a)") " done! Read ",ml_read," predictors from files"

      ! Allocate minimization cache
      if (.not. allocated(kmc%minim_cache)) then
        allocate(kmc%minim_cache(nZ,0:NN_SHELLS,MINIM_CACHE_SIZE))
        kmc%minim_cache=0
      endif

      open(newunit=outfileunit,file=xyz_out,action="write",iostat=ios,RECL=10000)
      if (ios/=0) then
        write(0,"(a)") "*** Error opening outfile "//xyz_out
        error=ios
      endif
    endif
  endif
  call mpi_bcast(error,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierror)
  if (error/=0) then
    call lmp%close()
    call mpi_finalize(mpierror)
    stop
  endif
  
  if (myrank==0) then
    t1=mpi_wtime()
  endif
  min_height=max(min_height,0.0_rk)
  max_height=min(max_height,kmc%system%lattice(3,3))
  call run_kmc(kmc,lmp,temp,dep_rate,max_steps,max_time,max_wall_time,min_height,max_height,&
               start_time,outfileunit,print_interval,write_interval,myrank,error=error)
  if (myrank==0) then
    t2=mpi_wtime()
    kmc%kmc_run_time=t2-t1
  endif

  if (myrank==0) then
    close(outfileunit)
    if (error/=0) then
      write(0,"(a)") "*** Error in running KMC"
    else
      if (verbosity>=VERBOSITY_QUIET) then
        write(6,"(a)") "# All done!"
      endif
      write(6,"(a)") "# Writing restart files..."
      call dump_data_final(kmc,error=internal_error)
      if (internal_error/=0) then
        write(0,"(a)") "*** Error writing KMC restart info"
      endif
    endif

    t2=mpi_wtime()
    write(6,"(a)") "# Total wall time: "//print_time(t2-start_time)
    write(6,"(a)") "#   Total KMC runtime: "//print_time(kmc%kmc_run_time)
    write(6,"(a)") "#     Getting rates:     "//print_time(kmc%get_rates_time)
    write(6,"(a)") "#       Overhead:          "//print_time(kmc%get_rates_overhead_time)
    write(6,"(a)") "#       Finding events:    "//print_time(kmc%event_finding_time)
    write(6,"(a)") "#       Create fingerprint:"//print_time(kmc%rate_fingerprint_time)
    write(6,"(a)") "#       Check cache:       "//print_time(kmc%rate_cache_time)
    write(6,"(a)") "#       Substr. LAE check: "//print_time(kmc%lae_sub_time)
    write(6,"(a)") "#       Getting LAEs:      "//print_time(kmc%lae_time)
    write(6,"(a)") "#         Init LAE:          "//print_time(kmc%lae_init_time)
    write(6,"(a)") "#         Add atoms:         "//print_time(kmc%lae_add_time)
    write(6,"(a)") "#         Calc connect:      "//print_time(kmc%lae_connect_time)
    write(6,"(a)") "#       Descriptors:       "//print_time(kmc%desc_time)
    write(6,"(a)") "#       Getting barriers:  "//print_time(kmc%get_barrier_time)
    write(6,"(a)") "#         ML predictions:    "//print_time(kmc%predict_time)
    write(6,"(a)") "#         ML fitting:        "//print_time(kmc%fit_time)
    write(6,"(a)") "#         Minimizations:     "//print_time(kmc%minimize_time)
    write(6,"(a)") "#         NEB calculations:  "//print_time(kmc%neb_barrier_time)
    write(6,"(a)") "#       Assinging rates:   "//print_time(kmc%assigning_time)
    write(6,"(a)") "#     Picking events:    "//print_time(kmc%pick_event_time)
    write(6,"(a)") "#     Executing events:    "//print_time(kmc%execute_jump_event_time)
    write(6,"(a)") "#       Create fingerprint:"//print_time(kmc%event_fingerprint_time)
    write(6,"(a)") "#       Check cache:       "//print_time(kmc%event_cache_time)
    write(6,"(a)") "#     Depositing atoms:  "//print_time(kmc%deposit_atom_time)
    write(6,"(a)") "#     Writing to files:  "//print_time(kmc%write_time)
  endif

  call lmp%close()
  call mpi_finalize(mpierror)

contains

  subroutine read_input()
    integer :: i,c,ios,n_initial_atoms_commands
    character(len=80) :: param_file,key,arg,initial_atoms_arg,periodic_args(3),termination_arg
    character(len=200) :: line

    ! Sensible defaults
    xmax=27
    ymax=47
    zmax=4
    periodic=[.true.,.true.,.false.]
    substrate_layers=2
    substrate_type=SUBSTRATE_TYPE_REGULAR
    n_initial_atoms_commands=0
    allocate(initial_atoms(10))
    allocate(initial_atoms_layers(10))
    allocate(initial_atoms_termination(10))
    allocate(initial_atoms_coords(3,10))
    initial_atoms=INITIAL_ATOMS_NONE
    initial_atoms_layers=0
    initial_atoms_termination=TERMINATION_NONE
    initial_atoms_coords=0.0d0
    relax_iter=1000
    neb_images=11
    neb_steps=1000
    neb_spring=1.0
    neb_error=0.1d0
    relax_tol=1.0d-7
    ml_tol=1.0d-2
    xyz_out="out.xyz"
    restart_file="kmc.restart"
    seed=12345
    n_data_min=100
    n_data_max=10000
    n_sparse_max=1000
    l_max=3
    n_max=4
    atomic_gaussian_width=1.0d0
    delta=20.0d0
    zeta=20.0d0
    covariance_type=COVARIANCE_DOT_PRODUCT
    regularisation=0.01d0
    max_time=0.0d0
    max_wall_time=604000 ! 800 seconds less than one week
    min_height=-huge(min_height)
    max_height=huge(max_height)
    write_interval=0.0d0
    verbosity=VERBOSITY_NORMAL
    debug=.false.
    track_MSD=.false.
    allow_exchange=.true.

    ! Initial values for checking duplicates
    pot_name=""
    pot_style=""
    lattice_constant=-1.0d0
    substrate_species=-1
    nZ=-1
    temp=-1.0d0
    attempt_frequency=-1.0d0
    max_steps=-1
    print_interval=-1
    lae_cutoff_nth=-1
    max_jump_nth=-1

    param_file="in.mlkmc"

    open(unit=10,file=param_file,action="read",iostat=ios)
    if (ios/=0) then
      write(0,"(a)") "*** Error opening parameter file "//trim(param_file)
      stop
    endif

    i=0
    do
      read(10,"(a)",iostat=ios) line
      if (ios/=0) exit
      i=i+1
      if (line(1:1)=="#" .or. line(1:1)=="") then
        cycle
      endif
      read(line,*,iostat=ios) key,arg
      select case(key)
        case("xmax")
          read(arg,*,iostat=ios) xmax
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read xmax from file "//trim(param_file)//" line ",i
            stop
          endif
          if (xmax<=0) then
            write(0,"(a,i0)") "*** Error: xmax>1 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("ymax")
          read(arg,*,iostat=ios) ymax
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read ymax from file "//trim(param_file)//" line ",i
            stop
          endif
          if (ymax<=0) then
            write(0,"(a,i0)") "*** Error: ymax>1 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("zmax")
          read(arg,*,iostat=ios) zmax
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read zmax from file "//trim(param_file)//" line ",i
            stop
          endif
          if (zmax<=0) then
            write(0,"(a,i0)") "*** Error: zmax>1 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("periodic")
          read(line,"(a8,a)",iostat=ios) key,arg
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read periodic from file "//trim(param_file)//" line ",i
            stop
          endif
          read(arg,*,iostat=ios) periodic_args
          do c=1,3
            arg=periodic_args(c)
            if (any(arg(1:1)==["T","t","Y","y","1"])) then
              periodic(c)=.true.
            elseif (any(arg(1:1)==["F","f","N","n","0"])) then
              periodic(c)=.false.
            else
              write(0,"(a,i0,a,i0)") "*** Error: unrecognised/unimplemented periodic flag ",c," "//trim(arg)//" on line ",i
            endif
          end do
        case("substrate_layers")
          read(arg,*,iostat=ios) substrate_layers
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read substrate_layers from file "//trim(param_file)//" line ",i
            stop
          endif
          if (substrate_layers<=0) then
            write(0,"(a,i0)") "*** Error: substrate_layers>1 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("substrate_species")
          read(arg,*,iostat=ios) substrate_species
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read substrate_species from file "//trim(param_file)//" line ",i
            stop
          endif
          if (substrate_species<=0) then
            write(0,"(a,i0)") "*** Error: substrate_species>1 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("substrate_type")
          select case(arg)
            case("regular")
              substrate_type=SUBSTRATE_TYPE_REGULAR
            case("weak")
              substrate_type=SUBSTRATE_TYPE_WEAK
            case default
              write(0,"(a,i0)") "*** Error: unimplemented substrate_type "//trim(arg)//" on line ",i
              stop
          end select
        case("relax_iter")
          read(arg,*,iostat=ios) relax_iter
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read relax_iter from file "//trim(param_file)//" line ",i
            stop
          endif
          if (relax_iter<=0) then
            write(0,"(a,i0)") "*** Error: relax_iter>1 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("neb_images")
          read(arg,*,iostat=ios) neb_images
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read neb_images from file "//trim(param_file)//" line ",i
            stop
          endif
          if (relax_iter<3) then
            write(0,"(a,i0)") "*** Error: relax_iter>=3 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("neb_steps")
          read(arg,*,iostat=ios) neb_steps
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read neb_steps from file "//trim(param_file)//" line ",i
            stop
          endif
          if (neb_steps<=0) then
            write(0,"(a,i0)") "*** Error: neb_steps>0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("neb_spring")
          read(arg,*,iostat=ios) neb_spring
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read neb_spring from file "//trim(param_file)//" line ",i
            stop
          endif
          if (neb_spring<=0) then
            write(0,"(a,i0)") "*** Error: neb_spring>0.0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("relax_tol")
          read(arg,*,iostat=ios) relax_tol
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read relax_tol from file "//trim(param_file)//" line ",i
            stop
          endif
          if (relax_tol<0) then
            write(0,"(a,i0)") "*** Error: ml_tol>=0.0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("ml_tol")
          read(arg,*,iostat=ios) ml_tol
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read ml_tol from file "//trim(param_file)//" line ",i
            stop
          endif
          if (ml_tol<0) then
            write(0,"(a,i0)") "*** Error: ml_tol>=0.0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("xyz_out")
          xyz_out=arg
        case("restart_file")
          restart_file=arg
        case("seed")
          read(arg,*,iostat=ios) seed
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read seed from file "//trim(param_file)//" line ",i
            stop
          endif
          if (seed<0) then
            write(0,"(a,i0)") "*** Error: seed>=0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("n_data_min")
          read(arg,*,iostat=ios) n_data_min
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read n_data_min from file "//trim(param_file)//" line ",i
            stop
          endif
          if (n_data_min<=0) then
            write(0,"(a,i0)") "*** Error: n_data_min>0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("n_data_max")
          read(arg,*,iostat=ios) n_data_max
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read n_data_max from file "//trim(param_file)//" line ",i
            stop
          endif
          if (n_data_max<=0) then
            write(0,"(a,i0)") "*** Error: n_data_max>0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("n_sparse_max")
          read(arg,*,iostat=ios) n_sparse_max
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read n_sparse_max from file "//trim(param_file)//" line ",i
            stop
          endif
          if (n_sparse_max<=0) then
            write(0,"(a,i0)") "*** Error: n_sparse_max>0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("neb_error")
          read(arg,*,iostat=ios) neb_error
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read neb_error from file "//trim(param_file)//" line ",i
            stop
          endif
          if (neb_error<=0) then
            write(0,"(a,i0)") "*** Error: neb_error>0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("l_max")
          read(arg,*,iostat=ios) l_max
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read l_max from file "//trim(param_file)//" line ",i
            stop
          endif
          if (l_max<=0) then
            write(0,"(a,i0)") "*** Error: l_max>0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("n_max")
          read(arg,*,iostat=ios) n_max
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read n_max from file "//trim(param_file)//" line ",i
            stop
          endif
          if (n_max<=0) then
            write(0,"(a,i0)") "*** Error: n_max>0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("delta")
          read(arg,*,iostat=ios) delta
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read delta from file "//trim(param_file)//" line ",i
            stop
          endif
          if (delta<=0) then
            write(0,"(a,i0)") "*** Error: delta>0.0 shouldn't be "//trim(arg)//" on line ",i
            stop
          endif
        case("atomic_gaussian_width")
          read(arg,*,iostat=ios) atomic_gaussian_width
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read atomic_gaussian_width from file "//trim(param_file)//" line ",i
            stop
          endif
          if (delta<=0) then
            write(0,"(a,i0)") "*** Error: atomic_gaussian_width>0.0 shouldn't be "//trim(arg)//" on line ",i
            stop
          endif
        case("zeta")
          read(arg,*,iostat=ios) zeta
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read zeta from file "//trim(param_file)//" line ",i
            stop
          endif
          if (zeta<0) then
            write(0,"(a,i0)") "*** Error: zeta>=0.0 shouldn't be "//trim(arg)//" on line ",i
            stop
          endif
        case("covariance_type")
          select case(arg)
            case("dot_product")
              covariance_type=COVARIANCE_DOT_PRODUCT
            case("ard_se")
              covariance_type=COVARIANCE_ARD_SE
            case default
              write(0,"(a,i0)") "*** Error: unknown/unimplemented covariance type "//trim(arg)//" on line ",i
              stop
          end select
        case("max_time")
          read(arg,*,iostat=ios) max_time
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read max_time from file "//trim(param_file)//" line ",i
            stop
          endif
          if (max_time<0) then
            write(0,"(a,i0)") "*** Error: max_time>=0.0 shouldn't be "//trim(arg)//" on line ",i
            stop
          endif
        case("max_wall_time")
          read(arg,*,iostat=ios) max_wall_time
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read max_wall_time from file "//trim(param_file)//" line ",i
            stop
          endif
          if (max_wall_time<0) then
            write(0,"(a,i0)") "*** Error: max_wall_time>=0.0 shouldn't be "//trim(arg)//" on line ",i
            stop
          endif
        case("min_height")
          read(arg,*,iostat=ios) min_height
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read min_height from file "//trim(param_file)//" line ",i
            stop
          endif
        case("max_height")
          read(arg,*,iostat=ios) max_height
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read max_height from file "//trim(param_file)//" line ",i
            stop
          endif
        case("verbosity")
          select case(arg)
            case("silent")
              verbosity=VERBOSITY_SILENT
            case("quiet")
              verbosity=VERBOSITY_QUIET
            case("normal")
              verbosity=VERBOSITY_NORMAL
            case("loud")
              verbosity=VERBOSITY_LOUD
            case("absurd")
              verbosity=VERBOSITY_ABSURD
            case default
              write(0,"(a,i0)") "*** Error: unrecognised/unimplemented verbosity "//trim(arg)//" on line ",i
              stop
          end select
        case("debug")
          if (any(arg(1:1)==["T","t","Y","y","1"])) then
            debug=.true.
          elseif (any(arg(1:1)==["F","f","N","n","0"])) then
            debug=.false.
          else
              write(0,"(a,i0)") "*** Error: unrecognised/unimplemented debug flag "//trim(arg)//" on line ",i
            stop
          endif
        case("track_MSD")
          if (any(arg(1:1)==["T","t","Y","y","1"])) then
            track_MSD=.true.
          elseif (any(arg(1:1)==["F","f","N","n","0"])) then
            track_MSD=.false.
          else
              write(0,"(a,i0)") "*** Error: unrecognised/unimplemented track_MSD flag "//trim(arg)//" on line ",i
            stop
          endif
        case("allow_exchange")
          if (any(arg(1:1)==["T","t","Y","y","1"])) then
            allow_exchange=.true.
          elseif (any(arg(1:1)==["F","f","N","n","0"])) then
            allow_exchange=.false.
          else
              write(0,"(a,i0)") "*** Error: unrecognised/unimplemented allow_exchange flag "//trim(arg)//" on line ",i
            stop
          endif
        case("regularisation")
          read(arg,*,iostat=ios) regularisation
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read regularisation from file "//trim(param_file)//" line ",i
            stop
          endif
          if (regularisation<0) then
            write(0,"(a,i0)") "*** Error: regularisation>=0.0 shouldn't be "//trim(arg)//" on line ",i
            stop
          endif
        case("pot_name")
          if (pot_name/="") then
            write(0,"(a,i0)") "*** Error: duplicate parameter pot_file on line ",i
            stop
          endif
          pot_name=arg
        case("pot_style")
          if (pot_style/="") then
            write(0,"(a,i0)") "*** Error: duplicate parameter pot_style on line ",i
            stop
          endif
          pot_style=arg
        case("lattice_constant")
          if (lattice_constant>=0.0d0) then
            write(0,"(a,i0)") "*** Error: duplicate parameter lattice_constant on line ",i
            stop
          endif
          read(arg,*,iostat=ios) lattice_constant
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read lattice_constant from file "//trim(param_file)//" line ",i
            stop
          endif
          if (lattice_constant<=0) then
            write(0,"(a,i0)") "*** Error: lattice_constant>0.0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("nZ")
          if (nZ>0) then
            write(0,"(a,i0)") "*** Error: duplicate parameter nZ on line ",i
            stop
          endif
          read(arg,*,iostat=ios) nZ
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read nZ from file "//trim(param_file)//" line ",i
            stop
          endif
          if (nZ<1) then
            write(0,"(a,i0)") "*** Error: nZ>=1 can't be "//trim(arg)//" on line ",i
            stop
          endif
          allocate(Z(nZ))
          allocate(dep_rate(nZ))
          Z=-1
          dep_rate=0.0d0 ! Default deposition rate
        case("initial_atoms")
          n_initial_atoms_commands=n_initial_atoms_commands+1
          if (n_initial_atoms_commands>size(initial_atoms)) then
            call size_up(initial_atoms,default_value=INITIAL_ATOMS_NONE,custom_buffer_size=10,error=error)
            if (error/=0) then
              write(0,"(a)") "*** Error sizing up initial_atoms array in read_input()"
              stop
            endif
            call size_up(initial_atoms_layers,default_value=0,custom_buffer_size=10,error=error)
            if (error/=0) then
              write(0,"(a)") "*** Error sizing up initial_atoms_layers array in read_input()"
              stop
            endif
            call size_up(initial_atoms_termination,default_value=TERMINATION_NONE,custom_buffer_size=10,error=error)
            if (error/=0) then
              write(0,"(a)") "*** Error sizing up initial_atoms_termination array in read_input()"
              stop
            endif
            call size_up(initial_atoms_coords,rank=2,default_value=0.0d0,custom_buffer_size=10,error=error)
            if (error/=0) then
              write(0,"(a)") "*** Error sizing up initial_atoms_coords array in read_input()"
              stop
            endif
          endif
          read(line,"(a13,a)",iostat=ios) key,arg
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read initial_atoms from file "//trim(param_file)//" line ",i
            stop
          endif
          read(arg,*,iostat=ios) key,initial_atoms_arg
          select case(key)
            case("random")
              initial_atoms(n_initial_atoms_commands)=INITIAL_ATOMS_RANDOM
              read(arg,*,iostat=ios) key,initial_atoms_layers(n_initial_atoms_commands)
              if (ios/=0) then
                write(0,"(a,i0)") "*** Error parsing the number of initial atoms from line ",i
                stop
              endif
            case("hexagon")
              initial_atoms(n_initial_atoms_commands)=INITIAL_ATOMS_HEXAGON
              read(arg,*,iostat=ios) key,initial_atoms_layers(n_initial_atoms_commands)
              if (ios/=0) then
                write(0,"(a,i0)") "*** Error parsing the number of hexagon layers from line ",i
                stop
              endif
            case("hexagon_plus_adatom")
              initial_atoms(n_initial_atoms_commands)=INITIAL_ATOMS_HEXAGON_PLUS_ADATOM
              read(arg,*,iostat=ios) key,initial_atoms_layers(n_initial_atoms_commands)
              if (ios/=0) then
                write(0,"(a,i0)") "*** Error parsing the number of hexagon layers from line ",i
                stop
              endif
            case("triangle")
              initial_atoms(n_initial_atoms_commands)=INITIAL_ATOMS_TRIANGLE
              read(arg,*,iostat=ios) key,initial_atoms_layers(n_initial_atoms_commands),termination_arg
              if (ios/=0) then
                write(0,"(a,i0)") "*** Error parsing the number of triangle layers and termination from line ",i
                stop
              endif
            case("triangle_plus_adatom")
              initial_atoms(n_initial_atoms_commands)=INITIAL_ATOMS_TRIANGLE_PLUS_ADATOM
              read(arg,*,iostat=ios) key,initial_atoms_layers(n_initial_atoms_commands),termination_arg
              if (ios/=0) then
                write(0,"(a,i0)") "*** Error parsing the number of triangle layers and termination from line ",i
                stop
              endif
            case("single")
              initial_atoms(n_initial_atoms_commands)=INITIAL_ATOMS_SINGLE
              read(arg,*,iostat=ios) key,initial_atoms_coords(:,n_initial_atoms_commands)
              if (ios/=0) then
                write(0,"(a,i0)") "*** Error parsing the initial atom coordinates from line ",i
                stop
              endif
            case default
              write(0,"(a,i0,a)") "*** Error: unsupported initial_atoms style "//trim(arg)//" on line ",i,&
                             "; supported styles are random, hexagon, hexagon_plus_adatom, triangle and triangle_plus_adatom"
              stop
          end select
          if (initial_atoms(n_initial_atoms_commands)==INITIAL_ATOMS_TRIANGLE &
              .or. initial_atoms(n_initial_atoms_commands)==INITIAL_ATOMS_TRIANGLE_PLUS_ADATOM) then
            if (any(termination_arg(1:1)==["A","a"])) then
              initial_atoms_termination(n_initial_atoms_commands)=TERMINATION_A
            elseif (any(termination_arg(1:1)==["B","b"])) then
              initial_atoms_termination(n_initial_atoms_commands)=TERMINATION_B
            else
              write(0,"(a,i0,a,i0)") "*** Error: initial_atoms_termination==[A,B] can't be ",trim(termination_arg)," on line ",i
              stop
            endif
          endif
        case("Z")
          if (Z(1)>0) then
            write(0,"(a,i0)") "*** Error: duplicate parameter Z on line ",i
            stop
          endif
          read(line,"(a1,a)",iostat=ios) key,arg
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read Z from file "//trim(param_file)//" line ",i
            stop
          endif
          read(arg,*,iostat=ios) Z
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't interpret Z "//trim(arg)//" line ",i
            stop
          endif
          if (any(Z<=0)) then
            write(0,"(a,i0)") "*** Error: any Z>0.0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("dep_rate")
          read(line,"(a8,a)",iostat=ios) key,arg
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read dep_rate from file "//trim(param_file)//" line ",i
            stop
          endif
          read(arg,*,iostat=ios) dep_rate
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't interpret dep_rate "//trim(arg)//" line ",i
            stop
          endif
          if (any(dep_rate<0)) then
            write(0,"(a,i0)") "*** Error: any dep_rate>=0.0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("T")
          if (temp>=0.0d0) then
            write(0,"(a,i0)") "*** Error: duplicate parameter T on line ",i
            stop
          endif
          read(arg,*,iostat=ios) temp
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read T from file "//trim(param_file)//" line ",i
            stop
          endif
          if (temp<=0) then
            write(0,"(a,i0)") "*** Error: T>0.0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("attempt_frequency")
          if (attempt_frequency>=0.0d0) then
            write(0,"(a,i0)") "*** Error: duplicate parameter attempt_frequency on line ",i
            stop
          endif
          read(arg,*,iostat=ios) attempt_frequency
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read attempt_frequency from file "//trim(param_file)//" line ",i
            stop
          endif
          if (attempt_frequency<=0) then
            write(0,"(a,i0)") "*** Error: attempt_frequency>0.0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("max_steps")
          if (max_steps>=0) then
            write(0,"(a,i0)") "*** Error: duplicate parameter max_steps on line ",i
            stop
          endif
          read(arg,*,iostat=ios) max_steps
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read max_steps from file "//trim(param_file)//" line ",i
            stop
          endif
          if (max_steps<0) then
            write(0,"(a,i0)") "*** Error: max_steps>=0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("print_interval")
          if (print_interval>=0) then
            write(0,"(a,i0)") "*** Error: duplicate parameter print_interval on line ",i
            stop
          endif
          read(arg,*,iostat=ios) print_interval
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read print_interval from file "//trim(param_file)//" line ",i
            stop
          endif
          if (print_interval<0) then
            write(0,"(a,i0)") "*** Error: print_interval>=0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case("write_interval")
          read(arg,*,iostat=ios) write_interval
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read write_interval from file "//trim(param_file)//" line ",i
            stop
          endif
       case("lae_cutoff_nth")
          if (lae_cutoff_nth>=0) then
            write(0,"(a,i0)") "*** Error: duplicate parameter lae_cutoff_nth on line ",i
            stop
          endif
          read(arg,*,iostat=ios) lae_cutoff_nth
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read lae_cutoff_nth from file "//trim(param_file)//" line ",i
            stop
          endif
          if (lae_cutoff_nth<0) then
            write(0,"(a,i0)") "*** Error: lae_cutoff_nth>=0 can't be "//trim(arg)//" on line ",i
            stop
          endif
       case("max_jump_nth")
          if (max_jump_nth>=0) then
            write(0,"(a,i0)") "*** Error: duplicate parameter max_jump_nth on line ",i
            stop
          endif
          read(arg,*,iostat=ios) max_jump_nth
          if (ios/=0) then
            write(0,"(a,i0)") "*** Error: couldn't read max_jump_nth from file "//trim(param_file)//" line ",i
            stop
          endif
          if (max_jump_nth<0) then
            write(0,"(a,i0)") "*** Error: max_jump_nth>=0 can't be "//trim(arg)//" on line ",i
            stop
          endif
        case default
          write(0,"(a,i0)") "*** Error: unknown keyword "//trim(key)//" on line ",i
          stop
      end select
    end do

    ! Check we got all the mandatory parameters
    if (pot_name=="") then
      write(0,"(a)") "*** Error: didn't find mandatory parameter pot_name in input file."
      stop
    endif
    if (pot_style=="") then
      write(0,"(a)") "*** Error: didn't find mandatory parameter pot_style in input file."
      stop
    endif
    if (lattice_constant<=0.0d0) then
      write(0,"(a)") "*** Error: didn't find mandatory parameter lattice_constant in input file."
      stop
    endif
    if (nZ<1) then
      write(0,"(a)") "*** Error: didn't find mandatory parameter nZ in input file."
      stop
    endif
    if (Z(1)<1) then
      write(0,"(a)") "*** Error: didn't find mandatory parameter Z in input file."
      stop
    endif
    if (temp<=0.0d0) then
      write(0,"(a)") "*** Error: didn't find mandatory parameter T in input file."
      stop
    endif
    if (attempt_frequency<=0.0d0) then
      write(0,"(a)") "*** Error: didn't find mandatory parameter attempt_frequency in input file."
      stop
    endif
    if (max_steps<0) then
      write(0,"(a)") "*** Error: didn't find mandatory parameter max_steps in input file."
      stop
    endif
    if (print_interval<0) then
      write(0,"(a)") "*** Error: didn't find mandatory parameter print_interval in input file."
      stop
    endif
    if (lae_cutoff_nth<=0) then
      write(0,"(a)") "*** Error: didn't find mandatory parameter lae_cutoff_nth in input file."
      stop
    endif
    if (max_jump_nth<=0) then
      write(0,"(a)") "*** Error: didn't find mandatory parameter max_jump_nth in input file."
      stop
    endif
  end subroutine read_input

  subroutine print_params()
    integer :: i

    write(6,"(a)") "# Starting program with these parameters:"
    write(6,"(a)") "# pot_name """//trim(pot_name)//""""
    write(6,"(a)") "# pot_style """//trim(pot_style)//""""
    write(6,"(a,f10.6)") "# lattice_constant ",lattice_constant
    write(6,"(a,i0)") "# xmax ",xmax
    write(6,"(a,i0)") "# ymax ",ymax
    write(6,"(a,i0)") "# zmax ",zmax
    write(6,"(a,3(l2))") "# periodic ",periodic
    write(6,"(a,i0)") "# substrate_layers ",substrate_layers
    write(6,"(a,i0)") "# substrate_species ",substrate_species
    write(6,"(a)") "# substrate_type "//trim(print_substrate_type(substrate_type))
    write(6,"(a,i0)") "# nZ ",nZ
    write(6,"(a,"//zstring//")") "# Z ",Z
    write(6,"(a,f10.6)") "# T ",temp
    write(6,"(a,g10.3)") "# attempt_frequency ",attempt_frequency
    write(6,"(a,i0)") "# max_steps ",max_steps
    write(6,"(a,i0)") "# print_interval ",print_interval
    write(6,"(a,g10.3,a)") "# write_interval ",write_interval," s"
    write(6,"(a,i0)") "# lae_cutoff_nth ",lae_cutoff_nth
    write(6,"(a,i0)") "# max_jump_nth ",max_jump_nth
    if (max_time>0.0d0) then
      write(6,"(a,g10.3,a)") "# max_time ",max_time," s"
    endif
    write(6,"(a,g10.3,a)") "# max_wall_time ",max_wall_time," s"
    if (min_height>-huge(min_height)) then
      write(6,"(a,g10.3,a)") "# min_height ",min_height," Å"
    endif
    if (max_height<huge(max_height)) then
      write(6,"(a,g10.3,a)") "# max_height ",max_height," Å"
    endif
    write(6,"(a)") "# verbosity "//trim(print_verbosity(verbosity))
    write(6,"(a)") "# debug "//debug
    write(6,"(a)") "# track_MSD "//track_MSD
    write(6,"(a)") "# allow_exchange "//allow_exchange
    write(6,"(a)") "# xyz_out "//trim(xyz_out)
    write(6,"(a,i0)") "# seed ",seed
    write(6,"(a,i0)") "# n_data_min ",n_data_min
    write(6,"(a,i0)") "# n_data_max ",n_data_max
    write(6,"(a,i0)") "# n_sparse_max ",n_sparse_max
    write(6,"(a,i0)") "# relax_iter ",relax_iter
    write(6,"(a,g10.1)") "# relax_tol ",relax_tol
    write(6,"(a,i0)") "# neb_images ",neb_images
    write(6,"(a,i0)") "# neb_steps ",neb_steps
    write(6,"(a,f10.6)") "# neb_spring ",neb_spring
    write(6,"(a,f10.6)") "# neb_error ",neb_error
    write(6,"(a,f10.6)") "# ml_tol ",ml_tol
    write(6,"(a,"//depstring//")") "# dep_rate ",dep_rate
    do i=1,size(initial_atoms)
      if (initial_atoms(i)/=INITIAL_ATOMS_NONE) then
        if (initial_atoms(i)==INITIAL_ATOMS_SINGLE) then
          write(6,"(a,3(f10.6))") "# initial_atoms "//trim(print_initial_atoms(initial_atoms(i)))//" ",&
                                  initial_atoms_coords(:,i)
        else
          write(6,"(a,i0)",advance="no") "# initial_atoms "//trim(print_initial_atoms(initial_atoms(i)))//" ",&
                                         initial_atoms_layers(i)
          if (initial_atoms(i)==INITIAL_ATOMS_TRIANGLE .or. initial_atoms(i)==INITIAL_ATOMS_TRIANGLE_PLUS_ADATOM) then
            write(6,"(a)") " "//trim(print_termination(initial_atoms_termination(i)))
          else
            write(6,*)
          endif
        endif
      endif
    end do
    write(6,"(a,i0)") "# l_max ",l_max
    write(6,"(a,i0)") "# n_max ",n_max
    write(6,"(a,f10.6)") "# atomic_gaussian_width ",atomic_gaussian_width
    write(6,"(a,f10.6)") "# delta ",delta
    write(6,"(a,f10.6)") "# zeta ",zeta
    write(6,"(a)") "# covariance_type "//trim(print_covariance_type(covariance_type))
    write(6,"(a,f10.6)") "# regularisation ",regularisation
  end subroutine print_params

  function print_verbosity(verbosity) result(verb_str)
    integer,intent(in) :: verbosity
    character(len=80) :: verb_str

    select case(verbosity)
      case(VERBOSITY_SILENT)
        verb_str="silent"
      case(VERBOSITY_QUIET)
        verb_str="quiet"
      case(VERBOSITY_NORMAL)
        verb_str="normal"
      case(VERBOSITY_LOUD)
        verb_str="loud"
      case(VERBOSITY_ABSURD)
        verb_str="absurd"
      case default
        verb_str=""
        write(0,"(a)") "*** Unimplemented verbosity in print_verbosity()"
    end select
  end function print_verbosity

  function print_initial_atoms(initial_atoms) result(init_atoms_str)
    integer,intent(in) :: initial_atoms
    character(len=80) :: init_atoms_str

    select case(initial_atoms)
      case(INITIAL_ATOMS_NONE)
        init_atoms_str="none"
      case(INITIAL_ATOMS_RANDOM)
        init_atoms_str="random"
      case(INITIAL_ATOMS_HEXAGON)
        init_atoms_str="hexagon"
      case(INITIAL_ATOMS_HEXAGON_PLUS_ADATOM)
        init_atoms_str="hexagon_plus_adatom"
      case(INITIAL_ATOMS_TRIANGLE)
        init_atoms_str="triangle"
      case(INITIAL_ATOMS_TRIANGLE_PLUS_ADATOM)
        init_atoms_str="triangle_plus_adatom"
      case(INITIAL_ATOMS_SINGLE)
        init_atoms_str="single"
      case default
        init_atoms_str=""
        write(0,"(a)") "*** Unimplemented initial_atoms in print_initial_atoms()"
    end select
  end function print_initial_atoms

  character(80) function print_termination(termination)
    integer,intent(in) :: termination

    select case(termination)
      case(TERMINATION_A)
        print_termination="A"
      case(TERMINATION_B)
        print_termination="B"
      case default
        write(0,"(a)") "*** Warning: no termination assigned to initial triangle"
    end select

  end function print_termination

  function print_substrate_type(substrate_type) result(sub_type_str)
    integer,intent(in) :: substrate_type
    character(len=80) :: sub_type_str

    select case(substrate_type)
      case(SUBSTRATE_TYPE_REGULAR)
        sub_type_str="regular"
      case(SUBSTRATE_TYPE_WEAK)
        sub_type_str="weak"
      case default
        sub_type_str=""
        write(0,"(a)") "*** Unimplemented substrate_type in print_substrate_type()"
    end select
  end function print_substrate_type

end program mlkmc
