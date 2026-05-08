! MLKMC - program for reading an atomic configuration and writing
! the corresponding descriptor in descriptor.dat
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

program make_descriptor
  use libAtoms_module
  use constants_module
  use descriptors_module
  implicit none

  integer :: l_max,n_max,nZ,error,ios,n_soap_species,weak_substrate,i
  integer,allocatable :: Z(:),soap_species_Z(:)
  double precision :: soap_cutoff,atomic_gaussian_width,cutoff_transition_width
  type(descriptor) :: desc
  type(descriptor_data) :: desc_data
  type(extendable_str) :: desc_str
  type(Atoms) :: lae
  type(CInOutput) :: infile
  character(len=80) :: filename

  call read_cmd()

  call system_initialise()
  !call verbosity_push(PRINT_SILENT)

  call initialise(infile,filename)
  call read(lae,infile,error=error)
  if (error/=0) then
    write(0,"(a,i0)") "*** Error: reading atoms from "//trim(filename)//" failed with error ",error
    stop
  endif
  write(6,"(i0,a)") lae%N," atoms read from file "//trim(filename)
  call calc_connect(lae)
  soap_cutoff=min(soap_cutoff,lae%cutoff)

  n_soap_species=nZ
  allocate(soap_species_Z(n_soap_species))
  soap_species_Z(:)=Z(:)
  if (weak_substrate==1) then
    deallocate(soap_species_Z)
    n_soap_species=n_soap_species+1
    allocate(soap_species_Z(n_soap_species))
    soap_species_Z(2:)=Z(:)
    soap_species_Z(1)=SUBSTRATE_Z
  endif

  desc_str="soap cutoff="//soap_cutoff//&
           " l_max="//l_max//" n_max="//n_max//&
           " cutoff_transition_width="//cutoff_transition_width//&
           " atom_gaussian_width="//atomic_gaussian_width//" n_Z="//nZ//" n_species="//n_soap_species//&
           " species_Z={"//soap_species_Z//"} Z={"//Z//"}"
  print *,string(desc_str)

  call initialise(desc,string(desc_str),error=error)
  if (error/=0) then
    write(0,"(a,i0)") "*** Error: couldn't initialise descriptor"
    stop
  endif
  write(6,"(a,i0)") "Descriptor dimensions: ",descriptor_dimensions(desc)

  call calc(desc,lae,desc_data,do_descriptor=.true.,args_str="atom_mask_name=jumping",error=error)
  if (error/=0) then
    write(0,"(a,i0)") "*** Error: couldn't form descriptor_data"
    stop
  endif
  open(unit=10,file="descriptor.dat",action="write",iostat=ios)
  if (ios/=0) then
    write(0,"(a,i0)") "*** Error opening file descriptor.dat"
    stop
  endif
  do i=1,size(desc_data%x(1)%data(:))
    write(10,*,iostat=ios) desc_data%x(1)%data(i)
  end do
  if (ios/=0) then
    write(0,"(a,i0)") "*** Error writing descriptor to descriptor.dat"
    stop
  endif
  close(10)
  write(6,"(a)") "Done writing descriptor to descriptor.dat!"


contains

  subroutine read_cmd()
    integer :: i,iarg,ios
    character(len=80) :: arg

    iarg=command_argument_count()
    if (iarg<9) then
      call get_command_argument(0,arg)
      write(0,"(a)") "Usage: "//trim(arg)//" file soap_cutoff(max) l_max n_max &
                      atomic_gaussian_width cutoff_transition_width nZ Z(:) weak_substrate"
      stop
    endif

    call get_command_argument(1,filename)
    call get_command_argument(2,arg)
    read(arg,*,iostat=ios) soap_cutoff
    if (ios/=0) then
      write(0,"(a)") "*** Error: invalid real "//trim(arg)//" for soap_cutoff"
      stop
    endif

    call get_command_argument(3,arg)
    read(arg,*,iostat=ios) l_max
    if (ios/=0) then
      write(0,"(a)") "*** Error: invalid integer "//trim(arg)//" for l_max"
      stop
    endif

    call get_command_argument(4,arg)
    read(arg,*,iostat=ios) n_max
    if (ios/=0) then
      write(0,"(a)") "*** Error: invalid integer "//trim(arg)//" for n_max"
      stop
    endif

    call get_command_argument(5,arg)
    read(arg,*,iostat=ios) atomic_gaussian_width
    if (ios/=0) then
      write(0,"(a)") "*** Error: invalid real "//trim(arg)//" for atomic_gaussian_width"
      stop
    endif

    call get_command_argument(6,arg)
    read(arg,*,iostat=ios) cutoff_transition_width
    if (ios/=0) then
      write(0,"(a)") "*** Error: invalid real "//trim(arg)//" for cutoff_transition_width"
      stop
    endif

    call get_command_argument(7,arg)
    read(arg,*,iostat=ios) nZ
    if (ios/=0) then
      write(0,"(a)") "*** Error: invalid integer "//trim(arg)//" for nZ"
      stop
    endif

    if (iarg<7+nZ) then
      write(0,"(a)") "*** Error: too few arguments for Z"
      call get_command_argument(0,arg)
      write(0,"(a)") "Usage: "//trim(arg)//" file soap_cutoff l_max n_max atomic_gaussian_width cutoff_transition_width nZ Z(:)"
      stop
    endif

    allocate(Z(nZ))
    do i=1,nZ
      call get_command_argument(7+i,arg)
      read(arg,*,iostat=ios) Z(i)
      if (ios/=0) then
        write(0,"(a)") "*** Error: invalid integer "//trim(arg)//" for Z("//i//")"
        stop
      endif
    end do

    call get_command_argument(7+nZ+1,arg)
    read(arg,*,iostat=ios) weak_substrate
    if (ios/=0) then
      write(0,"(a)") "*** Error: invalid integer "//trim(arg)//" for weak_substrate"
      stop
    endif

  end subroutine read_cmd
  
end program make_descriptor
