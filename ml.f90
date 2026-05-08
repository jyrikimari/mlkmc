! MLKMC - module for machine learning functions and subroutines
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

#include "QUIP/src/libAtoms/error.inc"

module ml_module
  use libAtoms_module, only : Inoutput
  use error_module
  use system_module, only : optional_default,system_abort,string_to_numerical
  use dictionary_module
  use extendable_str_module
  use gp_predict_module
  use gp_fit_module
  use constants_module
  use utils_module
  use geometry_module
  use descriptors_module
  use fox_wxml
  use FoX_sax, only : xml_t,dictionary_t,haskey,getvalue,parse,open_xml_file,open_xml_string,close_xml_t
  use task_manager_module
  implicit none

  type t_predictor
    type(gpSparse) :: gp_sp
    integer :: n_y,d,n_y_min,n_y_max,n_sparse_max
    double precision,allocatable :: x(:,:),y(:),y_error(:),covariance_cutoff(:)
    double precision :: delta,zeta,regularisation
    type(extendable_str) :: desc_str
    integer :: covariance_type
    logical :: initialised=.false.
    logical :: changed=.false.
    logical :: fitted=.false.
  end type t_predictor

  logical,save,private :: parse_matched_label,parse_in_predictor,parse_in_data_point
  integer,save,private :: parse_i_x
  type(t_predictor),pointer,private :: parse_predictor
  type(extendable_str),save,private :: parse_cur_data
  character(len=1024),save,private :: parse_predictor_label

contains

  subroutine initialise_gp(gp,predictor,error)
    type(gpFull),intent(out) :: gp
    type(t_predictor),intent(in) :: predictor
    integer,optional,intent(out) :: error

    integer :: internal_error
    type(descriptor) :: desc
    integer,dimension(:,:),allocatable :: permutations
    double precision,allocatable :: theta(:)

    if (present(error)) then
      error=0
    endif

    call gp_setParameters(gp,1,predictor%n_y,0,1d-5,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error: could not set gp parameters in initialise_gp()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

    call gp_setParameters(gp,1,predictor%d,predictor%n_y,0,predictor%delta,f0=0.0d0,&
                          covariance_type=predictor%covariance_type,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error: could not set gp coordinate parameters in initialise_gp()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    call gp_addDescriptor(gp,1,string(predictor%desc_str))

    call initialise(desc,string(predictor%desc_str))
    allocate(permutations(predictor%d,descriptor_n_permutations(desc)))
    call descriptor_permutations(desc,permutations)
    call gp_setPermutations(gp,1,permutations)

    if (predictor%covariance_type==COVARIANCE_DOT_PRODUCT) then
      call gp_setTheta(gp,1,zeta=predictor%zeta)
    else
      allocate(theta(predictor%d))
      theta=predictor%zeta
      call gp_setTheta(gp,1,zeta=predictor%zeta,theta=theta)
      deallocate(theta)
    endif

    call finalise(desc)
    deallocate(permutations)
  end subroutine initialise_gp

  subroutine initialise_predictor(predictor,n_y_min,n_y_max,n_sparse_max,desc_str,delta,zeta,covariance_type,regularisation,error)
    type(t_predictor),intent(inout) :: predictor
    integer,intent(in) :: n_y_min,n_y_max,n_sparse_max
    type(extendable_str),intent(in) :: desc_str
    double precision,optional,intent(in) :: delta,zeta,regularisation
    integer,optional,intent(in) :: covariance_type
    integer,optional,intent(out) :: error

    integer :: internal_error
    type(descriptor) :: desc

    if (present(error)) then
      error=0
    endif

    if (predictor%initialised) call finalise_predictor(predictor)
    predictor%fitted=.false.

    if (present(error)) then
      call initialise(desc,string(desc_str),error=internal_error)
      if (internal_error/=0) then
        write(0,"(a)") "*** Error initialising descriptor in initialise_predictor()"
        if (present(error)) then
          error=internal_error
        endif
        return
      endif
    else
      call initialise(desc,string(desc_str))
    endif

    predictor%n_y=0
    predictor%d=descriptor_dimensions(desc)
    predictor%n_y_min=n_y_min
    predictor%n_y_max=n_y_max
    predictor%n_sparse_max=n_sparse_max
    predictor%desc_str=desc_str
    allocate(predictor%x(predictor%d,ML_BUFFER_SIZE))
    allocate(predictor%y(ML_BUFFER_SIZE))
    allocate(predictor%y_error(ML_BUFFER_SIZE))
    allocate(predictor%covariance_cutoff(ML_BUFFER_SIZE))
    predictor%x=0.0d0
    predictor%y=0.0d0
    predictor%y_error=0.0d0
    predictor%covariance_cutoff=0.0d0

    if (present(delta)) then
      predictor%delta=delta
    else
      predictor%delta=2.0d0
    endif

    if (present(zeta)) then
      predictor%zeta=zeta
    else
      predictor%zeta=4.0d0
    endif

    if (present(regularisation)) then
      predictor%regularisation=regularisation
    else
      predictor%regularisation=1.0d-2
    endif

    if (present(covariance_type)) then
      predictor%covariance_type=covariance_type
    else
      predictor%covariance_type=COVARIANCE_DOT_PRODUCT
    endif

    predictor%initialised=.true.
    predictor%changed=.false.
  end subroutine initialise_predictor

  subroutine finalise_predictor(predictor)
    type(t_predictor),intent(inout) :: predictor

    if (.not. predictor%initialised) then
      return
    endif

    call finalise(predictor%gp_sp)
    predictor%n_y=0
    deallocate(predictor%x)
    deallocate(predictor%y)
    deallocate(predictor%y_error)
    deallocate(predictor%covariance_cutoff)
    predictor%initialised=.false.
    predictor%fitted=.false.
    predictor%changed=.false.

  end subroutine finalise_predictor

  subroutine add_data_point(predictor,x,y,y_error,covariance_cutoff,error)
    type(t_predictor),intent(inout) :: predictor
    double precision,intent(in) :: x(predictor%d)
    double precision,intent(in) :: y,y_error,covariance_cutoff
    integer,optional,intent(out) :: error

    integer :: n_y,internal_error

    if (present(error)) then
      error=0
    endif
    internal_error=0

    if (.not. predictor%initialised) then
      write(0,"(a)") "*** Error: predictor not initialised before add_data_point()"
      if (present(error)) then
        error=1
      endif
      return
    endif

    n_y=predictor%n_y
    if (n_y>=predictor%n_y_max) then
      write(0,"(a)") "*** Warning: no more room for data points in predictor"
      if (present(error)) then
        error=2
      endif
      return
    endif

    if (n_y>=size(predictor%y)) then
      call size_up(predictor%x,rank=RANK_SECOND,custom_buffer_size=ML_BUFFER_SIZE,error=internal_error)
      if (internal_error/=0) then
        write(0,"(a)") "*** Error sizing up predictor%x in add_data_point()"
        if (present(error)) then
          error=internal_error
        endif
        return
      endif
      call size_up(predictor%y,custom_buffer_size=ML_BUFFER_SIZE,error=internal_error)
      if (internal_error/=0) then
        write(0,"(a)") "*** Error sizing up predictor%y in add_data_point()"
        if (present(error)) then
          error=internal_error
        endif
        return
      endif
      call size_up(predictor%y_error,custom_buffer_size=ML_BUFFER_SIZE,error=internal_error)
      if (internal_error/=0) then
        write(0,"(a)") "*** Error sizing up predictor%y_error in add_data_point()"
        if (present(error)) then
          error=internal_error
        endif
        return
      endif
      call size_up(predictor%covariance_cutoff,custom_buffer_size=ML_BUFFER_SIZE,error=internal_error)
      if (internal_error/=0) then
        write(0,"(a)") "*** Error sizing up predictor%covariance_cutoff in add_data_point()"
        if (present(error)) then
          error=internal_error
        endif
        return
      endif
    endif

    n_y=n_y+1
    predictor%n_y=n_y
    predictor%x(:,n_y)=x
    predictor%y(n_y)=y
    predictor%y_error(n_y)=y_error
    predictor%covariance_cutoff(n_y)=covariance_cutoff
    predictor%changed=.true.

  end subroutine add_data_point

  logical function data_point_exists(predictor,x,y)
    type(t_predictor),intent(inout) :: predictor
    double precision,intent(in) :: x(predictor%d)
    double precision,optional,intent(out) :: y

    integer :: i
    double precision :: dist2

    data_point_exists=.false.
    if (present(y)) then
      y=0.0d0
    endif

    do i=1,predictor%n_y
      dist2=normsq(predictor%x(:,i)-x(:))
      if (dist2<DUPLICATE_TOL) then
        data_point_exists=.true.
        if (present(y)) then
          y=predictor%y(i)
        endif
        return
      endif
    end do

  end function data_point_exists

  subroutine fit_predictor(predictor,error)
    type(t_predictor),intent(inout) :: predictor
    integer,optional,intent(out) :: error

    type(gpFull) :: gp_full
    type(task_manager_type) :: task_manager
    integer :: i,current_x,current_y,n_sparse(1,1),sparseMethod,internal_error
    character(len=STRING_LENGTH) :: print_sparse_index

    if (present(error)) then
      error=0
    endif

    if (.not. predictor%initialised) then
      write(0,"(a)") "*** Error: predictor not initialised before fit_predictor()"
      if (present(error)) then
        error=1
      endif
      return
    endif

    if (predictor%n_y<=0) then
      write(0,"(a)") "*** Error: trying to fit predictor with no data in fit_predictor()"
      if (present(error)) then
        error=2
      endif
      return
    endif

    call initialise_gp(gp_full,predictor,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error initialising GP in fit_predictor"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    do i=1,predictor%n_y
      current_y=gp_addFunctionValue(gp_full,predictor%y(i),predictor%y_error(i),error=internal_error)
      if (internal_error/=0) then
        write(0,"(a,i0,a)") "*** Error adding y value ",i," in fit_predictor()"
        if (present(error)) then
          error=internal_error
        endif
        return
      endif
      current_x=gp_addCoordinates(gp_full,predictor%x(:,i),1,&
                                  cutoff_in=predictor%covariance_cutoff(i),current_y=current_y,error=internal_error)
      if (internal_error/=0) then
        write(0,"(a,i0,a)") "*** Error adding x value ",i," in fit_predictor()"
        if (present(error)) then
          error=internal_error
        endif
        return
      endif
    end do

    if (predictor%n_y<=predictor%n_sparse_max) then
      n_sparse=predictor%n_y
      sparseMethod=GP_SPARSE_NONE
    else
      n_sparse=predictor%n_sparse_max
      sparseMethod=GP_SPARSE_CUR_POINTS
    endif
    print_sparse_index=""
    call gp_sparsify(gp_full,n_sparseX=n_sparse,default_all=[.true.],&
                     sparseMethod=[sparseMethod],print_sparse_index=[print_sparse_index],&
                     use_actual_gpcov=.false.,unique_hash_tolerance=[1.0d-10],&
                     unique_descriptor_tolerance=[1.0d-10],error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error: couldn't sparsify gp_full in fit_predictor()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    call gp_covariance_sparse(gp_full,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error: couldn't solve the sparse covariance matrix in fit_predictor()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    call gpSparse_fit(predictor%gp_sp,gp_full,task_manager,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error: couldn't fit sparse gp in fit_predictor()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    call gpCoordinates_initialise_variance_estimate(predictor%gp_sp%coordinate(1),&
                                                    regularisation=predictor%regularisation,&
                                                    error=internal_error)
    if (internal_error/=0) then
      write(0,"(a,i0,a)") "*** Error: couldn't initialise variance estimate in fit_predictor()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif

    predictor%fitted=.true.
    predictor%changed=.true.

    call finalise(gp_full)

  end subroutine fit_predictor

  subroutine predict(predictor,x,mean,variance,error)
    type(t_predictor),intent(inout) :: predictor
    double precision,intent(in) :: x(predictor%d)
    double precision,intent(out) :: mean,variance
    integer,optional,intent(out) :: error

    integer :: internal_error

    if (present(error)) then
      error=0
    endif

    if (.not. predictor%fitted) then
      write(0,"(a)") "*** Error: predictor not fitted before predict()"
      if (present(error)) then
        error=1
      endif
      return
    endif
    mean=gp_predict(predictor%gp_sp%coordinate(1),x,variance_estimate=variance,&
                    do_variance_estimate=.true.,error=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error: couldn't get barrier prediction in predict()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
  end subroutine predict

  subroutine print_predictor_xml(predictor,filename,error)
    type(t_predictor),intent(inout) :: predictor
    character(len=*),intent(in) :: filename
    integer,optional,intent(out) :: error

    type(xmlf_t) :: xf
    integer :: i,internal_error

    if (present(error)) then
      error=0
    endif

    if (.not. predictor%initialised) then
      write(0,"(a)") "*** Error: predictor not initialised before print_predictor_xml()"
      if (present(error)) then
        error=1
      endif
      return
    endif

    if (.not. predictor%changed) then
      write(0,"(a)") "# Predictor hasn't changed since the last write, returning"
      return
    endif

    call xml_OpenFile(filename,xf,addDecl=.false.,iostat=internal_error)
    if (internal_error/=0) then
      write(0,"(a)") "*** Error: couldn't open file "//trim(filename)//" in print_predictor_xml()"
      if (present(error)) then
        error=internal_error
      endif
      return
    endif
    call xml_NewElement(xf,"Predictor")
    call xml_AddAttribute(xf,"n_y",""//predictor%n_y)
    call xml_AddAttribute(xf,"d",""//predictor%d)
    call xml_AddAttribute(xf,"n_y_min",""//predictor%n_y_min)
    call xml_AddAttribute(xf,"n_y_max",""//predictor%n_y_max)
    call xml_AddAttribute(xf,"n_sparse_max",""//predictor%n_sparse_max)
    call xml_AddAttribute(xf,"delta",""//predictor%delta)
    call xml_AddAttribute(xf,"zeta",""//predictor%zeta)
    call xml_AddAttribute(xf,"covariance_type",""//predictor%covariance_type)
    call xml_AddAttribute(xf,"regularisation",""//predictor%regularisation)
    call xml_AddAttribute(xf,"fitted",""//predictor%fitted)
    call xml_NewElement(xf,"desc_str")
    call xml_AddCharacters(xf,string(predictor%desc_str))
    call xml_EndElement(xf,"desc_str")
    do i=1,predictor%n_y
      call xml_NewElement(xf,"data_point")
      call xml_AddAttribute(xf,"i",""//i)
      call xml_AddAttribute(xf,"y",""//predictor%y(i))
      call xml_AddAttribute(xf,"y_error",""//predictor%y_error(i))
      call xml_AddAttribute(xf,"covariance_cutoff",""//predictor%covariance_cutoff(i))
      call xml_NewElement(xf,"x")
      call xml_AddCharacters(xf,""//predictor%x(:,i)//"")
      call xml_EndElement(xf,"x")
      call xml_EndElement(xf,"data_point")
    end do
    if (predictor%fitted) then
      call gp_printXML(predictor%gp_sp,xf)
    endif
    call xml_EndElement(xf,"Predictor")
    call xml_Close(xf)
    predictor%changed=.false.
  end subroutine print_predictor_xml

  subroutine read_predictor_xml(predictor,filename,error)
    type(t_predictor),intent(inout),target :: predictor
    character(len=*),intent(in) :: filename
    integer,optional,intent(out) :: error

    type(Inoutput) :: xmlf
    type(extendable_str) :: params_str

    integer :: internal_error

    if (present(error)) then
      error=0
    endif

    call initialise(xmlf,filename)
    call initialise(params_str)
    call read(params_str,xmlf%unit,convert_to_string=.true.)
    call finalise(xmlf)
    call predictor_readXML_string(predictor,string(params_str))
    if (predictor%fitted) then
      call gp_readXML(predictor%gp_sp,string(params_str))
      call gpCoordinates_initialise_variance_estimate(predictor%gp_sp%coordinate(1),&
                                                      regularisation=predictor%regularisation,&
                                                      error=internal_error)
      if (internal_error/=0) then
        write(0,"(a)") "*** Error initialising variance estimate in read_predictor_xml()"
        if (present(error)) then
          error=internal_error
        endif
      endif
    endif
    call finalise(params_str)
    predictor%changed=.false.

  end subroutine read_predictor_xml

  subroutine predictor_readXML_string(predictor,params_str,label,error)
    type(t_predictor),intent(inout),target :: predictor
    character(len=*),intent(in) :: params_str
    character(len=*),optional,intent(in) :: label
    integer,optional,intent(out) :: error

    type(xml_t) :: xp

    INIT_ERROR(error)

    call open_xml_string(xp,params_str)
    call predictor_readXML(predictor,xp,label,error)
    call close_xml_t(xp)

  end subroutine predictor_readXML_string

  subroutine predictor_readXML(predictor,xp,label,error)
    type(t_predictor),intent(inout), target :: predictor
    type(xml_t),intent(inout) :: xp
    character(len=*),intent(in),optional :: label
    integer,optional,intent(out) :: error

    INIT_ERROR(error)

    if (predictor%initialised) call finalise_predictor(predictor)

    parse_in_predictor=.false.
    parse_matched_label=.false.
    parse_predictor => predictor
    parse_predictor_label=optional_default("",label)

    call initialise(parse_cur_data)
    call parse(xp,&
         characters_handler=predictor_characters_handler,&
         startElement_handler=predictor_startElement_handler,&
         endElement_handler=predictor_endElement_handler)

    call finalise(parse_cur_data)

    predictor%initialised=.true.

  end subroutine predictor_readXML

  subroutine predictor_characters_handler(in)
    character(len=*),intent(in) :: in

    if(parse_in_predictor) then
       call concat(parse_cur_data,in,keep_lf=.false.,lf_to_whitespace=.true.)
    endif
  end subroutine predictor_characters_handler

  subroutine predictor_startElement_handler(URI,localname,name,attributes)
    character(len=*),intent(in)   :: URI
    character(len=*),intent(in)   :: localname
    character(len=*),intent(in)   :: name
    type(dictionary_t),intent(in) :: attributes

    integer :: stat,i

    character(len=1024) :: value

    if (name=="Predictor") then
      if (parse_in_predictor) then
        call system_abort("*** Error: predictor_startElement_handler entered predictor with parse_in_predictor true. Bug?")
      endif

      if (parse_matched_label) return

      call FoX_get_value(attributes,"label",value,stat)
      if (stat/=0) value=""

      if (len(trim(parse_predictor_label))>0) then
        if (trim(value)==trim(parse_predictor_label)) then
          parse_matched_label=.true.
          parse_in_predictor=.true.
        else
          parse_in_predictor=.false.
        endif
      else
        parse_in_predictor=.true.
      endif

      if (parse_in_predictor) then
        if (parse_predictor%initialised) call finalise_predictor(parse_predictor)

        call FoX_get_value(attributes,"n_y",value,stat)
        if (stat==0) then
          read(value,*) parse_predictor%n_y
        else
          call system_abort("*** Error: could not find the number of data points attribute 'n_y'")
        endif

        call FoX_get_value(attributes,"d",value,stat)
        if (stat==0) then
          read(value,*) parse_predictor%d
        else
          call system_abort("*** Error: could not find the dimensions attribute 'd'")
        endif

        allocate (parse_predictor%x(parse_predictor%d,parse_predictor%n_y))
        allocate (parse_predictor%y(parse_predictor%n_y))
        allocate (parse_predictor%y_error(parse_predictor%n_y))
        allocate (parse_predictor%covariance_cutoff(parse_predictor%n_y))

        call FoX_get_value(attributes,"n_y_min",value,stat)
        if (stat==0) then
          read(value,*) parse_predictor%n_y_min
        else
          call system_abort("*** Error: could not find the n_y_min attribute")
        endif

        call FoX_get_value(attributes,"n_y_max",value,stat)
        if (stat==0) then
          read(value,*) parse_predictor%n_y_max
        else
          call system_abort("*** Error: could not find the n_y_max attribute")
        endif

        call FoX_get_value(attributes,"n_sparse_max",value,stat)
        if (stat==0) then
          read(value,*) parse_predictor%n_sparse_max
        else
          call system_abort("*** Error: could not find the n_sparse_max attribute")
        endif

        call FoX_get_value(attributes,"delta",value,stat)
        if (stat==0) then
          read(value,*) parse_predictor%delta
        else
          call system_abort("*** Error: could not find the delta attribute")
        endif

        call FoX_get_value(attributes,"zeta",value,stat)
        if (stat==0) then
          read(value,*) parse_predictor%zeta
        else
          call system_abort("*** Error: could not find the zeta attribute")
        endif

        call FoX_get_value(attributes,"covariance_type",value,stat)
        if (stat==0) then
          read(value,*) parse_predictor%covariance_type
        else
          call system_abort("*** Error: could not find the covariance_type attribute")
        endif

        call FoX_get_value(attributes,"regularisation",value,stat)
        if (stat==0) then
          read(value,*) parse_predictor%regularisation
        else
          call system_abort("*** Error: could not find the regularisation attribute")
        endif

        call FoX_get_value(attributes,"fitted",value,stat)
        if (stat==0) then
          read(value,*) parse_predictor%fitted
        else
          call system_abort("*** Error: could not find the fitted attribute")
        endif

      endif

    elseif (parse_in_predictor .and. name=="desc_str") then
      call zero(parse_cur_data)

    elseif (parse_in_predictor .and. name=="data_point") then

      parse_in_data_point=.true.

      call FoX_get_value(attributes,"i",value,stat)
      if (stat==0) then
        read(value,*) i
      else
        call system_abort("*** Error: could not find the i attribute")
      endif

      call FoX_get_value(attributes,"y",value,stat)
      if (stat==0) then
        read(value,*) parse_predictor%y(i)
      else
        call system_abort("*** Error: could not find the y attribute")
      endif

      call FoX_get_value(attributes,"y_error",value,stat)
      if (stat==0) then
        read(value,*) parse_predictor%y_error(i)
      else
        call system_abort("*** Error: could not find the y_error attribute")
      endif

      call FoX_get_value(attributes,"covariance_cutoff",value,stat)
      if (stat==0) then
        read(value,*) parse_predictor%covariance_cutoff(i)
      else
        call system_abort("*** Error: could not find the covariance_cutoff attribute")
      endif

      parse_i_x=i

      call zero(parse_cur_data)

    elseif (parse_in_predictor .and. parse_in_data_point .and. name=="x") then

      call zero(parse_cur_data)

    endif

  end subroutine predictor_startElement_handler

  subroutine predictor_endElement_handler(URI,localname,name)
    character(len=*),intent(in)   :: URI
    character(len=*),intent(in)   :: localname
    character(len=*),intent(in)   :: name

    if (parse_in_predictor) then
      if (name=="Predictor") then
        parse_in_predictor=.false.
      elseif (name=="desc_str") then
        parse_predictor%desc_str=parse_cur_data
      elseif (name=="data_point") then
        if (.not. allocated(parse_predictor%x)) then
         call system_abort("*** Error: predictor_endElement_handler: x not allocated")
        endif
        if (.not. allocated(parse_predictor%y)) then
          call system_abort("*** Error: predictor_endElement_handler: y not allocated")
        endif
        if (.not. allocated(parse_predictor%y_error)) then
          call system_abort("*** Error: predictor_endElement_handler: y_error not allocated")
        endif
        if (.not. allocated(parse_predictor%covariance_cutoff)) then
          call system_abort("*** Error: predictor_endElement_handler: covariance_cutoff not allocated")
        endif
        parse_in_data_point=.false.
      elseif (name=="x") then
        call string_to_numerical(string(parse_cur_data),parse_predictor%x(:,parse_i_x))
      endif
    endif

  end subroutine predictor_endElement_handler

  subroutine FoX_get_value(attributes,key,val,stat)
     type(dictionary_t),intent(in) :: attributes
     character(len=*),intent(in) :: key
     character(len=*),intent(inout) :: val
     integer, intent(out),optional :: stat

     if (HasKey(attributes,key)) then
       val = GetValue(attributes, trim(key))
       if (present(stat)) stat = 0
     else
       val = ""
       if (present(stat)) stat = 1
     endif
  end subroutine FoX_get_value

  function print_covariance_type(covariance_type) result(cov_str)
    integer,intent(in) :: covariance_type
    character(len=80) :: cov_str

    select case(covariance_type)
      case(COVARIANCE_DOT_PRODUCT)
        cov_str="dot_product"
      case(COVARIANCE_ARD_SE)
        cov_str="ard_se"
      case default
        cov_str="Unknown! Bug?"
    end select
  end function print_covariance_type

end module ml_module
