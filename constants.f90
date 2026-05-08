! MLKMC - module defining constants used elsewhere
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

module constants_module
  use mpi
  implicit none

  ! Double precision
  integer,parameter :: my_dp=selected_real_kind(15,307)
  ! Single precision
  integer,parameter :: sp=selected_real_kind(6,37)
  ! Real kind 
  integer,parameter :: rk=my_dp
  integer,parameter :: MPI_FLOAT_TYPE=MPI_DOUBLE

  double precision,parameter :: kB=8.617333262d-5

  integer,parameter :: EVENT_VALID=0
  integer,parameter :: EVENT_EMPTY=-1
  integer,parameter :: EVENT_LAE_CHANGED=-2
  integer,parameter :: EVENT_FORBIDDEN=-3

  integer,parameter :: RANK_FIRST=1
  integer,parameter :: RANK_SECOND=2
  integer,parameter :: RANK_THIRD=3
  integer,parameter :: RANK_LAST_TWO=4
  integer,parameter :: RANK_ALL=5
  
  integer,parameter :: VERBOSITY_SILENT=0
  integer,parameter :: VERBOSITY_QUIET=1
  integer,parameter :: VERBOSITY_NORMAL=2
  integer,parameter :: VERBOSITY_LOUD=3
  integer,parameter :: VERBOSITY_ABSURD=4

  integer,parameter :: POS_INITIAL=0
  integer,parameter :: POS_FINAL=1

  integer,parameter :: EVENT_TYPE_HOP=0
  integer,parameter :: EVENT_TYPE_EXCHANGE=1
  integer,parameter :: EVENT_TYPE_DEPOSITION=2
  integer,parameter :: EVENT_TYPE_NONE=-1

  integer,parameter :: SUBSTRATE_NOT=0
  integer,parameter :: SUBSTRATE_FCC_HCP=1
  integer,parameter :: SUBSTRATE_HCP_FCC=2

  integer,parameter :: INITIAL_ATOMS_NONE=0
  integer,parameter :: INITIAL_ATOMS_RANDOM=1
  integer,parameter :: INITIAL_ATOMS_HEXAGON=2
  integer,parameter :: INITIAL_ATOMS_HEXAGON_PLUS_ADATOM=3
  integer,parameter :: INITIAL_ATOMS_TRIANGLE=4
  integer,parameter :: INITIAL_ATOMS_TRIANGLE_PLUS_ADATOM=5
  integer,parameter :: INITIAL_ATOMS_SINGLE=6
  
  integer,parameter :: SUBSTRATE_TYPE_REGULAR=0
  integer,parameter :: SUBSTRATE_TYPE_WEAK=1

  integer,parameter :: SUBSTRATE_Z=-1

  integer,parameter :: TERMINATION_NONE=0
  integer,parameter :: TERMINATION_A=1
  integer,parameter :: TERMINATION_B=2

  integer,parameter :: INACTIVE_DELETION_THRESHOLD=1000
  integer,parameter :: INACTIVE_DELETION_INTERVAL=10000

  integer,parameter :: MINIM_CACHE_SIZE=1000

  integer,parameter :: MINIM_RESULT_UNINIT=0
  integer,parameter :: MINIM_RESULT_SUCCESS=1
  integer,parameter :: MINIM_RESULT_FAIL=2

  real(kind=rk),parameter :: RATE_LOWERING_FACTOR=0.5_rk

  real(kind=rk),parameter,dimension(3) :: UNIT_CELL=[sqrt(6.0_rk)/2.0_rk,&
                                                     1.0_rk/sqrt(2.0_rk),&
                                                     3.0_rk/sqrt(3.0_rk)]

  integer :: N_NN_MAX=12
  integer :: LOCK_THRESHOLD=9

  real(kind=rk),parameter,dimension(6) :: NEIGH_DIST_FACT=[1.0_rk/sqrt(2.0_rk),1.0_rk,      sqrt(1.5_rk),&
                                                           sqrt(2.0_rk),       sqrt(2.5_rk),sqrt(3.0_rk)]

  integer,parameter :: NN_SHELLS=16
  real(kind=rk),parameter :: JUMP_DIST2_INDEX(NN_SHELLS)=[ 1.5_rk, 3.0_rk, 6.0_rk, 9.0_rk,10.5_rk,12.0_rk,&
                                                          15.0_rk,18.0_rk,19.5_rk,21.0_rk,&
                                                          22.5_rk,24.0_rk,25.5_rk,27.0_rk,31.5_rk,36.0_rk]/18.0_rk

  integer,parameter :: TRUE_1NN_SHELL=4
  integer,parameter :: TRUE_2NN_SHELL=8
  integer,parameter :: CLOSE_NEIGH_SHELL=12
  integer,parameter :: TRUE_3NN_SHELL=14
  integer,parameter :: TRUE_4NN_SHELL=NN_SHELLS

  integer,parameter :: JUMP_SHELLS(6)=[2,4,6,8,10,14]
  integer,parameter :: N_JUMP_SHELLS=size(JUMP_SHELLS)
  integer,parameter :: TRUE_1NN_JUMP_SHELL=2
  integer,parameter :: TRUE_2NN_JUMP_SHELL=4
  integer,parameter :: TRUE_3NN_JUMP_SHELL=6
  integer,parameter :: NN_SHELL_TO_JUMP_SHELL(NN_SHELLS)=[0,1,0,2,0,3,0,4,0,&
                                                          5,0,0,0,6,0,0]

  real(kind=rk),parameter,dimension(3,6) :: P_BASE=reshape([0.0_rk,          0.0_rk,                      0.0_rk,&
                                                  sqrt(6.0_rk)/4.0_rk,       1.0_rk/(2.0_rk*sqrt(2.0_rk)),0.0_rk,&
                                                  2.0_rk/3.0_rk*sqrt(1.5_rk),0.0_rk,                      1.0_rk/sqrt(3.0_rk),&
                                                  1.0_rk/6.0_rk*sqrt(1.5_rk),1.0_rk/(2.0_rk*sqrt(2.0_rk)),1.0_rk/sqrt(3.0_rk),&
                                                  sqrt(6.0_rk)/6.0_rk,       0.0_rk,                      2.0_rk/sqrt(3.0_rk),&
                                                  5.0_rk/6.0_rk*sqrt(1.5_rk),1.0_rk/(2.0_rk*sqrt(2.0_rk)),2.0_rk/sqrt(3.0_rk)],&
                                                  shape(P_BASE))

  double precision,parameter :: DUPLICATE_TOL=1.0d-10
  real(kind=rk),parameter :: DIST2_TOL=1.0e-5_rk
  real(kind=rk),parameter :: COLLISION2_TOL=1e-9_rk

  integer,parameter :: BUFFER_SIZE=100
  integer,parameter :: ML_BUFFER_SIZE=10
  integer,parameter :: NN_BUFFER_SIZE=3

end module constants_module
