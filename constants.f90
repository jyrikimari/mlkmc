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

  integer,parameter :: ORIENTATION_UNDEFINED=0
  integer,parameter :: ORIENTATION_111=1
  integer,parameter :: ORIENTATION_110=2
  integer,parameter :: ORIENTATION_100=3

  integer,parameter :: CRYSTAL_STRUCTURE_UNDEFINED=0
  integer,parameter :: CRYSTAL_STRUCTURE_CLOSE_PACKED=1 ! fcc any orientation or hcp (0001)
  integer,parameter :: CRYSTAL_STRUCTURE_BCC=2

  real(kind=rk),parameter :: RATE_LOWERING_FACTOR=0.5_rk

  ! Unit cells for different orientations and crystal structures
  real(kind=rk),parameter,dimension(3,3,2) :: UNIT_CELL=reshape([sqrt(6.0_rk)/2.0_rk,& ! fcc {111} or hcp (0001)
                                                                 1.0_rk/sqrt(2.0_rk),&
                                                                 3.0_rk/sqrt(3.0_rk),&
                                                                 1.0_rk,& ! fcc {110}
                                                                 1.0_rk/sqrt(2.0_rk),&
                                                                 1.0_rk/sqrt(2.0_rk),&
                                                                 1.0_rk,& ! fcc {100}
                                                                 1.0_rk,&
                                                                 1.0_rk,&
                                                                 sqrt(6.0_rk),& ! bcc {111}
                                                                 sqrt(2.0_rk),&
                                                                 sqrt(3.0_rk),&
                                                                 sqrt(2.0_rk),& ! bcc {110}
                                                                 1.0_rk,&
                                                                 sqrt(2.0_rk),&
                                                                 1.0_rk,& !bcc {100}
                                                                 1.0_rk,&
                                                                 1.0_rk],shape(UNIT_CELL))


  integer,dimension(2) :: N_NN_MAX=[12,8]
  integer,dimension(2) :: N_NN_MIN=[3,2]
  integer,dimension(2) :: LOCK_THRESHOLD=[9,6]

  real(kind=rk),parameter,dimension(6,2) :: NEIGH_DIST_FACT=reshape([1.0_rk/sqrt(2.0_rk), 1.0_rk,      sqrt(1.5_rk),& ! fcc/hcp
                                                                sqrt(2.0_rk),        sqrt(2.5_rk),sqrt(3.0_rk),&
                                                                0.5_rk*sqrt(3.0_rk), 1.0_rk,      sqrt(2.0_rk),& ! bcc
                                                                0.5_rk*sqrt(11.0_rk),sqrt(3.0_rk),2.0_rk],shape(NEIGH_DIST_FACT))

  integer,parameter :: NN_SHELLS=16
  real(kind=rk),parameter,dimension(2) :: JUMP_DIST2_INDEX(NN_SHELLS,2)=reshape([ 1.5_rk, 3.0_rk, 6.0_rk, 9.0_rk,10.5_rk,& ! fcc/hcp
                                                                                 12.0_rk,15.0_rk,18.0_rk,19.5_rk,21.0_rk,&
                                                                                 22.5_rk,24.0_rk,25.5_rk,27.0_rk,31.5_rk,36.0_rk,&
                                                                                 13.5_rk,18.0_rk,36.0_rk,49.5_rk,54.0_rk,& ! bcc
                                                                                 72.0_rk,& ! Unused from here on for bcc
                                                                                         72.0_rk,72.0_rk,72.0_rk,72.0_rk,&
                                                                                 72.0_rk,72.0_rk,72.0_rk,72.0_rk,72.0_rk,72.0_rk&
                                                                                 ]/18.0_rk,shape(JUMP_DIST2_INDEX))

  integer,parameter,dimension(2) :: TRUE_1NN_SHELL=[4,1]
  integer,parameter,dimension(2) :: TRUE_2NN_SHELL=[8,2]
  integer,parameter,dimension(2) :: CLOSE_NEIGH_SHELL=[12,3]
  integer,parameter,dimension(2) :: TRUE_3NN_SHELL=[14,3]
  integer,parameter,dimension(2) :: TRUE_4NN_SHELL=[NN_SHELLS,4]

  integer,parameter :: N_JUMP_SHELLS=6
  integer,parameter,dimension(N_JUMP_SHELLS,2) :: JUMP_SHELLS=reshape([2,4,6,8,10,14,&
                                                                       1,2,3,4, 5, 6],shape(JUMP_SHELLS))
  integer,parameter,dimension(2) :: TRUE_1NN_JUMP_SHELL=[2,1]
  integer,parameter,dimension(2) :: TRUE_2NN_JUMP_SHELL=[4,2]
  integer,parameter,dimension(2) :: TRUE_3NN_JUMP_SHELL=[6,3]
  integer,parameter,dimension(2) :: TRUE_4NN_JUMP_SHELL=[-1,4]
  integer,parameter,dimension(NN_SHELLS,2) :: NN_SHELL_TO_JUMP_SHELL=reshape([0,1,0,2,0,3,0,4,0,& ! fcc/hcp
                                                                              5,0,0,0,6,0,0,&
                                                                              1,2,3,4,5,6,0,0,0,& ! bcc
                                                                              0,0,0,0,0,0,0],shape(NN_SHELL_TO_JUMP_SHELL))

  integer,parameter,dimension(3,2) :: P_MAX=reshape([6,2,4,12,4,2],shape(P_MAX))
  integer,parameter,dimension(3,2) :: ATOMS_IN_LAYER=reshape([2,1,2,2,2,1],shape(ATOMS_IN_LAYER))
  real(kind=rk),parameter,dimension(3,12,3,2) :: P_BASE=reshape([0.0_rk,           0.0_rk,                      0.0_rk,& ! fcc {111}
                                             sqrt(6.0_rk)/4.0_rk,       1.0_rk/(2.0_rk*sqrt(2.0_rk)),0.0_rk,&
                                             2.0_rk/3.0_rk*sqrt(1.5_rk),0.0_rk,                      1.0_rk/sqrt(3.0_rk),&
                                             1.0_rk/6.0_rk*sqrt(1.5_rk),1.0_rk/(2.0_rk*sqrt(2.0_rk)),1.0_rk/sqrt(3.0_rk),&
                                             sqrt(6.0_rk)/6.0_rk,       0.0_rk,                      2.0_rk/sqrt(3.0_rk),&
                                             5.0_rk/6.0_rk*sqrt(1.5_rk),1.0_rk/(2.0_rk*sqrt(2.0_rk)),2.0_rk/sqrt(3.0_rk),& ! Unused from here on
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,& ! fcc {110}
                                             0.5_rk,                    0.5_rk/sqrt(2.0_rk),         0.5_rk/sqrt(2.0_rk),& ! Unused from here on
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,& ! fcc {100}
                                             0.5_rk,                    0.5_rk,                      0.0_rk,&
                                             0.5_rk,                    0.0_rk,                      0.5_rk,&
                                             0.0_rk,                    0.5_rk,                      0.5_rk,& ! Unused from here on for fcc {100}
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,& ! bcc {111}
                                             1.0_rk/2.0_rk*sqrt(6.0_rk),1.0_rk/2.0_rk*sqrt(2.0_rk),  0.0_rk,&
                                             1.0_rk/3.0_rk*sqrt(6.0_rk),0.0_rk,                      1.0_rk/6.0_rk*sqrt(3.0_rk),&
                                             5.0_rk/6.0_rk*sqrt(6.0_rk),1.0_rk/2.0_rk*sqrt(2.0_rk),  1.0_rk/6.0_rk*sqrt(3.0_rk),&
                                             2.0_rk/3.0_rk*sqrt(6.0_rk),0.0_rk,                      1.0_rk/3.0_rk*sqrt(3.0_rk),&
                                             1.0_rk/6.0_rk*sqrt(6.0_rk),1.0_rk/2.0_rk*sqrt(2.0_rk),  1.0_rk/3.0_rk*sqrt(3.0_rk),&
                                             0.0_rk,                    0.0_rk,                      1.0_rk/2.0_rk*sqrt(3.0_rk),&
                                             1.0_rk/2.0_rk*sqrt(6.0_rk),1.0_rk/2.0_rk*sqrt(2.0_rk),  1.0_rk/2.0_rk*sqrt(3.0_rk),&
                                             1.0_rk/3.0_rk*sqrt(6.0_rk),0.0_rk,                      2.0_rk/3.0_rk*sqrt(3.0_rk),&
                                             5.0_rk/6.0_rk*sqrt(6.0_rk),1.0_rk/2.0_rk*sqrt(2.0_rk),  2.0_rk/3.0_rk*sqrt(3.0_rk),&
                                             2.0_rk/3.0_rk*sqrt(6.0_rk),0.0_rk,                      5.0_rk/6.0_rk*sqrt(3.0_rk),&
                                             1.0_rk/6.0_rk*sqrt(6.0_rk),1.0_rk/2.0_rk*sqrt(2.0_rk),  5.0_rk/6.0_rk*sqrt(3.0_rk),&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,& ! bcc {110}
                                             sqrt(2.0_rk)/2.0_rk,       0.5_rk,                      0.0_rk,&
                                             0.0_rk,                    0.5_rk,                      sqrt(2.0_rk)/2.0_rk,&
                                             sqrt(2.0_rk)/2.0_rk,       0.0_rk,                      sqrt(2.0_rk)/2.0_rk,& ! Unused from here on
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,& ! bcc {100}
                                             0.5_rk,                    0.5_rk,                      0.5_rk,& ! Unused from here on
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk,&
                                             0.0_rk,                    0.0_rk,                      0.0_rk],shape(P_BASE))

  double precision,parameter :: DUPLICATE_TOL=1.0d-10
  real(kind=rk),parameter :: DIST2_TOL=1.0e-5_rk
  real(kind=rk),parameter :: COLLISION2_TOL=1e-9_rk

  integer,parameter :: BUFFER_SIZE=100
  integer,parameter :: ML_BUFFER_SIZE=10
  integer,parameter :: NN_BUFFER_SIZE=3

end module constants_module
