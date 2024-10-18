#!/usr/bin/csh

#  This program developed in FORTRAN is a Finite Element solver for linear-
#  static analyses as well as linearized stability analyses. It is inherently
#  coupled to the open-source panel method APAME for providing fluid-structure-
#  interaction capabilites.
#    
#  Copyright (C) 2024 TUD Dresden University of Technology
# 
#  This file is part of FiPPS².
# 
#  FiPPS² is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FiPPS² is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#machinename=$(hostname)

  @ anz_rechner = 1					# Zu verknuepfende Rechneranzahl
  
  @ anz_prozesse = $1					# Zu startende Prezessanzahl in GEOpS
  
  set rechner = (mlr038u:$1)       		        # Rechner die in MPI-Umgebung verknuepft 
  							# werden
  @ ii = 1
  
  echo $rechner[$ii] > machines.file
  
   @ ii++

  while ($ii <= $anz_rechner)

    echo $rechner[$ii] >> machines.file
    @ ii++
  
  end

  #mpiexec -n $anz_prozesse ./FiPPS -mat_mumps_sym 1 $*
#/usr/bin/time  mpiexec -machinefile machines.file -n $anz_prozesse ./FiPPS/FiPPS -mat_mumps_sym 1 $*
/usr/bin/time  mpiexec -f machines.file -n $anz_prozesse ./FiPPS $*
  #mpiexec -n $anz_prozesse valgrind --tool=massif --log-file=valgrind.log --massif-out-file=massif.out.fipps ./FiPPS -mat_mumps_sym 1 $*
  #mpiexec -machinefile machines.file -n $anz_prozesse valgrind --tool=memcheck --leak-check=full --show-reachable=yes --log-file=valgrind.log --track-origins=yes ./FiPPS -mat_mumps_sym 1
  #mpiexec -n $anz_prozesse valgrind --tool=memcheck --leak-check=full --show-reachable=yes --log-file=valgrind_1.log ./FiPPS #-ksp_type preonly -pc_type cholesky -pc_factor_mat_solver_package mumps



  #mpiexec -n $anz_prozesse ./FiPPS -mat_mumps_sym 1 $*
  #/usr/bin/time  mpiexec -n $anz_prozesse ./FiPPS -mat_mumps_sym 1 $*
  #/usr/bin/time  mpiexec -n $anz_prozesse ./FiPPS $*
  #/usr/bin/time  mpiexec -n $anz_prozesse ./ex10 -f0 /mnt/autofs/data1/lft/btmpg/fe-solver/FiPPS_mpi/example_small/mat.bin -matload_type sbaij -pc_type cholesky -pc_factor_mat_solver_package mumps -ksp_monitor_true_residual $*
# -rhs 0
#/usr/bin/time mpirun ./FiPPS -eps_largest_magnitude -eps_nev 10 -mat_type sbaij -st_ksp_type preonly -st_pc_type cholesky -st_pc_factor_mat_solver_package mumps -ksp_type preonly -pc_type cholesky -pc_factor_mat_solver_package mumps -mat_mumps_sym 1 -mat_mumps_icntl_4 0
#mpirun valgrind --tool=memcheck -q --num-callers=20 --log-file=valgrind.log --track-origins=yes ./FiPPS -mat_mumps_sym 1 -malloc off
#mpirun valgrind --tool=massif --log-file=valgrind.log ./FiPPS -mat_mumps_sym 1
#mpirun -info -mat_view_info -memory_info -malloc_log
#-log_trace start and end of all events: good for hanging code
#-draw_pause -1
