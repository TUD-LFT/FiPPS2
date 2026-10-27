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

#===============================================================================
# makeconfig.mk:  common configuration file for FiPPS
#===============================================================================

VERSION=2.1.1
REVISION=1

#===============================================================================

APAME_LIB   = ../ApameSolver/apame.a
PANEL2D_LIB = ../PANEL2D/panel2d.a
XFOILWRAPPER_LIB = ../XfoilWrapper/xfoilwrapper.a
FIPPS2_LIB  = libfipps.a

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

FIPPS := FiPPS
FIPPS_VERSION := include/version.include

# Fortran compiler
FC  = mpif90
# FC  = mpiifort

# Options for default and debug
FC_DEF_FLAGS = -Dmod_slepc -O3 -fPIC
FC_DBG_FLAGS = -Dmod_slepc -g -traceback -fPIC -O0 -check all -warn unused -warn uninitialized -warn interfaces -check noarg_temp_created
