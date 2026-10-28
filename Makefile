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

#-------------------------------------------------------------------------------
# Makefile for FIPPS
#-------------------------------------------------------------------------------

include makeconfig.mk

.PHONY: default debug
.DEFAULT_GOAL = default

default : FCFLAGS = $(FC_DEF_FLAGS)
debug : FCFLAGS = $(FC_DBG_FLAGS)

# Compile the default rules for each package

default debug:
	( cd ApameSolver ; $(MAKE) $@)
	( cd PANEL2D ; $(MAKE) $@)
	( cd XfoilWrapper ; $(MAKE) $@)
	( echo -n "CHARACTER(len=10), PARAMETER :: version=\"$(VERSION)-R$(REVISION)" > $(FIPPS_VERSION) )
	( echo -n "\"" >> $(FIPPS_VERSION) ) 
	( cd include ; $(MAKE) $@ )
	( cd source ; $(MAKE) flib_$@ )
	( cd outputLibrary ; $(MAKE) $@ )
	( cd source ; $(MAKE) $@ )

# Remove all files
clean::
	( cd ApameSolver ; $(MAKE) clean)
	( cd PANEL2D ; $(MAKE) clean)
	( cd XfoilWrapper ; $(MAKE) clean)
	( cd include ; $(MAKE) clean2 )
	( cd source ; $(MAKE) clean2 )
	( cd outputLibrary ; $(MAKE) clean2 )
	rm -rf $(FIPPS) $(FIPPS_VERSION)
